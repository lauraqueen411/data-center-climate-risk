"""Facility-point (109 Oregon sites) ensemble PUE/WUE -- fast-path stopgap.

Same physics as ``run_pue_wue_gridded.py`` (the three Lei & Masanet 2022
archetypes in ``cooling_archetypes.py``, SSP3-7.0, three periods), but
sampled at the 109 IM3/PNNL Oregon facility points instead of every LOCA2
grid cell in the Oregon bbox. Exists because the full-grid run is currently
compute-infeasible: ``CoolProp.HAPropsSI``'s array support is a per-element
loop internally (measured ~15.6us/point, not true vectorization), so cost
scales linearly with point count -- ~9,100 grid cells means ~5 days
single-threaded per GCM; 109 facility points means ~87 minutes. See the
conversation/task-doc discussion in ``docs/tasks/ensemble_computation.md``'s
Phase 0 report for the full arithmetic.

This is a deliberate, temporary scope reduction (facility points, not full
Oregon coverage), not the task's actual deliverable -- it exists to (1) get
real ensemble PUE/WUE numbers now rather than waiting on a rewrite, and
(2) provide real, CoolProp-computed ground truth to validate a
vectorized-numpy psychrometrics replacement against (see
``cooling_archetypes_fast.py``, once built), before trusting that
replacement at full-grid scale.

Output: one CSV per (GCM, archetype, period), written and checkpointed
(skip-if-exists) as each combination finishes, so progress is visible
incrementally rather than only at the end -- mirrors
``run_pue_wue_gridded.py``'s per-file checkpointing philosophy at facility
scale. Columns match ``run_pue_processing.py``'s (1.2's) daily-timeseries
convention: per facility, per day, both raw weather inputs and both daily
pairings' PUE/WUE plus their mean.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from climate_risk_dc.climate.cooling_archetypes import (
    Case1Config,
    Case2Config,
    Case5Config,
    check_pue_wue_plausible,
    elevation_adjusted_pressure_pa,
    load_cop_water_cooled_gp,
)
from climate_risk_dc.climate.heat_indices import clip_to_bbox, load_tasmax_tasmin_c, select_year_range
from climate_risk_dc.climate.loca2_io import (
    resolve_loca2_file,
    sample_datacube_at_points,
    temperature_path_config,
)
from climate_risk_dc.geo import load_oregon_data_centers
from run_pue_wue_gridded import (
    COP_2_RELATIVE_PATH,
    PERIODS,
    GriddedRunConfig,
    _scenario_for_period,
    get_cached_humidity,
    load_gcm_manifest,
    resolve_gcm_member,
)

def _archetypes_for_backend(backend: str) -> dict:
    """CoolProp (slow, ground-truth) vs. numpy-vectorized (fast) archetype functions.

    Same ``Case1Config``/``Case2Config``/``Case5Config`` either way -- only
    the psychrometrics implementation differs (see
    ``cooling_archetypes_fast.py``'s Assumptions log for the validated
    accuracy tradeoff, ~50x+ faster, PUE within ~0.1-0.35% and WUE within
    ~0.03 L/kWh of the CoolProp path on random test conditions).
    """
    if backend == "fast":
        from climate_risk_dc.climate.cooling_archetypes_fast import pue_wue_case1, pue_wue_case2, pue_wue_case5
    else:
        from climate_risk_dc.climate.cooling_archetypes import pue_wue_case1, pue_wue_case2, pue_wue_case5
    return {
        "ae-chiller": (pue_wue_case1, Case1Config),
        "we-chiller": (pue_wue_case2, Case2Config),
        "chiller-only": (pue_wue_case5, Case5Config),
    }


def sample_elevation_at_points(elevation_path: Path, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Nearest-neighbor elevation (m) at facility points from the static LOCA elevation field.

    Not ``geo.get_facility_elevation_m`` -- that wraps ``rasterio`` for
    GDAL-readable rasters (USGS 3DEP GeoTIFF, per its own docstring); the
    LOCA elevation field is a NetCDF ``(lat, lon)`` DataArray, so this
    mirrors ``loca2_io.sample_datacube_at_points``'s nearest-neighbor
    approach instead, for a field with no time dimension.
    """
    elev = xr.open_dataset(elevation_path)["Elevation"]
    lon_vals = np.asarray(elev["lon"].values)
    if lon_vals.max() > 180.0:
        lon_query = np.where(lon < 0, lon + 360.0, lon)
    else:
        lon_query = lon
    lat_vals = np.asarray(elev["lat"].values)
    lat_idx = np.abs(lat_vals[:, None] - lat[None, :]).argmin(axis=0)
    lon_idx = np.abs(lon_vals[:, None] - lon_query[None, :]).argmin(axis=0)
    return elev.values[lat_idx, lon_idx]


def run_one_combination(
    gcm: str,
    member: str,
    archetype_name: str,
    period_name: str,
    root: Path,
    args: argparse.Namespace,
    cfg: GriddedRunConfig,
    cop_gp,
    facilities: pd.DataFrame,
    p_atm_pa: np.ndarray,
    out_dir: Path,
    archetypes: dict,
) -> None:
    output_path = out_dir / f"facility_pue_wue_{gcm}_{cfg.scenario_future}_{archetype_name}_{period_name}_{args.backend}.csv"
    if output_path.exists():
        print(f"  [skip, exists] {output_path.name}")
        return

    t0 = time.time()
    start_year, end_year = PERIODS[period_name]
    scenario = _scenario_for_period(period_name, cfg)
    lon = facilities["lon"].to_numpy()
    lat = facilities["lat"].to_numpy()

    temp_cfg = temperature_path_config(Path(args.loca2_root), model=gcm, member=member)
    tasmax_path = resolve_loca2_file(temp_cfg, scenario, "tasmax", start_year, end_year)
    tasmin_path = resolve_loca2_file(temp_cfg, scenario, "tasmin", start_year, end_year)
    tasmax_c, tasmin_c = load_tasmax_tasmin_c(tasmax_path, tasmin_path)
    tasmax_c = clip_to_bbox(tasmax_c, cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
    tasmin_c = clip_to_bbox(tasmin_c, cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
    tasmax_c = select_year_range(tasmax_c, start_year, end_year)
    tasmin_c = select_year_range(tasmin_c, start_year, end_year)

    # Shared cache dir with run_pue_wue_gridded.py: same Oregon(+buffer)-clipped
    # humidity subset serves both scripts, so whichever runs first builds it.
    hursmax_pct, hursmin_pct = get_cached_humidity(
        gcm, member, scenario, period_name, start_year, end_year,
        Path(args.loca2_humidity_root), root / "outputs" / "ensemble_pue_wue" / "humidity_cache", cfg,
    )

    tasmax_fac = sample_datacube_at_points(tasmax_c, lon, lat)
    tasmin_fac = sample_datacube_at_points(tasmin_c, lon, lat)
    hursmax_fac = sample_datacube_at_points(hursmax_pct, lon, lat)
    hursmin_fac = sample_datacube_at_points(hursmin_pct, lon, lat)
    dates = tasmax_c["time"].values
    n_facility, n_time = tasmax_fac.shape

    p_atm_flat = np.broadcast_to(p_atm_pa[:, np.newaxis], (n_facility, n_time)).ravel()

    archetype_fn, config_cls = archetypes[archetype_name]
    config = config_cls()
    pue_hot, wue_hot = archetype_fn(tasmax_fac.ravel(), hursmin_fac.ravel(), p_atm_flat, cop_gp, config)
    pue_cool, wue_cool = archetype_fn(tasmin_fac.ravel(), hursmax_fac.ravel(), p_atm_flat, cop_gp, config)
    pue_daily = ((pue_hot + pue_cool) / 2.0).reshape(n_facility, n_time)
    wue_daily = ((wue_hot + wue_cool) / 2.0).reshape(n_facility, n_time)

    try:
        check_pue_wue_plausible(pue_daily, wue_daily)
    except ValueError as exc:
        print(f"  Sanity check WARNING for {gcm}/{archetype_name}/{period_name}: {exc}")

    df = pd.DataFrame(
        {
            "facility_id": np.repeat(facilities["id"].to_numpy(), n_time),
            "facility_name": np.repeat(facilities["name"].to_numpy(), n_time),
            "date": np.tile(dates, n_facility),
            "tasmax_c": tasmax_fac.ravel(),
            "tasmin_c": tasmin_fac.ravel(),
            "hursmax_pct": hursmax_fac.ravel(),
            "hursmin_pct": hursmin_fac.ravel(),
            "pue_hot_pairing": pue_hot,
            "pue_cool_pairing": pue_cool,
            "pue_daily": pue_daily.ravel(),
            "wue_hot_pairing": wue_hot,
            "wue_cool_pairing": wue_cool,
            "wue_daily": wue_daily.ravel(),
        }
    )
    tmp_path = output_path.with_suffix(".csv.tmp")
    df.to_csv(tmp_path, index=False)
    tmp_path.rename(output_path)

    elapsed = time.time() - t0
    print(
        f"  [{elapsed:6.1f}s] {output_path.name}  "
        f"(PUE mean={pue_daily.mean():.4f}, WUE mean={wue_daily.mean():.4f} L/kWh)"
    )


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    cfg = GriddedRunConfig()
    out_dir = root / "outputs" / "ensemble_pue_wue_facilities"
    out_dir.mkdir(parents=True, exist_ok=True)
    elevation_path = root / "data" / "elevation.LOCA_2016-04-02.nc"
    dc_csv = root / "data" / "data_centers_im3_pnnl" / "im3_open_source_data_center_atlas_v2026.02.09.csv"

    manifest = load_gcm_manifest(root)
    cop_gp = load_cop_water_cooled_gp(root / COP_2_RELATIVE_PATH)
    facilities = load_oregon_data_centers(dc_csv)
    print(f"Loaded {len(facilities)} Oregon facilities")

    lon = facilities["lon"].to_numpy()
    lat = facilities["lat"].to_numpy()
    elevation_m = sample_elevation_at_points(elevation_path, lon, lat)
    p_atm_pa = elevation_adjusted_pressure_pa(elevation_m)
    print(
        f"Facility elevation: min={elevation_m.min():.0f}m max={elevation_m.max():.0f}m "
        f"-> pressure min={p_atm_pa.min():.0f}Pa max={p_atm_pa.max():.0f}Pa"
    )

    archetypes = _archetypes_for_backend(args.backend)
    print(f"Backend: {args.backend}")

    for gcm in args.gcms:
        member = resolve_gcm_member(manifest, gcm)
        print(f"=== {gcm} (member {member}) ===")
        gcm_t0 = time.time()
        for archetype_name in archetypes:
            for period_name in PERIODS:
                run_one_combination(
                    gcm, member, archetype_name, period_name, root, args, cfg, cop_gp, facilities, p_atm_pa, out_dir,
                    archetypes,
                )
        print(f"=== {gcm} total: {time.time() - gcm_t0:.1f}s ===")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Facility-point (109 sites) ensemble PUE/WUE.")
    parser.add_argument("--gcms", nargs="+", required=True, help="One or more GCM names to process sequentially.")
    parser.add_argument(
        "--backend",
        choices=["fast", "coolprop"],
        default="fast",
        help="'fast' = vectorized numpy psychrometrics (default, ~50x+ faster, validated against CoolProp to "
        "within ~0.1-0.35%% PUE / ~0.03 L/kWh WUE). 'coolprop' = the original reference-repo psychrometrics -- "
        "confirmed impractically slow at this scale (HAPropsSI's array support degrades badly, worse than "
        "linear, above ~100K-point arrays), kept only for point comparisons, not full runs.",
    )
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    parser.add_argument(
        "--loca2-root",
        default="/home/redmond/data/pub/scripps_downscaled_CMIP6/LOCA2/CONUS_regions_split",
        help="Root path for LOCA2 CONUS_regions_split temperature files",
    )
    parser.add_argument(
        "--loca2-humidity-root",
        default="/home/redmond/data/pub/scripps_downscaled_CMIP6/LOCA2_humidity",
        help="Root path for the LOCA CMIP6 humidity companion files",
    )
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
