"""Gridded, full-ensemble PUE/WUE projections for Oregon (ensemble_computation.md).

Extends 1.2's facility-level, single-GCM PUE pipeline to every LOCA2 grid
cell in Oregon, the SSP3-7.0 GCM ensemble, and three cooling archetypes
(Lei & Masanet 2022's Cases 1, 2, 5 -- see
``climate_risk_dc.climate.cooling_archetypes`` for the physics and its
Assumptions log for departures from the published reference code). This
script's scope, per the task doc, stops at daily gridded PUE/WUE NetCDFs:
GCMs in, files out -- ensemble aggregation, climatology, and figures are a
separate downstream step that reads these files.

CLI: one required argument, ``--gcms``, e.g.::

    python run_pue_wue_gridded.py --gcms ACCESS-CM2 CESM2-LENS

Within one invocation, the given GCM(s) are processed **sequentially**, each
one fully vectorized internally (no per-cell/per-day Python loops -- see
``cooling_archetypes.py``). Running multiple GCMs concurrently is the
operator's job (Phase 3): start one ``--gcms <model>`` process per
concurrent slot, and set ``OMP_NUM_THREADS``/``MKL_NUM_THREADS``/
``OPENBLAS_NUM_THREADS`` before launching so concurrent processes don't
oversubscribe this shared 64-core machine -- this script does not set those
itself, since it has no way to know how many sibling processes will run
alongside it.

Assumptions log
----------------
- **Scenario**: fixed SSP3-7.0 (not separately selectable), per the task
  doc -- a different forcing pathway than Section 1.1's SSP5-8.5, carried
  into Limitations, not reconciled here.
- **Domain**: a simple Oregon bounding box (lon -124.6 to -116.5, lat 41.95
  to 46.3), per explicit instruction -- not the state polygon the task doc
  originally specified, and not the dissolved-watershed boundary
  ``run_heat_processing.py`` builds for that purpose. This pulls in slivers
  of Washington, Idaho, California, and Nevada near the corners; carry that
  into Limitations if this run's grid-cell count is reported as "Oregon."
- **GCM/member selection**: read from ``outputs/ensemble_pue_wue/gcm_manifest.csv``
  (built during this task's Phase 1 inventory) rather than re-scanning the
  archive per invocation or hardcoding ``r1i1p1f1`` -- that manifest already
  records the auto-detected member per model and which of the 27 candidate
  GCMs actually have SSP3-7.0 for both temperature and humidity (20 of 27).
- **Daily pairing**: Tmax+RHmin / Tmin+RHmax, averaged -- same documented
  approximation as ``run_pue_processing.py`` (1.2), for the same reason
  (LOCA2 and its humidity companion are daily-only; the source models were
  validated at hourly resolution). Applies identically to all three cases
  here.
- **Elevation-adjusted pressure** computed fresh per (GCM, period) from the
  static ``elevation.LOCA_2016-04-02.nc`` field, nearest-neighbor-aligned to
  that call's actual tasmax grid rather than assumed identical across GCMs
  -- cheap relative to the rest of the pipeline, and avoids assuming grid
  identity across models that in principle could differ slightly.
- **Humidity Oregon-bbox extraction is cached per (GCM, scenario, period)**
  under ``outputs/ensemble_pue_wue/humidity_cache/`` -- the humidity
  companion archive is CONUS-wide and unclipped (~40GB/file); this avoids
  re-reading the full CONUS grid once per archetype for the same
  (GCM, period). A small buffer (0.25 deg) beyond the analysis bbox is kept
  in the cache in case of minor grid-alignment differences from temperature
  or elevation; the final bbox clip to the exact analysis extent happens
  after loading the cache.
- **Why the humidity read is slow, and why it's read plain rather than
  dask-chunked.** Both the humidity and (pre-clipped) temperature archives
  store data with internal NetCDF chunks of ``(1, n_lat, n_lon)`` -- one
  time step, the *entire* spatial extent, zlib level 9 (confirmed via each
  file's own ``encoding["chunksizes"]``). A lat/lon bounding-box selection
  therefore gives **no** decompression savings on its own: every day's
  full-CONUS chunk must be decompressed regardless of how small the
  requested spatial slice is, since a compressed chunk can't be partially
  decompressed. A year-range selection *does* help (whole out-of-range time
  chunks are skipped at the HDF5 level without decompression -- this
  already happens implicitly for the historical file, which spans
  1950-2014 on disk but only 1985-2014 is requested). Since each time
  step's chunk is independent, decompression is in principle embarrassingly
  parallel across the time axis -- **tried and empirically rejected**:
  ``xr.open_dataset(..., chunks={"time": ...})`` under a dask
  ``scheduler="threads"`` pool measured *slower* than plain single-threaded
  ``xr.open_dataset`` on a 2-year test slice (10.5s vs. 6-7.5s), because the
  underlying HDF5/netCDF4 C library serializes access via a global lock --
  multiple Python threads calling into it don't get real concurrency, just
  overhead (confirmed via ``top -H``: one thread at 99.9% CPU, all others
  at 0%). A ``scheduler="processes"`` pool would sidestep that lock (each
  process gets its own HDF5 handle) but adds real engineering complexity
  for a read that, measured plainly, costs on the order of 100s per
  variable per 30-year period -- a few minutes, once per (GCM, period),
  then cached; not worth it given that cost. If this ever needs to be
  faster, Phase 4's own suggestion -- converting the Oregon-clipped
  extracts to Zarr with spatially-aware chunking for repeated access -- is
  the real fix, not thread/process parallelism over the existing NetCDF
  chunk layout.
- **Checkpointing**: the (GCM, archetype, period) output file is the
  resume/skip unit (Phase 3) -- an existing output file for a given
  combination is skipped entirely, so a crashed or incrementally-extended
  ``--gcms`` invocation does not recompute finished work.
- **Ocean/no-data cell masking.** The Oregon bbox above (not the state
  polygon) includes Pacific Ocean cells that LOCA2 leaves NaN (confirmed:
  ~5% of the bbox on a sample day) -- CoolProp's ``HAPropsSI`` raises
  ``ValueError`` on NaN input rather than propagating it, which crashed
  every (GCM, archetype, period) combination in the first overnight run
  before any output was written. Fixed by masking non-finite cells out
  before calling into ``archetype_fn`` and NaN-filling them back in
  (``_archetype_fn_nan_safe``) -- this is the ocean/no-data masking Phase 1
  of the task doc calls for, applied at the point of flattening rather than
  as a separate land-mask pass. Chosen over switching to the state polygon
  (which would shrink but not eliminate the no-data cell count) to keep the
  bbox-vs-polygon deviation as the only domain change from the task doc.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from climate_risk_dc.climate.heat_indices import clip_to_bbox, load_tasmax_tasmin_c, select_year_range
from climate_risk_dc.climate.loca2_io import (
    humidity_path_config,
    load_hursmax_hursmin_pct,
    log_humidity_sanity,
    resolve_loca2_file,
    temperature_path_config,
)
from climate_risk_dc.climate.cooling_archetypes import (
    Case1Config,
    Case2Config,
    Case5Config,
    check_pue_wue_plausible,
    elevation_adjusted_pressure_pa,
    load_cop_water_cooled_gp,
    pue_wue_case1,
    pue_wue_case2,
    pue_wue_case5,
)


@dataclass(frozen=True)
class GriddedRunConfig:
    """Fixed parameters for the ensemble gridded run -- see Assumptions log."""

    scenario_future: str = "ssp370"
    scenario_hist: str = "historical"
    lon_min: float = -124.6
    lon_max: float = -116.5
    lat_min: float = 41.95
    lat_max: float = 46.3
    humidity_cache_buffer_deg: float = 0.25


PERIODS: dict[str, tuple[int, int]] = {
    "historical": (1985, 2014),
    "midcentury": (2045, 2074),
    "endcentury": (2075, 2100),  # 26 years, not 30 -- see methods doc Time Horizons
}

ARCHETYPES = {
    "ae-chiller": (pue_wue_case1, Case1Config),
    "we-chiller": (pue_wue_case2, Case2Config),
    "chiller-only": (pue_wue_case5, Case5Config),
}

COP_2_RELATIVE_PATH = Path("data") / "lei_masanet_2022" / "COP_2.pkl"
MANIFEST_RELATIVE_PATH = Path("outputs") / "ensemble_pue_wue" / "gcm_manifest.csv"


def _scenario_for_period(period_name: str, cfg: GriddedRunConfig) -> str:
    return cfg.scenario_hist if period_name == "historical" else cfg.scenario_future


def load_gcm_manifest(root: Path) -> pd.DataFrame:
    path = root / MANIFEST_RELATIVE_PATH
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found -- run the Phase 1 GCM inventory scan first "
            "(records auto-detected member IDs and SSP3-7.0 availability per model)."
        )
    return pd.read_csv(path)


def resolve_gcm_member(manifest: pd.DataFrame, gcm: str) -> str:
    row = manifest[manifest["model"] == gcm]
    if row.empty:
        raise ValueError(f"'{gcm}' is not in the GCM manifest -- not one of the 27 candidate LOCA2 models.")
    row = row.iloc[0]
    if row["included_in_ensemble"] != "yes":
        raise ValueError(
            f"'{gcm}' is excluded from the ensemble ({row['exclusion_reason']}) -- "
            "not eligible for --gcms."
        )
    return str(row["member"])


def elevation_pressure_on_grid(elevation_path: Path, target: xr.DataArray, cfg: GriddedRunConfig) -> xr.DataArray:
    """Elevation-adjusted pressure (Pa), aligned to ``target``'s (lat, lon) grid.

    Nearest-neighbor-aligned rather than assumed grid-identical to
    ``target`` -- the elevation file is a different LOCA product ("LOCA",
    not "LOCA2"; see CLAUDE.md) from the temperature/humidity archives.
    """
    elev = xr.open_dataset(elevation_path)["Elevation"]
    elev = clip_to_bbox(
        elev,
        cfg.lon_min - cfg.humidity_cache_buffer_deg,
        cfg.lon_max + cfg.humidity_cache_buffer_deg,
        cfg.lat_min - cfg.humidity_cache_buffer_deg,
        cfg.lat_max + cfg.humidity_cache_buffer_deg,
    )
    elev_on_grid = elev.reindex(lat=target["lat"], lon=target["lon"], method="nearest")
    return xr.apply_ufunc(elevation_adjusted_pressure_pa, elev_on_grid)


def get_cached_humidity(
    gcm: str,
    member: str,
    scenario: str,
    period_name: str,
    start_year: int,
    end_year: int,
    humidity_root: Path,
    cache_dir: Path,
    cfg: GriddedRunConfig,
) -> tuple[xr.DataArray, xr.DataArray]:
    """Oregon(+buffer)-clipped hursmax/hursmin for one (GCM, period), cached to disk.

    The humidity companion archive is CONUS-wide and unclipped (Phase 4) --
    this is the one-time subsetting pass, reused across all three
    archetypes for the same (GCM, period) rather than re-reading the full
    CONUS grid per archetype. Read plain (not dask-chunked): measured
    slower under a dask thread pool (the HDF5/netCDF4 backend serializes
    access via a global lock, so threads add overhead without real
    concurrency -- see this script's Assumptions log) and a process pool is
    not worth the added complexity for a read that plainly costs on the
    order of 100s per variable per 30-year period.

    Written atomically (to a ``.tmp`` path, then renamed) so a killed or
    crashed run never leaves a truncated file that a later invocation's
    checkpoint check would mistake for a complete cache.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"hursmax_hursmin_{gcm}_{scenario}_{period_name}.nc"

    if cache_path.exists():
        cached = xr.open_dataset(cache_path)
        return cached["hursmax"], cached["hursmin"]

    humidity_cfg = humidity_path_config(humidity_root, model=gcm, member=member)
    hursmax_path = resolve_loca2_file(humidity_cfg, scenario, "hursmax", start_year, end_year)
    hursmin_path = resolve_loca2_file(humidity_cfg, scenario, "hursmin", start_year, end_year)
    log_humidity_sanity(hursmax_path, hursmin_path)

    hursmax_pct, hursmin_pct = load_hursmax_hursmin_pct(hursmax_path, hursmin_path)
    hursmax_pct = clip_to_bbox(
        hursmax_pct,
        cfg.lon_min - cfg.humidity_cache_buffer_deg,
        cfg.lon_max + cfg.humidity_cache_buffer_deg,
        cfg.lat_min - cfg.humidity_cache_buffer_deg,
        cfg.lat_max + cfg.humidity_cache_buffer_deg,
    )
    hursmin_pct = clip_to_bbox(
        hursmin_pct,
        cfg.lon_min - cfg.humidity_cache_buffer_deg,
        cfg.lon_max + cfg.humidity_cache_buffer_deg,
        cfg.lat_min - cfg.humidity_cache_buffer_deg,
        cfg.lat_max + cfg.humidity_cache_buffer_deg,
    )
    hursmax_pct = select_year_range(hursmax_pct, start_year, end_year)
    hursmin_pct = select_year_range(hursmin_pct, start_year, end_year)

    out = xr.Dataset({"hursmax": hursmax_pct, "hursmin": hursmin_pct})
    encoding = {v: {"zlib": True, "complevel": 4, "dtype": "float32"} for v in out.data_vars}
    print("  Reading + clipping humidity (single-threaded, this is the slow one-time pass)...")
    t0 = time.time()
    tmp_path = cache_path.with_suffix(".nc.tmp")
    out.to_netcdf(tmp_path, encoding=encoding)
    tmp_path.rename(cache_path)
    print(f"  Cached humidity subset in {time.time() - t0:.1f}s: {cache_path}")
    cached = xr.open_dataset(cache_path)
    return cached["hursmax"], cached["hursmin"]


def _archetype_fn_nan_safe(
    archetype_fn,
    t_c_flat: np.ndarray,
    rh_pct_flat: np.ndarray,
    p_atm_flat: np.ndarray,
    cop_gp,
    archetype_config,
) -> tuple[np.ndarray, np.ndarray]:
    """Run ``archetype_fn`` only on finite inputs, NaN-filling the rest.

    The Oregon analysis bbox (not the state polygon -- see this script's
    Assumptions log) includes ocean/no-data LOCA2 cells (~5% of the bbox,
    confirmed empirically), which are NaN in tasmax/tasmin/humidity. CoolProp's
    ``HAPropsSI`` raises ``ValueError`` on NaN input rather than propagating
    it, so every element must be masked out *before* calling into
    ``archetype_fn`` rather than relying on NaN to flow through -- this is
    the ocean/no-data masking Phase 1 of the task doc calls for, applied at
    the point of flattening rather than as a separate land-mask pass, since
    that's the only place common to all three archetypes.
    """
    valid = np.isfinite(t_c_flat) & np.isfinite(rh_pct_flat) & np.isfinite(p_atm_flat)
    pue_flat = np.full(t_c_flat.shape, np.nan, dtype=float)
    wue_flat = np.full(t_c_flat.shape, np.nan, dtype=float)
    if valid.any():
        pue_valid, wue_valid = archetype_fn(
            t_c_flat[valid], rh_pct_flat[valid], p_atm_flat[valid], cop_gp, archetype_config
        )
        pue_flat[valid] = pue_valid
        wue_flat[valid] = wue_valid
    return pue_flat, wue_flat


def compute_gridded_pue_wue(
    tasmax_c: xr.DataArray,
    tasmin_c: xr.DataArray,
    hursmax_pct: xr.DataArray,
    hursmin_pct: xr.DataArray,
    p_atm_pa: xr.DataArray,
    archetype_fn,
    archetype_config,
    cop_gp,
) -> tuple[np.ndarray, np.ndarray]:
    """Daily gridded PUE/WUE via the Tmax+RHmin / Tmin+RHmax pairing (1.2's convention).

    Broadcasts the static per-grid-cell ``p_atm_pa`` across the time axis by
    flattening; returns arrays shaped ``(time, lat, lon)``, NaN at any
    ocean/no-data cell (see ``_archetype_fn_nan_safe``).
    """
    shape = tasmax_c.shape  # (time, lat, lon)
    p_atm_flat = np.broadcast_to(p_atm_pa.values[np.newaxis, :, :], shape).ravel()

    pue_hot, wue_hot = _archetype_fn_nan_safe(
        archetype_fn, tasmax_c.values.ravel(), hursmin_pct.values.ravel(), p_atm_flat, cop_gp, archetype_config
    )
    pue_cool, wue_cool = _archetype_fn_nan_safe(
        archetype_fn, tasmin_c.values.ravel(), hursmax_pct.values.ravel(), p_atm_flat, cop_gp, archetype_config
    )
    pue = ((pue_hot + pue_cool) / 2.0).reshape(shape)
    wue = ((wue_hot + wue_cool) / 2.0).reshape(shape)
    return pue, wue


def run_one_combination(
    gcm: str,
    member: str,
    archetype_name: str,
    period_name: str,
    root: Path,
    args: argparse.Namespace,
    cfg: GriddedRunConfig,
    cop_gp,
    elevation_path: Path,
    out_dir: Path,
    log_path: Path,
) -> float | None:
    """Run one (GCM, archetype, period) combination; return elapsed seconds, or None if skipped."""
    output_path = out_dir / f"pue_wue_{gcm}_{cfg.scenario_future}_{archetype_name}_{period_name}.nc"
    if output_path.exists():
        print(f"  [skip, exists] {output_path.name}")
        return None

    t0 = time.time()
    start_year, end_year = PERIODS[period_name]
    scenario = _scenario_for_period(period_name, cfg)

    temp_cfg = temperature_path_config(Path(args.loca2_root), model=gcm, member=member)
    tasmax_path = resolve_loca2_file(temp_cfg, scenario, "tasmax", start_year, end_year)
    tasmin_path = resolve_loca2_file(temp_cfg, scenario, "tasmin", start_year, end_year)
    tasmax_c, tasmin_c = load_tasmax_tasmin_c(tasmax_path, tasmin_path)
    tasmax_c = clip_to_bbox(tasmax_c, cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
    tasmin_c = clip_to_bbox(tasmin_c, cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
    tasmax_c = select_year_range(tasmax_c, start_year, end_year)
    tasmin_c = select_year_range(tasmin_c, start_year, end_year)

    hursmax_pct, hursmin_pct = get_cached_humidity(
        gcm,
        member,
        scenario,
        period_name,
        start_year,
        end_year,
        Path(args.loca2_humidity_root),
        out_dir / "humidity_cache",
        cfg,
    )
    hursmax_pct = clip_to_bbox(hursmax_pct, cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
    hursmin_pct = clip_to_bbox(hursmin_pct, cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
    hursmax_pct = hursmax_pct.reindex(lat=tasmax_c["lat"], lon=tasmax_c["lon"], method="nearest")
    hursmin_pct = hursmin_pct.reindex(lat=tasmax_c["lat"], lon=tasmax_c["lon"], method="nearest")

    if tasmax_c.sizes["time"] != hursmax_pct.sizes["time"]:
        raise ValueError(
            f"{gcm}/{period_name}: temperature has {tasmax_c.sizes['time']} days, humidity has "
            f"{hursmax_pct.sizes['time']} -- expected matching calendars (humidity is computed "
            "directly from this GCM's own downscaled Tmax/Tmin)."
        )

    p_atm_pa = elevation_pressure_on_grid(elevation_path, tasmax_c, cfg)

    archetype_fn, config_cls = ARCHETYPES[archetype_name]
    pue, wue = compute_gridded_pue_wue(
        tasmax_c, tasmin_c, hursmax_pct, hursmin_pct, p_atm_pa, archetype_fn, config_cls(), cop_gp
    )

    try:
        check_pue_wue_plausible(pue, wue)
    except ValueError as exc:
        print(f"  Sanity check WARNING for {gcm}/{archetype_name}/{period_name}: {exc}")

    out = xr.Dataset(
        {
            "pue": (("time", "lat", "lon"), pue.astype("float32")),
            "wue": (("time", "lat", "lon"), wue.astype("float32")),
        },
        coords={"time": tasmax_c["time"], "lat": tasmax_c["lat"], "lon": tasmax_c["lon"]},
        attrs={
            "gcm": gcm,
            "member": member,
            "scenario": scenario,
            "archetype": archetype_name,
            "period": period_name,
            "period_years": f"{start_year}-{end_year}",
            "period_length_note": (
                "26 years, not 30 -- see methods doc Time Horizons" if period_name == "endcentury" else "30 years"
            ),
            "daily_pairing": "Tmax+RHmin and Tmin+RHmax, averaged",
            "wue_units": "L/kWh",
        },
    )
    encoding = {v: {"zlib": True, "complevel": 4, "dtype": "float32"} for v in out.data_vars}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_output_path = output_path.with_suffix(".nc.tmp")
    out.to_netcdf(tmp_output_path, encoding=encoding)
    tmp_output_path.rename(output_path)  # atomic-ish: a killed run never leaves a file the checkpoint check treats as done

    elapsed = time.time() - t0
    _append_log(
        log_path,
        {
            "gcm": gcm,
            "member": member,
            "archetype": archetype_name,
            "period": period_name,
            "elapsed_s": round(elapsed, 2),
            "output_path": str(output_path),
            "timestamp": pd.Timestamp.utcnow().isoformat(),
        },
    )
    print(f"  [{elapsed:6.1f}s] {output_path.name}")
    return elapsed


def _append_log(log_path: Path, record: dict) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a") as f:
        f.write(json.dumps(record) + "\n")


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    cfg = GriddedRunConfig()
    out_dir = root / "outputs" / "ensemble_pue_wue"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "run_log.jsonl"
    elevation_path = root / "data" / "elevation.LOCA_2016-04-02.nc"

    manifest = load_gcm_manifest(root)
    cop_gp = load_cop_water_cooled_gp(root / COP_2_RELATIVE_PATH)

    for gcm in args.gcms:
        member = resolve_gcm_member(manifest, gcm)
        print(f"=== {gcm} (member {member}) ===")
        gcm_t0 = time.time()
        for archetype_name in ARCHETYPES:
            for period_name in PERIODS:
                run_one_combination(
                    gcm, member, archetype_name, period_name, root, args, cfg, cop_gp, elevation_path, out_dir, log_path
                )
        gcm_elapsed = time.time() - gcm_t0
        _append_log(log_path, {"gcm": gcm, "gcm_total_elapsed_s": round(gcm_elapsed, 2), "timestamp": pd.Timestamp.utcnow().isoformat()})
        print(f"=== {gcm} total: {gcm_elapsed:.1f}s ===")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Gridded Oregon PUE/WUE for one or more LOCA2 GCMs (SSP3-7.0).")
    parser.add_argument(
        "--gcms",
        nargs="+",
        required=True,
        help="One or more GCM names to process sequentially (e.g. --gcms ACCESS-CM2 CESM2-LENS). "
        "Must be in outputs/ensemble_pue_wue/gcm_manifest.csv with included_in_ensemble=yes.",
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
