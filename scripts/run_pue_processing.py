"""Deliverable 1.2: Chronic heat risk (PUE) processing for Oregon data centers.

Assumptions log
----------------
- Three periods, matching 1.1's corrected structure (see
  ``docs/tasks/0.2-fix-loca2-window-resolver-and-time-horizons.md``):
  historical (1985-2014), mid-century (2045-2074, primary), end-of-century
  (2075-2100, conservative secondary). Each maps to exactly one native
  LOCA2/humidity file via ``climate_risk_dc.climate.loca2_io.resolve_loca2_file``.
- **Daily pairing is a genuine step down from the model's validated temporal
  resolution -- flagged prominently, not a minor detail.** Lei & Masanet's
  model (Section 3.4 of the paper) was built and validated against *hourly*
  (T, RH) pairs. LOCA2 and its humidity companion only provide daily
  tasmax/tasmin + hursmax/hursmin, not hourly series, for either historical
  or projected periods -- switching the humidity source to the LOCA CMIP6
  companion (see ``loca2_io.py``) fixed the spatial/GCM-lineage consistency
  problem (temperature and humidity now come from the same downscaling
  lineage) but does **not** fix this temporal-resolution gap. This script
  approximates a representative daily PUE by pairing Tmax with RHmin and
  Tmin with RHmax (typical diurnal anti-correlation: hot afternoons are
  drier, cool nights are more humid), computing PUE for both pairings, and
  averaging the two for a per-day estimate, then averaging across all days
  in the period for the facility-level result. This is an explicit,
  documented approximation, not a validated substitute for hourly
  simulation -- treat resulting PUE trends as directionally informative,
  not hourly-accurate.
- Humidity companion files are CONUS-wide and unclipped, unlike the
  pre-clipped LOCA2 temperature archive (a ~40GB historical file vs. ~3GB
  for temperature) -- confirmed via inspection, 2026-08-21. Materializing
  the full grid was measured to hang; this script clips to the same
  Oregon-ish bounding box as ``run_heat_processing.py`` immediately after
  opening (via ``heat_indices.clip_to_bbox``), before extracting per-facility
  time series with ``loca2_io.sample_datacube_at_points`` (see that
  function's docstring for why plain numpy nearest-neighbor indexing is used
  instead of xarray's vectorized ``.sel()`` -- the latter was also measured
  to be impractically slow against this backend).
- ``PUEModelConfig`` defaults (Table A.1 midpoints, see ``climate/pue.py``)
  represent a "typical hyperscale AE+adiabatic system," not any specific
  Oregon facility's real equipment. Same generic-equipment caveat as the
  1.1 elevation-adjustment placeholder: ``p_atm_pa`` is fixed at standard
  atmosphere (101325 Pa) for every facility -- no DEM exists in the repo yet
  to apply a real per-facility elevation adjustment.
- No raster/map output, per the methods doc -- facility-level table only.
- In addition to the facility-level period-mean summary table, this script
  writes the full daily per-facility time series
  (``outputs/month2_pue/facility_pue_timeseries.csv``: date, both raw
  Tmax/Tmin/RHmax/RHmin inputs, both daily pairings' PUE, and their daily
  mean) so the expensive per-day model evaluation (the dominant cost of this
  script -- an implicit wetbulb solve per facility-day) doesn't need to be
  rerun to answer a question the period-mean table can't, e.g. inspecting a
  specific hot day or recomputing the period mean with a different
  aggregation.
"""

from __future__ import annotations

import argparse
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
    sample_datacube_at_points,
    temperature_path_config,
)
from climate_risk_dc.climate.pue import PUEModelConfig, check_pue_plausible, pue_ae_chiller
from climate_risk_dc.geo import load_oregon_data_centers


@dataclass(frozen=True)
class RunConfig:
    """Configuration for PUE processing run."""

    model: str = "ACCESS-CM2"
    member: str = "r1i1p1f1"
    scenario_hist: str = "historical"
    scenario_fut: str = "ssp585"
    hist_start: int = 1985
    hist_end: int = 2014
    mid_start: int = 2045
    mid_end: int = 2074
    eoc_start: int = 2075
    eoc_end: int = 2100
    lon_min: float = -125.0
    lon_max: float = -116.0
    lat_min: float = 42.0
    lat_max: float = 46.5


def _first_finite_kelvin(da: xr.DataArray) -> float:
    if "time" in da.dims:
        da = da.isel(time=0)
    arr = np.asarray(da.values, dtype=float).ravel()
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        raise ValueError("No finite values found for temperature sanity check.")
    return float(finite[0])


def _log_temperature_sanity(tasmax_path: Path, tasmin_path: Path) -> None:
    """Kelvin-bounds sanity check, mirroring run_heat_processing.py's ``_log_unit_sanity``."""
    dsx = xr.open_dataset(tasmax_path)
    dsn = xr.open_dataset(tasmin_path)
    tasmax_k = dsx["tasmax"] if "tasmax" in dsx.data_vars else next(iter(dsx.data_vars.values()))
    tasmin_k = dsn["tasmin"] if "tasmin" in dsn.data_vars else next(iter(dsn.data_vars.values()))

    sample_max = _first_finite_kelvin(tasmax_k)
    sample_min = _first_finite_kelvin(tasmin_k)
    print(f"Kelvin sample tasmax={sample_max:.3f}, tasmin={sample_min:.3f}")
    for v in (sample_max, sample_min):
        if v < 180.0 or v > 330.0:
            raise ValueError(f"Temperature sample {v:.3f} is outside plausible Kelvin bounds [180, 330].")


def _daily_facility_pue(
    tasmax_c: xr.DataArray,
    tasmin_c: xr.DataArray,
    hursmax_pct: xr.DataArray,
    hursmin_pct: xr.DataArray,
    lon: np.ndarray,
    lat: np.ndarray,
    pue_cfg: PUEModelConfig,
) -> dict[str, np.ndarray]:
    """Per-facility daily PUE (both pairings + their mean) plus the underlying weather inputs.

    Returns arrays of shape ``(n_facility, n_time)`` (except ``dates``, shape
    ``(n_time,)``) so the caller can both persist the full daily series and
    reduce it to a period mean without recomputing the model.
    """
    tasmax_fac = sample_datacube_at_points(tasmax_c, lon, lat)
    tasmin_fac = sample_datacube_at_points(tasmin_c, lon, lat)
    hursmax_fac = sample_datacube_at_points(hursmax_pct, lon, lat)
    hursmin_fac = sample_datacube_at_points(hursmin_pct, lon, lat)

    if not (tasmax_fac.shape == tasmin_fac.shape == hursmax_fac.shape == hursmin_fac.shape):
        raise ValueError(
            "Temperature and humidity facility time series have mismatched shapes "
            f"(tasmax={tasmax_fac.shape}, tasmin={tasmin_fac.shape}, "
            f"hursmax={hursmax_fac.shape}, hursmin={hursmin_fac.shape}) -- expected the same "
            "calendar for both, since the humidity companion is computed directly from LOCA2's "
            "own downscaled Tmax/Tmin."
        )

    shape = tasmax_fac.shape
    pue_hot_pairing = pue_ae_chiller(
        tasmax_fac.ravel(), hursmin_fac.ravel(), pue_cfg.p_atm_pa, pue_cfg.power_it_kw, pue_cfg
    ).reshape(shape)
    pue_cool_pairing = pue_ae_chiller(
        tasmin_fac.ravel(), hursmax_fac.ravel(), pue_cfg.p_atm_pa, pue_cfg.power_it_kw, pue_cfg
    ).reshape(shape)
    pue_daily = (pue_hot_pairing + pue_cool_pairing) / 2.0

    return {
        "dates": tasmax_c["time"].values,
        "tasmax_c": tasmax_fac,
        "tasmin_c": tasmin_fac,
        "hursmax_pct": hursmax_fac,
        "hursmin_pct": hursmin_fac,
        "pue_hot_pairing": pue_hot_pairing,
        "pue_cool_pairing": pue_cool_pairing,
        "pue_daily": pue_daily,
    }


def compute_facility_pue(root: Path, args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (facility-level period-mean summary, full daily per-facility time series)."""
    cfg = RunConfig()
    pue_cfg = PUEModelConfig()

    dc_csv = root / "data" / "data_centers_im3_pnnl" / "im3_open_source_data_center_atlas_v2026.02.09.csv"
    dc = load_oregon_data_centers(dc_csv)
    lon = dc["lon"].to_numpy()
    lat = dc["lat"].to_numpy()
    facility_id = dc["id"].to_numpy()
    facility_name = dc["name"].to_numpy()

    temp_cfg = temperature_path_config(Path(args.loca2_root), model=cfg.model, member=cfg.member)
    humidity_cfg = humidity_path_config(Path(args.loca2_humidity_root), model=cfg.model, member=cfg.member)

    periods = [
        ("historical_1985_2014", cfg.scenario_hist, cfg.hist_start, cfg.hist_end),
        ("ssp585_2045_2074", cfg.scenario_fut, cfg.mid_start, cfg.mid_end),
        ("ssp585_2075_2100", cfg.scenario_fut, cfg.eoc_start, cfg.eoc_end),
    ]

    period_pue: dict[str, np.ndarray] = {}
    timeseries_frames: list[pd.DataFrame] = []
    for label, scenario, start_year, end_year in periods:
        tasmax_path = resolve_loca2_file(temp_cfg, scenario, "tasmax", start_year, end_year)
        tasmin_path = resolve_loca2_file(temp_cfg, scenario, "tasmin", start_year, end_year)
        hursmax_path = resolve_loca2_file(humidity_cfg, scenario, "hursmax", start_year, end_year)
        hursmin_path = resolve_loca2_file(humidity_cfg, scenario, "hursmin", start_year, end_year)
        print(
            f"Using {label}:\n  tasmax={tasmax_path}\n  tasmin={tasmin_path}\n"
            f"  hursmax={hursmax_path}\n  hursmin={hursmin_path}"
        )
        _log_temperature_sanity(tasmax_path, tasmin_path)
        log_humidity_sanity(hursmax_path, hursmin_path)

        tasmax_c, tasmin_c = load_tasmax_tasmin_c(tasmax_path, tasmin_path)
        tasmax_c = clip_to_bbox(tasmax_c, cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
        tasmin_c = clip_to_bbox(tasmin_c, cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
        tasmax_c = select_year_range(tasmax_c, start_year, end_year)
        tasmin_c = select_year_range(tasmin_c, start_year, end_year)

        hursmax_pct, hursmin_pct = load_hursmax_hursmin_pct(hursmax_path, hursmin_path)
        hursmax_pct = clip_to_bbox(hursmax_pct, cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
        hursmin_pct = clip_to_bbox(hursmin_pct, cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
        hursmax_pct = select_year_range(hursmax_pct, start_year, end_year)
        hursmin_pct = select_year_range(hursmin_pct, start_year, end_year)

        daily = _daily_facility_pue(tasmax_c, tasmin_c, hursmax_pct, hursmin_pct, lon, lat, pue_cfg)
        period_pue[label] = daily["pue_daily"].mean(axis=1)
        print(f"{label}: facility PUE mean={period_pue[label].mean():.4f}, max={period_pue[label].max():.4f}")

        n_facility, n_time = daily["pue_daily"].shape
        timeseries_frames.append(
            pd.DataFrame(
                {
                    "facility_id": np.repeat(facility_id, n_time),
                    "facility_name": np.repeat(facility_name, n_time),
                    "period": label,
                    "date": np.tile(daily["dates"], n_facility),
                    "tasmax_c": daily["tasmax_c"].ravel(),
                    "tasmin_c": daily["tasmin_c"].ravel(),
                    "hursmax_pct": daily["hursmax_pct"].ravel(),
                    "hursmin_pct": daily["hursmin_pct"].ravel(),
                    "pue_hot_pairing": daily["pue_hot_pairing"].ravel(),
                    "pue_cool_pairing": daily["pue_cool_pairing"].ravel(),
                    "pue_daily": daily["pue_daily"].ravel(),
                }
            )
        )

    result = pd.DataFrame(
        {
            "facility_id": facility_id,
            "facility_name": facility_name,
            "facility_lon": lon,
            "facility_lat": lat,
            "pue_historical_1985_2014": period_pue["historical_1985_2014"],
            "pue_ssp585_2045_2074": period_pue["ssp585_2045_2074"],
            "pue_ssp585_2075_2100": period_pue["ssp585_2075_2100"],
        }
    )
    result["pue_mid_century_delta"] = result["pue_ssp585_2045_2074"] - result["pue_historical_1985_2014"]
    result["pue_end_of_century_delta"] = result["pue_ssp585_2075_2100"] - result["pue_historical_1985_2014"]

    timeseries = pd.concat(timeseries_frames, ignore_index=True)
    return result, timeseries


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out_dir = root / "outputs" / "month2_pue"
    out_dir.mkdir(parents=True, exist_ok=True)

    result, timeseries = compute_facility_pue(root, args)

    out_path = out_dir / "facility_pue.csv"
    result.to_csv(out_path, index=False)
    print(f"Wrote {len(result)} facility rows to {out_path}")

    timeseries_path = out_dir / "facility_pue_timeseries.csv"
    timeseries.to_csv(timeseries_path, index=False)
    print(f"Wrote {len(timeseries)} daily facility-period rows to {timeseries_path}")

    neg_mid = int((result["pue_mid_century_delta"] < 0).sum())
    neg_eoc = int((result["pue_end_of_century_delta"] < 0).sum())
    print(
        f"Sanity check: {neg_mid} facilities with negative mid-century PUE delta, "
        f"{neg_eoc} with negative end-of-century PUE delta (expect close to 0 for both under SSP585, "
        "though PUE deltas -- unlike heat exceedance-day counts -- can plausibly be very small or "
        "slightly negative for a few facilities given daily-pairing noise)."
    )
    for col in ("pue_historical_1985_2014", "pue_ssp585_2045_2074", "pue_ssp585_2075_2100"):
        try:
            check_pue_plausible(result[col].to_numpy())
            print(f"Sanity check: all {col} values within plausible [1.0, 3.0] range.")
        except ValueError as exc:
            print(f"Sanity check WARNING for {col}: {exc}")

    cfg = RunConfig()
    pue_cfg = PUEModelConfig()
    manifest = pd.DataFrame(
        [
            {
                "model": cfg.model,
                "member": cfg.member,
                "historical_window": f"{cfg.hist_start}-{cfg.hist_end}",
                "mid_century_window": f"{cfg.mid_start}-{cfg.mid_end}",
                "end_of_century_window": f"{cfg.eoc_start}-{cfg.eoc_end}",
                "pue_variant": "AE_Chiller (airside economizer + adiabatic cooling)",
                "p_atm_pa": pue_cfg.p_atm_pa,
                "power_it_kw": pue_cfg.power_it_kw,
                "daily_pairing": "Tmax+RHmin and Tmin+RHmax, averaged per day",
            }
        ]
    )
    manifest.to_csv(out_dir / "run_manifest.csv", index=False)
    print(f"Wrote outputs to {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run LOCA2/humidity PUE processing for Oregon.")
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root path (default: current directory)",
    )
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
