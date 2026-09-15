"""Extreme-heat index processing across the LOCA2 SSP3-7.0 ensemble (replaces 1.1's single-GCM SSP5-8.5 figures).

Extends ``run_heat_processing.py`` (locked to ACCESS-CM2/SSP5-8.5, per that
script's own documented limitation) to the same 20-GCM SSP3-7.0 ensemble
already used for the facility PUE/WUE run (``docs/tasks/pue_wue_analysis.md``),
so Section 1.1's extreme-heat figures and Figure 7's exceedance-day/PUE
comparison are on the same ensemble and scenario footing rather than mixing
a single-GCM SSP5-8.5 result with a 20-GCM SSP3-7.0 one.

Unlike the PUE/WUE gridded attempt (shelved as computationally infeasible --
CoolProp's per-point cost), heat-index computation is plain threshold/degree-day
arithmetic on tasmax/tasmin with no CoolProp or GP-model step: benchmarked at
~55s per (GCM, period) for the full Oregon grid, so a full 20-GCM x 3-period
gridded ensemble (~55 min sequential) is tractable, unlike PUE/WUE.

This script computes and writes, per (GCM, period): the four full-grid indices
(cdd, tmax_annual_max, tmax_threshold_days, heatwave_events) as GeoTIFFs, and
the same four indices sampled at the 109 Oregon facility points as one CSV.
``scripts/aggregate_heat_ensemble.py`` collapses these per-GCM outputs into the
ensemble mean/median/5th/95th-percentile rasters, maps, and facility summary
that replace 1.1's figures.

CLI: one required argument, ``--gcms``, matching ``run_pue_wue_gridded.py``'s
own CLI contract exactly (sequential within one invocation; concurrency across
GCMs is the operator's job, launching one process per GCM/small batch).

Assumptions log
----------------
- **Scenario**: SSP3-7.0 for the two future periods (not SSP5-8.5, which the
  original single-GCM ``run_heat_processing.py`` used) -- this is the whole
  point of this script, per the user's explicit "replace all 1.1 figures to
  use ensemble data" direction. The original SSP5-8.5/ACCESS-CM2 outputs
  under ``outputs/month2_heat/`` are left untouched (not deleted or
  overwritten in place) as the previously-reported single-GCM benchmark --
  this ensemble run is written to a new ``outputs/month2_heat_ensemble/``
  tree, which supersedes the old figures for the chapter without destroying
  the prior result.
- **Period labels**: ``historical``/``midcentury``/``endcentury``, matching
  the PUE/WUE facility ensemble's convention (not ``run_heat_processing.py``'s
  own ``historical_1985_2014``/``ssp585_2045_2074``/``ssp585_2075_2100``
  labels), so the two ensembles' output files line up by period name for
  Figure 7's direct pairing.
- **GCM/member selection**: read from the same
  ``outputs/ensemble_pue_wue/gcm_manifest.csv`` the PUE/WUE run built --
  reused as-is (not re-scanned) so both ensembles share the exact same GCM
  list and member IDs; a mismatch there would silently break Figure 7's
  premise that both metrics come from "the same ensemble."
- **Domain and boundary**: same Oregon-polygon boundary (WBD watershed
  dissolve) and bbox as the original ``run_heat_processing.py`` -- unlike the
  PUE/WUE gridded script's plain bbox (a documented, separate deviation for
  that pipeline), this one already masks to the actual state polygon, so no
  ocean/no-data NaN-crash risk analogous to the PUE/WUE gridded bug applies
  here (confirmed: no CoolProp-style calls in this codepath to crash on NaN
  in the first place -- threshold/degree-day arithmetic propagates NaN
  silently, which is the correct behavior for masked-out cells).
- **Facility sampling**: sampled from each period's *written* GeoTIFF via
  ``geo.sample_raster_at_points`` (nearest-cell), reusing that existing,
  tested helper rather than writing a new in-memory DataArray sampler --
  consistent with the project's convention of reusing raster point-sampling
  code rather than duplicating it.
- **Checkpointing**: per (GCM, period, key) raster and per (GCM, period)
  facility CSV are the resume/skip units, mirroring
  ``run_pue_wue_gridded.py``'s own checkpointing philosophy.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from climate_risk_dc.climate.heat_indices import (
    HeatIndexConfig,
    clip_to_bbox,
    compute_period_heat_indices,
    mask_to_polygon,
)
from climate_risk_dc.climate.loca2_io import resolve_loca2_file, temperature_path_config
from climate_risk_dc.geo import load_oregon_data_centers, sample_raster_at_points
from run_heat_processing import _da_to_geotiff, _log_unit_sanity, _oregon_boundary
from run_pue_wue_gridded import load_gcm_manifest, resolve_gcm_member

PERIODS: dict[str, tuple[int, int]] = {
    "historical": (1985, 2014),
    "midcentury": (2045, 2074),
    "endcentury": (2075, 2100),  # 26 years, not 30 -- see methods doc Time Horizons
}
INDEX_KEYS = ("cdd", "tmax_annual_max", "tmax_threshold_days", "heatwave_events")
MANIFEST_RELATIVE_PATH = Path("outputs") / "ensemble_pue_wue" / "gcm_manifest.csv"


@dataclass(frozen=True)
class HeatEnsembleConfig:
    scenario_hist: str = "historical"
    scenario_fut: str = "ssp370"
    threshold_c: float = 35.0
    cdd_base_c: float = 18.3
    min_duration_days: int = 3
    lon_min: float = -125.0
    lon_max: float = -116.0
    lat_min: float = 42.0
    lat_max: float = 46.5


def _scenario_for_period(period_name: str, cfg: HeatEnsembleConfig) -> str:
    return cfg.scenario_hist if period_name == "historical" else cfg.scenario_fut


def run_one_gcm_period(
    gcm: str,
    member: str,
    period_name: str,
    root: Path,
    args: argparse.Namespace,
    cfg: HeatEnsembleConfig,
    hi_cfg: HeatIndexConfig,
    boundary_poly,
    boundary_gdf,
    facilities,
    per_gcm_dir: Path,
    facilities_dir: Path,
) -> None:
    facility_csv_path = facilities_dir / f"facility_heat_indices_{gcm}_{cfg.scenario_fut}_{period_name}.csv"
    start_year, end_year = PERIODS[period_name]
    scenario = _scenario_for_period(period_name, cfg)

    raster_paths = {key: per_gcm_dir / gcm / f"{key}_{period_name}.tif" for key in INDEX_KEYS}
    need_rasters = not all(p.exists() for p in raster_paths.values())
    if not need_rasters and facility_csv_path.exists():
        print(f"  [skip, exists] {gcm}/{period_name}")
        return

    path_cfg = temperature_path_config(Path(args.loca2_root), model=gcm, member=member)
    tasmax_file = resolve_loca2_file(path_cfg, scenario, "tasmax", start_year, end_year)
    tasmin_file = resolve_loca2_file(path_cfg, scenario, "tasmin", start_year, end_year)
    _log_unit_sanity(tasmax_file, tasmin_file)

    idx = compute_period_heat_indices(
        tasmax_path=tasmax_file, tasmin_path=tasmin_file, start_year=start_year, end_year=end_year, config=hi_cfg
    )

    sampled = {"facility_id": facilities["id"].to_numpy(), "facility_name": facilities["name"].to_numpy()}
    for key in INDEX_KEYS:
        da = clip_to_bbox(idx[key], cfg.lon_min, cfg.lon_max, cfg.lat_min, cfg.lat_max)
        da = mask_to_polygon(da, boundary_poly)
        tif_path = raster_paths[key]
        tif_path.parent.mkdir(parents=True, exist_ok=True)
        if not tif_path.exists():
            _da_to_geotiff(da, tif_path)
        sampled[key] = sample_raster_at_points(tif_path, facilities["lon"].to_numpy(), facilities["lat"].to_numpy())

    facilities_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(sampled).to_csv(facility_csv_path, index=False)
    print(f"  [done] {gcm}/{period_name} -> {facility_csv_path.name}")


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    cfg = HeatEnsembleConfig()
    hi_cfg = HeatIndexConfig(cdd_base_c=cfg.cdd_base_c, threshold_c=cfg.threshold_c, min_duration_days=cfg.min_duration_days)

    out_dir = root / "outputs" / "month2_heat_ensemble"
    per_gcm_dir = out_dir / "per_gcm"
    facilities_dir = out_dir / "facilities"

    watershed_gdb = root / "data" / "watershed_boundaries" / "WBD_National_GDB.gdb"
    dc_csv = root / "data" / "data_centers_im3_pnnl" / "im3_open_source_data_center_atlas_v2026.02.09.csv"
    boundary_gdf = _oregon_boundary(watershed_gdb)
    boundary_poly = boundary_gdf.geometry.iloc[0]
    dc_gdf = load_oregon_data_centers(dc_csv)
    facilities = pd.DataFrame({"id": dc_gdf["id"], "name": dc_gdf["name"], "lon": dc_gdf["lon"], "lat": dc_gdf["lat"]})

    manifest = load_gcm_manifest(root)

    for gcm in args.gcms:
        member = resolve_gcm_member(manifest, gcm)
        print(f"=== {gcm} (member {member}) ===")
        for period_name in PERIODS:
            run_one_gcm_period(
                gcm, member, period_name, root, args, cfg, hi_cfg, boundary_poly, boundary_gdf, facilities,
                per_gcm_dir, facilities_dir,
            )

    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Wrote per-GCM outputs to {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extreme-heat index ensemble (SSP3-7.0) for one or more LOCA2 GCMs.")
    parser.add_argument("--gcms", nargs="+", required=True, help="One or more GCM names, processed sequentially.")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    parser.add_argument(
        "--loca2-root",
        default="/home/redmond/data/pub/scripps_downscaled_CMIP6/LOCA2/CONUS_regions_split",
        help="Root path for LOCA2 CONUS_regions_split temperature files",
    )
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
