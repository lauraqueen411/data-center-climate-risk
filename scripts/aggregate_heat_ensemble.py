"""Collapse scripts/run_heat_indices_ensemble.py's per-GCM output into the ensemble figures that replace 1.1's original single-GCM/SSP5-8.5 maps.

Reads every GCM's per-(period, index) GeoTIFF under
``outputs/month2_heat_ensemble/per_gcm/<gcm>/`` and every GCM's per-period
facility CSV under ``outputs/month2_heat_ensemble/facilities/``, and writes:
- ``outputs/month2_heat_ensemble/rasters/{key}_{period}_{stat}.tif`` for
  stat in (mean, median, p05, p95) -- mean is the new headline map,
  matching OCA7's "ensemble mean" convention; p05/p95 are the supplementary
  spread rasters.
- ``outputs/month2_heat_ensemble/maps/{key}_{period}_mean.png`` -- quicklook,
  reusing ``run_heat_processing.py``'s own ``_plot_quicklook`` unchanged.
- ``outputs/month2_heat_ensemble/facility_heat_indices_ensemble_summary.csv``
  -- ensemble mean/median/5th/95th-percentile per facility, per period, per
  index, using the same shared ``ensemble_summary`` module the PUE/WUE
  analysis uses.

Assumptions log
----------------
- **Grid alignment across GCMs assumed identical** (same LOCA2
  0.0625-degree grid, same bbox/polygon mask applied identically per GCM) --
  checked, not assumed: every per-GCM raster's shape is compared before
  stacking, and a mismatch raises rather than silently broadcasting/cropping.
- **Ensemble stats via nan-aware per-cell reduction**
  (``ensemble_summary.ensemble_raster_stats``) -- a cell that's ocean/no-data
  in some GCMs' rasters (all four indices are polygon-masked identically, so
  in practice this is an inland-boundary rounding case, not open ocean) does
  not null out the ensemble statistic at that cell for the GCMs where it is
  valid.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import rioxarray  # noqa: F401 -- registers the "rasterio" xarray engine used below

from climate_risk_dc.climate.ensemble_summary import ensemble_raster_stats, ensemble_summary_stats
from climate_risk_dc.geo import load_oregon_data_centers
from run_heat_processing import _oregon_boundary, _plot_quicklook
from run_heat_indices_ensemble import INDEX_KEYS, PERIODS
from run_pue_wue_gridded import load_gcm_manifest


def _write_geotiff_like(reference_path: Path, data: np.ndarray, out_path: Path) -> None:
    with rasterio.open(reference_path) as ref:
        profile = ref.profile
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(data.astype("float32"), 1)


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out_dir = root / "outputs" / "month2_heat_ensemble"
    per_gcm_dir = out_dir / "per_gcm"
    facilities_dir = out_dir / "facilities"
    raster_dir = out_dir / "rasters"
    map_dir = out_dir / "maps"

    manifest = load_gcm_manifest(root)
    gcms = manifest[manifest["included_in_ensemble"] == "yes"]["model"].tolist()
    print(f"Aggregating across {len(gcms)} GCMs: {gcms}")

    watershed_gdb = root / "data" / "watershed_boundaries" / "WBD_National_GDB.gdb"
    dc_csv = root / "data" / "data_centers_im3_pnnl" / "im3_open_source_data_center_atlas_v2026.02.09.csv"
    boundary_gdf = _oregon_boundary(watershed_gdb)
    dc_gdf = load_oregon_data_centers(dc_csv)

    for period in PERIODS:
        for key in INDEX_KEYS:
            paths = [per_gcm_dir / gcm / f"{key}_{period}.tif" for gcm in gcms]
            missing = [p for p in paths if not p.exists()]
            if missing:
                raise FileNotFoundError(f"Missing per-GCM rasters for {key}/{period}: {missing[:3]}...")

            arrays = []
            for p in paths:
                with rasterio.open(p) as src:
                    arrays.append(src.read(1))
            shapes = {a.shape for a in arrays}
            if len(shapes) != 1:
                raise ValueError(f"{key}/{period}: inconsistent raster shapes across GCMs: {shapes}")
            stack = np.stack(arrays, axis=0)

            stats = ensemble_raster_stats(stack)
            for stat_name, stat_array in stats.items():
                out_path = raster_dir / f"{key}_{period}_{stat_name}.tif"
                _write_geotiff_like(paths[0], stat_array, out_path)

            mean_da = rioxarray.open_rasterio(raster_dir / f"{key}_{period}_mean.tif", masked=True).squeeze("band", drop=True)
            mean_da = mean_da.rename({"x": "lon", "y": "lat"})
            _plot_quicklook(
                mean_da, boundary_gdf, dc_gdf,
                title=f"{key} ensemble mean ({period}, SSP3-7.0, n={len(gcms)} GCMs)",
                out_path=map_dir / f"{key}_{period}_mean.png",
            )
            print(f"  {key}/{period}: wrote mean/median/p05/p95 rasters + mean quicklook")

    facility_parts = []
    for gcm in gcms:
        for period in PERIODS:
            csv_path = facilities_dir / f"facility_heat_indices_{gcm}_ssp370_{period}.csv"
            df = pd.read_csv(csv_path)
            df["gcm"] = gcm
            df["period"] = period
            facility_parts.append(df)
    facility_all = pd.concat(facility_parts, ignore_index=True)

    # Group by facility_id ONLY -- pandas groupby drops NaN keys by default,
    # and many facilities have a NaN facility_name in the raw atlas (same
    # issue already fixed once in build_pue_wue_ensemble_summary.py).
    summary_parts = []
    for key in INDEX_KEYS:
        summary_parts.append(
            ensemble_summary_stats(facility_all, key, ["facility_id", "period"]).assign(metric=key)
        )
    facility_summary = pd.concat(summary_parts, ignore_index=True)
    facility_names = facility_all[["facility_id", "facility_name"]].drop_duplicates("facility_id")
    facility_summary = facility_summary.merge(facility_names, on="facility_id", how="left")
    facility_summary.to_csv(out_dir / "facility_heat_indices_ensemble_summary.csv", index=False)
    print(f"Wrote facility_heat_indices_ensemble_summary.csv ({len(facility_summary)} rows)")

    manifest_out = pd.DataFrame([{"n_gcms": len(gcms), "gcms": ",".join(gcms), "scenario_fut": "ssp370", "periods": ",".join(PERIODS)}])
    manifest_out.to_csv(out_dir / "run_manifest.csv", index=False)
    print(f"Wrote outputs to {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate per-GCM heat-index ensemble outputs into ensemble stats.")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
