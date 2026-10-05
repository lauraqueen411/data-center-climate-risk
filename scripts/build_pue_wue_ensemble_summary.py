"""Phase 1 of docs/tasks/pue_wue_analysis.md: ensemble summary statistics for the facility PUE/WUE run.

Reads the 180 raw daily facility CSVs
(``outputs/ensemble_pue_wue_facilities/facility_pue_wue_{gcm}_ssp370_{archetype}_{period}_fast.csv``,
20 GCMs x 3 archetypes x 3 periods, ~1.19M daily rows each) and aggregates
them into the summary tables the chapter's figures/tables are built from,
using the shared ``climate_risk_dc.climate.ensemble_summary`` module (also
used by the extreme-heat ensemble) rather than one-off per-figure code.

Archetype note: the task doc's Phase 0 review found "Case 3" in
``docs/tasks/pue_wue_analysis.md``/``docs/data-center-chapter-outline.md`` was
a typo for **Case 2** (``we-chiller``, large-scale waterside economizer +
water-cooled chiller) -- Case 2 is what was actually run (confirmed: "Case 3"
appears nowhere in the codebase). This script and everything downstream uses
the real archetype set: Case 1 (``ae-chiller``), Case 2 (``we-chiller``),
Case 5 (``chiller-only``).

Outputs (``outputs/ensemble_pue_wue_facilities/analysis/``):
- ``facility_period_summary.csv`` -- one row per (facility, gcm, archetype,
  period): period-mean and period-95th-percentile PUE/WUE. The per-file
  raw-daily-row aggregation step (~214M rows total across all 180 files if
  loaded raw) -- each file is aggregated and discarded immediately, so peak
  memory stays at one file's ~1.19M rows.
- ``facility_month_summary.csv`` -- one row per (facility, gcm, archetype,
  period, month): monthly-mean PUE/WUE, for the seasonal-climatology figure.
- ``ensemble_summary_statewide.csv`` -- per (archetype, period, metric):
  ensemble mean/median/5th/95th-percentile across all 109 sites x 20 GCMs.
- ``ensemble_summary_by_site.csv`` -- per (facility, archetype, period,
  metric): ensemble mean/median/5th/95th-percentile across GCMs only.
- ``ensemble_delta.csv`` -- per (facility, gcm, archetype): historical value
  plus mid-century/end-century deltas, for both PUE and WUE.
- ``ecoregion_representative_facilities.csv`` -- one row per EPA Level III
  ecoregion (Omernik 1987) that contains at least one facility: the single
  facility nearest that ecoregion's own facility centroid, plus how many
  facilities the ecoregion contains and how far the pick sits from the
  centroid. Feeds Figure 9 / Table 5 (the sub-region breakdown) in
  ``scripts/build_pue_wue_ensemble_figures.py``.

Assumptions log
----------------
- **One representative facility per ecoregion, not an average of the
  ecoregion's climate inputs** -- follows Lei & Masanet (2022)'s own
  precedent of representing a climate region with one representative point
  (their per-IECC/ASHRAE-zone representative city) rather than spatially
  averaging temperature/humidity across many points first, which would run
  into Jensen's-inequality problems given the model's thresholded
  economizer/chiller switching. The representative facility is the one
  nearest the ecoregion's own facility-cluster centroid (computed in
  ``CRS_EQUAL_AREA`` metres via ``geo.select_representative_facility``), not
  the ecoregion polygon's geographic centroid -- this ties the pick to where
  exposure actually clusters. No new model run is needed: the representative
  facility's daily PUE/WUE was already computed as part of the full 109-site
  run, so this step is a selection, not a computation.
- **95th percentile, not max**, for the per-facility-per-gcm period summary
  -- a single extreme day's outlier shouldn't set a facility's whole-period
  "extreme" characterization; matches the same p95-over-max choice already
  documented in the shared ``ensemble_summary`` module for the *ensemble*
  level, applied here at the *daily* level too for consistency.
- **Facility lon/lat** joined from the raw IM3/PNNL atlas CSV directly
  (state_abb == "OR"), not from ``climate_risk_dc.geo`` -- this script lives
  in ``scripts/`` alongside the run scripts (not a decoupled notebook), so
  importing the package is fine here; the analysis notebook downstream stays
  decoupled per the project's existing notebook convention.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import geopandas as gpd
import pandas as pd

from climate_risk_dc.climate.ensemble_summary import compute_period_deltas, ensemble_summary_stats
from climate_risk_dc.config import load_dataset_paths
from climate_risk_dc.geo import assign_us_l3_ecoregion, points_from_lonlat, select_representative_facility

FNAME_RE = re.compile(
    r"facility_pue_wue_(?P<gcm>.+)_ssp370_(?P<archetype>ae-chiller|we-chiller|chiller-only)_"
    r"(?P<period>historical|midcentury|endcentury)_fast\.csv"
)
ARCHETYPE_ORDER = ["ae-chiller", "we-chiller", "chiller-only"]
PERIOD_ORDER = ["historical", "midcentury", "endcentury"]


def aggregate_raw_files(fac_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate all raw daily facility CSVs to period-level and month-level summaries."""
    files = sorted(fac_dir.glob("facility_pue_wue_*_fast.csv"))
    if not files:
        raise FileNotFoundError(f"No facility_pue_wue_*_fast.csv files found under {fac_dir}")

    period_parts = []
    month_parts = []
    for f in files:
        m = FNAME_RE.match(f.name)
        if not m:
            raise ValueError(f"Unparseable filename: {f.name}")
        meta = m.groupdict()
        df = pd.read_csv(f, usecols=["facility_id", "date", "pue_daily", "wue_daily"])

        # Group by facility_id ONLY -- pandas groupby drops NaN keys by default,
        # and many facilities have a NaN facility_name in the raw atlas, so
        # including it as a second group key silently loses those facilities.
        period_agg = df.groupby("facility_id", as_index=False).agg(
            pue_mean=("pue_daily", "mean"),
            pue_p95=("pue_daily", lambda s: s.quantile(0.95)),
            wue_mean=("wue_daily", "mean"),
            wue_p95=("wue_daily", lambda s: s.quantile(0.95)),
        )
        period_agg["gcm"] = meta["gcm"]
        period_agg["archetype"] = meta["archetype"]
        period_agg["period"] = meta["period"]
        period_parts.append(period_agg)

        df["month"] = pd.to_datetime(df["date"]).dt.month
        month_agg = df.groupby(["facility_id", "month"], as_index=False).agg(
            pue_mean=("pue_daily", "mean"), wue_mean=("wue_daily", "mean")
        )
        month_agg["gcm"] = meta["gcm"]
        month_agg["archetype"] = meta["archetype"]
        month_agg["period"] = meta["period"]
        month_parts.append(month_agg)

    period_summary = pd.concat(period_parts, ignore_index=True)
    month_summary = pd.concat(month_parts, ignore_index=True)
    return period_summary, month_summary


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    fac_dir = root / "outputs" / "ensemble_pue_wue_facilities"
    out_dir = fac_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    dc_csv = root / "data" / "data_centers_im3_pnnl" / "im3_open_source_data_center_atlas_v2026.02.09.csv"
    dc_atlas = pd.read_csv(dc_csv)
    facilities_geo = dc_atlas[dc_atlas["state_abb"] == "OR"][["id", "name", "lon", "lat"]].rename(
        columns={"id": "facility_id", "name": "facility_name"}
    )

    ecoregions = gpd.read_file(load_dataset_paths().ecoregion_l3_shp)
    facilities_points = points_from_lonlat(facilities_geo)
    facilities_labeled = assign_us_l3_ecoregion(facilities_points, ecoregions)
    assert facilities_labeled["US_L3NAME"].notna().all(), "Facilities outside every Oregon Level III ecoregion polygon"
    facilities_geo = pd.DataFrame(facilities_labeled.drop(columns="geometry"))

    representative = select_representative_facility(facilities_labeled, "US_L3NAME")
    representative = representative.rename(
        columns={"US_L3NAME": "ecoregion", "n_facilities_in_group": "n_facilities_in_ecoregion"}
    )[["ecoregion", "facility_id", "facility_name", "lon", "lat", "n_facilities_in_ecoregion", "distance_to_centroid_m"]]
    representative.to_csv(out_dir / "ecoregion_representative_facilities.csv", index=False)
    print("Ecoregion representative facilities (nearest to each ecoregion's own facility centroid):")
    print(representative.to_string(index=False))

    print("Aggregating 180 raw daily facility CSVs (this reads ~23GB total)...")
    period_summary, month_summary = aggregate_raw_files(fac_dir)
    n_expected = len(ARCHETYPE_ORDER) * len(PERIOD_ORDER) * 20  # 20 GCMs confirmed in gcm_manifest.csv
    n_actual = period_summary[["gcm", "archetype", "period"]].drop_duplicates().shape[0]
    print(f"Aggregated {n_actual} (gcm, archetype, period) combinations (expected {n_expected})")

    period_summary = period_summary.merge(facilities_geo, on="facility_id", how="left")
    assert period_summary["lon"].notna().all(), "Some facility_ids failed to join to the facility inventory"

    period_summary.to_csv(out_dir / "facility_period_summary.csv", index=False)
    month_summary.to_csv(out_dir / "facility_month_summary.csv", index=False)
    print(f"Wrote facility_period_summary.csv ({len(period_summary)} rows), facility_month_summary.csv ({len(month_summary)} rows)")

    statewide_parts, by_site_parts, delta_parts = [], [], []
    for metric in ("pue_mean", "wue_mean"):
        statewide_parts.append(
            ensemble_summary_stats(period_summary, metric, ["archetype", "period"]).assign(metric=metric)
        )
        by_site_parts.append(
            ensemble_summary_stats(period_summary, metric, ["facility_id", "archetype", "period"]).assign(metric=metric)
        )
        # id_cols includes archetype so one call handles all three archetypes at once
        # (pivoting only on period) -- looping per archetype and outer-merging the
        # per-archetype frames back together on a key that includes archetype itself
        # doesn't work: disjoint archetype values across frames never match, so an
        # "outer" merge just concatenates rows with _x/_y-suffixed NaN columns instead
        # of widening them.
        delta_parts.append(compute_period_deltas(period_summary, metric, ["facility_id", "gcm", "archetype"]))

    ensemble_summary_statewide = pd.concat(statewide_parts, ignore_index=True)
    ensemble_summary_by_site = pd.concat(by_site_parts, ignore_index=True)
    ensemble_summary_by_site = ensemble_summary_by_site.merge(facilities_geo, on="facility_id", how="left")

    ensemble_delta = delta_parts[0].merge(delta_parts[1], on=["facility_id", "gcm", "archetype"], how="outer")

    ensemble_summary_statewide.to_csv(out_dir / "ensemble_summary_statewide.csv", index=False)
    ensemble_summary_by_site.to_csv(out_dir / "ensemble_summary_by_site.csv", index=False)
    ensemble_delta.to_csv(out_dir / "ensemble_delta.csv", index=False)
    print("Wrote ensemble_summary_statewide.csv, ensemble_summary_by_site.csv, ensemble_delta.csv")

    manifest = pd.DataFrame(
        [{"n_files": 180, "n_gcms": 20, "archetypes": ",".join(ARCHETYPE_ORDER), "periods": ",".join(PERIOD_ORDER), "backend": "fast"}]
    )
    manifest.to_csv(out_dir / "run_manifest.csv", index=False)
    print(f"Wrote outputs to {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build PUE/WUE ensemble summary statistics (Phase 1 of pue_wue_analysis.md).")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
