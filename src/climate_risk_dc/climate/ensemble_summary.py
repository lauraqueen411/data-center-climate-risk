"""Shared ensemble summary statistics for the LOCA2 multi-GCM ensemble.

Used by both the PUE/WUE facility ensemble (``docs/tasks/pue_wue_analysis.md``)
and the extreme-heat facility/gridded ensemble
(``scripts/run_heat_indices_ensemble.py``) -- one reusable module rather than
duplicating the same mean/median/p5/p95/delta logic per metric, per the
project's convention of building a summary module, not one-off per-figure
code (see ``pue_wue_analysis.md`` Phase 1).

Assumptions log
----------------
- **5th/95th percentile, not min/max**, for ensemble spread -- matches
  OCA7's own reporting convention (``docs/oca7-structure-and-style-reference.md``:
  "ensemble mean plus 5th-95th percentile range across models, reported in
  parentheses"), not a project-specific choice.
- **Percentile interpolation**: pandas/numpy default (linear interpolation
  between the two nearest ranks) -- with a 20-27 member ensemble this is a
  minor, unremarked convention choice, not a documented sensitivity.
- ``ensemble_summary_stats`` is deliberately generic over ``group_cols`` --
  the same function computes both the statewide-fleet summary (grouped by
  archetype/period only, collapsing every facility x GCM row together) and
  the per-site summary (grouped by facility_id/archetype/period, collapsing
  GCMs only) depending on what's passed, rather than two near-duplicate
  functions.
- ``compute_period_deltas`` assumes exactly one historical value per
  ``id_cols`` group (one row after pivoting) -- raises via a clear KeyError
  from the pivot itself if a group is missing the historical period, rather
  than silently producing NaN deltas.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def ensemble_summary_stats(df: pd.DataFrame, value_col: str, group_cols: list[str]) -> pd.DataFrame:
    """Mean, median, 5th/95th percentile, and sample count of ``value_col`` within each group.

    Pass ``group_cols`` including ``facility_id`` for a per-site summary
    (collapsing across GCMs only), or excluding it for a statewide-fleet
    summary (collapsing across every facility x GCM row in the group).
    """
    grouped = df.groupby(group_cols, observed=True)[value_col]
    out = grouped.agg(
        mean="mean",
        median="median",
        p05=lambda s: s.quantile(0.05),
        p95=lambda s: s.quantile(0.95),
        n="size",
    )
    return out.reset_index()


def compute_period_deltas(
    df: pd.DataFrame,
    value_col: str,
    id_cols: list[str],
    period_col: str = "period",
    historical_label: str = "historical",
    future_labels: tuple[str, ...] = ("midcentury", "endcentury"),
) -> pd.DataFrame:
    """Future-minus-historical delta of ``value_col``, per ``id_cols`` group (e.g. facility x GCM x archetype).

    Returns one row per ``id_cols`` combination with the historical value and
    a ``delta_{future_label}`` column per future period.
    """
    wide = df.pivot_table(index=id_cols, columns=period_col, values=value_col)
    missing = [historical_label, *future_labels]
    absent = [p for p in missing if p not in wide.columns]
    if absent:
        raise KeyError(f"Period(s) {absent} not present in {period_col!r} column -- got {list(wide.columns)}")
    out = wide[[historical_label]].rename(columns={historical_label: f"{value_col}_{historical_label}"})
    for fut in future_labels:
        out[f"{value_col}_delta_{fut}"] = wide[fut] - wide[historical_label]
    return out.reset_index()


def ensemble_raster_stats(stack: np.ndarray) -> dict[str, np.ndarray]:
    """Per-cell mean/median/p5/p95 across a stack of per-GCM 2D rasters, shape ``(n_gcm, lat, lon)``.

    ``nanmean``/``nanpercentile`` propagate a cell's own NaN pattern (e.g. a
    coastline cell that's ocean/no-data for some GCMs' grids but not others)
    without letting one GCM's missing cell null out the whole ensemble
    statistic at that cell -- a cell with zero finite values across all GCMs
    still returns NaN, not a fabricated number.
    """
    return {
        "mean": np.nanmean(stack, axis=0),
        "median": np.nanmedian(stack, axis=0),
        "p05": np.nanpercentile(stack, 5, axis=0),
        "p95": np.nanpercentile(stack, 95, axis=0),
    }
