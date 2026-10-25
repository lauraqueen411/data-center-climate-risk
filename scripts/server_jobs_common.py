"""Shared helpers for the October 23 server jobs (draft 4 of the data center chapter).

Not used by the pipeline. The ``server_jobs_*.py`` scripts import the draft 3
conventions from ``build_draft3_values.py`` (ecoregion abbreviations, order,
rounding, and ``eco``: equal-weight ecoregion mean of facilities within each
model, then mean and p25/p75 across the 20 models) so the numbers are
computed exactly as draft 3 computed them.

Assumptions log
----------------
- **Facility pressure** uses the nearest LOCA elevation cell and
  ``elevation_adjusted_pressure_pa``, the same rule as
  ``run_pue_wue_facilities.sample_elevation_at_points`` and
  ``fig_show_everything._facility_table``. Reimplemented here (five lines)
  rather than imported, because importing either script pulls in the full
  run or plotting stack.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from build_draft3_values import ABBR, MONTHS, N_MODELS, NAME, ORDER, eco, t  # noqa: F401 (re-exported)
from climate_risk_dc.climate.cooling_archetypes import elevation_adjusted_pressure_pa

PERIODS = {"historical": (1985, 2014), "midcentury": (2045, 2074), "endcentury": (2075, 2100)}
PERIOD_LABEL = {"historical": "1985–2014", "midcentury": "2045–2074", "endcentury": "2075–2100"}
FAC_DIR = Path("outputs/ensemble_pue_wue_facilities")
OUT_DIR = Path("docs/results-figures/server_jobs_oct23")
DRAFT3 = Path("docs/results-figures/draft3")
FNAME_RE = re.compile(
    r"facility_pue_wue_(?P<gcm>.+)_ssp370_(?P<archetype>ae-chiller|we-chiller|chiller-only)_"
    r"(?P<period>historical|midcentury|endcentury)_(?P<backend>fast|coolprop)\.csv"
)


def facilities(root: Path, with_pressure: bool = False) -> pd.DataFrame:
    """facility_id, ecoregion, abbr, grid_cell, lon, lat (+ pressure_pa)."""
    fac = pd.read_csv(root / FAC_DIR / "analysis/fig_show_everything/facility_grid_cells.csv")
    fac["abbr"] = fac["ecoregion"].map(ABBR)
    assert len(fac) == 109 and fac["abbr"].notna().all() and fac["grid_cell"].nunique() == 22
    if with_pressure:
        elev = xr.open_dataset(root / "data/elevation.LOCA_2016-04-02.nc")["Elevation"]
        e_lat, e_lon = np.asarray(elev["lat"].values), np.asarray(elev["lon"].values)
        lon = np.where(fac["lon"] < 0, fac["lon"] + 360.0, fac["lon"]) if e_lon.max() > 180 else fac["lon"].to_numpy()
        ilat = np.abs(e_lat[:, None] - fac["lat"].to_numpy()[None, :]).argmin(axis=0)
        ilon = np.abs(e_lon[:, None] - lon[None, :]).argmin(axis=0)
        fac["pressure_pa"] = elevation_adjusted_pressure_pa(elev.values[ilat, ilon])
        assert fac["pressure_pa"].between(80_000, 102_000).all(), "facility pressure outside 80-102 kPa"
    return fac


def models(root: Path) -> list[str]:
    m = pd.read_csv(root / "outputs/ensemble_pue_wue/gcm_manifest.csv")
    out = m.loc[m["included_in_ensemble"] == "yes", "model"].tolist()
    assert len(out) == N_MODELS
    return out


def daily_path(root: Path, gcm: str, archetype: str, period: str, backend: str = "fast") -> Path:
    return root / FAC_DIR / f"facility_pue_wue_{gcm}_ssp370_{archetype}_{period}_{backend}.csv"


def eco_rows(per_model: pd.DataFrame, col: str, kind: str, diff: bool, **keys) -> list[dict]:
    """Ecoregion rows (mean, p25, p75 across models), formatted with draft 3 decimals, plus raw values."""
    e = eco(per_model, col)
    return [{**keys, "ecoregion": a, "mean": t(e.loc[a, "mean"], kind, diff), "p25": t(e.loc[a, "p25"], kind, diff),
             "p75": t(e.loc[a, "p75"], kind, diff), "_mean": e.loc[a, "mean"]} for a in ORDER]


def compare(label: str, ours: pd.Series, ref: pd.Series) -> list[str]:
    """Compare formatted strings aligned on index; return mismatch lines (empty if all match)."""
    ours, ref = ours.astype(str), ref.astype(str)
    ref = ref.reindex(ours.index)
    bad = ours[ours != ref]
    print(f"[validate] {label}: {len(ours) - len(bad)}/{len(ours)} match")
    return [f"{label} {i}: ours {ours[i]} vs draft3 {ref[i]}" for i in bad.index]
