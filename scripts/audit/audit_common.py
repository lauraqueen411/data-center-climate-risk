"""Shared helpers for the Part A server checks (docs/tasks/methodology-audit.md).

Every check writes CSVs to ``outputs/audit_part_a/`` and prints a short
summary. Nothing here modifies model outputs.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

FAC_DIR = ROOT / "outputs/ensemble_pue_wue_facilities"
OLD_DIR = FAC_DIR / "superseded_rh_up60_2026-10-08"
OUT_DIR = ROOT / "outputs/audit_part_a"
OUT_DIR.mkdir(parents=True, exist_ok=True)

LOCA2_ROOT = Path("/home/redmond/data/pub/scripps_downscaled_CMIP6/LOCA2/CONUS_regions_split")
HUMIDITY_ROOT = Path("/home/redmond/data/pub/scripps_downscaled_CMIP6/LOCA2_humidity")
ELEVATION = ROOT / "data/elevation.LOCA_2016-04-02.nc"
ATLAS = ROOT / "data/data_centers_im3_pnnl/im3_open_source_data_center_atlas_v2026.02.09.csv"

PERIODS = {"historical": (1985, 2014), "midcentury": (2045, 2074), "endcentury": (2075, 2100)}
ARCHETYPES = ("ae-chiller", "we-chiller", "chiller-only")


def manifest() -> pd.DataFrame:
    m = pd.read_csv(ROOT / "outputs/ensemble_pue_wue/gcm_manifest.csv")
    return m[m["included_in_ensemble"] == "yes"].reset_index(drop=True)


def gcms() -> list[str]:
    return list(manifest()["model"])


def member(gcm: str) -> str:
    m = manifest()
    return str(m.loc[m["model"] == gcm, "member"].iloc[0])


def scenario(period: str) -> str:
    return "historical" if period == "historical" else "ssp370"


def daily_path(gcm: str, archetype: str, period: str, folder: Path = FAC_DIR, backend: str = "fast") -> Path:
    return folder / f"facility_pue_wue_{gcm}_ssp370_{archetype}_{period}_{backend}.csv"


def facilities() -> pd.DataFrame:
    """Facility id, lon, lat (as the runner loads them), ecoregion and grid cell where available."""
    from climate_risk_dc.geo import load_oregon_data_centers

    fac = load_oregon_data_centers(ATLAS)
    out = pd.DataFrame({"facility_id": fac["id"].to_numpy(), "lon": fac["lon"].to_numpy(), "lat": fac["lat"].to_numpy()})
    cells = ROOT / "outputs/ensemble_pue_wue_facilities/analysis/fig_show_everything/facility_grid_cells.csv"
    if cells.exists():
        out = out.merge(pd.read_csv(cells, usecols=["facility_id", "ecoregion", "grid_cell", "grid_lat", "grid_lon"]),
                        on="facility_id", how="left")
    return out


def pressure_by_facility(fac: pd.DataFrame) -> dict[int, float]:
    """Elevation-adjusted pressure exactly as run_pue_wue_facilities.py computes it."""
    from climate_risk_dc.climate.cooling_archetypes import elevation_adjusted_pressure_pa
    from run_pue_wue_facilities import sample_elevation_at_points

    elev = sample_elevation_at_points(ELEVATION, fac["lon"].to_numpy(), fac["lat"].to_numpy())
    return dict(zip(fac["facility_id"], elevation_adjusted_pressure_pa(elev)))


def write(df: pd.DataFrame, name: str) -> Path:
    path = OUT_DIR / name
    df.to_csv(path, index=False)
    print(f"  wrote {path.relative_to(ROOT)} ({len(df)} rows)")
    return path


def lon_query(lon: float, lon_values: np.ndarray) -> float:
    """Facility longitude in the file's own convention (0-360 or +/-180)."""
    return lon + 360.0 if np.nanmax(lon_values) > 180.0 and lon < 0 else lon
