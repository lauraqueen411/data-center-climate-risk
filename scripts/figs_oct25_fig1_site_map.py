"""Figure 1 redraw (reviewer comment, October 2026): darker ecoregion
boundaries and four labelled cities named in the text.

Same map as ``build_pue_wue_ensemble_figures.figure1_site_map`` (study-area
map; copied to ``docs/methodology-figures/study_area_map.png``): same
projection, extent, ecoregion fills, facility dots, labels, figure size
(8 x 4.8 in) and dpi (140). Its helpers and constants are imported, not
copied. Changes:

- Ecoregion boundaries are drawn as a separate edge-only layer above the
  state lines, dark grey (``BOUNDARY_COLOR``) and slightly thicker.
- The Dalles, Portland, Prineville and Eugene are marked with white stars
  and labelled with a white halo. Stars sit beneath the facility dots so
  the dots at Prineville and The Dalles stay visible.
- Ecoregion labels get a thin white halo so the darker boundaries do not
  cut through them.

Inputs: ``outputs/ensemble_pue_wue_facilities/analysis/ensemble_summary_by_site.csv``
(facility coordinates and ecoregion, as in the original) and the Level III
ecoregion shapefile.

Outputs (``docs/results-figures/figures_oct25/``): ``fig1_site_map.{png,pdf}``,
``preview_fig1_old_vs_new.png``.

Assumptions log
----------------
- **City coordinates** are city-centre points (WGS84, from the US Census
  Gazetteer / GNIS): Portland 45.5152 N, 122.6784 W; The Dalles 45.5946 N,
  121.1787 W; Prineville 44.2999 N, 120.8345 W; Eugene 44.0521 N,
  123.0868 W.
- **Label offsets** (``CITIES``, degrees) were set by eye so labels sit off
  the facility dots and the ecoregion labels.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import pandas as pd

from build_pue_wue_ensemble_figures import (
    ECOREGION_FILL_COLORS,
    ECOREGION_FILL_NODATA,
    INK_PRIMARY,
    INK_SECONDARY,
    _apply_style,
    _oregon_axes,
)
from climate_risk_dc.config import load_dataset_paths
from climate_risk_dc.geo import to_equal_area, to_geographic
from figs_oct25_common import OUT_DIR, side_by_side

BOUNDARY_COLOR = "#4a4a4a"
BOUNDARY_WIDTH = 0.9
# name: (lon, lat, label dx, label dy, ha)
CITIES = {
    "Portland": (-122.6784, 45.5152, 0.12, 0.13, "left"),
    "The Dalles": (-121.1787, 45.5946, -0.1, 0.2, "center"),
    "Prineville": (-120.8345, 44.2999, 0.17, -0.13, "left"),
    "Eugene": (-123.0868, 44.0521, -0.14, 0.0, "right"),
}
HALO = [pe.withStroke(linewidth=2.5, foreground="white")]


def plot(by_site: pd.DataFrame, path: Path) -> None:
    sites = by_site[["facility_id", "lon", "lat", "US_L3NAME"]].drop_duplicates("facility_id")
    assert len(sites) == 109
    fig = plt.figure(figsize=(8, 4.8))
    ax = _oregon_axes(fig)

    ecoregions = gpd.read_file(load_dataset_paths().ecoregion_l3_shp).to_crs("EPSG:4326")
    for row in ecoregions.itertuples():
        ax.add_geometries([row.geometry], crs=ccrs.PlateCarree(),
                          facecolor=ECOREGION_FILL_COLORS.get(row.US_L3NAME, ECOREGION_FILL_NODATA),
                          edgecolor="none", zorder=0)
    ax.add_geometries(list(ecoregions.geometry), crs=ccrs.PlateCarree(), facecolor="none",
                      edgecolor=BOUNDARY_COLOR, linewidth=BOUNDARY_WIDTH, zorder=1.5)

    ax.scatter(sites["lon"], sites["lat"], color=INK_PRIMARY, s=18, edgecolor="white", linewidth=0.3, zorder=2,
               transform=ccrs.PlateCarree())

    # Ecoregion labels exactly as in figure1_site_map
    n_by_eco = sites.groupby("US_L3NAME")["facility_id"].count()
    label_overrides = {
        "Eastern Cascades Slopes and Foothills": "Eastern Cascades\nSlopes and Foothills",
        "Klamath Mountains/California High North Coast Range": "Klamath Mountains",
    }
    label_y_offsets = {
        "Willamette Valley": -0.35,
        "Eastern Cascades Slopes and Foothills": -0.6,
        "Klamath Mountains/California High North Coast Range": 0.15,
    }
    centroids = to_geographic(to_equal_area(gpd.read_file(load_dataset_paths().ecoregion_l3_shp)).assign(
        geometry=lambda g: g.centroid))
    for eco_row in centroids.itertuples():
        name = eco_row.US_L3NAME
        if name not in n_by_eco.index:
            continue
        ax.annotate(f"{label_overrides.get(name, name)}\n(n={n_by_eco[name]})",
                    xy=(eco_row.geometry.x, eco_row.geometry.y + label_y_offsets.get(name, 0)),
                    transform=ccrs.PlateCarree(), fontsize=12, color=INK_SECONDARY, ha="center", va="center", zorder=3,
                    path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])

    for name, (lon, lat, dx, dy, ha) in CITIES.items():
        ax.scatter(lon, lat, marker="*", s=75, facecolor="white", edgecolor=INK_PRIMARY, linewidth=0.8, zorder=1.8,
                   transform=ccrs.PlateCarree())
        ax.text(lon + dx, lat + dy, name, transform=ccrs.PlateCarree(), fontsize=10, color=INK_PRIMARY, ha=ha,
                va="center", zorder=5, path_effects=HALO)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(path.with_suffix(f".{ext}"), bbox_inches="tight")
    plt.close(fig)


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out = root / OUT_DIR
    out.mkdir(parents=True, exist_ok=True)
    _apply_style()
    by_site = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/ensemble_summary_by_site.csv")
    plot(by_site, out / "fig1_site_map.png")
    side_by_side(root / "docs/methodology-figures/study_area_map.png", out / "fig1_site_map.png",
                 out / "preview_fig1_old_vs_new.png")
    print(f"Wrote {out / 'fig1_site_map.png'} (+ .pdf, preview)")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Figure 1 redraw: darker ecoregion boundaries, labelled cities.")
    p.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return p.parse_args()


if __name__ == "__main__":
    run(parse_args())
