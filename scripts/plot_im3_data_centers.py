#!/usr/bin/env python3
"""Plot the IM3/PNNL data-center atlas as a point map.

The script reads the repository's IM3 PNNL data-center dataset, builds a
GeoDataFrame from the geospatial layer when available (falling back to the CSV
lon/lat columns), and saves a static map image to the repository outputs
folder.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import geopandas as gpd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from climate_risk_dc.config import load_dataset_paths


def load_dataset() -> gpd.GeoDataFrame:
    """Load the IM3 PNNL data-center dataset as a GeoDataFrame."""
    datasets = load_dataset_paths()
    gpkg_path = datasets.data_centers_gpkg
    csv_path = datasets.data_centers_csv

    if gpkg_path.exists():
        gdf = gpd.read_file(gpkg_path, layer="point")
        if gdf.crs is None:
            gdf = gdf.set_crs("EPSG:4326")
        return gdf

    if not csv_path.exists():
        raise FileNotFoundError(f"Could not find dataset files: {gpkg_path} or {csv_path}")

    df = pd.read_csv(csv_path)
    geometry = gpd.points_from_xy(df["lon"], df["lat"])
    return gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")


def plot_dataset(gdf: gpd.GeoDataFrame, output_path: Path) -> None:
    """Create a static map and save it to disk."""
    gdf = gdf.copy()
    gdf["type"] = gdf["type"].fillna("unknown")
    gdf["sqft"] = pd.to_numeric(gdf.get("sqft"), errors="coerce").fillna(0)

    # Make marker sizes comparable across records.
    gdf["marker_size"] = np.clip(np.sqrt(gdf["sqft"] / 1_000) + 20, 20, 250)

    fig, ax = plt.subplots(figsize=(14, 8), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#f7f7f7")

    categories = sorted(gdf["type"].dropna().unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(1, len(categories))))

    for category, color in zip(categories, colors):
        subset = gdf[gdf["type"] == category]
        if subset.empty:
            continue
        ax.scatter(
            subset.geometry.x,
            subset.geometry.y,
            s=subset["marker_size"],
            c=[color],
            alpha=0.75,
            edgecolors="none",
            label=category,
        )

    ax.set_xlim(-170, -60)
    ax.set_ylim(10, 75)
    ax.set_title("IM3 PNNL data centers (all records)", fontsize=14, fontweight="bold")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True, alpha=0.2)
    ax.legend(title="Type", loc="upper left", frameon=True, fontsize=9)
    ax.set_aspect("equal", adjustable="box")
    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    mpl.use("Agg")
    gdf = load_dataset()
    output_path = REPO_ROOT / "outputs" / "im3_pnnl_data_centers_map.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plot_dataset(gdf, output_path)
    print(f"Saved map to {output_path}")
    print(f"Loaded {len(gdf)} records")
    print("Columns:", ", ".join(gdf.columns.tolist()))


if __name__ == "__main__":
    main()
