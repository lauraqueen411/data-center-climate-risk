"""Helpers for loading and plotting fiber-provider density in Oregon."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import h3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Polygon

from ..config import CRS_GEOGRAPHIC, load_dataset_paths


def load_fiber_oregon() -> pd.DataFrame:
    """Load the fiber dataset and filter to Oregon rows."""
    datasets = load_dataset_paths()
    path = datasets.fiber_csv
    if not path.exists():
        raise FileNotFoundError(f"Fiber CSV not found: {path}")

    df = pd.read_csv(path)
    return df[df["state_usps"] == "OR"].copy()


def plot_fiber_density_oregon(
    df: pd.DataFrame,
    output_path: str | Path | None = None,
) -> None:
    """Create a geographic map of fiber-provider counts in Oregon using H3 hexes."""
    if df.empty:
        raise ValueError("Fiber dataframe is empty.")

    density = (
        df.groupby("h3_res8_id", dropna=False)
        .agg(provider_count=("provider_id", "nunique"), location_count=("location_id", "count"))
        .reset_index()
    )
    density = density.dropna(subset=["h3_res8_id"])
    if density.empty:
        raise ValueError("No H3 hex IDs available for plotting.")

    geoms = []
    values = []
    for _, row in density.iterrows():
        try:
            boundary = h3.cell_to_boundary(row["h3_res8_id"])
            coords = [(lng, lat) for lat, lng in boundary]
            coords.append(coords[0])
            geoms.append(Polygon(coords))
            values.append(row["provider_count"])
        except Exception:
            continue

    if not geoms:
        raise ValueError("Could not convert any H3 cells to polygons.")

    hex_gdf = gpd.GeoDataFrame({"provider_count": values}, geometry=geoms, crs=CRS_GEOGRAPHIC)
    hex_gdf = hex_gdf.to_crs(CRS_GEOGRAPHIC)

    fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
    hex_gdf.plot(
        ax=ax,
        column="provider_count",
        cmap="viridis",
        edgecolor="white",
        linewidth=0.3,
        legend=True,
    )
    ax.set_title("Fiber provider density in Oregon (H3 res8)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_axis_off()
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
