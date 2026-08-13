"""Helpers for loading and plotting public water-source geometries."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..config import CRS_GEOGRAPHIC, load_dataset_paths


def load_public_water_sources() -> gpd.GeoDataFrame:
    """Load the public water-service-area shapefile as a GeoDataFrame."""
    datasets = load_dataset_paths()
    path = datasets.water_service_areas_shp
    if not path.exists():
        raise FileNotFoundError(f"Water-source shapefile not found: {path}")

    gdf = gpd.read_file(path)
    if gdf.crs is None:
        gdf = gdf.set_crs(CRS_GEOGRAPHIC)
    return gdf.to_crs(CRS_GEOGRAPHIC)


def plot_public_water_sources(
    gdf: gpd.GeoDataFrame,
    output_path: str | Path | None = None,
) -> None:
    """Create a simple plot of public water-source geometries."""
    fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
    if gdf.crs is None:
        gdf = gdf.set_crs(CRS_GEOGRAPHIC)
    gdf.to_crs(CRS_GEOGRAPHIC).plot(
        ax=ax,
        color="lightsteelblue",
        edgecolor="black",
        linewidth=0.6,
    )
    ax.set_title("Public water sources")
    ax.set_axis_off()
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
