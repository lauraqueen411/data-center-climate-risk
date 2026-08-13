"""Helpers for loading and plotting watershed boundary data."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..config import CRS_GEOGRAPHIC, load_dataset_paths


def load_watersheds(level: str = "WBDHU4") -> gpd.GeoDataFrame:
    """Load a coarser watershed boundary layer from the unpacked geodatabase."""
    datasets = load_dataset_paths()
    path = datasets.watershed_gdb_zip
    gdb_path = path.parent / "unpacked" / "WBD_National_GDB.gdb"
    if not gdb_path.exists():
        raise FileNotFoundError(f"Unpacked watershed geodatabase not found: {gdb_path}")

    layer_table = gpd.list_layers(gdb_path)
    available_layers = set(layer_table["name"])
    if level not in available_layers:
        raise ValueError(f"Unsupported watershed layer '{level}'. Available layers: {sorted(available_layers)}")

    gdf = gpd.read_file(gdb_path, layer=level)
    if gdf.crs is None:
        gdf = gdf.set_crs(CRS_GEOGRAPHIC)
    return gdf.to_crs(CRS_GEOGRAPHIC)


def plot_watersheds(
    gdf: gpd.GeoDataFrame,
    output_path: str | Path | None = None,
) -> None:
    """Create a simple plot of watershed boundary geometries."""
    fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
    if gdf.crs is None:
        gdf = gdf.set_crs(CRS_GEOGRAPHIC)
    gdf.to_crs(CRS_GEOGRAPHIC).plot(
        ax=ax,
        color="lightgray",
        edgecolor="black",
        linewidth=0.5,
    )
    ax.set_title("Watershed boundaries")
    ax.set_axis_off()
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
