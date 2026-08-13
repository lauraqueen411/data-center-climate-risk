"""Exposure dataset loaders and the unified asset-library builder.

Each loader function reads a raw dataset, filters to the western-US region of
interest, and returns a clean :class:`~geopandas.GeoDataFrame` with a
consistent CRS (WGS84).

The :func:`build_exposure_library` function in :mod:`library` combines all
loaders into a single GeoDataFrame where each row is a data-center asset
enriched with exposure attributes (nearest fiber distance, water-service-area
membership, watershed HUC codes, etc.).
"""

from __future__ import annotations

try:
    from .data_centers import load_data_centers
except ModuleNotFoundError:  # pragma: no cover - optional loader may be absent
    load_data_centers = None

from .fiber import load_fiber_oregon, plot_fiber_density_oregon

try:
    from .library import build_exposure_library
except ModuleNotFoundError:  # pragma: no cover - optional loader may be absent
    build_exposure_library = None

from .water import load_public_water_sources, plot_public_water_sources

try:
    from .watersheds import load_watersheds
except ModuleNotFoundError:  # pragma: no cover - optional loader may be absent
    load_watersheds = None

__all__ = [
    "load_public_water_sources",
    "plot_public_water_sources",
    "load_fiber_oregon",
    "plot_fiber_density_oregon",
]
if load_data_centers is not None:
    __all__.append("load_data_centers")
if build_exposure_library is not None:
    __all__.append("build_exposure_library")
if load_watersheds is not None:
    __all__.append("load_watersheds")
