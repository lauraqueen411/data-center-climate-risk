"""Geospatial helper utilities shared across exposure loaders.

These helpers keep CRS handling explicit and consistent. Distance and area
calculations are always performed in an equal-area projected CRS
(:data:`climate_risk_dc.config.CRS_EQUAL_AREA`) rather than in degrees.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio

from .config import CRS_EQUAL_AREA, CRS_GEOGRAPHIC


def points_from_lonlat(
    df: pd.DataFrame,
    lon_col: str = "lon",
    lat_col: str = "lat",
    crs: str = CRS_GEOGRAPHIC,
) -> gpd.GeoDataFrame:
    """Build a point :class:`~geopandas.GeoDataFrame` from lon/lat columns.

    Parameters
    ----------
    df:
        Input table containing longitude and latitude columns.
    lon_col, lat_col:
        Names of the longitude and latitude columns (decimal degrees).
    crs:
        CRS of the input coordinates. Defaults to WGS84.

    Returns
    -------
    geopandas.GeoDataFrame
        Copy of ``df`` with a ``geometry`` column of points. Rows with missing
        coordinates are dropped.
    """
    valid = df.dropna(subset=[lon_col, lat_col]).copy()
    geometry = gpd.points_from_xy(valid[lon_col], valid[lat_col])
    return gpd.GeoDataFrame(valid, geometry=geometry, crs=crs)


def load_oregon_data_centers(csv_path: Path) -> gpd.GeoDataFrame:
    """Load the IM3/PNNL data-center atlas CSV, filtered to Oregon facilities.

    Parameters
    ----------
    csv_path:
        Path to the IM3/PNNL open-source data-center atlas CSV.

    Returns
    -------
    geopandas.GeoDataFrame
        Oregon rows (``state_abb == "OR"``) with a point ``geometry`` column
        built from ``lon``/``lat``, in :data:`CRS_GEOGRAPHIC`.
    """
    df = pd.read_csv(csv_path)
    df = df[df["state_abb"] == "OR"].copy()
    return points_from_lonlat(df)


def sample_raster_at_points(raster_path: Path, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Sample a single-band raster's values at a set of lon/lat points.

    Parameters
    ----------
    raster_path:
        Path to a single-band raster readable by ``rasterio``.
    lon, lat:
        Point coordinates (decimal degrees) in the raster's CRS.

    Returns
    -------
    numpy.ndarray
        Sampled band-1 values, one per input point. Points that fall outside
        the raster extent, or land on a nodata cell, are returned as ``NaN``
        rather than the raw nodata sentinel.
    """
    with rasterio.open(raster_path) as src:
        points = list(zip(lon, lat, strict=False))
        values = np.array([v[0] for v in src.sample(points)], dtype=float)
        if src.nodata is not None:
            values = np.where(np.isclose(values, src.nodata), np.nan, values)
    return values


def get_facility_elevation_m(lon: np.ndarray, lat: np.ndarray, dem_path: Path) -> np.ndarray:
    """Sample ground elevation (metres) at facility coordinates from a DEM.

    Thin wrapper around :func:`sample_raster_at_points` so callers don't need
    to know it's raster sampling under the hood. Kept generic (not heat-index
    specific) so Phase 2.2's flood-elevation work can reuse it against the
    same or a different DEM instead of duplicating a sampler.

    Parameters
    ----------
    lon, lat:
        Facility coordinates (decimal degrees), matching the DEM's CRS.
    dem_path:
        Path to a single-band elevation raster (e.g. USGS 3DEP).

    Returns
    -------
    numpy.ndarray
        Elevation in metres per point; ``NaN`` where a point falls outside
        the DEM extent or on a nodata cell.
    """
    return sample_raster_at_points(dem_path, lon, lat)


def to_equal_area(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Reproject to the project equal-area CRS for metric computations."""
    return gdf.to_crs(CRS_EQUAL_AREA)


def to_geographic(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Reproject to WGS84 for web mapping / lon-lat output."""
    return gdf.to_crs(CRS_GEOGRAPHIC)


def nearest_distance_m(
    assets: gpd.GeoDataFrame,
    features: gpd.GeoDataFrame,
    result_col: str = "nearest_distance_m",
    feature_id_col: str | None = None,
) -> gpd.GeoDataFrame:
    """Compute distance from each asset to the nearest feature, in metres.

    Both inputs are reprojected to the equal-area CRS internally so distances
    are in metres regardless of the incoming CRS.

    Parameters
    ----------
    assets:
        Asset geometries (e.g. data centers). The returned frame preserves the
        asset CRS.
    features:
        Feature geometries to measure distance to (e.g. fiber locations).
    result_col:
        Name of the output distance column (metres).
    feature_id_col:
        Optional column in ``features`` whose value for the nearest feature is
        attached to each asset (e.g. a fiber ``location_id``).

    Returns
    -------
    geopandas.GeoDataFrame
        Copy of ``assets`` with the nearest-distance column (and optionally the
        nearest feature id) appended.
    """
    if assets.crs is None or features.crs is None:
        raise ValueError("Both inputs must have a defined CRS.")

    assets_ea = to_equal_area(assets)
    features_ea = to_equal_area(features)[
        [c for c in ([feature_id_col] if feature_id_col else []) + ["geometry"]]
    ]

    joined = gpd.sjoin_nearest(
        assets_ea,
        features_ea,
        how="left",
        distance_col=result_col,
    )
    # sjoin_nearest can emit duplicate rows on ties; keep the first per asset.
    joined = joined[~joined.index.duplicated(keep="first")]

    out = assets.copy()
    out[result_col] = joined[result_col].to_numpy()
    if feature_id_col is not None:
        out[feature_id_col] = joined[feature_id_col].to_numpy()
    return out


def spatial_label(
    assets: gpd.GeoDataFrame,
    polygons: gpd.GeoDataFrame,
    label_cols: list[str],
) -> gpd.GeoDataFrame:
    """Attach polygon attribute(s) to each asset via point-in-polygon join.

    Parameters
    ----------
    assets:
        Point (or other) geometries to label.
    polygons:
        Polygon layer providing the labels.
    label_cols:
        Columns from ``polygons`` to attach to each asset.

    Returns
    -------
    geopandas.GeoDataFrame
        Copy of ``assets`` with ``label_cols`` appended (NaN where an asset
        falls outside all polygons).
    """
    polys = polygons.to_crs(assets.crs)[label_cols + ["geometry"]]
    joined = gpd.sjoin(assets, polys, how="left", predicate="within")
    joined = joined[~joined.index.duplicated(keep="first")]

    out = assets.copy()
    for col in label_cols:
        out[col] = joined[col].to_numpy()
    return out
