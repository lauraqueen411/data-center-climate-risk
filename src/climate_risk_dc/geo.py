"""Geospatial helper utilities shared across exposure loaders.

These helpers keep CRS handling explicit and consistent. Distance and area
calculations are always performed in an equal-area projected CRS
(:data:`climate_risk_dc.config.CRS_EQUAL_AREA`) rather than in degrees.
"""

from __future__ import annotations

import geopandas as gpd
import pandas as pd

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
