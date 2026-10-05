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


# Approximate bounding boxes for informal Oregon sub-regions, used only for
# exploratory grouping in notebooks -- not authoritative geographic or
# administrative boundaries (no county/region shapefile exists in this repo).
# Derived by inspecting the actual facility_lon/facility_lat clusters in
# outputs/month2_heat/facility_heat_indices.csv (109 Oregon data centers).
# Checked in order; first match wins. "Columbia Gorge" and "Columbia Basin /
# Umatilla" are kept separate even though both sit along the Columbia River
# corridor -- the Gorge box covers the river canyon proper (Troutdale-to-
# Dalles), while Columbia Basin covers the geographically distinct
# Boardman/Umatilla high-desert cluster ~100mi further east (the bulk of the
# AWS us-west-2 footprint, ~59 of 109 facilities). Verified zero-unclassified,
# zero-ambiguous against all 109 facilities as of 2026-08-24.
OREGON_REGION_BOXES: list[tuple[str, float, float, float, float]] = [
    # (region, lon_min, lon_max, lat_min, lat_max)
    ("Portland area", -123.1, -122.55, 45.4, 45.65),
    ("Willamette Valley / Southwest", -123.5, -122.3, 42.0, 45.4),
    ("Columbia Gorge", -121.5, -120.9, 45.5, 45.75),
    ("Central Oregon", -121.5, -120.3, 43.8, 44.6),
    ("Columbia Basin / Umatilla", -119.8, -119.0, 45.75, 46.0),
]

OREGON_REGION_ORDER: list[str] = [region for region, *_ in OREGON_REGION_BOXES]


def assign_oregon_region(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Label each point with an approximate Oregon sub-region.

    Coarse lon/lat bounding-box classification against
    :data:`OREGON_REGION_BOXES` -- an informal grouping for exploratory
    plots, not a spatial join against real administrative or physiographic
    boundaries. A point outside all five boxes is labelled
    ``"Unclassified"`` rather than silently dropped or mis-assigned into the
    nearest box.

    Parameters
    ----------
    lon, lat:
        Point coordinates (decimal degrees).

    Returns
    -------
    numpy.ndarray
        Region label per point (dtype object), one of
        :data:`OREGON_REGION_ORDER` or ``"Unclassified"``.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    labels = np.full(lon.shape, "Unclassified", dtype=object)
    for region, lon_min, lon_max, lat_min, lat_max in OREGON_REGION_BOXES:
        in_box = (
            (labels == "Unclassified")
            & (lon >= lon_min)
            & (lon <= lon_max)
            & (lat >= lat_min)
            & (lat <= lat_max)
        )
        labels[in_box] = region
    return labels


def assign_us_l3_ecoregion(
    assets: gpd.GeoDataFrame, ecoregions: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """Label each asset with its EPA Level III ecoregion (Omernik 1987 scheme).

    Thin wrapper around :func:`spatial_label` for the specific
    ``US_L3CODE``/``US_L3NAME`` fields in EPA's ecoregion shapefiles (see
    :data:`climate_risk_dc.config.DatasetPaths.ecoregion_l3_shp`) -- kept
    generic (point-in-polygon join, not ecoregion-specific logic) so it stays
    a thin pass-through rather than duplicating ``spatial_label``.

    An asset outside every polygon gets ``NaN`` for both fields (from
    ``spatial_label``'s left join), unlike the informal
    :func:`assign_oregon_region`, which labels such points
    ``"Unclassified"`` -- the two functions differ deliberately: the EPA
    shapefile is expected to fully cover the state, so an unmatched point
    signals a real data problem worth surfacing as NaN, not a named category.

    Parameters
    ----------
    assets:
        Point geometries to label (e.g. data center facilities).
    ecoregions:
        EPA Level III ecoregion polygons with ``US_L3CODE``/``US_L3NAME``
        columns (loaded from ``DatasetPaths.ecoregion_l3_shp``).

    Returns
    -------
    geopandas.GeoDataFrame
        Copy of ``assets`` with ``US_L3CODE`` and ``US_L3NAME`` columns
        appended.
    """
    return spatial_label(assets, ecoregions, ["US_L3CODE", "US_L3NAME"])


def select_representative_facility(
    assets: gpd.GeoDataFrame, group_col: str
) -> pd.DataFrame:
    """Pick the single facility nearest each group's own facility centroid.

    Used to choose one representative data center per ecoregion for the
    sub-regional PUE/WUE breakdown, following Lei & Masanet (2022)'s own
    convention of representing a climate region with one representative
    point run through the model in full, rather than averaging climate
    inputs across many points first (which would run into Jensen's
    inequality problems given the model's thresholded economizer/chiller
    switching behavior -- see ``docs/data-center-chapter-outline.md``).

    Deliberately simple: the "centroid" is the group's own mean facility
    position (not the ecoregion polygon's geographic centroid), computed in
    the equal-area CRS so distances are in metres, not degrees -- consistent
    with the project's existing distance-computation convention
    (:func:`nearest_distance_m`). This ties the representative point to
    where facilities actually cluster (exposure-weighted) rather than an
    arbitrary polygon shape. A group with one member trivially returns that
    member (distance 0).

    Parameters
    ----------
    assets:
        Point geometries with a ``group_col`` column to group by (e.g.
        facilities already labelled with ``US_L3NAME``).
    group_col:
        Column name to group by (e.g. ``"US_L3NAME"``).

    Returns
    -------
    pandas.DataFrame
        One row per group: all of ``assets``' original columns for the
        selected representative facility, plus ``n_facilities_in_group`` and
        ``distance_to_centroid_m``. Rows with a missing ``group_col`` value
        are dropped rather than grouped together.
    """
    ea = to_equal_area(assets).dropna(subset=[group_col]).copy()
    ea["_x"] = ea.geometry.x
    ea["_y"] = ea.geometry.y

    rows = []
    for _, sub in ea.groupby(group_col, observed=True):
        centroid_x, centroid_y = sub["_x"].mean(), sub["_y"].mean()
        distance_m = np.hypot(sub["_x"] - centroid_x, sub["_y"] - centroid_y)
        rep_idx = distance_m.idxmin()
        rep = sub.loc[rep_idx].drop(labels=["_x", "_y"]).to_dict()
        rep["n_facilities_in_group"] = len(sub)
        rep["distance_to_centroid_m"] = float(distance_m.loc[rep_idx])
        rows.append(rep)

    return pd.DataFrame(rows).drop(columns=["geometry"])


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
