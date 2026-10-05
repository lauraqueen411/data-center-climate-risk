"""Unit tests for the ecoregion-labelling and representative-facility helpers.

Uses small synthetic polygons/points, not the real EPA ecoregion shapefile
-- these are illustrative fixtures only, exercising the join and
nearest-to-centroid logic, not real Oregon ecoregion boundaries or facility
locations.
"""

from __future__ import annotations

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, Polygon

from climate_risk_dc.geo import assign_us_l3_ecoregion, select_representative_facility


def _synthetic_ecoregions() -> gpd.GeoDataFrame:
    """Two adjacent 1x1-degree square ecoregions sharing the lon=0 edge."""
    west = Polygon([(-1, -1), (0, -1), (0, 1), (-1, 1)])
    east = Polygon([(0, -1), (1, -1), (1, 1), (0, 1)])
    return gpd.GeoDataFrame(
        {"US_L3CODE": ["1", "2"], "US_L3NAME": ["West Region", "East Region"]},
        geometry=[west, east],
        crs="EPSG:4326",
    )


def test_assign_us_l3_ecoregion_labels_interior_points() -> None:
    assets = gpd.GeoDataFrame(
        {"facility_id": [1, 2]},
        geometry=[Point(-0.5, 0.5), Point(0.5, 0.5)],
        crs="EPSG:4326",
    )
    labeled = assign_us_l3_ecoregion(assets, _synthetic_ecoregions())
    assert list(labeled["US_L3NAME"]) == ["West Region", "East Region"]


def test_assign_us_l3_ecoregion_outside_all_polygons_is_nan() -> None:
    assets = gpd.GeoDataFrame(
        {"facility_id": [1]},
        geometry=[Point(10.0, 10.0)],
        crs="EPSG:4326",
    )
    labeled = assign_us_l3_ecoregion(assets, _synthetic_ecoregions())
    assert pd.isna(labeled["US_L3NAME"].iloc[0])


def test_select_representative_facility_picks_nearest_to_group_centroid() -> None:
    # Group "A": three points on a line at x=0,1,2m (after equal-area
    # reprojection) -- centroid sits at x=1m, so the middle point should win.
    assets = gpd.GeoDataFrame(
        {"facility_id": [1, 2, 3], "region": ["A", "A", "A"]},
        geometry=[Point(-123.0, 45.0), Point(-122.99991, 45.0), Point(-122.99982, 45.0)],
        crs="EPSG:4326",
    )
    rep = select_representative_facility(assets, "region")
    assert len(rep) == 1
    assert rep["facility_id"].iloc[0] == 2


def test_select_representative_facility_single_member_group_returns_that_member() -> None:
    assets = gpd.GeoDataFrame(
        {"facility_id": [1], "region": ["B"]},
        geometry=[Point(-120.0, 44.0)],
        crs="EPSG:4326",
    )
    rep = select_representative_facility(assets, "region")
    assert rep["facility_id"].iloc[0] == 1
    assert rep["n_facilities_in_group"].iloc[0] == 1
    assert rep["distance_to_centroid_m"].iloc[0] == 0.0


def test_select_representative_facility_drops_missing_group_rows() -> None:
    assets = gpd.GeoDataFrame(
        {"facility_id": [1, 2], "region": ["A", None]},
        geometry=[Point(-120.0, 44.0), Point(-119.0, 44.0)],
        crs="EPSG:4326",
    )
    rep = select_representative_facility(assets, "region")
    assert list(rep["facility_id"]) == [1]
