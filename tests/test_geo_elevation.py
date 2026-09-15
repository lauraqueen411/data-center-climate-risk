"""Unit tests for the elevation-adjustment helpers added for the 1.1 close-out.

These tests use a small synthetic in-memory raster, not a real DEM -- no
DEM/elevation raster exists anywhere in this repo yet (see
scripts/extract_facility_heat_indices.py's Assumptions log), so the
elevation and threshold values below are illustrative fixtures only, not
real facility numbers.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

from climate_risk_dc.climate.heat_indices import HeatIndexConfig, elevation_adjusted_threshold
from climate_risk_dc.geo import (
    OREGON_REGION_BOXES,
    assign_oregon_region,
    get_facility_elevation_m,
    sample_raster_at_points,
)


def _write_synthetic_dem(path: Path) -> None:
    """Write a 3x3 synthetic elevation grid with one nodata cell.

    Grid (row=lat descending, col=lon ascending), 1-degree cells centered on
    integer degrees:

              lon=-122  lon=-121  lon=-120
    lat=45         0.0     500.0    1000.0
    lat=44       100.0     600.0   nodata
    lat=43       200.0     700.0    1200.0
    """
    data = np.array(
        [
            [0.0, 500.0, 1000.0],
            [100.0, 600.0, -9999.0],
            [200.0, 700.0, 1200.0],
        ],
        dtype="float32",
    )
    transform = from_origin(-122.5, 45.5, 1.0, 1.0)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=3,
        width=3,
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=transform,
        nodata=-9999.0,
    ) as dst:
        dst.write(data, 1)


def test_sample_raster_at_points_matches_known_grid(tmp_path: Path) -> None:
    dem_path = tmp_path / "synthetic_dem.tif"
    _write_synthetic_dem(dem_path)

    lon = np.array([-122.0, -121.0, -120.0])
    lat = np.array([45.0, 44.0, 43.0])
    values = sample_raster_at_points(dem_path, lon, lat)

    np.testing.assert_allclose(values, [0.0, 600.0, 1200.0])


def test_sample_raster_at_points_nodata_returns_nan(tmp_path: Path) -> None:
    dem_path = tmp_path / "synthetic_dem.tif"
    _write_synthetic_dem(dem_path)

    values = sample_raster_at_points(dem_path, np.array([-120.0]), np.array([44.0]))
    assert np.isnan(values[0])


def test_get_facility_elevation_m_is_thin_wrapper(tmp_path: Path) -> None:
    dem_path = tmp_path / "synthetic_dem.tif"
    _write_synthetic_dem(dem_path)

    values = get_facility_elevation_m(np.array([-122.0]), np.array([45.0]), dem_path)
    assert values[0] == 0.0


def test_elevation_adjusted_threshold_sea_level_unchanged() -> None:
    config = HeatIndexConfig()
    assert elevation_adjusted_threshold(35.0, 0.0, config) == 35.0


def test_elevation_adjusted_threshold_applies_isa_lapse_rate() -> None:
    config = HeatIndexConfig(lapse_rate_c_per_1000m=6.5)
    # Illustrative fixture only: 1000m is a hypothetical elevation, not a
    # real facility's -- no DEM exists yet to pull real elevations from.
    adjusted = elevation_adjusted_threshold(35.0, 1000.0, config)
    assert adjusted == 28.5


def test_assign_oregon_region_matches_known_clusters() -> None:
    lon = np.array([-122.9, -121.2, -120.87, -119.6, -123.0])
    lat = np.array([45.55, 45.63, 44.29, 45.85, 44.93])
    labels = assign_oregon_region(lon, lat)
    assert list(labels) == [
        "Portland area",
        "Columbia Gorge",
        "Central Oregon",
        "Columbia Basin / Umatilla",
        "Willamette Valley / Southwest",
    ]


def test_assign_oregon_region_point_outside_all_boxes_is_unclassified() -> None:
    labels = assign_oregon_region(np.array([-116.0]), np.array([42.0]))
    assert labels[0] == "Unclassified"


def test_oregon_region_boxes_do_not_overlap_in_their_interiors() -> None:
    # Regression guard: if a future edit widens a box until two overlap, a
    # point could silently land in whichever region happens to be checked
    # first instead of raising. Sample points just inside each box's
    # boundary (not the exact edge -- adjacent boxes are allowed to share a
    # touching edge, e.g. Portland/Willamette both include lat=45.4; the
    # function's first-match-wins order resolves that deterministically, and
    # no real facility falls on it -- see assign_oregon_region's docstring
    # and the zero-ambiguous verification against the real 109 facilities).
    eps = 1e-6
    for region, lon_min, lon_max, lat_min, lat_max in OREGON_REGION_BOXES:
        lon_pts = [lon_min + eps, lon_max - eps, (lon_min + lon_max) / 2]
        lat_pts = [lat_min + eps, lat_max - eps, (lat_min + lat_max) / 2]
        for lo in lon_pts:
            for la in lat_pts:
                matches = [
                    r
                    for r, blo, bhi, bla, bha in OREGON_REGION_BOXES
                    if blo <= lo <= bhi and bla <= la <= bha
                ]
                assert matches == [region], (region, lo, la, matches)
