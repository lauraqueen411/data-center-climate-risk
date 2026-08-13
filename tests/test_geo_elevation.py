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
from climate_risk_dc.geo import get_facility_elevation_m, sample_raster_at_points


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
