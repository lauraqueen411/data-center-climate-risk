"""Heat-index calculations from LOCA2 tasmin/tasmax files.

Assumptions
-----------
- Inputs are LOCA2 daily near-surface air temperature fields in Kelvin.
- Output indices are computed in Celsius after explicit unit conversion.
- Heatwave events are defined as runs of at least ``min_duration_days`` where
  daily Tmax exceeds an absolute threshold.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr
import shapely
from shapely.geometry.base import BaseGeometry


@dataclass(frozen=True)
class HeatIndexConfig:
    """Configuration for heat index calculations."""

    cdd_base_c: float = 18.3
    threshold_c: float = 35.0
    min_duration_days: int = 3
    lapse_rate_c_per_1000m: float = 6.5


def _open_temperature_dataarray(path: Path | str, expected_var: str) -> xr.DataArray:
    """Open temperature variable as DataArray from either DataArray/Dataset files."""
    ds = xr.open_dataset(path)
    if expected_var in ds.data_vars:
        da = ds[expected_var]
    elif len(ds.data_vars) == 1:
        da = next(iter(ds.data_vars.values()))
    else:
        raise KeyError(
            f"Could not find expected variable '{expected_var}' in {path}. "
            f"Available variables: {list(ds.data_vars)}"
        )
    return da


def kelvin_to_celsius(da_k: xr.DataArray) -> xr.DataArray:
    """Convert Kelvin to Celsius with explicit arithmetic."""
    da_c = da_k - 273.15
    da_c.attrs.update(da_k.attrs)
    da_c.attrs["units"] = "degC"
    return da_c


def load_tasmax_tasmin_c(
    tasmax_path: Path | str,
    tasmin_path: Path | str,
) -> tuple[xr.DataArray, xr.DataArray]:
    """Load tasmax and tasmin and convert from Kelvin to Celsius."""
    tasmax_k = _open_temperature_dataarray(tasmax_path, "tasmax")
    tasmin_k = _open_temperature_dataarray(tasmin_path, "tasmin")
    return kelvin_to_celsius(tasmax_k), kelvin_to_celsius(tasmin_k)


def select_year_range(da: xr.DataArray, start_year: int, end_year: int) -> xr.DataArray:
    """Subset DataArray to inclusive year range."""
    year = da["time"].dt.year
    return da.sel(time=(year >= start_year) & (year <= end_year))


def annual_cdd(tasmax_c: xr.DataArray, tasmin_c: xr.DataArray, base_c: float) -> xr.DataArray:
    """Annual Cooling Degree Days from daily mean temperature in Celsius."""
    tmean_c = (tasmax_c + tasmin_c) / 2.0
    cdd_daily = xr.where(tmean_c > base_c, tmean_c - base_c, 0.0)
    return cdd_daily.groupby("time.year").sum("time")


def annual_tmax_max(tasmax_c: xr.DataArray) -> xr.DataArray:
    """Annual maximum daily Tmax."""
    return tasmax_c.groupby("time.year").max("time")


def annual_tmax_threshold_days(tasmax_c: xr.DataArray, threshold_c: float) -> xr.DataArray:
    """Annual count of days with Tmax above absolute threshold."""
    above = tasmax_c > threshold_c
    return above.groupby("time.year").sum("time")


def annual_heatwave_event_count(
    tasmax_c: xr.DataArray,
    threshold_c: float,
    min_duration_days: int,
) -> xr.DataArray:
    """Annual heatwave event count for threshold exceedance runs."""
    hot = tasmax_c > threshold_c

    # A run becomes "qualifying" at the first day where a rolling window of
    # length min_duration_days is all True. Counting transitions into that
    # state yields one event per qualifying run.
    def _events_for_year(year_hot: xr.DataArray) -> xr.DataArray:
        qualifying = (
            year_hot.rolling(time=min_duration_days, min_periods=min_duration_days).sum()
            >= min_duration_days
        )
        starts = qualifying & (~qualifying.shift(time=1, fill_value=False))
        return starts.sum("time").astype(np.int16)

    return hot.groupby("time.year").map(_events_for_year)


def period_mean_annual_index(annual_da: xr.DataArray) -> xr.DataArray:
    """Mean across annual index stack for a single period map output."""
    out = annual_da.mean("year")
    out.attrs.update(annual_da.attrs)
    return out


def normalize_longitude_180(da: xr.DataArray) -> xr.DataArray:
    """Normalize longitude coordinate to [-180, 180] if needed."""
    if "lon" not in da.coords:
        raise KeyError("Expected longitude coordinate named 'lon'.")

    lon = da["lon"]
    if float(lon.max()) > 180.0:
        new_lon = ((lon + 180.0) % 360.0) - 180.0
        da = da.assign_coords(lon=new_lon).sortby("lon")
    return da


def clip_to_bbox(
    da: xr.DataArray,
    lon_min: float,
    lon_max: float,
    lat_min: float,
    lat_max: float,
) -> xr.DataArray:
    """Clip DataArray to lat/lon bounding box."""
    da = normalize_longitude_180(da)
    lat_asc = bool(da["lat"][0] < da["lat"][-1])
    lon_asc = bool(da["lon"][0] < da["lon"][-1])

    lat_slice = slice(lat_min, lat_max) if lat_asc else slice(lat_max, lat_min)
    lon_slice = slice(lon_min, lon_max) if lon_asc else slice(lon_max, lon_min)
    return da.sel(lat=lat_slice, lon=lon_slice)


def mask_to_polygon(da_2d: xr.DataArray, polygon: BaseGeometry) -> xr.DataArray:
    """Mask 2-D (lat, lon) DataArray by polygon containment at cell centers."""
    if tuple(da_2d.dims) != ("lat", "lon"):
        raise ValueError(f"Expected dims ('lat', 'lon'), got {da_2d.dims}")

    lon2d, lat2d = np.meshgrid(da_2d["lon"].values, da_2d["lat"].values)
    if hasattr(shapely, "contains_xy"):
        inside = shapely.contains_xy(polygon, lon2d, lat2d)
    else:
        from shapely import vectorized

        inside = vectorized.contains(polygon, lon2d, lat2d)
    return da_2d.where(inside)


def compute_period_heat_indices(
    tasmax_path: Path | str,
    tasmin_path: Path | str,
    start_year: int,
    end_year: int,
    config: HeatIndexConfig,
) -> dict[str, xr.DataArray]:
    """Compute required period-mean annual heat indices from file paths."""
    tasmax_c, tasmin_c = load_tasmax_tasmin_c(tasmax_path, tasmin_path)
    tasmax_c = select_year_range(tasmax_c, start_year, end_year)
    tasmin_c = select_year_range(tasmin_c, start_year, end_year)

    cdd_annual = annual_cdd(tasmax_c, tasmin_c, config.cdd_base_c)
    tmax_max_annual = annual_tmax_max(tasmax_c)
    tmax_days_annual = annual_tmax_threshold_days(tasmax_c, config.threshold_c)
    hw_events_annual = annual_heatwave_event_count(
        tasmax_c,
        threshold_c=config.threshold_c,
        min_duration_days=config.min_duration_days,
    )

    cdd = period_mean_annual_index(cdd_annual)
    cdd.attrs["long_name"] = "Mean annual Cooling Degree Days"
    cdd.attrs["units"] = "degC day"

    tmax_max = period_mean_annual_index(tmax_max_annual)
    tmax_max.attrs["long_name"] = "Mean annual Tmax maximum"
    tmax_max.attrs["units"] = "degC"

    tmax_days = period_mean_annual_index(tmax_days_annual)
    tmax_days.attrs["long_name"] = f"Mean annual days with Tmax > {config.threshold_c}C"
    tmax_days.attrs["units"] = "days/year"

    heatwave_events = period_mean_annual_index(hw_events_annual)
    heatwave_events.attrs["long_name"] = (
        f"Mean annual heatwave event count (Tmax > {config.threshold_c}C for >= {config.min_duration_days} days)"
    )
    heatwave_events.attrs["units"] = "events/year"

    return {
        "cdd": cdd,
        "tmax_annual_max": tmax_max,
        "tmax_threshold_days": tmax_days,
        "heatwave_events": heatwave_events,
        "tmax_c": tasmax_c,
        "tasmin_c": tasmin_c,
    }


def elevation_adjusted_threshold(
    base_threshold_c: float,
    elevation_m: float,
    config: HeatIndexConfig,
) -> float:
    """Lower an absolute Tmax threshold for a facility's elevation.

    **Placeholder rate — not yet confirmed against ASHRAE guidance.** The
    methods doc calls for adjusting the exceedance threshold downward at
    high-elevation sites ("cooling performance drops") but does not specify
    an exact formula or rate. This implementation stands in with the
    standard-atmosphere (ISA) environmental lapse rate,
    ``config.lapse_rate_c_per_1000m`` (default 6.5 degC/km) — a generic
    free-air temperature-vs-altitude relationship, not a validated
    psychrometric or cooling-equipment derating curve. Treat the output as
    illustrative until a real ASHRAE-referenced rate (or facility-specific
    cooling-performance curve) replaces this default.

    Parameters
    ----------
    base_threshold_c:
        Flat (sea-level) exceedance threshold in Celsius.
    elevation_m:
        Facility elevation in metres above sea level.
    config:
        :class:`HeatIndexConfig`; only ``lapse_rate_c_per_1000m`` is used.

    Returns
    -------
    float
        ``base_threshold_c`` reduced by the lapse rate applied over
        ``elevation_m``.
    """
    return base_threshold_c - (config.lapse_rate_c_per_1000m / 1000.0) * elevation_m


def threshold_frequency_summary(
    tmax_c: xr.DataArray,
    threshold_c: float,
) -> dict[str, float]:
    """Summarize annual threshold-day frequency for simple plausibility checks."""
    annual_days = annual_tmax_threshold_days(tmax_c, threshold_c)
    pooled = annual_days.values.ravel()
    pooled = pooled[np.isfinite(pooled)]
    if pooled.size == 0:
        return {"p05": np.nan, "p50": np.nan, "p95": np.nan, "max": np.nan}

    return {
        "p05": float(np.percentile(pooled, 5)),
        "p50": float(np.percentile(pooled, 50)),
        "p95": float(np.percentile(pooled, 95)),
        "max": float(np.max(pooled)),
    }
