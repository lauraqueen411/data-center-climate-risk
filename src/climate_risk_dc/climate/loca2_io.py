"""Shared LOCA2-family archive file resolution and loading helpers.

Used by both ``scripts/run_heat_processing.py`` (LOCA2 temperature,
tasmax/tasmin) and ``scripts/run_pue_processing.py`` (LOCA CMIP6 humidity
companion, hursmax/hursmin) so both scripts resolve native-window files the
same, correct way rather than duplicating (or re-diverging on) the logic.

Assumptions log
----------------
- The temperature and humidity trees on the HPC share the same ssp585
  native-window boundaries (2015-2044, 2045-2074, 2075-2100) and the same
  single historical file (1950-2014), but differ structurally: the
  temperature tree has an extra ``n_west`` regional-split directory level and
  filename suffix, and a different LOCA2 version stamp
  (``v20220413`` vs. humidity's ``v20250218``). ``Loca2PathConfig`` captures
  both differences so one resolver serves both trees.
- ``resolve_loca2_file`` requires the requested period to be fully contained
  in exactly one native-window file and raises otherwise, rather than
  silently picking the "closest" file or truncating -- this replaces the
  ``end >= 2099`` heuristic that silently dropped 2070-2074 from every
  ssp585 heat output. See
  ``docs/tasks/0.2-fix-loca2-window-resolver-and-time-horizons.md``.
- Every period actually used in this project (historical 1985-2014, a subset
  of the single 1950-2014 historical file; mid-century 2045-2074 and
  end-of-century 2075-2100, each equal to a native ssp585 window) is
  satisfiable by exactly one file. A period spanning multiple native windows
  would need multi-file concatenation, which this module intentionally does
  not implement -- ``resolve_loca2_file`` raises rather than guessing.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr

YearSpan = tuple[int, int]


@dataclass(frozen=True)
class Loca2PathConfig:
    """Path-template parameters for one LOCA2-family archive tree.

    ``has_region_split`` and ``version_stamp`` capture the structural
    differences between the temperature and humidity companion trees
    confirmed via ``ls`` on the HPC (see
    ``docs/tasks/1.2-pue-analysis.md``, Task 2): the humidity tree has no
    ``n_west`` regional-split directory level (or filename suffix), and its
    files carry a different LOCA2 version stamp than the temperature files.
    """

    root: Path
    model: str = "ACCESS-CM2"
    member: str = "r1i1p1f1"
    has_region_split: bool = True
    region_label: str = "n_west"
    version_stamp: str = "v20220413"


def temperature_path_config(
    root: Path | str,
    model: str = "ACCESS-CM2",
    member: str = "r1i1p1f1",
) -> Loca2PathConfig:
    """Path config for the LOCA2 temperature tree (tasmax/tasmin, n_west split)."""
    return Loca2PathConfig(
        root=Path(root),
        model=model,
        member=member,
        has_region_split=True,
        region_label="n_west",
        version_stamp="v20220413",
    )


def humidity_path_config(
    root: Path | str,
    model: str = "ACCESS-CM2",
    member: str = "r1i1p1f1",
) -> Loca2PathConfig:
    """Path config for the LOCA CMIP6 humidity companion tree (hursmax/hursmin, no region split)."""
    return Loca2PathConfig(
        root=Path(root),
        model=model,
        member=member,
        has_region_split=False,
        version_stamp="v20250218",
    )


def _variable_dir(cfg: Loca2PathConfig, scenario: str, variable: str) -> Path:
    if cfg.has_region_split:
        return (
            cfg.root / cfg.model / cfg.region_label / "0p0625deg" / cfg.member / scenario / variable
        )
    return cfg.root / cfg.model / "0p0625deg" / cfg.member / scenario / variable


def _year_span(path: Path) -> YearSpan | None:
    m = re.search(r"\.(\d{4})-(\d{4})\.LOCA", path.name)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def _daily_files(cfg: Loca2PathConfig, scenario: str, variable: str) -> list[Path]:
    """All non-aggregate (not monthly/yearly) files for one variable/scenario."""
    var_dir = _variable_dir(cfg, scenario, variable)
    all_matches = sorted(var_dir.glob(f"{variable}.{cfg.model}.{scenario}.{cfg.member}*.nc"))
    daily = [p for p in all_matches if ".monthly." not in p.name and ".yearly." not in p.name]
    if not daily:
        raise FileNotFoundError(f"No daily LOCA2 files found in {var_dir}")
    return daily


def resolve_loca2_windows(
    cfg: Loca2PathConfig,
    scenario: str,
    variable: str,
    start_year: int,
    end_year: int,
) -> list[Path]:
    """Return native-window files overlapping ``[start_year, end_year]``, sorted by window start."""
    daily = _daily_files(cfg, scenario, variable)
    overlapping: list[tuple[YearSpan, Path]] = []
    for p in daily:
        span = _year_span(p)
        if span is None:
            continue
        if span[0] <= end_year and span[1] >= start_year:
            overlapping.append((span, p))

    if not overlapping:
        available = sorted({s for s in (_year_span(p) for p in daily) if s is not None})
        raise FileNotFoundError(
            f"No native LOCA2 window overlaps requested period {start_year}-{end_year} for "
            f"{variable}/{scenario}. Available native windows: {available}."
        )
    overlapping.sort(key=lambda item: item[0])
    return [p for _, p in overlapping]


def resolve_loca2_file(
    cfg: Loca2PathConfig,
    scenario: str,
    variable: str,
    start_year: int,
    end_year: int,
) -> Path:
    """Resolve exactly one native-window file that fully contains ``[start_year, end_year]``.

    Raises ``ValueError`` rather than picking the "closest" file or silently
    truncating if the requested period isn't fully covered by exactly one
    native window -- see this module's Assumptions log.
    """
    candidates = resolve_loca2_windows(cfg, scenario, variable, start_year, end_year)
    containing = [
        p
        for p in candidates
        if (span := _year_span(p)) is not None and span[0] <= start_year and span[1] >= end_year
    ]
    if len(containing) != 1:
        raise ValueError(
            f"Expected exactly one native LOCA2 file window containing {start_year}-{end_year} "
            f"for {variable}/{scenario}; found {len(containing)} (of {len(candidates)} overlapping "
            f"candidates: {[p.name for p in candidates]}). A period spanning multiple native "
            "windows needs multi-file concatenation, which this function intentionally does not "
            "perform."
        )
    return containing[0]


def _open_variable(path: Path, expected_var: str) -> xr.DataArray:
    ds = xr.open_dataset(path)
    if expected_var in ds.data_vars:
        return ds[expected_var]
    if len(ds.data_vars) == 1:
        return next(iter(ds.data_vars.values()))
    raise KeyError(
        f"Could not find expected variable '{expected_var}' in {path}. "
        f"Available variables: {list(ds.data_vars)}"
    )


def load_hursmax_hursmin_pct(
    hursmax_path: Path,
    hursmin_path: Path,
) -> tuple[xr.DataArray, xr.DataArray]:
    """Load LOCA CMIP6 humidity companion hursmax/hursmin (relative humidity, percent).

    No unit conversion needed -- already in percent, unlike LOCA2
    temperature's Kelvin-to-Celsius conversion in ``heat_indices.py``. The
    on-disk variable names are actually ``rel_humid_max``/``rel_humid_min``,
    not ``hursmax``/``hursmin`` (confirmed via inspection, 2026-08-21) --
    ``_open_variable``'s single-data-var fallback resolves this since each
    file has exactly one data variable, so the ``hursmax``/``hursmin``
    argument here is documentation of intent rather than a literal lookup key.
    """
    hursmax = _open_variable(hursmax_path, "hursmax")
    hursmin = _open_variable(hursmin_path, "hursmin")
    return hursmax, hursmin


def sample_datacube_at_points(da: xr.DataArray, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Nearest-neighbor facility time series from an in-memory ``(time, lat, lon)`` DataArray.

    Indexes by nearest grid cell using plain numpy on an already-materialized
    (or cheap-to-materialize -- e.g. already spatially clipped via
    ``heat_indices.clip_to_bbox``) array, rather than xarray's vectorized
    ``.sel(..., method="nearest")``. The latter was measured to be
    impractically slow against the netCDF4 backend for LOCA2-scale grids
    (multi-minute timeout even after spatial clipping, 2026-08-21) -- the
    backend does not support efficient combined multi-dimension fancy
    indexing, so this trades a slightly less idiomatic API for a method that
    actually finishes. Caller is responsible for materializing/clipping
    ``da`` to a tractable size first; the humidity companion files in
    particular are CONUS-wide (unlike the pre-clipped LOCA2 temperature
    files) and must be clipped before this is called.

    Parameters
    ----------
    da:
        DataArray with dims ``(time, lat, lon)``.
    lon, lat:
        Facility coordinates, decimal degrees. If ``da``'s longitude is in
        0-360 convention, normalize both (e.g. via
        ``heat_indices.normalize_longitude_180`` on ``da``) before calling.

    Returns
    -------
    numpy.ndarray
        Shape ``(len(lon), da.sizes["time"])``.
    """
    if tuple(da.dims) != ("time", "lat", "lon"):
        raise ValueError(f"Expected dims ('time', 'lat', 'lon'), got {da.dims}")

    lat_vals = np.asarray(da["lat"].values)
    lon_vals = np.asarray(da["lon"].values)
    lat_idx = np.abs(lat_vals[:, None] - np.asarray(lat)[None, :]).argmin(axis=0)
    lon_idx = np.abs(lon_vals[:, None] - np.asarray(lon)[None, :]).argmin(axis=0)

    arr = da.values  # (time, lat, lon)
    return arr[:, lat_idx, lon_idx].T  # (facility, time)


def _first_finite_value(da: xr.DataArray) -> float:
    """First finite value, read from a small slice rather than the full array.

    The humidity companion files are CONUS-wide and unclipped (unlike the
    n_west-clipped LOCA2 temperature files) -- a full historical file is
    ~40GB. Materializing the whole array just to find one sample value would
    be needlessly slow (and was observed to hang in practice, 2026-08-21), so
    this restricts to a single time step before pulling values into memory.
    """
    if "time" in da.dims:
        da = da.isel(time=0)
    arr = np.asarray(da.values, dtype=float).ravel()
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        raise ValueError("No finite values found for unit sanity check.")
    return float(finite[0])


def log_humidity_sanity(hursmax_path: Path, hursmin_path: Path) -> None:
    """Sanity-check sampled humidity values fall within plausible [0, 100]% bounds.

    Mirrors ``run_heat_processing.py``'s Kelvin-bounds check (``_log_unit_sanity``).
    """
    hursmax, hursmin = load_hursmax_hursmin_pct(hursmax_path, hursmin_path)
    sample_max = _first_finite_value(hursmax)
    sample_min = _first_finite_value(hursmin)
    print(f"Humidity sample hursmax={sample_max:.3f}%, hursmin={sample_min:.3f}%")
    for v in (sample_max, sample_min):
        if v < 0.0 or v > 100.0:
            raise ValueError(f"Humidity sample {v:.3f}% is outside plausible bounds [0, 100].")
