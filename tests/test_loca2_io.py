"""Unit tests for the LOCA2 native-window file resolver.

Builds a synthetic directory tree that mirrors the real HPC layout (three
non-overlapping ssp585 windows, one historical file) so these tests would
have caught the original ``end >= 2099`` bug: requesting the 2045-2074
period must resolve to the 2045-2074 file, not silently fall through to
2075-2100.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from climate_risk_dc.climate.loca2_io import (
    humidity_path_config,
    resolve_loca2_file,
    resolve_loca2_windows,
    temperature_path_config,
)

FUTURE_WINDOWS = ["2015-2044", "2045-2074", "2075-2100"]


def _touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()


def _build_temperature_tree(root: Path, model: str, member: str, variable: str) -> None:
    var_dir = root / model / "n_west" / "0p0625deg" / member
    hist_dir = var_dir / "historical" / variable
    _touch(hist_dir / f"{variable}.{model}.historical.{member}.1950-2014.LOCA_16thdeg_v20220413.n_west.nc")

    fut_dir = var_dir / "ssp585" / variable
    for window in FUTURE_WINDOWS:
        _touch(fut_dir / f"{variable}.{model}.ssp585.{member}.{window}.LOCA_16thdeg_v20220413.n_west.nc")
        # Aggregate siblings that must be excluded from resolution.
        _touch(fut_dir / f"{variable}.{model}.ssp585.{member}.{window}.LOCA_16thdeg_v20220413.monthly.n_west.nc")
        _touch(fut_dir / f"{variable}.{model}.ssp585.{member}.{window}.LOCA_16thdeg_v20220413.yearly.n_west.nc")


def _build_humidity_tree(root: Path, model: str, member: str, variable: str) -> None:
    var_dir = root / model / "0p0625deg" / member
    hist_dir = var_dir / "historical" / variable
    _touch(hist_dir / f"{variable}.{model}.historical.{member}.1950-2014.LOCA_16thdeg_v20250218.nc")

    fut_dir = var_dir / "ssp585" / variable
    for window in FUTURE_WINDOWS:
        _touch(fut_dir / f"{variable}.{model}.ssp585.{member}.{window}.LOCA_16thdeg_v20250218.nc")
        _touch(fut_dir / f"{variable}.{model}.ssp585.{member}.{window}.LOCA_16thdeg_v20250218.monthly.nc")


def test_resolve_mid_century_picks_2045_2074_not_2075_2100(tmp_path: Path) -> None:
    _build_temperature_tree(tmp_path, "ACCESS-CM2", "r1i1p1f1", "tasmax")
    cfg = temperature_path_config(tmp_path, "ACCESS-CM2", "r1i1p1f1")

    resolved = resolve_loca2_file(cfg, "ssp585", "tasmax", 2045, 2074)

    assert resolved.name == "tasmax.ACCESS-CM2.ssp585.r1i1p1f1.2045-2074.LOCA_16thdeg_v20220413.n_west.nc"


def test_resolve_end_of_century_includes_year_2100(tmp_path: Path) -> None:
    _build_temperature_tree(tmp_path, "ACCESS-CM2", "r1i1p1f1", "tasmax")
    cfg = temperature_path_config(tmp_path, "ACCESS-CM2", "r1i1p1f1")

    resolved = resolve_loca2_file(cfg, "ssp585", "tasmax", 2075, 2100)

    assert resolved.name == "tasmax.ACCESS-CM2.ssp585.r1i1p1f1.2075-2100.LOCA_16thdeg_v20220413.n_west.nc"


def test_resolve_rejects_period_not_matching_a_single_window(tmp_path: Path) -> None:
    """Regression guard: the old 2070-2099 period must now fail loudly, not silently truncate."""
    _build_temperature_tree(tmp_path, "ACCESS-CM2", "r1i1p1f1", "tasmax")
    cfg = temperature_path_config(tmp_path, "ACCESS-CM2", "r1i1p1f1")

    with pytest.raises(ValueError):
        resolve_loca2_file(cfg, "ssp585", "tasmax", 2070, 2099)


def test_resolve_historical_subset_of_single_file(tmp_path: Path) -> None:
    _build_temperature_tree(tmp_path, "ACCESS-CM2", "r1i1p1f1", "tasmax")
    cfg = temperature_path_config(tmp_path, "ACCESS-CM2", "r1i1p1f1")

    resolved = resolve_loca2_file(cfg, "historical", "tasmax", 1985, 2014)

    assert resolved.name == "tasmax.ACCESS-CM2.historical.r1i1p1f1.1950-2014.LOCA_16thdeg_v20220413.n_west.nc"


def test_resolve_excludes_monthly_and_yearly_aggregates(tmp_path: Path) -> None:
    _build_temperature_tree(tmp_path, "ACCESS-CM2", "r1i1p1f1", "tasmax")
    cfg = temperature_path_config(tmp_path, "ACCESS-CM2", "r1i1p1f1")

    windows = resolve_loca2_windows(cfg, "ssp585", "tasmax", 2045, 2074)

    assert len(windows) == 1
    assert ".monthly." not in windows[0].name
    assert ".yearly." not in windows[0].name


def test_humidity_resolver_has_no_region_split_and_different_version_stamp(tmp_path: Path) -> None:
    _build_humidity_tree(tmp_path, "ACCESS-CM2", "r1i1p1f1", "hursmin")
    cfg = humidity_path_config(tmp_path, "ACCESS-CM2", "r1i1p1f1")

    resolved = resolve_loca2_file(cfg, "ssp585", "hursmin", 2045, 2074)

    assert resolved.name == "hursmin.ACCESS-CM2.ssp585.r1i1p1f1.2045-2074.LOCA_16thdeg_v20250218.nc"
    assert "n_west" not in str(resolved)


def test_resolve_raises_clear_error_when_no_window_overlaps(tmp_path: Path) -> None:
    _build_temperature_tree(tmp_path, "ACCESS-CM2", "r1i1p1f1", "tasmax")
    cfg = temperature_path_config(tmp_path, "ACCESS-CM2", "r1i1p1f1")

    with pytest.raises(FileNotFoundError):
        resolve_loca2_file(cfg, "ssp585", "tasmax", 2200, 2210)
