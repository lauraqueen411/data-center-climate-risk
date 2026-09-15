"""Unit tests for the ported Lei & Masanet (2022) gridded cooling archetypes.

Mirrors ``tests/test_pue.py``'s style. Regression values below were computed
from this module itself after validating it against the (index-corrected,
per this module's Assumptions log) scalar reference implementation in
``nuoaleon/Data-Center-Water-footprint``'s ``simulation_funs_DC.py`` on 60
random (T_oa, RH_oa, P_oa) points spanning -15 to 45 degC, 5-100% RH, and
85-101 kPa -- all three cases matched to floating-point precision
(max error ~1e-13). That cross-check script is not part of this repo (it
imports the cloned reference repo directly), so the values here are frozen
regression points rather than a live comparison, same pattern as
``test_pue.py``'s ``_SOURCE_REPO_EXPECTED_PUE``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from climate_risk_dc.climate.cooling_archetypes import (
    Case1Config,
    Case2Config,
    Case5Config,
    check_pue_wue_plausible,
    elevation_adjusted_pressure_pa,
    load_cop_water_cooled_gp,
    pue_wue_case1,
    pue_wue_case2,
    pue_wue_case5,
)

# data/ is gitignored (present on disk, not in the repo) per CLAUDE.md -- same
# pattern as elevation.LOCA_2016-04-02.nc; skip rather than fail if a fresh
# checkout hasn't fetched it yet.
_COP_2_PATH = Path(__file__).resolve().parent.parent / "data" / "lei_masanet_2022" / "COP_2.pkl"


@pytest.fixture(scope="module")
def cop_gp():
    if not _COP_2_PATH.exists():
        pytest.skip(f"{_COP_2_PATH} not present (gitignored data file, see CLAUDE.md)")
    return load_cop_water_cooled_gp(_COP_2_PATH)


def test_case1_regression_at_mild_conditions(cop_gp) -> None:
    pue, wue = pue_wue_case1(np.array([20.0]), np.array([50.0]), 101325.0, cop_gp, Case1Config())
    np.testing.assert_allclose(pue[0], 1.11376797483027, rtol=1e-9)
    np.testing.assert_allclose(wue[0], 0.0, atol=1e-9)


def test_case1_wue_nonzero_when_economizer_shuts_off(cop_gp) -> None:
    """Regression guard for the Case 1 indexing bug this module corrects.

    The published PUE_WUE_AE_Chiller (uncorrected) returns WUE = 0.0
    unconditionally -- for the archetype the paper itself says has the
    *strongest* climate-driven WUE signal. Hot, humid conditions should
    shut off the airside economizer and produce real cooling-tower water use.
    """
    pue, wue = pue_wue_case1(np.array([38.0]), np.array([70.0]), 101325.0, cop_gp, Case1Config())
    assert wue[0] > 0.5


def test_case2_regression_at_mild_conditions(cop_gp) -> None:
    pue, wue = pue_wue_case2(np.array([20.0]), np.array([50.0]), 101325.0, cop_gp, Case2Config())
    np.testing.assert_allclose(pue[0], 1.1339562609734686, rtol=1e-9)
    np.testing.assert_allclose(wue[0], 2.0875127293876075, rtol=1e-9)


def test_case5_regression_at_mild_conditions(cop_gp) -> None:
    pue, wue = pue_wue_case5(np.array([20.0]), np.array([50.0]), 101325.0, cop_gp, Case5Config())
    np.testing.assert_allclose(pue[0], 1.6634576841003643, rtol=1e-9)
    np.testing.assert_allclose(wue[0], 2.9496357304550784, rtol=1e-9)


@pytest.mark.parametrize(
    "fn,config",
    [
        (pue_wue_case1, Case1Config()),
        (pue_wue_case2, Case2Config()),
        (pue_wue_case5, Case5Config()),
    ],
)
def test_array_input_matches_scalar_calls(fn, config, cop_gp) -> None:
    t = np.array([-10.0, 5.0, 25.0, 40.0])
    rh = np.array([80.0, 70.0, 45.0, 20.0])

    batched_pue, batched_wue = fn(t, rh, 101325.0, cop_gp, config)
    scalar_pue = np.array([fn(np.array([ti]), np.array([rhi]), 101325.0, cop_gp, config)[0][0] for ti, rhi in zip(t, rh)])
    scalar_wue = np.array([fn(np.array([ti]), np.array([rhi]), 101325.0, cop_gp, config)[1][0] for ti, rhi in zip(t, rh)])

    np.testing.assert_allclose(batched_pue, scalar_pue)
    np.testing.assert_allclose(batched_wue, scalar_wue)


def test_check_pue_wue_plausible_accepts_typical_values() -> None:
    check_pue_wue_plausible(np.array([1.1, 1.3, 1.5]), np.array([0.0, 2.0, 4.0]))


def test_check_pue_wue_plausible_rejects_out_of_range() -> None:
    with pytest.raises(ValueError):
        check_pue_wue_plausible(np.array([1.1, 0.5]), np.array([1.0, 1.0]))
    with pytest.raises(ValueError):
        check_pue_wue_plausible(np.array([1.1, 1.2]), np.array([1.0, 50.0]))


def test_elevation_adjusted_pressure_matches_sea_level_standard_atmosphere() -> None:
    np.testing.assert_allclose(elevation_adjusted_pressure_pa(0.0), 101325.0)


def test_elevation_adjusted_pressure_decreases_with_elevation() -> None:
    p = elevation_adjusted_pressure_pa(np.array([0.0, 1000.0, 3400.0]))
    assert np.all(np.diff(p) < 0)
