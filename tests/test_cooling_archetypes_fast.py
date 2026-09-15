"""Unit tests for the vectorized-numpy psychrometrics fast path.

Mirrors ``tests/test_cooling_archetypes.py``'s style. Regression values were
computed from this module itself after validating it end-to-end against
the CoolProp-based ``cooling_archetypes.py`` on 60 random (T_oa, RH_oa, P_oa)
points spanning -15 to 45 degC, 5-100% RH, 85-101 kPa: PUE agreed within
0.09-0.35% relative error and WUE within 0.001-0.034 L/kWh absolute error
across all three cases -- see ``cooling_archetypes_fast.py``'s Assumptions
log for the per-primitive validation against CoolProp itself.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from climate_risk_dc.climate.cooling_archetypes import Case1Config, Case2Config, Case5Config
from climate_risk_dc.climate.cooling_archetypes_fast import (
    load_cop_water_cooled_gp,
    pue_wue_case1,
    pue_wue_case2,
    pue_wue_case5,
)

_COP_2_PATH = Path(__file__).resolve().parent.parent / "data" / "lei_masanet_2022" / "COP_2.pkl"


@pytest.fixture(scope="module")
def cop_gp():
    if not _COP_2_PATH.exists():
        pytest.skip(f"{_COP_2_PATH} not present (gitignored data file, see CLAUDE.md)")
    return load_cop_water_cooled_gp(_COP_2_PATH)


def test_case1_regression_at_mild_conditions(cop_gp) -> None:
    pue, wue = pue_wue_case1(np.array([20.0]), np.array([50.0]), 101325.0, cop_gp, Case1Config())
    np.testing.assert_allclose(pue[0], 1.114110035726075, rtol=1e-9)
    np.testing.assert_allclose(wue[0], 0.0, atol=1e-9)


def test_case2_regression_at_mild_conditions(cop_gp) -> None:
    pue, wue = pue_wue_case2(np.array([20.0]), np.array([50.0]), 101325.0, cop_gp, Case2Config())
    np.testing.assert_allclose(pue[0], 1.1343657029660112, rtol=1e-9)
    np.testing.assert_allclose(wue[0], 2.087439076646512, rtol=1e-9)


def test_case5_regression_at_mild_conditions(cop_gp) -> None:
    pue, wue = pue_wue_case5(np.array([20.0]), np.array([50.0]), 101325.0, cop_gp, Case5Config())
    np.testing.assert_allclose(pue[0], 1.6647927042647905, rtol=1e-9)
    np.testing.assert_allclose(wue[0], 2.949159772100321, rtol=1e-9)


@pytest.mark.parametrize(
    "fn,config",
    [
        (pue_wue_case1, Case1Config()),
        (pue_wue_case2, Case2Config()),
        (pue_wue_case5, Case5Config()),
    ],
)
def test_agrees_with_coolprop_reference(fn, config, cop_gp) -> None:
    """Cross-check against the CoolProp-based implementation at a fixed set of points.

    Loose tolerances (this module trades a small, validated accuracy loss --
    see the module docstring -- for ~50x+ speed at scale); this guards
    against a future edit silently introducing a much larger discrepancy,
    not bit-for-bit agreement.
    """
    from climate_risk_dc.climate import cooling_archetypes as ca

    t = np.array([-10.0, 5.0, 25.0, 40.0])
    rh = np.array([80.0, 70.0, 45.0, 20.0])
    slow_fn = {pue_wue_case1: ca.pue_wue_case1, pue_wue_case2: ca.pue_wue_case2, pue_wue_case5: ca.pue_wue_case5}[fn]

    pue_fast, wue_fast = fn(t, rh, 95000.0, cop_gp, config)
    pue_slow, wue_slow = slow_fn(t, rh, 95000.0, cop_gp, config)

    np.testing.assert_allclose(pue_fast, pue_slow, rtol=0.01)
    np.testing.assert_allclose(wue_fast, wue_slow, atol=0.1)


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
