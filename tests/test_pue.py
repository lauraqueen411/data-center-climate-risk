"""Unit tests for the ported Lei & Masanet AE-chiller PUE model.

Mirrors ``tests/test_geo_elevation.py``'s synthetic-fixture style: pure-math
checks on known inputs, not real facility data.
"""

from __future__ import annotations

import numpy as np
import pytest

from climate_risk_dc.climate.pue import (
    PUEModelConfig,
    air_density,
    chiller_cop,
    check_pue_plausible,
    dewpoint_temperature_c,
    pue_ae_chiller,
    supply_air_temperature_c,
)

# Exact worked-example parameters_AE vector from the source repo's own
# "Simple examples for using the PUE simulation tool.ipynb"
# (https://github.com/nuoaleon/Data-center-PUE-prediction-tool), in the
# w_aech positional order documented in PUE_simulation_tool.ipynb's
# PUE_AE_Chiller.
_SOURCE_REPO_PARAM_NAMES = [
    "ups_e",
    "pd_lr",
    "l_percentage",
    "t_sa_lower_c",
    "t_sa_upper_c",
    "dp_sa_lower_c",
    "dp_sa_upper_c",
    "rh_sa_lower",
    "rh_sa_upper",
    "delta_t_air_c",
    "delta_t_water_c",
    "delta_t_ct_c",
    "fan_pressure_crac_pa",
    "fan_e_crac",
    "fan_pressure_ct_pa",
    "fan_e_ct",
    "pump_pressure_hd_pa",
    "pump_e_hd",
    "pump_pressure_cw_pa",
    "pump_e_cw",
    "pump_pressure_ct_pa",
    "pump_e_ct",
    "at_ct_c",
    "chiller_load",
    "p_cop",
    "l_to_g",
]

_SOURCE_REPO_PARAM_VALUES = [
    9.80302734e-01,
    5.60546875e-03,
    1.81445313e-03,
    1.03671875e01,
    3.41953125e01,
    -1.04970703e01,
    1.58320312e01,
    9.01953125e-02,
    6.76464844e-01,
    1.47840712e01,
    6.91894531e00,
    4.28710938e00,
    6.85058594e02,
    8.29931641e-01,
    2.03417969e02,
    8.97802734e-01,
    7.31582031e06,
    7.04492187e-01,
    1.46519937e05,
    7.90039062e-01,
    2.12001820e05,
    7.16601562e-01,
    6.33818359e00,
    2.24023438e-01,
    1.05488281e-01,
    1.11914062e00,
]

# Reported by the source repo's own worked example:
# PUE_AE_Chiller(parameters_AE, 20, 40, 101325, 1)
_SOURCE_REPO_EXPECTED_PUE = 1.0746922969797756


def test_regression_matches_source_repo_worked_example() -> None:
    """The port must reproduce the source repo's own worked-example PUE exactly."""
    config = PUEModelConfig(**dict(zip(_SOURCE_REPO_PARAM_NAMES, _SOURCE_REPO_PARAM_VALUES)))

    pue = pue_ae_chiller(20.0, 40.0, 101325.0, 1.0, config)

    assert pue.shape == (1,)
    np.testing.assert_allclose(pue[0], _SOURCE_REPO_EXPECTED_PUE, rtol=1e-9)


def test_default_config_plausible_at_mild_conditions() -> None:
    config = PUEModelConfig()
    pue = pue_ae_chiller(20.0, 40.0, config.p_atm_pa, config.power_it_kw, config)
    assert 1.0 <= pue[0] <= 1.3


def test_pue_increases_with_outdoor_temperature() -> None:
    """Hotter, more humid outdoor air should never produce a lower PUE."""
    config = PUEModelConfig()
    t = np.array([0.0, 15.0, 30.0, 40.0])
    rh = np.full_like(t, 50.0)

    pue = pue_ae_chiller(t, rh, config.p_atm_pa, config.power_it_kw, config)

    assert np.all(np.diff(pue) >= 0)


def test_pue_array_input_matches_scalar_calls() -> None:
    config = PUEModelConfig()
    t = np.array([5.0, 25.0, 35.0])
    rh = np.array([70.0, 45.0, 20.0])

    batched = pue_ae_chiller(t, rh, config.p_atm_pa, config.power_it_kw, config)
    scalars = np.array(
        [pue_ae_chiller(float(ti), float(rhi), config.p_atm_pa, config.power_it_kw, config)[0] for ti, rhi in zip(t, rh)]
    )

    np.testing.assert_allclose(batched, scalars)


def test_check_pue_plausible_accepts_typical_values() -> None:
    check_pue_plausible(np.array([1.05, 1.12, 1.30]))


def test_check_pue_plausible_rejects_out_of_range() -> None:
    with pytest.raises(ValueError):
        check_pue_plausible(np.array([1.1, 0.5, 1.2]))
    with pytest.raises(ValueError):
        check_pue_plausible(np.array([1.1, 5.0, 1.2]))


def test_dewpoint_temperature_below_dry_bulb() -> None:
    """Dewpoint must be at or below dry-bulb temperature for any RH <= 100%."""
    t_dp = dewpoint_temperature_c(25.0, 60.0)
    assert t_dp < 25.0


def test_supply_air_temperature_clips_to_bounds() -> None:
    assert supply_air_temperature_c(10.0, 30.0, 5.0) == 10.0
    assert supply_air_temperature_c(10.0, 30.0, 35.0) == 30.0
    assert supply_air_temperature_c(10.0, 30.0, 20.0) == 20.0


def test_air_density_decreases_with_temperature() -> None:
    assert air_density(101325.0, 30.0) < air_density(101325.0, 0.0)


def test_chiller_cop_is_deterministic_regression() -> None:
    cop = chiller_cop(0.5, 20.0)
    assert cop == pytest.approx(16.798751 + 0.008840 * 400 - 8.528931 * 0.25 + 9.885294 * 0.5 - 0.760621 * 20 + 0.084615 * 0.5 * 20)
