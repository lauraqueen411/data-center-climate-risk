"""Power Usage Effectiveness (PUE) model: airside economizer + adiabatic cooling.

Ported from Nuoa Lei's own reference implementation for Lei & Masanet (2020),
*Energy* 201, 117556, https://doi.org/10.1016/j.energy.2020.117556
(https://github.com/nuoaleon/Data-center-PUE-prediction-tool,
``PUE_internalfunctions.ipynb`` + ``PUE_simulation_tool.ipynb``'s
``PUE_AE_Chiller``). This is the "open-source engineering model" the methods
doc cites in ref [8]/[11], and the AE (airside economizer + adiabatic
cooling) variant is confirmed as the model variant used throughout this
project -- see ``docs/tasks/1.2-pue-analysis.md``.

Assumptions log
----------------
- Ported directly from the source notebooks; formulas, variable names, and
  computation order are kept close to the original rather than re-derived,
  so the port can be checked line-by-line against the source. The one
  deliberate departure: the source notebooks do ``from numpy import *``,
  which on current numpy (2.x, where ``max``/``min`` were added to
  ``numpy.__all__``) silently shadows Python's builtin ``max``/``min`` and
  breaks the ``Airside_economizer_use_advanced`` branch that calls
  ``max(0, ...)`` on a scalar. This module uses ``import numpy as np``
  and Python's builtin ``max``/``min`` instead -- a compatibility fix, not a
  modeling change. Verified against the source repo's own worked example
  (``tests/test_pue.py``): ``PUE_AE_Chiller(parameters_AE, 20, 40, 101325, 1)``
  reproduces the repo's reported PUE ~= 1.0746922969797756 exactly.
- ``PUEModelConfig`` defaults to the **midpoint** of each parameter's bounds
  in **Table A.1 of the paper itself** (not just the GitHub sensitivity
  notebook -- see the ``l_to_g`` note below), representing a "typical
  hyperscale AE+adiabatic system," not any specific Oregon facility's real
  equipment spec. 28 of the 29 bounds were checked directly against Table
  A.1 and match the GitHub notebook (``F.1.1_AE_sensitivity analysis.ipynb``,
  ``problem_AE_Chiller``) exactly. The one exception is ``l_to_g``, the
  cooling-tower liquid-gas ratio: the GitHub notebook narrows it to
  ``[1, 3]`` for its own sensitivity run, but Table A.1 states the general
  model's range as ``[0.2, 4]`` (citing refs [83-85] in the paper). This
  module uses the paper's wider range (midpoint 2.1), not the notebook's.
- ``p_atm_pa`` (default 101325 Pa, standard atmosphere) is a marked
  extension point for a future per-facility elevation adjustment once a DEM
  exists in the repo (``geo.get_facility_elevation_m``).
- ``power_it_kw`` defaults to 1 kW, matching the source repo's convention:
  PUE is a dimensionless ratio (total facility power / IT power), so the
  absolute IT power scale does not change the result unless it interacts
  with load-dependent terms (it doesn't, in this model).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import fsolve


@dataclass(frozen=True)
class PUEModelConfig:
    """AE-variant (airside economizer + adiabatic cooling) equipment parameters.

    All 26 core parameters default to the midpoint of the bounds published in
    Table A.1 of Lei & Masanet (2020); see this module's Assumptions log for
    the one correction to the GitHub notebook's narrower ``l_to_g`` bound.
    """

    ups_e: float = 0.945
    pd_lr: float = 0.01
    l_percentage: float = 0.001
    t_sa_lower_c: float = 14.0
    t_sa_upper_c: float = 31.0
    dp_sa_lower_c: float = -10.5
    dp_sa_upper_c: float = 21.0
    rh_sa_lower: float = 0.14
    rh_sa_upper: float = 0.75
    delta_t_air_c: float = 16.67
    delta_t_water_c: float = 7.5
    delta_t_ct_c: float = 5.0
    fan_pressure_crac_pa: float = 750.0
    fan_e_crac: float = 0.775
    fan_pressure_ct_pa: float = 250.0
    fan_e_ct: float = 0.775
    pump_pressure_hd_pa: float = 7_000_000.0
    pump_e_hd: float = 0.7
    pump_pressure_cw_pa: float = 143_658.0
    pump_e_cw: float = 0.7
    pump_pressure_ct_pa: float = 208_660.0
    pump_e_ct: float = 0.7
    at_ct_c: float = 4.75
    chiller_load: float = 0.5
    p_cop: float = 0.0
    l_to_g: float = 2.1

    #: Marked extension point -- see Assumptions log.
    p_atm_pa: float = 101325.0
    #: Reference IT load; PUE is dimensionless so this is a convention, not a
    #: real facility's IT capacity.
    power_it_kw: float = 1.0


def dewpoint_temperature_c(t_c: float, rh_pct: float) -> float:
    """Dewpoint temperature (C) from dry-bulb temperature and relative humidity (%)."""
    log_term = np.log(rh_pct / 100 * np.exp((18.678 - t_c / 234.5) * (t_c / (257.14 + t_c))))
    return 257.14 * log_term / (18.678 - log_term)


def _wetbulb_enthalpy_kj_per_kg(t_wb_c: float, p_atm_pa: float) -> float:
    """Moist-air enthalpy (kJ/kg dry air) at saturation, evaluated at a candidate wetbulb temperature."""
    p_sat_twb = 611.2 * np.exp(((18.678 - t_wb_c / 234.5) * t_wb_c) / (t_wb_c + 257.14))
    d_sat_twb = 0.62198 * p_sat_twb / (p_atm_pa - p_sat_twb)
    return 1.01 * t_wb_c + (2500 + 1.84 * t_wb_c) * d_sat_twb


def wetbulb_temperature_c(t_oa_c: float, rh_oa_pct: float, p_atm_pa: float, p_sat_pa: float) -> float:
    """Outdoor-air wetbulb temperature (C), solved implicitly via enthalpy matching.

    Scalar by construction (uses ``scipy.optimize.fsolve``); ``pue_ae_chiller``
    loops over facilities rather than vectorizing this solve directly.
    """
    p_vapor = p_sat_pa * rh_oa_pct / 100
    d_oa = 0.62198 * p_vapor / (p_atm_pa - p_vapor)
    e_oa = 1.01 * t_oa_c + (2500 + 1.84 * t_oa_c) * d_oa

    def _residual(t_wb: np.ndarray) -> np.ndarray:
        return _wetbulb_enthalpy_kj_per_kg(t_wb[0], p_atm_pa) - e_oa

    return float(fsolve(_residual, x0=[-20.0])[0])


def supply_air_temperature_c(t_sa_lower_c: float, t_sa_upper_c: float, t_oa_c: float) -> float:
    """CRAC supply-air drybulb temperature (C), clipped to the equipment's supply range."""
    if t_oa_c < t_sa_lower_c:
        return t_sa_lower_c
    if t_oa_c > t_sa_upper_c:
        return t_sa_upper_c
    return t_oa_c


def supply_humidity_boundary_g_per_kg(
    t_sa_c: float,
    p_atm_pa: float,
    dp_sa_lower_c: float,
    dp_sa_upper_c: float,
    rh_sa_lower: float,
    rh_sa_upper: float,
) -> tuple[float, float]:
    """Lower/upper bound of supply-air absolute humidity (g/kg) given dewpoint and RH limits."""
    p_sat = 611.21 * np.exp(((18.678 - t_sa_c / 234.5) * t_sa_c) / (t_sa_c + 257.14))

    p_act_l = 611.21 * np.exp((dp_sa_lower_c * 18.678) / (257.14 + dp_sa_lower_c))
    rh_l = p_act_l / p_sat
    if rh_sa_lower <= rh_l <= rh_sa_upper:
        d_lower = 621.98 * p_act_l / (p_atm_pa - p_act_l)
    elif rh_l > rh_sa_upper:
        d_lower = 621.98 * (p_sat * rh_sa_upper) / (p_atm_pa - p_sat * rh_sa_upper)
    else:
        d_lower = 621.98 * (p_sat * rh_sa_lower) / (p_atm_pa - p_sat * rh_sa_lower)

    p_act_u = 611.21 * np.exp((dp_sa_upper_c * 18.678) / (257.14 + dp_sa_upper_c))
    rh_u = p_act_u / p_sat
    if rh_sa_lower <= rh_u <= rh_sa_upper:
        d_upper = 621.98 * p_act_u / (p_atm_pa - p_act_u)
    elif rh_u > rh_sa_upper:
        d_upper = 621.98 * (p_sat * rh_sa_upper) / (p_atm_pa - p_sat * rh_sa_upper)
    else:
        d_upper = 621.98 * (p_sat * rh_sa_lower) / (p_atm_pa - p_sat * rh_sa_lower)

    return d_lower, d_upper


def airside_economizer_state(
    e_oa: float,
    p_atm_pa: float,
    t_sa_c: float,
    d_oa: float,
    t_oa_c: float,
    dp_sa_lower_c: float,
    dp_sa_upper_c: float,
    rh_sa_lower: float,
    rh_sa_upper: float,
) -> tuple[int, int, float, float, float, float, float]:
    """Airside-economizer usage decision and conditioned outdoor-air state.

    Returns ``(ae_use, hd_use, hd_g_per_kg, t_oa_cd_c, e_oa_cd, e_sa, d_sa)``:
    whether the economizer is used, whether adiabatic humidification is used,
    the humidification amount, the conditioned outdoor-air temperature and
    enthalpy, and the supply-air enthalpy/humidity.
    """
    d_threshold_1 = supply_humidity_boundary_g_per_kg(27, p_atm_pa, -9, 15, 0, 15)[1]
    e_threshold_1 = 1.01 * 27 + (2500 + 1.84 * 27) * d_threshold_1 / 1000

    if e_oa < e_threshold_1:
        ae_use = 1
        e_oa_cd = e_oa
        e_sa = e_oa
        if t_oa_c <= t_sa_c:
            p_sat_threshold_11 = 611.21 * np.exp(((18.678 - t_sa_c / 234.5) * t_sa_c) / (t_sa_c + 257.14))
            d_threshold_11 = 621.98 * (p_sat_threshold_11 * 20 / 100) / (p_atm_pa - p_sat_threshold_11 * 20 / 100)
            if d_oa >= d_threshold_11:
                hd_use = 0
                hd_g_per_kg = 0.0
                t_oa_cd = t_oa_c
                d_oa_cd = d_oa
                d_sa = d_oa_cd
            else:
                hd_use = 1
                d_oa_cd = d_threshold_11
                d_sa = d_threshold_11
                hd_g_per_kg = max(d_threshold_11 - d_oa, 0.0)
                t_oa_cd = (e_sa - 2.5 * d_oa_cd) / (1.01 + 1.84 * d_oa_cd / 1000)
        else:
            hd_use = 1
            t_oa_cd = t_sa_c
            d_oa_cd = (e_oa_cd - 1.01 * t_oa_cd) * 1000 / (2500 + 1.84 * t_oa_cd)
            d_sa = d_oa_cd
            hd_g_per_kg = max(d_oa_cd - d_oa, 0.0)
    else:
        p_sat_threshold_2 = 611.21 * np.exp(((18.678 - t_sa_c / 234.5) * t_sa_c) / (t_sa_c + 257.14))
        d_threshold_2 = 621.98 * (p_sat_threshold_2 * 80 / 100) / (p_atm_pa - p_sat_threshold_2 * 80 / 100)
        e_threshold_2 = 1.01 * t_sa_c + (2500 + 1.84 * t_sa_c) * d_threshold_2 / 1000

        if e_oa <= e_threshold_2:
            e_oa_cd = e_oa
            ae_use = 1
            hd_use = 1
            e_sa = e_oa
            t_oa_cd = t_sa_c
            d_oa_cd = (e_oa_cd - 1.01 * t_oa_cd) * 1000 / (2500 + 1.84 * t_oa_cd)
            d_sa = d_oa_cd
            hd_g_per_kg = max(d_oa_cd - d_oa, 0.0)
        else:
            ae_use = 0
            hd_use = 0
            hd_g_per_kg = 0.0
            d_sa_lower, d_sa_upper = supply_humidity_boundary_g_per_kg(
                t_sa_c, p_atm_pa, dp_sa_lower_c, dp_sa_upper_c, rh_sa_lower, rh_sa_upper
            )
            d_sa = (d_sa_lower + d_sa_upper) / 2
            e_sa = 1.01 * t_sa_c + (2500 + 1.84 * t_sa_c) * d_sa / 1000
            # Not used downstream when ae_use == 0 (multiplied through by AE_use in
            # the source model); kept for return-shape parity with the source.
            t_oa_cd = t_sa_c
            e_oa_cd = e_sa

    return ae_use, hd_use, hd_g_per_kg, t_oa_cd, e_oa_cd, e_sa, d_sa


def air_density(p_atm_pa: float, t_air_c: float) -> float:
    """Air density (kg/m^3) from atmospheric pressure and dry-bulb temperature."""
    return 1.293 * (p_atm_pa / 101325) * (273.15 / (273.15 + t_air_c))


def fan_power_kw(mass_flow_rate_kg_s: float, air_density_kg_m3: float, fan_pressure_pa: float, fan_e: float) -> float:
    """Fan shaft power (kW) from volumetric flow, static pressure, and efficiency."""
    return mass_flow_rate_kg_s / air_density_kg_m3 * fan_pressure_pa / fan_e / 1000


def pump_power_kw(
    mass_flow_rate_kg_s: float, pump_pressure_pa: float, pump_e: float, liquid_density_kg_m3: float
) -> float:
    """Pump shaft power (kW) from flow, static pressure, and efficiency."""
    return pump_pressure_pa * mass_flow_rate_kg_s / (1000 * pump_e * liquid_density_kg_m3)


def cooling_requirement_kw(
    power_it_kw: float, ups_e: float, pd_lr: float, l_percentage: float
) -> tuple[float, float, float, float, float, float, float]:
    """Facility heat load (kW) broken out by source.

    Returns ``(powerloss_ups, powerloss_pd, power_lighting, heat_people,
    dc_heat_ac, dc_heat_lc, dc_heat)``. This project's AE-variant use only
    ``dc_heat_ac`` (100% air-cooled IT load, per the source model's AE
    configuration).
    """
    ac_percentage = 1.0
    lc_percentage = 1 - ac_percentage

    powerloss_ups = power_it_kw / ups_e - power_it_kw
    powerloss_pd = power_it_kw / (1 - pd_lr) - power_it_kw
    power_lighting = power_it_kw * l_percentage
    heat_people = 0.0

    dc_heat = power_it_kw + powerloss_ups + powerloss_pd + power_lighting + heat_people
    dc_heat_ac = dc_heat * ac_percentage
    dc_heat_lc = dc_heat * lc_percentage
    return powerloss_ups, powerloss_pd, power_lighting, heat_people, dc_heat_ac, dc_heat_lc, dc_heat


def chiller_cop(load: float, t_c: float) -> float:
    """Chiller coefficient of performance (kW cooling / kW electric), regression fit."""
    return 16.798751 + 0.008840 * t_c**2 - 8.528931 * load**2 + 9.885294 * load - 0.760621 * t_c + 0.084615 * load * t_c


def _pue_ae_chiller_single(t_oa_c: float, rh_oa_pct: float, p_atm_pa: float, power_it_kw: float, config: PUEModelConfig) -> float:
    """Single-facility PUE for the AE (airside economizer + adiabatic cooling) variant."""
    p_sat = 611.21 * np.exp(((18.678 - t_oa_c / 234.5) * t_oa_c) / (t_oa_c + 257.14))
    d_oa = 622 * (rh_oa_pct / 100 * p_sat) / (p_atm_pa - rh_oa_pct / 100 * p_sat)
    t_wb_oa = wetbulb_temperature_c(t_oa_c, rh_oa_pct, p_atm_pa, p_sat)
    density_oa = air_density(p_atm_pa, t_oa_c)
    e_oa = 1.01 * t_oa_c + (2500 + 1.84 * t_oa_c) * d_oa / 1000

    _, _, _, _, dc_heat_ac, _, _ = cooling_requirement_kw(power_it_kw, config.ups_e, config.pd_lr, config.l_percentage)

    t_sa = supply_air_temperature_c(config.t_sa_lower_c, config.t_sa_upper_c, t_oa_c)

    ae_use, _hd_use, hd_g_per_kg, _t_oa_cd, e_oa_cd, e_sa, _d_sa = airside_economizer_state(
        e_oa,
        p_atm_pa,
        t_sa,
        d_oa,
        t_oa_c,
        config.dp_sa_lower_c,
        config.dp_sa_upper_c,
        config.rh_sa_lower,
        config.rh_sa_upper,
    )

    e_ra = e_sa + 1.01 * config.delta_t_air_c
    m_sa = dc_heat_ac / (e_ra - e_sa)
    if ae_use == 1:
        m_ra = max(0.0, m_sa * (e_sa - e_oa_cd) / (e_ra - e_oa_cd))
    else:
        m_ra = m_sa
    m_oa_cd = m_sa - m_ra
    m_oa = m_oa_cd / (1 + hd_g_per_kg / 1000)

    ae_heat_removed = ae_use * m_oa_cd * (e_ra - e_oa_cd)
    density_sa = air_density(p_atm_pa, t_sa)
    power_fan_crac = fan_power_kw(m_sa, density_sa, config.fan_pressure_crac_pa, config.fan_e_crac)

    hd_amount = m_oa_cd - m_oa
    power_pump_hd = pump_power_kw(hd_amount, config.pump_pressure_hd_pa, config.pump_e_hd, 1000)

    cop_chiller = chiller_cop(config.chiller_load, t_wb_oa + config.at_ct_c) * (1 + config.p_cop)

    chiller_heat_removed = dc_heat_ac - ae_heat_removed
    power_chiller = chiller_heat_removed / cop_chiller

    m_sw = chiller_heat_removed / (4.2 * config.delta_t_water_c)
    power_pump_cw = pump_power_kw(m_sw, config.pump_pressure_cw_pa, config.pump_e_cw, 1000)

    ct_heat_removed = chiller_heat_removed + power_chiller
    m_ct = ct_heat_removed / (4.2 * config.delta_t_ct_c)
    power_pump_ct = pump_power_kw(m_ct, config.pump_pressure_ct_pa, config.pump_e_ct, 1000)
    power_fan_ct = fan_power_kw(m_ct / config.l_to_g, density_oa, config.fan_pressure_ct_pa, config.fan_e_ct)

    powerloss_ups, powerloss_pd, power_lighting, _, _, _, _ = cooling_requirement_kw(
        power_it_kw, config.ups_e, config.pd_lr, config.l_percentage
    )

    total_power_kw = (
        power_it_kw
        + powerloss_ups
        + powerloss_pd
        + power_lighting
        + power_fan_crac
        + power_pump_hd
        + power_chiller
        + power_pump_cw
        + power_pump_ct
        + power_fan_ct
    )
    return total_power_kw / power_it_kw


def pue_ae_chiller(
    t_oa_c: np.ndarray,
    rh_oa_pct: np.ndarray,
    p_atm_pa: np.ndarray | float,
    power_it_kw: np.ndarray | float,
    config: PUEModelConfig = PUEModelConfig(),
) -> np.ndarray:
    """PUE for the AE (airside economizer + adiabatic cooling) chiller variant.

    Vectorized over facilities via an explicit loop around the per-facility
    scalar computation (the underlying model solves an implicit wetbulb
    equation per point via ``scipy.optimize.fsolve`` and branches on
    economizer usage, so it isn't a natural fit for elementwise numpy ops).

    Parameters
    ----------
    t_oa_c, rh_oa_pct:
        Outdoor-air dry-bulb temperature (C) and relative humidity (%),
        array-like, one value per facility/timestep.
    p_atm_pa, power_it_kw:
        Scalar or array-like (broadcast to ``t_oa_c``'s shape); see
        ``PUEModelConfig`` for their defaults and the elevation-adjustment
        extension point.
    config:
        Equipment parameters; defaults to the Table A.1 midpoints.

    Returns
    -------
    numpy.ndarray
        PUE (dimensionless), same shape as ``t_oa_c``.
    """
    t_oa_c = np.atleast_1d(np.asarray(t_oa_c, dtype=float))
    rh_oa_pct = np.atleast_1d(np.asarray(rh_oa_pct, dtype=float))
    p_atm_pa = np.broadcast_to(np.asarray(p_atm_pa, dtype=float), t_oa_c.shape)
    power_it_kw = np.broadcast_to(np.asarray(power_it_kw, dtype=float), t_oa_c.shape)

    pue = np.empty_like(t_oa_c)
    for i in range(t_oa_c.size):
        pue[i] = _pue_ae_chiller_single(
            float(t_oa_c[i]), float(rh_oa_pct[i]), float(p_atm_pa[i]), float(power_it_kw[i]), config
        )
    return pue


def check_pue_plausible(pue: np.ndarray, lower: float = 1.0, upper: float = 3.0) -> None:
    """Sanity-check PUE values fall within a plausible range, mirroring ``_log_unit_sanity``.

    PUE is bounded below by 1.0 by definition (total facility power can't be
    less than IT power); the paper's own AE contour plot (Fig. 6) puts
    typical hyperscale AE+adiabatic PUE around 1.0-1.2, so an upper bound of
    3.0 is deliberately generous -- wide enough not to flag legitimate
    extreme-heat days, but still low enough to catch a broken unit
    conversion or a diverged ``fsolve`` call.
    """
    finite = pue[np.isfinite(pue)]
    if finite.size == 0:
        raise ValueError("No finite PUE values found for plausibility check.")
    if finite.min() < lower or finite.max() > upper:
        raise ValueError(
            f"PUE values outside plausible range [{lower}, {upper}]: "
            f"min={finite.min():.4f}, max={finite.max():.4f}."
        )
