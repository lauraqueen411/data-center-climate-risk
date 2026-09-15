"""Vectorized-numpy psychrometrics replacement for ``cooling_archetypes.py``.

Same physics, same equipment configs (imported, not duplicated), same
function signatures as ``cooling_archetypes.py`` -- the only difference is
that every ``CoolProp.HumidAirProp.HAPropsSI`` call is replaced with a
standard ASHRAE moist-air formula evaluated as a plain numpy array
expression. Exists because ``HAPropsSI``'s array support is a per-element
internal loop, not true vectorization (measured ~15.6us/point, the same
order as calling it in an explicit Python loop) -- at Oregon's full ~9,100
grid cells this makes the CoolProp path cost ~5 days/GCM single-threaded;
these formulas measured ~50x faster for the wet-bulb solve alone (3.2M
pts/s vs. CoolProp's 64K pts/s), all pure numpy elementwise arithmetic.

Assumptions log
----------------
- **Validated against CoolProp itself, not re-derived from a textbook and
  trusted on faith.** Each primitive below was checked against
  ``HAPropsSI`` on 200 random (T, RH, P) points spanning -15 to 45 degC,
  5-100% RH, 85-101 kPa. Typical (median) agreement: enthalpy within
  ~140 J/kg, dewpoint within ~0.21 degC, specific volume within ~0.5%.
  Worst-case agreement is concentrated at extreme conditions rather than
  spread evenly: enthalpy's largest absolute error (~663 J/kg, ~0.35%
  relative) occurred at a hot/humid/high-altitude point (43 degC, 81% RH,
  85 kPa); dewpoint's largest error (~3.15 degC) occurred at a cold, very
  dry point (-3 degC, 6% RH) -- a known limitation of the Arden Buck
  dewpoint inversion (same formula ``pue.py`` already uses for the 2020
  model) far from typical conditions, not a bug introduced here. These are
  engineering approximations, not exact -- treat this module's output as a
  fast path to be cross-checked against the CoolProp-based
  ``cooling_archetypes.py`` on real conditions (see
  ``scripts/run_pue_wue_facilities.py``'s CoolProp output) before trusting
  it at full-grid scale, not as automatically equivalent.
- **Specific volume/density** uses the ideal-gas humid-air formula
  (``R_da * T_K * (1 + 1.6078*W) / P``) rather than CoolProp's real-gas
  equation of state -- median ~0.5% error, worst-case ~5.6% at the same
  hot/humid/high-altitude extreme as enthalpy. Affects fan/pump power terms
  (which scale with 1/density), not the chiller COP or the core
  economizer/humidity decision logic.
- **Wet-bulb temperature** is solved via 8 iterations of vectorized Newton's
  method (matching outside-air enthalpy to saturated enthalpy at a
  candidate wet-bulb temperature) rather than ``scipy.optimize.fsolve``
  per-point (as ``pue.py``'s own ``wetbulb_temperature_c`` already does for
  the 2020 model) -- 8 iterations converges to sub-0.001 degC residual
  well within CoolProp's own ~0.5 degC agreement envelope, so is not the
  limiting source of error here.
- Everything else (equipment configs, branch logic, GP chiller COP,
  elevation-adjusted pressure) is identical to ``cooling_archetypes.py`` --
  imported, not reimplemented, so a fix or config change there does not
  need to be duplicated here.
"""

from __future__ import annotations

import numpy as np

from climate_risk_dc.climate.cooling_archetypes import (  # noqa: F401 (re-exported for callers)
    Case1Config,
    Case2Config,
    Case5Config,
    check_pue_wue_plausible,
    chiller_cop_water_cooled,
    elevation_adjusted_pressure_pa,
    load_cop_water_cooled_gp,
)

_R_DA = 287.055  # J/(kg K), specific gas constant for dry air
_FAN_POWER = lambda m, rho, dp, e: m / rho * dp / e / 1000.0  # noqa: E731
_PUMP_POWER = lambda m, dp, e, rho: dp * m / (1000.0 * e * rho)  # noqa: E731


# --------------------------------------------------------------------------
# Vectorized numpy psychrometric primitives (validated against CoolProp --
# see this module's Assumptions log)
# --------------------------------------------------------------------------


def _p_sat_pa(t_c: np.ndarray) -> np.ndarray:
    """Saturation vapor pressure (Pa), Arden Buck equation -- same formula pue.py uses."""
    return 611.2 * np.exp(((18.678 - t_c / 234.5) * t_c) / (t_c + 257.14))


def _w_from_rh(t_c: np.ndarray, rh_frac: np.ndarray | float, p_pa: np.ndarray) -> np.ndarray:
    p_v = rh_frac * _p_sat_pa(t_c)
    return 0.62198 * p_v / (p_pa - p_v)


def _w_from_tdp(tdp_c: np.ndarray, p_pa: np.ndarray) -> np.ndarray:
    p_v = _p_sat_pa(tdp_c)
    return 0.62198 * p_v / (p_pa - p_v)


def _tdp_from_rh(t_c: np.ndarray, rh_frac: np.ndarray) -> np.ndarray:
    log_term = np.log(rh_frac * np.exp((18.678 - t_c / 234.5) * (t_c / (257.14 + t_c))))
    return 257.14 * log_term / (18.678 - log_term)


def _h_from_t_w(t_c: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Moist-air specific enthalpy (J/kg dry air), ASHRAE formula."""
    return 1006.0 * t_c + w * (2_501_000.0 + 1860.0 * t_c)


def _t_from_h_w(h: np.ndarray, w: np.ndarray) -> np.ndarray:
    return (h - w * 2_501_000.0) / (1006.0 + 1860.0 * w)


def _vha_from_t_w(t_c: np.ndarray, w: np.ndarray, p_pa: np.ndarray) -> np.ndarray:
    """Specific volume of moist air per kg dry air (m^3/kg), ideal-gas based."""
    return _R_DA * (t_c + 273.15) * (1.0 + 1.6078 * w) / p_pa


def _wetbulb_from_rh(t_c: np.ndarray, rh_frac: np.ndarray, p_pa: np.ndarray, n_iter: int = 8) -> np.ndarray:
    """Vectorized Newton iteration for wet-bulb temperature (degC)."""
    h_oa = _h_from_t_w(t_c, _w_from_rh(t_c, rh_frac, p_pa))
    twb = np.asarray(t_c, dtype=float) * 0.7  # initial guess, below dry-bulb
    for _ in range(n_iter):
        h_sat = _h_from_t_w(twb, _w_from_rh(twb, 1.0, p_pa))
        eps = 1e-4
        h_sat2 = _h_from_t_w(twb + eps, _w_from_rh(twb + eps, 1.0, p_pa))
        df = (h_sat2 - h_sat) / eps
        twb = twb - (h_sat - h_oa) / df
    return twb


# --------------------------------------------------------------------------
# Case 1 economizer logic (vectorized), numpy-primitive version of
# cooling_archetypes._air_side_economizer
# --------------------------------------------------------------------------


def _air_side_economizer_fast(
    t_up_c: float,
    t_lw_c: float,
    dp_up_c: float,
    dp_lw_c: float,
    rh_up_pct: float,
    rh_lw_pct: float,
    t_oa_c: np.ndarray,
    rh_oa_pct: np.ndarray,
    p_oa_pa: np.ndarray,
    delta_t_air_c: float,
) -> dict[str, np.ndarray]:
    t_oa_c = np.asarray(t_oa_c, dtype=float)
    rh_oa_pct = np.asarray(rh_oa_pct, dtype=float)
    p_oa_pa = np.broadcast_to(np.asarray(p_oa_pa, dtype=float), t_oa_c.shape).copy()
    rh_oa_frac = rh_oa_pct / 100.0
    rh_up_frac, rh_lw_frac = rh_up_pct / 100.0, rh_lw_pct / 100.0

    h_oa = _h_from_t_w(t_oa_c, _w_from_rh(t_oa_c, rh_oa_frac, p_oa_pa))
    d_oa = _w_from_rh(t_oa_c, rh_oa_frac, p_oa_pa)

    mask_a = (t_oa_c >= t_lw_c) & (t_oa_c < t_up_c)
    mask_b = t_oa_c < t_lw_c
    mask_c = t_oa_c >= t_up_c

    # ---- Branch A: T_lw <= T_oa < T_up, thresholds at T_oa ----
    d_up_a = np.minimum(_w_from_rh(t_oa_c, rh_up_frac, p_oa_pa), _w_from_tdp(np.full_like(t_oa_c, dp_up_c), p_oa_pa))
    d_lw_a = np.maximum(_w_from_rh(t_oa_c, rh_lw_frac, p_oa_pa), _w_from_tdp(np.full_like(t_oa_c, dp_lw_c), p_oa_pa))
    cond_a1 = (d_oa >= d_lw_a) & (d_oa <= d_up_a)
    cond_a2 = d_oa > d_up_a

    h_cd_a2 = _h_from_t_w(np.full_like(t_oa_c, dp_up_c), d_up_a)
    t_sa_a2 = np.full_like(t_oa_c, max(t_lw_c, dp_up_c))
    h_sa_a2 = _h_from_t_w(t_sa_a2, d_up_a)

    t_cd_a3 = _t_from_h_w(h_oa, d_lw_a)
    t_sa_a3 = np.maximum(t_lw_c, t_cd_a3)
    h_sa_a3 = _h_from_t_w(t_sa_a3, d_lw_a)

    ae_use_a = np.where(cond_a1, 1, np.where(cond_a2, 0, 1))
    hd_use_a = np.where(cond_a1, 0, np.where(cond_a2, 0, 1))
    dhd_use_a = np.where(cond_a1, 0, np.where(cond_a2, 1, 0))
    delta_d_a = np.where(cond_a1, 0.0, np.where(cond_a2, d_oa - d_up_a, d_lw_a - d_oa))
    d_sa_a = np.where(cond_a1, d_oa, np.where(cond_a2, d_up_a, d_lw_a))
    t_cd_a = np.where(cond_a1, t_oa_c, np.where(cond_a2, dp_up_c, t_cd_a3))
    h_cd_a = np.where(cond_a1, h_oa, np.where(cond_a2, h_cd_a2, h_oa))
    t_sa_a = np.where(cond_a1, t_oa_c, np.where(cond_a2, t_sa_a2, t_sa_a3))
    h_sa_a = np.where(cond_a1, h_oa, np.where(cond_a2, h_sa_a2, h_sa_a3))

    # ---- Branch B: T_oa < T_lw, thresholds at T_lw ----
    t_lw_arr = np.full_like(t_oa_c, t_lw_c)
    d_up_b = np.minimum(_w_from_rh(t_lw_arr, rh_up_frac, p_oa_pa), _w_from_tdp(np.full_like(t_oa_c, dp_up_c), p_oa_pa))
    d_lw_b = np.maximum(_w_from_rh(t_lw_arr, rh_lw_frac, p_oa_pa), _w_from_tdp(np.full_like(t_oa_c, dp_lw_c), p_oa_pa))
    cond_b1 = d_oa < d_lw_b
    cond_b2 = (d_oa >= d_lw_b) & (d_oa <= d_up_b)
    cond_b3 = d_oa > d_up_b

    t_cd_b1 = _t_from_h_w(h_oa, d_lw_b)
    t_sa_b1 = np.maximum(t_lw_c, t_oa_c)
    h_sa_b1 = _h_from_t_w(t_sa_b1, d_lw_b)

    t_sa_b2 = np.maximum(t_lw_c, t_oa_c)
    h_sa_b2 = _h_from_t_w(t_sa_b2, d_oa)

    d_up_at_upsasp = _w_from_tdp(np.full_like(t_oa_c, dp_up_c), p_oa_pa)  # T_up + Tdp(dp_up) -- Tdp-derived W is T-independent
    cond_b3a = cond_b3 & (d_oa <= d_up_at_upsasp)
    cond_b3b = cond_b3 & (d_oa > d_up_at_upsasp)
    h_sa_b3a = _h_from_t_w(np.full_like(t_oa_c, t_lw_c), d_oa)
    t_sa_b3b_val = max(t_lw_c, dp_up_c)
    h_cd_b3b = _h_from_t_w(np.full_like(t_oa_c, dp_up_c), d_up_at_upsasp)
    h_sa_b3b = _h_from_t_w(np.full_like(t_oa_c, t_sa_b3b_val), d_up_at_upsasp)

    hd_use_b = np.where(cond_b1, 1, 0)
    dhd_use_b = np.where(cond_b3b, 1, 0)
    delta_d_b = np.where(cond_b1, d_lw_b - d_oa, np.where(cond_b3b, d_oa - d_up_at_upsasp, 0.0))
    d_sa_b = np.where(cond_b1, d_lw_b, np.where(cond_b2, d_oa, np.where(cond_b3a, d_oa, d_up_at_upsasp)))
    h_cd_b = np.where(cond_b1, h_oa, np.where(cond_b2, h_oa, np.where(cond_b3a, h_oa, h_cd_b3b)))
    t_sa_b = np.where(cond_b1, t_sa_b1, np.where(cond_b2, t_sa_b2, np.where(cond_b3a, t_lw_c, t_sa_b3b_val)))
    h_sa_b = np.where(cond_b1, h_sa_b1, np.where(cond_b2, h_sa_b2, np.where(cond_b3a, h_sa_b3a, h_sa_b3b)))

    # ---- Branch C: T_oa >= T_up (unconditional simplified control policy) ----
    t_up_arr = np.full_like(t_oa_c, t_up_c)
    d_up_c_ = np.minimum(_w_from_rh(t_up_arr, rh_up_frac, p_oa_pa), _w_from_tdp(np.full_like(t_oa_c, dp_up_c), p_oa_pa))
    h_cd_c = _h_from_t_w(t_up_arr, d_up_c_)

    ae_use_c = np.zeros_like(t_oa_c)
    hd_use_c = np.zeros_like(t_oa_c)
    dhd_use_c = np.zeros_like(t_oa_c)
    delta_d_c = np.zeros_like(t_oa_c)
    d_sa_c = d_up_c_
    t_sa_c = t_up_arr
    h_sa_c = h_cd_c

    def _select(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
        return np.select([mask_a, mask_b, mask_c], [a, b, c])

    return {
        "ae_use": _select(ae_use_a, np.ones_like(t_oa_c), ae_use_c),
        "hd_use": _select(hd_use_a, hd_use_b, hd_use_c),
        "d_sa": _select(d_sa_a, d_sa_b, d_sa_c),
        "delta_d": _select(delta_d_a, delta_d_b, delta_d_c),
        "t_sa": _select(t_sa_a, t_sa_b, t_sa_c),
        "h_sa": _select(h_sa_a, h_sa_b, h_sa_c),
        "h_cd": _select(h_cd_a, h_cd_b, h_cd_c),
    }


def _chiller_system_fast(
    t_up_c: float, t_lw_c: float, dp_up_c: float, dp_lw_c: float, rh_up_pct: float, rh_lw_pct: float, p_oa_pa: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Numpy-primitive version of ``cooling_archetypes._chiller_system`` (same midpoint deviation)."""
    t_sa = np.full_like(np.broadcast_to(p_oa_pa, np.shape(p_oa_pa)), max(t_up_c, t_lw_c), dtype=float)
    d_up = np.minimum(_w_from_rh(t_sa, rh_up_pct / 100.0, p_oa_pa), _w_from_tdp(np.full_like(t_sa, dp_up_c), p_oa_pa))
    d_lw = np.maximum(_w_from_rh(t_sa, rh_lw_pct / 100.0, p_oa_pa), _w_from_tdp(np.full_like(t_sa, dp_lw_c), p_oa_pa))
    d_sa = (d_up + d_lw) / 2.0
    h_sa = _h_from_t_w(t_sa, d_sa)
    return t_sa, d_sa, h_sa


def _cooling_tower_fast(
    t_oa_c: np.ndarray,
    rh_oa_pct: np.ndarray,
    p_oa_pa: np.ndarray,
    at_ct_c: float,
    q_kw: np.ndarray,
    delta_t_ct_c: float,
    windage_p: float,
    cc: float,
    lg_ratio: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    twb_oa = _wetbulb_from_rh(t_oa_c, rh_oa_pct / 100.0, p_oa_pa)
    latent_heat_vap = -0.0013 * (twb_oa + at_ct_c + delta_t_ct_c / 2) ** 2 - 2.3097 * (
        twb_oa + at_ct_c + delta_t_ct_c / 2
    ) + 2500.5
    m_ct = q_kw / (4.184 * delta_t_ct_c)
    water_evaporated = q_kw / latent_heat_vap
    m_air = m_ct / lg_ratio
    water_windage = m_ct * windage_p
    water_drainoff = np.maximum(water_evaporated / (cc - 1) - water_windage, 0.0)
    return m_air, water_evaporated, water_windage, water_drainoff


def _waterside_economizer_fast(
    t_sfw_c: np.ndarray, t_rfw_c: np.ndarray, twb_oa_c: np.ndarray, at_ct_c: float, at_he_c: float, cooling_required_kw: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    approach = twb_oa_c + at_ct_c + at_he_c
    full_use = approach <= t_sfw_c
    partial_use = (approach > t_sfw_c) & (approach < t_rfw_c)
    return (full_use | partial_use).astype(int), np.where(
        full_use, cooling_required_kw, np.where(partial_use, cooling_required_kw * ((t_rfw_c - approach) / (t_rfw_c - t_sfw_c)), 0.0)
    )


def pue_wue_case1(
    t_oa_c: np.ndarray, rh_oa_pct: np.ndarray, p_atm_pa: np.ndarray | float, cop_gp, config: Case1Config = Case1Config()
) -> tuple[np.ndarray, np.ndarray]:
    power_it = 1.0
    t_oa_c = np.asarray(t_oa_c, dtype=float)
    rh_oa_pct = np.asarray(rh_oa_pct, dtype=float)
    p_atm_pa = np.broadcast_to(np.asarray(p_atm_pa, dtype=float), t_oa_c.shape).copy()

    q = power_it + (power_it / config.ups_e - power_it) + (power_it / (1 - config.pd_lr) - power_it) + power_it * config.l_percentage

    ae = _air_side_economizer_fast(
        config.t_up_c, config.t_lw_c, config.dp_up_c, config.dp_lw_c, config.rh_up_pct, config.rh_lw_pct,
        t_oa_c, rh_oa_pct, p_atm_pa, config.delta_t_air_c,
    )
    ae_use, hd_use, delta_d, t_sa, h_sa, h_cd, d_sa = (
        ae["ae_use"], ae["hd_use"], ae["delta_d"], ae["t_sa"], ae["h_sa"], ae["h_cd"], ae["d_sa"]
    )

    h_ra = h_sa + 1010.0 * config.delta_t_air_c  # ASHRAE H in J/kg (vs. CoolProp's own kJ/kg scaling in the source -- see note below
    # NOTE: cooling_archetypes.py works in kJ/kg (H/1000) throughout because the source repo's H is in J/kg
    # and the repo itself divides by 1000. Here H is already J/kg from _h_from_t_w, so convert once to kJ/kg
    # to reuse the exact same downstream arithmetic as cooling_archetypes.py.
    h_sa_kj, h_cd_kj, h_ra_kj = h_sa / 1000.0, h_cd / 1000.0, (h_sa + 1010.0 * config.delta_t_air_c) / 1000.0

    m_sa = q / (h_ra_kj - h_sa_kj)
    m_cd = m_sa * (h_sa_kj - h_ra_kj) / (h_cd_kj - h_ra_kj)
    m_cd_dry = m_cd / (1 + d_sa)
    heat = ae_use * m_cd * (h_ra_kj - h_cd_kj)
    density_sa = 1.0 / _vha_from_t_w(t_sa, d_sa, p_atm_pa)
    power_fan_crac = _FAN_POWER(m_sa, density_sa, config.fan_pressure_crac_pa, config.fan_e_crac)
    hd_amount_ae = np.maximum(hd_use * m_cd_dry * delta_d, 0.0)

    chiller_heat_removed_pre = q - heat
    q_heat_latent = np.where(
        chiller_heat_removed_pre != 0, chiller_heat_removed_pre / config.shr - chiller_heat_removed_pre, 0.0
    )
    chiller_heat_removed = chiller_heat_removed_pre + q_heat_latent
    hd_amount = np.maximum(q_heat_latent / 2266.0, 0.0)
    power_pump_hd = _PUMP_POWER(hd_amount, config.pump_pressure_hd_pa, config.pump_e_hd, 1000.0)

    twb_oa = _wetbulb_from_rh(t_oa_c, rh_oa_pct / 100.0, p_atm_pa)
    cop_chiller = chiller_cop_water_cooled(twb_oa, config.at_ct_c, config.chiller_load, config.pcop, cop_gp)
    power_chiller = chiller_heat_removed / cop_chiller

    m_sw = chiller_heat_removed / (4.184 * config.delta_t_water_c)
    power_pump_cw = _PUMP_POWER(m_sw, config.pump_pressure_cw_pa, config.pump_e_cw, 1000.0)

    ct_heat_removed = chiller_heat_removed + power_chiller
    m_ct = ct_heat_removed / (4.184 * config.delta_t_ct_c)
    power_pump_ct = _PUMP_POWER(m_ct, config.pump_pressure_ct_pa, config.pump_e_ct, 1000.0)

    m_ct_air, water_evap, water_windage, water_drainoff = _cooling_tower_fast(
        t_oa_c, rh_oa_pct, p_atm_pa, config.at_ct_c, ct_heat_removed, config.delta_t_ct_c, config.windage_p, config.cc, config.lg_ratio
    )
    density_oa = 1.0 / _vha_from_t_w(t_oa_c, _w_from_rh(t_oa_c, rh_oa_pct / 100.0, p_atm_pa), p_atm_pa)
    power_fan_ct = _FAN_POWER(m_ct_air, density_oa, config.fan_pressure_ct_pa, config.fan_e_ct)

    power_comp = (
        power_it + (power_it / config.ups_e - power_it) + (power_it / (1 - config.pd_lr) - power_it) + power_it * config.l_percentage
        + power_fan_crac + power_pump_hd + power_chiller + power_pump_cw + power_pump_ct + power_fan_ct
    )
    pue = power_comp / power_it
    water_comp = (
        (hd_amount + hd_amount_ae) + np.maximum(water_evap, 0.0) + np.maximum(water_windage, 0.0) + np.maximum(water_drainoff, 0.0)
    )
    wue = water_comp * 3600.0 / power_it
    return pue, wue


def pue_wue_case2(
    t_oa_c: np.ndarray, rh_oa_pct: np.ndarray, p_atm_pa: np.ndarray | float, cop_gp, config: Case2Config = Case2Config()
) -> tuple[np.ndarray, np.ndarray]:
    power_it = 1.0
    t_oa_c = np.asarray(t_oa_c, dtype=float)
    rh_oa_pct = np.asarray(rh_oa_pct, dtype=float)
    p_atm_pa = np.broadcast_to(np.asarray(p_atm_pa, dtype=float), t_oa_c.shape).copy()

    twb_oa = _wetbulb_from_rh(t_oa_c, rh_oa_pct / 100.0, p_atm_pa)
    q = power_it + (power_it / config.ups_e - power_it) + (power_it / (1 - config.pd_lr) - power_it) + power_it * config.l_percentage
    q_heat_latent = q / config.shr - q

    t_sa, d_sa, h_sa = _chiller_system_fast(
        config.t_up_c, config.t_lw_c, config.dp_up_c, config.dp_lw_c, config.rh_up_pct, config.rh_lw_pct, p_atm_pa
    )
    h_sa_kj = h_sa / 1000.0
    t_ra = t_sa + config.delta_t_air_c
    m_sa = np.full_like(t_oa_c, q / (1.01 * config.delta_t_air_c))
    density_sa = 1.0 / _vha_from_t_w(t_sa, d_sa, p_atm_pa)
    power_fan_crac = _FAN_POWER(m_sa, density_sa, config.fan_pressure_crac_pa, config.fan_e_crac)

    hd_amount = np.full_like(t_oa_c, max(q_heat_latent / 2266.0, 0.0))
    power_pump_hd = _PUMP_POWER(hd_amount, config.pump_pressure_hd_pa, config.pump_e_hd, 1000.0)

    cooling_required = q + q_heat_latent
    t_sfw = t_ra - (t_ra - t_sa) / config.hte
    t_rfw = t_sfw + config.delta_t_water_c
    _we_use, we_heat_removed = _waterside_economizer_fast(
        t_sfw, t_rfw, twb_oa, config.at_ct_c, config.at_he_c, np.full_like(t_oa_c, cooling_required)
    )
    m_sfw = we_heat_removed / (4.2 * config.delta_t_water_c)
    power_pump_we = _PUMP_POWER(m_sfw, config.pump_pressure_we_pa, config.pump_e_we, 1000.0)

    chiller_heat_removed = cooling_required - we_heat_removed
    m_sw = chiller_heat_removed / (4.2 * config.delta_t_water_c)
    power_pump_cw = _PUMP_POWER(m_sw, config.pump_pressure_cw_pa, config.pump_e_cw, 1000.0)

    cop_chiller = chiller_cop_water_cooled(twb_oa, config.at_ct_c, config.chiller_load, config.pcop, cop_gp)
    power_chiller = chiller_heat_removed / cop_chiller

    m_ct = np.full_like(t_oa_c, q / (4.184 * config.delta_t_ct_c))
    power_pump_ct = _PUMP_POWER(m_ct, config.pump_pressure_ct_pa, config.pump_e_ct, 1000.0)

    m_ct_air, water_evap, water_windage, water_drainoff = _cooling_tower_fast(
        t_oa_c, rh_oa_pct, p_atm_pa, config.at_ct_c, np.full_like(t_oa_c, q), config.delta_t_ct_c, config.windage_p, config.cc, config.lg_ratio
    )
    density_oa = 1.0 / _vha_from_t_w(t_oa_c, _w_from_rh(t_oa_c, rh_oa_pct / 100.0, p_atm_pa), p_atm_pa)
    power_fan_ct = _FAN_POWER(m_ct_air, density_oa, config.fan_pressure_ct_pa, config.fan_e_ct)

    power_comp = (
        power_it + (power_it / config.ups_e - power_it) + (power_it / (1 - config.pd_lr) - power_it) + power_it * config.l_percentage
        + power_fan_crac + power_pump_hd + power_pump_we + power_chiller + power_pump_cw + power_pump_ct + power_fan_ct
    )
    pue = power_comp / power_it
    water_comp = hd_amount + np.maximum(water_evap, 0.0) + np.maximum(water_windage, 0.0) + np.maximum(water_drainoff, 0.0)
    wue = water_comp * 3600.0 / power_it
    return pue, wue


def pue_wue_case5(
    t_oa_c: np.ndarray, rh_oa_pct: np.ndarray, p_atm_pa: np.ndarray | float, cop_gp, config: Case5Config = Case5Config()
) -> tuple[np.ndarray, np.ndarray]:
    power_it = 1.0
    t_oa_c = np.asarray(t_oa_c, dtype=float)
    rh_oa_pct = np.asarray(rh_oa_pct, dtype=float)
    p_atm_pa = np.broadcast_to(np.asarray(p_atm_pa, dtype=float), t_oa_c.shape).copy()

    q = power_it + (power_it / config.ups_e - power_it) + (power_it / (1 - config.pd_lr) - power_it) + power_it * config.l_percentage
    q_heat_latent = q / config.shr - q

    t_sa, d_sa, h_sa = _chiller_system_fast(
        config.t_up_c, config.t_lw_c, config.dp_up_c, config.dp_lw_c, config.rh_up_pct, config.rh_lw_pct, p_atm_pa
    )
    m_sa = np.full_like(t_oa_c, q / (1.01 * config.delta_t_air_c))
    density_sa = 1.0 / _vha_from_t_w(t_sa, d_sa, p_atm_pa)
    power_fan_crac = _FAN_POWER(m_sa, density_sa, config.fan_pressure_crac_pa, config.fan_e_crac)

    hd_amount = np.full_like(t_oa_c, max(q_heat_latent / 2266.0, 0.0))
    power_pump_hd = _PUMP_POWER(hd_amount, config.pump_pressure_hd_pa, config.pump_e_hd, 1000.0)

    chiller_heat_removed = np.full_like(t_oa_c, q_heat_latent + q)
    m_sw = chiller_heat_removed / (4.2 * config.delta_t_water_c)
    power_pump_cw = _PUMP_POWER(m_sw, config.pump_pressure_cw_pa, config.pump_e_cw, 1000.0)

    twb_oa = _wetbulb_from_rh(t_oa_c, rh_oa_pct / 100.0, p_atm_pa)
    cop_chiller = chiller_cop_water_cooled(twb_oa, config.at_ct_c, config.chiller_load, config.pcop, cop_gp)
    power_chiller = chiller_heat_removed / cop_chiller

    ct_heat_removed = chiller_heat_removed + power_chiller
    m_ct = ct_heat_removed / (4.184 * config.delta_t_ct_c)
    power_pump_ct = _PUMP_POWER(m_ct, config.pump_pressure_ct_pa, config.pump_e_ct, 1000.0)

    m_ct_air, water_evap, water_windage, water_drainoff = _cooling_tower_fast(
        t_oa_c, rh_oa_pct, p_atm_pa, config.at_ct_c, ct_heat_removed, config.delta_t_ct_c, config.windage_p, config.cc, config.lg_ratio
    )
    density_oa = 1.0 / _vha_from_t_w(t_oa_c, _w_from_rh(t_oa_c, rh_oa_pct / 100.0, p_atm_pa), p_atm_pa)
    power_fan_ct = _FAN_POWER(m_ct_air, density_oa, config.fan_pressure_ct_pa, config.fan_e_ct)

    power_comp = (
        power_it + (power_it / config.ups_e - power_it) + (power_it / (1 - config.pd_lr) - power_it) + power_it * config.l_percentage
        + power_fan_crac + power_pump_hd + power_chiller + power_pump_cw + power_pump_ct + power_fan_ct
    )
    pue = power_comp / power_it
    water_comp = hd_amount + np.maximum(water_evap, 0.0) + np.maximum(water_windage, 0.0) + np.maximum(water_drainoff, 0.0)
    wue = water_comp * 3600.0 / power_it
    return pue, wue
