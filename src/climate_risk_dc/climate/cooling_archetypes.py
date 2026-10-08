"""Gridded PUE/WUE cooling archetypes for the Oregon ensemble extension.

Ported from Nuoa Lei's reference implementation for Lei & Masanet (2022),
*Resources, Conservation & Recycling* 182, 106323,
https://doi.org/10.1016/j.resconrec.2022.106323
(https://github.com/nuoaleon/Data-Center-Water-footprint,
``simulation_funs_DC.py``). This is a *different* source repo/paper than
``pue.py`` (which ports Lei & Masanet 2020's ``Data-center-PUE-prediction-tool``
for the facility-level, single-GCM Section 1.2 analysis) -- see
``docs/tasks/ensemble_computation.md``'s Phase 0 report. The two are not
numerically comparable: this module uses CoolProp for psychrometrics and a
``GaussianProcessRegressor`` (fit on external chiller performance data) for
chiller COP, versus ``pue.py``'s hand-derived psychrometric formulas and a
polynomial COP regression, and this module additionally computes WUE, which
``pue.py`` does not.

Archetypes implemented (of the paper's 10 Cases, all sharing one GP model):

- **Case 1** (``pue_wue_case1``) -- large-scale, airside economizer +
  adiabatic cooling + water-cooled chiller. The paper's own validation
  facilities (Prineville, The Dalles) run this configuration.
- **Case 2** (``pue_wue_case2``) -- large-scale, waterside economizer +
  water-cooled chiller.
- **Case 5** (``pue_wue_case5``) -- midsize, water-cooled chiller, no
  economizer (the no-free-cooling counterfactual).

Cases using ``COP_DX.pkl``/``COP_AC.pkl`` (direct expansion, air-cooled
chiller) are deliberately not implemented here: neither has a cooling tower,
so their WUE is structurally just a small sensible-heat-ratio-driven
humidification term with no meaningful climate signal (confirmed against the
paper's own Sobol total-effect indices, Fig. 7), and the Oregon IM3/PNNL
facility inventory (109 facilities) has zero entries in the paper's "small
DC" size class those cases target. ``COP_AC.pkl`` is also currently broken
under every available scikit-learn version this project tested (``.predict()``
raises ``AttributeError: 'GaussianProcessRegressor' object has no attribute
'_y_train_std'``), which independently rules it out.

Assumptions log
----------------
- **Case 1 indexing bug, corrected.** The published
  ``PUE_WUE_AE_Chiller`` calls ``Air_side_economizer(...)[i]`` with indices
  that do not match that function's own documented return order
  (``AE_use, HD_use, DHD_use, Steam_use, delta_d, T_sa, H_sa, T_cd, H_cd,
  d_sa, H_oa, chiller_energy, heating_energy``) -- e.g. it reads index 2
  (``DHD_use``) expecting ``T_sa``, and index 3 (``Steam_use``, always 0)
  expecting ``H_sa``. Running the published code as-is (confirmed on this
  project's own hardware) produces ``WUE = 0.0`` unconditionally, with a
  ``RuntimeWarning`` for an invalid division inside ``Cooling_Tower`` --
  for the one archetype the paper itself identifies as having the
  *strongest* climate-driven WUE signal (Fig. 6, Section 4.2). The repo's
  own ``demo.ipynb`` calls a different case (``PUE_WUE_WE_Chiller_Colo``)
  and never exercises this path, which is presumably why the bug shipped
  unnoticed. This module uses the corrected indices
  (``T_sa`` <- index 5, ``H_sa`` <- index 6, ``H_cd`` <- index 8,
  ``d_sa`` <- index 9), verified to produce physically sensible PUE/WUE
  across a spread of test conditions (PUE 1.11-1.35, WUE 0-2.6 L/kWh over
  mild-to-hot-humid conditions, matching the paper's Fig. 5/6 ranges for
  Case 1) -- confirmed by direct assumption, not by re-deriving the physics.
  Cases 2 and 5 were checked against the same failure mode (their simpler
  ``Chiller_system`` 3-tuple calls) and are correctly indexed as published.
- **Vectorization approach.** ``Air_side_economizer``'s branching logic
  (three temperature regimes, each with 2-3 humidity sub-cases) is
  reproduced with ``numpy.where``/``numpy.select`` cascades over boolean
  masks, evaluating every branch over the full ``(N,)`` flattened
  grid-cell x time array and selecting per-element, rather than
  ``numpy.vectorize`` or a Python loop over cells/days -- a per-cell Python
  loop measured at ~50-100us/call (dominated by ``Air_side_economizer``'s
  several ``HAPropsSI`` sub-calls) would cost hours per (GCM, archetype)
  at this project's grid x daily-time scale, versus seconds for the
  array-vectorized form (``CoolProp.HumidAirProp.HAPropsSI`` accepts numpy
  array inputs directly, confirmed on this project's CoolProp 8.0.0).
  Validated against the (index-corrected) scalar reference implementation
  on a random grid of conditions -- see ``tests/test_cooling_archetypes.py``.
- **GP chiller COP model (``COP_2.pkl``).** Fit under scikit-learn
  0.22.2.post1 (per the pickle's own recorded ``_sklearn_version``, not the
  0.23.1 string found via a raw byte search, which was a different,
  unrelated occurrence in the pickle). Verified to predict identically (to
  floating-point precision) under this project's scikit-learn 1.9.0 as
  under 0.22.2.post1/0.23.1 -- no separate pinned environment needed for
  this specific model (unlike ``COP_AC.pkl``, which is broken under every
  version tested, and ``COP_DX.pkl``, which showed real extrapolation risk
  -- see the ensemble task's Phase 0 report; both are out of scope here).
  Training envelope: ``Twb_oa + AT_CT`` in [15.6, 29.4] degC, ``Chiller_load``
  in [0.25, 1.0] (42 training points). Oregon's SSP3-7.0 grid is expected to
  exceed this envelope at both tails (cold winter days, hot summer
  extremes) -- extrapolation behavior should be checked per-GCM against
  this range before trusting COP predictions on those days (see the task's
  Phase 0 report for the RBF+WhiteKernel reversion-to-prior-mean behavior
  observed just outside this range).
- **Parameter defaults.** Each ``CaseNConfig`` defaults every equipment
  parameter to the midpoint of the range published in Table B.1 of the
  paper itself (not re-derived from Latin-hypercube sample code), extracted
  from the PDF's text layer (not OCR'd from a rendered image) via
  ``pdftotext -layout``, cross-checked against the already-ported 2020-paper
  Table A.1 values in ``pue.py`` where the same parameter appears in both
  papers (``pcop`` bounds and ``l_to_g`` for the AE archetype match exactly
  between the two independently-extracted tables).
- **Midpoint corrections (2026-10-08).** A line-by-line check against
  Table B.1 found three defaults that were not midpoints: Case 1
  ``rh_up_pct`` was 60.0 (range 60-95%, now 77.5), Case 5
  ``delta_t_water_c`` was 5.0 (range 5-10 C, now 7.5) and Case 5
  ``fan_pressure_ct_pa`` was 250 (range 200-400 Pa, now 300). The
  facility ensemble was rerun for Cases 1 and 5 with the corrected values;
  the earlier outputs are kept under
  ``outputs/ensemble_pue_wue_facilities/superseded_rh_up60_2026-10-08/``.
- **Percent vs. fraction convention**, matching the reference repo's own
  mixed convention rather than normalizing it away (to keep the port
  line-by-line checkable): ``UPS_e``, ``PD_lr``, ``L_percentage``, ``SHR``,
  and ``Windage_p`` are used directly (undivided) in the source formulas,
  so their defaults are stored as fractions (e.g. ``ups_e=0.945`` for a
  published "90-99%" range). ``RH_up``/``RH_lw``/``RH_oa`` are divided by
  100 inside the source's own ``HAPropsSI`` calls, so their defaults are
  stored as raw percentages (e.g. ``rh_up=60.0``, not ``0.60``). ``pcop``
  ("COP relative error to regressed value", published as a percent range)
  is used as ``(1 + pcop)`` in the source, so it is stored as a fraction
  (e.g. Case 1's published "-11 - 11%" -> ``pcop=0.0`` default, bounds
  [-0.11, 0.11]).
- **Case 5's ``pcop`` is asymmetric and always non-positive**
  ([-0.40, 0.0], vs. Case 1/2's symmetric [-0.11, 0.11]): the GP was fit on
  chiller performance data more representative of well-optimized
  large-scale systems, so applying it to a midsize archetype needs a
  systematic downward correction, per the paper's own Table B.1 -- not a
  transcription oddity.
- **Elevation-adjusted pressure** (``elevation_adjusted_pressure_pa``) uses
  the standard NOAA/ICAO barometric formula for the troposphere. Replaces
  the standard-atmosphere assumption (101325 Pa) both papers' reference
  implementations use, per the ensemble task's Phase 2 instructions --
  matters far more here than for the facility-level pipeline, since
  Oregon's grid spans sea level to Cascades summits.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from CoolProp.HumidAirProp import HAPropsSI


# --------------------------------------------------------------------------
# Elevation-adjusted pressure
# --------------------------------------------------------------------------


def elevation_adjusted_pressure_pa(elevation_m: np.ndarray | float) -> np.ndarray:
    """Atmospheric pressure (Pa) from elevation via the standard barometric formula.

    Standard NOAA/ICAO troposphere approximation,
    ``P = P0 * (1 - 2.25577e-5 * z)**5.25588``, valid to ~11 km -- comfortably
    covers Oregon's sea-level-to-Cascades-summit range (~3400 m at Mt. Hood).
    Replaces the 101325 Pa standard-atmosphere default both the 2020- and
    2022-paper reference implementations use.
    """
    elevation_m = np.asarray(elevation_m, dtype=float)
    return 101325.0 * (1.0 - 2.25577e-5 * elevation_m) ** 5.25588


# --------------------------------------------------------------------------
# Equipment configs -- Table B.1 midpoints (Lei & Masanet, 2022)
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Case1Config:
    """Case 1: large-scale, airside economizer + adiabatic cooling + water-cooled chiller.

    Defaults are the midpoint of Table B.1's Case 1 column; see this
    module's Assumptions log for the percent/fraction storage convention.
    """

    ups_e: float = 0.945  # [0.90, 0.99]
    pd_lr: float = 0.01  # [0, 0.02]
    l_percentage: float = 0.001  # [0, 0.002]
    delta_t_air_c: float = 16.65  # [13.9, 19.4]
    fan_pressure_crac_pa: float = 650.0  # [300, 1000]
    fan_e_crac: float = 0.775  # [0.65, 0.90]
    pump_pressure_hd_pa: float = 7_000_000.0  # [6.3e6, 7.7e6]
    pump_e_hd: float = 0.7  # [0.6, 0.8]
    at_ct_c: float = 4.75  # [2.8, 6.7]
    chiller_load: float = 0.5  # [0.2, 0.8]
    delta_t_water_c: float = 7.5  # [5, 10]
    pump_pressure_cw_pa: float = 143_650.0  # [114900, 172400]
    pump_e_cw: float = 0.7  # [0.6, 0.8]
    delta_t_ct_c: float = 5.0  # [4, 6]
    pump_pressure_ct_pa: float = 208_650.0  # [166900, 250400]
    pump_e_ct: float = 0.7  # [0.6, 0.8]
    windage_p: float = 0.002525  # [0.00005, 0.005] (published as 0.005-0.5%)
    cc: float = 9.0  # [3, 15]
    fan_pressure_ct_pa: float = 250.0  # [100, 400]
    fan_e_ct: float = 0.775  # [0.65, 0.90]
    shr: float = 0.97  # [0.95, 0.99] -- drives a second (non-AE) latent-heat humidification term
    lg_ratio: float = 2.1  # [0.2, 4]
    t_up_c: float = 31.0  # [27, 35]
    t_lw_c: float = 14.0  # [10, 18]
    dp_up_c: float = 21.0  # [15, 27]
    dp_lw_c: float = -10.5  # [-12, -9]
    rh_up_pct: float = 77.5  # [60, 95]
    rh_lw_pct: float = 14.0  # raw percent, not fraction -- see Assumptions log
    pcop: float = 0.0  # [-0.11, 0.11]

    #: Marked extension point -- overridden per grid cell via
    #: ``elevation_adjusted_pressure_pa`` rather than left at this default.
    p_atm_pa: float = 101325.0


@dataclass(frozen=True)
class Case2Config:
    """Case 2: large-scale, waterside economizer + water-cooled chiller."""

    ups_e: float = 0.945  # [0.90, 0.99]
    pd_lr: float = 0.01  # [0, 0.02]
    l_percentage: float = 0.001  # [0, 0.002]
    shr: float = 0.97  # [0.95, 0.99]
    delta_t_air_c: float = 16.65  # [13.9, 19.4]
    fan_pressure_crac_pa: float = 500.0  # [300, 700]
    fan_e_crac: float = 0.775  # [0.65, 0.90]
    pump_pressure_hd_pa: float = 7_000_000.0  # [6.3e6, 7.7e6]
    pump_e_hd: float = 0.7  # [0.6, 0.8]
    hte: float = 0.8  # [0.7, 0.9]
    delta_t_water_c: float = 7.5  # [5, 10]
    at_ct_c: float = 4.75  # [2.8, 6.7]
    at_he_c: float = 2.25  # [1.7, 2.8]
    pump_pressure_we_pa: float = 143_650.0  # [114900, 172400]
    pump_e_we: float = 0.7  # [0.6, 0.8]
    pump_pressure_cw_pa: float = 143_650.0  # [114900, 172400]
    pump_e_cw: float = 0.7  # [0.6, 0.8]
    chiller_load: float = 0.5  # [0.2, 0.8]
    delta_t_ct_c: float = 5.0  # [4, 6]
    pump_pressure_ct_pa: float = 208_650.0  # [166900, 250400]
    pump_e_ct: float = 0.7  # [0.6, 0.8]
    windage_p: float = 0.002525  # [0.00005, 0.005]
    cc: float = 9.0  # [3, 15]
    lg_ratio: float = 2.1  # [0.2, 4]
    fan_pressure_ct_pa: float = 250.0  # [100, 400]
    fan_e_ct: float = 0.775  # [0.65, 0.90]
    t_up_c: float = 31.0  # [27, 35]
    t_lw_c: float = 14.0  # [10, 18]
    dp_up_c: float = 21.0  # [15, 27]
    dp_lw_c: float = -10.5  # [-12, -9]
    rh_up_pct: float = 75.0  # [60, 90]
    rh_lw_pct: float = 14.0  # [8, 20]
    pcop: float = 0.0  # [-0.11, 0.11]

    p_atm_pa: float = 101325.0


@dataclass(frozen=True)
class Case5Config:
    """Case 5: midsize, water-cooled chiller, no economizer."""

    ups_e: float = 0.87  # [0.80, 0.94]
    pd_lr: float = 0.035  # [0.02, 0.05]
    l_percentage: float = 0.035  # [0.02, 0.05]
    shr: float = 0.97  # [0.95, 0.99]
    delta_t_air_c: float = 7.5  # [5, 10]
    fan_pressure_crac_pa: float = 650.0  # [400, 900]
    fan_e_crac: float = 0.70  # [0.60, 0.80]
    pump_pressure_hd_pa: float = 7_000_000.0  # [6.3e6, 7.7e6]
    pump_e_hd: float = 0.7  # [0.6, 0.8]
    hte: float = 0.775  # [0.65, 0.90]
    delta_t_water_c: float = 7.5  # [5, 10]
    pump_pressure_cw_pa: float = 143_650.0  # [114900, 172400]
    pump_e_cw: float = 0.7  # [0.6, 0.8]
    at_ct_c: float = 4.75  # [2.8, 6.7]
    chiller_load: float = 0.3  # [0.1, 0.5]
    delta_t_ct_c: float = 5.0  # [4, 6]
    pump_pressure_ct_pa: float = 208_650.0  # [166900, 250400]
    pump_e_ct: float = 0.7  # [0.6, 0.8]
    windage_p: float = 0.002525  # [0.00005, 0.005]
    cc: float = 7.5  # [3, 12]
    fan_pressure_ct_pa: float = 300.0  # [200, 400]
    fan_e_ct: float = 0.70  # [0.60, 0.80]
    lg_ratio: float = 1.1  # [0.2, 2]
    t_up_c: float = 29.5  # [27, 32]
    t_lw_c: float = 16.5  # [15, 18]
    dp_up_c: float = 21.0  # [15, 27]
    dp_lw_c: float = -10.5  # [-12, -9]
    rh_up_pct: float = 70.0  # [60, 80]
    rh_lw_pct: float = 20.0  # [10, 30]
    pcop: float = -0.20  # [-0.40, 0.0] -- asymmetric, see Assumptions log

    p_atm_pa: float = 101325.0


# --------------------------------------------------------------------------
# Shared vectorized psychrometric/mechanical helpers
# --------------------------------------------------------------------------

_FAN_POWER = lambda m, rho, dp, e: m / rho * dp / e / 1000.0  # noqa: E731
_PUMP_POWER = lambda m, dp, e, rho: dp * m / (1000.0 * e * rho)  # noqa: E731


def _cooling_tower(
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
    """Vectorized port of ``Cooling_Tower``.

    Returns ``(m_air, water_evaporated, water_windage, water_drainoff)``
    (kg/s each, except ``m_air`` which is mass flow); WUE assembly happens
    in the caller per the source's own structure.
    """
    h_oa = HAPropsSI("H", "T", t_oa_c + 273.15, "P", p_oa_pa, "R", rh_oa_pct / 100.0)
    t_wb_oa = HAPropsSI("Twb", "H", h_oa, "T", t_oa_c + 273.15, "P", p_oa_pa) - 273.15
    latent_heat_vap = -0.0013 * (t_wb_oa + at_ct_c + delta_t_ct_c / 2) ** 2 - 2.3097 * (
        t_wb_oa + at_ct_c + delta_t_ct_c / 2
    ) + 2500.5
    m_ct = q_kw / (4.184 * delta_t_ct_c)
    t_ct = t_wb_oa + at_ct_c + delta_t_ct_c / 2  # noqa: F841 (documents the source's naming; folded into latent_heat_vap)
    water_evaporated = q_kw / latent_heat_vap
    m_air = m_ct / lg_ratio
    water_windage = m_ct * windage_p
    water_drainoff = np.maximum(water_evaporated / (cc - 1) - water_windage, 0.0)
    return m_air, water_evaporated, water_windage, water_drainoff


def _chiller_system(
    t_up_c: float,
    t_lw_c: float,
    dp_up_c: float,
    dp_lw_c: float,
    rh_up_pct: float,
    rh_lw_pct: float,
    p_oa_pa: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized port of ``Chiller_system`` (used by Cases 2 and 5).

    The source picks a random supply-air absolute humidity within
    ``[d_lw, d_up]`` via ``np.random.uniform`` -- reproduced here with the
    midpoint instead, since a random draw per grid cell/timestep would make
    the gridded output non-reproducible run-to-run for no physical reason
    (the source's own comment flags this exact line as a placeholder:
    "how to select absolute humidity"). Documented deviation, not a silent
    change: this shifts Case 2/5's PUE/WUE by a small, bounded amount
    relative to a literal port, in exchange for reproducibility.
    """
    t_sa = np.full_like(np.broadcast_to(p_oa_pa, np.shape(p_oa_pa)), max(t_up_c, t_lw_c), dtype=float)
    d_up = np.minimum(
        HAPropsSI("W", "T", t_sa + 273.15, "RH", rh_up_pct / 100.0, "P", p_oa_pa),
        HAPropsSI("W", "T", t_sa + 273.15, "Tdp", dp_up_c + 273.15, "P", p_oa_pa),
    )
    d_lw = np.maximum(
        HAPropsSI("W", "T", t_sa + 273.15, "RH", rh_lw_pct / 100.0, "P", p_oa_pa),
        HAPropsSI("W", "T", t_sa + 273.15, "Tdp", dp_lw_c + 273.15, "P", p_oa_pa),
    )
    d_sa = (d_up + d_lw) / 2.0
    h_sa = HAPropsSI("H", "T", t_sa + 273.15, "P", p_oa_pa, "W", d_sa)
    return t_sa, d_sa, h_sa


def _waterside_economizer(
    t_sfw_c: np.ndarray,
    t_rfw_c: np.ndarray,
    twb_oa_c: np.ndarray,
    at_ct_c: float,
    at_he_c: float,
    cooling_required_kw: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized port of ``waterside_economizer`` (Case 2)."""
    approach = twb_oa_c + at_ct_c + at_he_c
    full_use = approach <= t_sfw_c
    partial_use = (approach > t_sfw_c) & (approach < t_rfw_c)
    we_heat_removed = np.where(
        full_use,
        cooling_required_kw,
        np.where(
            partial_use,
            cooling_required_kw * ((t_rfw_c - approach) / (t_rfw_c - t_sfw_c)),
            0.0,
        ),
    )
    use = (full_use | partial_use).astype(int)
    return use, we_heat_removed


def _air_side_economizer(
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
    """Vectorized port of ``Air_side_economizer`` (Case 1 only).

    Reproduces all three temperature regimes (``T_oa < T_lw``,
    ``T_lw <= T_oa < T_up``, ``T_oa >= T_up``) and their humidity
    sub-branches via boolean masks over the full array, per this module's
    Assumptions log. Returns a dict keyed by the source's own return-tuple
    names, at the *correct* indices (see Assumptions log re: the published
    indexing bug).
    """
    t_oa_c = np.asarray(t_oa_c, dtype=float)
    rh_oa_pct = np.asarray(rh_oa_pct, dtype=float)
    p_oa_pa = np.broadcast_to(np.asarray(p_oa_pa, dtype=float), t_oa_c.shape).copy()

    h_oa = HAPropsSI("H", "T", t_oa_c + 273.15, "P", p_oa_pa, "RH", rh_oa_pct / 100.0)
    d_oa = HAPropsSI("W", "T", t_oa_c + 273.15, "RH", rh_oa_pct / 100.0, "P", p_oa_pa)

    mask_a = (t_oa_c >= t_lw_c) & (t_oa_c < t_up_c)
    mask_b = t_oa_c < t_lw_c
    mask_c = t_oa_c >= t_up_c

    # ---- Branch A: T_lw <= T_oa < T_up, thresholds evaluated at T_oa ----
    d_up_a = np.minimum(
        HAPropsSI("W", "T", t_oa_c + 273.15, "RH", rh_up_pct / 100.0, "P", p_oa_pa),
        HAPropsSI("W", "T", t_oa_c + 273.15, "Tdp", dp_up_c + 273.15, "P", p_oa_pa),
    )
    d_lw_a = np.maximum(
        HAPropsSI("W", "T", t_oa_c + 273.15, "RH", rh_lw_pct / 100.0, "P", p_oa_pa),
        HAPropsSI("W", "T", t_oa_c + 273.15, "Tdp", dp_lw_c + 273.15, "P", p_oa_pa),
    )
    cond_a1 = (d_oa >= d_lw_a) & (d_oa <= d_up_a)
    cond_a2 = d_oa > d_up_a
    # cond_a3 = d_oa < d_lw_a (implicit else)

    h_cd_a2 = HAPropsSI("H", "T", np.full_like(t_oa_c, dp_up_c) + 273.15, "P", p_oa_pa, "W", d_up_a)
    t_sa_a2 = np.full_like(t_oa_c, max(t_lw_c, dp_up_c))
    h_sa_a2 = HAPropsSI("H", "T", t_sa_a2 + 273.15, "P", p_oa_pa, "W", d_up_a)

    t_cd_a3 = HAPropsSI("T", "H", h_oa, "W", d_lw_a, "P", p_oa_pa) - 273.15
    t_sa_a3 = np.maximum(t_lw_c, t_cd_a3)
    h_sa_a3 = HAPropsSI("H", "T", t_sa_a3 + 273.15, "P", p_oa_pa, "W", d_lw_a)

    ae_use_a = np.where(cond_a1, 1, np.where(cond_a2, 0, 1))
    hd_use_a = np.where(cond_a1, 0, np.where(cond_a2, 0, 1))
    dhd_use_a = np.where(cond_a1, 0, np.where(cond_a2, 1, 0))
    delta_d_a = np.where(cond_a1, 0.0, np.where(cond_a2, d_oa - d_up_a, d_lw_a - d_oa))
    d_sa_a = np.where(cond_a1, d_oa, np.where(cond_a2, d_up_a, d_lw_a))
    t_cd_a = np.where(cond_a1, t_oa_c, np.where(cond_a2, dp_up_c, t_cd_a3))
    h_cd_a = np.where(cond_a1, h_oa, np.where(cond_a2, h_cd_a2, h_oa))
    t_sa_a = np.where(cond_a1, t_oa_c, np.where(cond_a2, t_sa_a2, t_sa_a3))
    h_sa_a = np.where(cond_a1, h_oa, np.where(cond_a2, h_sa_a2, h_sa_a3))
    chiller_energy_a = np.where(cond_a2, h_oa - h_cd_a2, 0.0)

    # ---- Branch B: T_oa < T_lw, thresholds evaluated at T_lw (constant) ----
    t_lw_arr = np.full_like(t_oa_c, t_lw_c)
    d_up_b = np.minimum(
        HAPropsSI("W", "T", t_lw_arr + 273.15, "RH", rh_up_pct / 100.0, "P", p_oa_pa),
        HAPropsSI("W", "T", t_lw_arr + 273.15, "Tdp", dp_up_c + 273.15, "P", p_oa_pa),
    )
    d_lw_b = np.maximum(
        HAPropsSI("W", "T", t_lw_arr + 273.15, "RH", rh_lw_pct / 100.0, "P", p_oa_pa),
        HAPropsSI("W", "T", t_lw_arr + 273.15, "Tdp", dp_lw_c + 273.15, "P", p_oa_pa),
    )
    cond_b1 = d_oa < d_lw_b
    cond_b2 = (d_oa >= d_lw_b) & (d_oa <= d_up_b)
    cond_b3 = d_oa > d_up_b

    t_cd_b1 = HAPropsSI("T", "H", h_oa, "W", d_lw_b, "P", p_oa_pa) - 273.15
    t_sa_b1 = np.maximum(t_lw_c, t_oa_c)
    h_sa_b1 = HAPropsSI("H", "T", t_sa_b1 + 273.15, "P", p_oa_pa, "W", d_lw_b)

    t_sa_b2 = np.maximum(t_lw_c, t_oa_c)
    h_sa_b2 = HAPropsSI("H", "T", t_sa_b2 + 273.15, "P", p_oa_pa, "W", d_oa)

    d_up_at_upsasp = HAPropsSI("W", "T", np.full_like(t_oa_c, t_up_c) + 273.15, "Tdp", dp_up_c + 273.15, "P", p_oa_pa)
    cond_b3a = cond_b3 & (d_oa <= d_up_at_upsasp)
    cond_b3b = cond_b3 & (d_oa > d_up_at_upsasp)
    h_sa_b3a = HAPropsSI("H", "T", np.full_like(t_oa_c, t_lw_c) + 273.15, "P", p_oa_pa, "W", d_oa)
    t_sa_b3b_val = max(t_lw_c, dp_up_c)
    h_cd_b3b = HAPropsSI("H", "T", np.full_like(t_oa_c, dp_up_c) + 273.15, "P", p_oa_pa, "W", d_up_at_upsasp)
    h_sa_b3b = HAPropsSI("H", "T", np.full_like(t_oa_c, t_sa_b3b_val) + 273.15, "P", p_oa_pa, "W", d_up_at_upsasp)
    chiller_energy_b3b = np.maximum(h_cd_b3b - h_oa, 0.0)

    ae_use_b = np.ones_like(t_oa_c)
    hd_use_b = np.where(cond_b1, 1, 0)
    dhd_use_b = np.where(cond_b3b, 1, 0)
    delta_d_b = np.where(cond_b1, d_lw_b - d_oa, np.where(cond_b3b, d_oa - d_up_at_upsasp, 0.0))
    d_sa_b = np.where(cond_b1, d_lw_b, np.where(cond_b2, d_oa, np.where(cond_b3a, d_oa, d_up_at_upsasp)))
    h_cd_b = np.where(cond_b1, h_oa, np.where(cond_b2, h_oa, np.where(cond_b3a, h_oa, h_cd_b3b)))
    t_cd_b = np.where(cond_b1, t_cd_b1, np.where(cond_b2, t_oa_c, np.where(cond_b3a, t_oa_c, dp_up_c)))
    t_sa_b = np.where(cond_b1, t_sa_b1, np.where(cond_b2, t_sa_b2, np.where(cond_b3a, t_lw_c, t_sa_b3b_val)))
    h_sa_b = np.where(cond_b1, h_sa_b1, np.where(cond_b2, h_sa_b2, np.where(cond_b3a, h_sa_b3a, h_sa_b3b)))
    chiller_energy_b = np.where(cond_b3b, chiller_energy_b3b, 0.0)

    # ---- Branch C: T_oa >= T_up (unconditional simplified control policy) ----
    t_up_arr = np.full_like(t_oa_c, t_up_c)
    d_up_c_ = np.minimum(
        HAPropsSI("W", "T", t_up_arr + 273.15, "RH", rh_up_pct / 100.0, "P", p_oa_pa),
        HAPropsSI("W", "T", t_up_arr + 273.15, "Tdp", dp_up_c + 273.15, "P", p_oa_pa),
    )
    h_cd_c = HAPropsSI("H", "T", t_up_arr + 273.15, "P", p_oa_pa, "W", d_up_c_)
    h_ra_c = h_cd_c + 1.01 * delta_t_air_c * 1000.0
    chiller_energy_c = h_ra_c - h_cd_c

    ae_use_c = np.zeros_like(t_oa_c)
    hd_use_c = np.zeros_like(t_oa_c)
    dhd_use_c = np.zeros_like(t_oa_c)
    delta_d_c = np.zeros_like(t_oa_c)
    d_sa_c = d_up_c_
    t_cd_c = t_up_arr
    t_sa_c = t_up_arr
    h_sa_c = h_cd_c

    # ---- Combine the three temperature regimes ----
    def _select(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
        return np.select([mask_a, mask_b, mask_c], [a, b, c])

    return {
        "ae_use": _select(ae_use_a, ae_use_b, ae_use_c),
        "hd_use": _select(hd_use_a, hd_use_b, hd_use_c),
        "dhd_use": _select(dhd_use_a, dhd_use_b, dhd_use_c),
        "delta_d": _select(delta_d_a, delta_d_b, delta_d_c),
        "t_sa": _select(t_sa_a, t_sa_b, t_sa_c),
        "h_sa": _select(h_sa_a, h_sa_b, h_sa_c),
        "t_cd": _select(t_cd_a, t_cd_b, t_cd_c),
        "h_cd": _select(h_cd_a, h_cd_b, h_cd_c),
        "d_sa": _select(d_sa_a, d_sa_b, d_sa_c),
        "chiller_energy": _select(chiller_energy_a, chiller_energy_b, chiller_energy_c),
        "h_oa": h_oa,
    }


# --------------------------------------------------------------------------
# GP chiller COP (water-cooled, shared across all three cases)
# --------------------------------------------------------------------------


def chiller_cop_water_cooled(
    twb_oa_c: np.ndarray,
    at_ct_c: float,
    chiller_load: float,
    pcop: float,
    cop_gp,
) -> np.ndarray:
    """Chiller COP via the water-cooled GP model (``COP_2.pkl``), batched.

    Per the ensemble task's Phase 2 instructions: build one ``(N, 2)``
    feature array spanning every grid cell x timestep and call
    ``cop_gp.predict()`` once, rather than looping. ``cop_gp`` is the
    already-loaded (``pickle.load``) ``GaussianProcessRegressor`` --
    load it once per process and pass it in, not per call.
    """
    twb_oa_c = np.asarray(twb_oa_c, dtype=float)
    features = np.column_stack(
        [twb_oa_c.ravel() + at_ct_c, np.full(twb_oa_c.size, chiller_load)]
    )
    cop = cop_gp.predict(features).reshape(twb_oa_c.shape)
    return cop * (1.0 + pcop)


def load_cop_water_cooled_gp(pickle_path):
    """Load ``COP_2.pkl`` (water-cooled chiller GP) once per process.

    Verified to predict identically under this project's scikit-learn
    1.9.0 as under the fit-time 0.22.2.post1 -- see this module's
    Assumptions log. An ``InconsistentVersionWarning`` from scikit-learn is
    expected and safe to ignore for this specific pickle (checked, not
    assumed); it is not safe to assume the same for a different pickle
    (``COP_AC.pkl`` raises on ``.predict()`` under every version tested).
    """
    import pickle

    with open(pickle_path, "rb") as f:
        return pickle.load(f)


def pue_wue_case1(
    t_oa_c: np.ndarray,
    rh_oa_pct: np.ndarray,
    p_atm_pa: np.ndarray | float,
    cop_gp,
    config: Case1Config = Case1Config(),
) -> tuple[np.ndarray, np.ndarray]:
    """PUE/WUE for Case 1 (large-scale, AE + adiabatic cooling + water-cooled chiller).

    ``t_oa_c``, ``rh_oa_pct`` are array-like of any shape (e.g. flattened
    ``(lat, lon, time)``); ``p_atm_pa`` broadcasts to that shape (typically
    per-grid-cell elevation-adjusted pressure, invariant across the time
    axis). Returns ``(pue, wue)`` arrays of the same shape, WUE in L/kWh.
    """
    power_it = 1.0
    t_oa_c = np.asarray(t_oa_c, dtype=float)
    rh_oa_pct = np.asarray(rh_oa_pct, dtype=float)
    p_atm_pa = np.broadcast_to(np.asarray(p_atm_pa, dtype=float), t_oa_c.shape).copy()

    q = (
        power_it
        + (power_it / config.ups_e - power_it)
        + (power_it / (1 - config.pd_lr) - power_it)
        + power_it * config.l_percentage
    )

    ae = _air_side_economizer(
        config.t_up_c,
        config.t_lw_c,
        config.dp_up_c,
        config.dp_lw_c,
        config.rh_up_pct,
        config.rh_lw_pct,
        t_oa_c,
        rh_oa_pct,
        p_atm_pa,
        config.delta_t_air_c,
    )
    ae_use = ae["ae_use"]
    hd_use = ae["hd_use"]
    delta_d = ae["delta_d"]
    t_sa = ae["t_sa"]
    h_sa = ae["h_sa"] / 1000.0
    h_cd = ae["h_cd"] / 1000.0
    d_sa = ae["d_sa"]

    h_ra = h_sa + 1.01 * config.delta_t_air_c
    m_sa = q / (h_ra - h_sa)
    m_cd = m_sa * (h_sa - h_ra) / (h_cd - h_ra)
    m_cd_dry = m_cd / (1 + d_sa)
    heat = ae_use * m_cd * (h_ra - h_cd)
    density_sa = 1.0 / HAPropsSI("Vha", "T", t_sa + 273.15, "W", d_sa, "P", p_atm_pa)
    power_fan_crac = _FAN_POWER(m_sa, density_sa, config.fan_pressure_crac_pa, config.fan_e_crac)
    hd_amount_ae = np.maximum(hd_use * m_cd_dry * delta_d, 0.0)

    chiller_heat_removed_pre = q - heat
    q_heat_latent = np.where(
        chiller_heat_removed_pre != 0,
        chiller_heat_removed_pre / config.shr - chiller_heat_removed_pre,
        0.0,
    )
    chiller_heat_removed = chiller_heat_removed_pre + q_heat_latent
    hd_amount = np.maximum(q_heat_latent / 2266.0, 0.0)
    power_pump_hd = _PUMP_POWER(hd_amount, config.pump_pressure_hd_pa, config.pump_e_hd, 1000.0)

    twb_oa = HAPropsSI("Twb", "T", t_oa_c + 273.15, "RH", rh_oa_pct / 100.0, "P", p_atm_pa) - 273.15
    cop_chiller = chiller_cop_water_cooled(twb_oa, config.at_ct_c, config.chiller_load, config.pcop, cop_gp)
    power_chiller = chiller_heat_removed / cop_chiller

    m_sw = chiller_heat_removed / (4.184 * config.delta_t_water_c)
    power_pump_cw = _PUMP_POWER(m_sw, config.pump_pressure_cw_pa, config.pump_e_cw, 1000.0)

    ct_heat_removed = chiller_heat_removed + power_chiller
    m_ct = ct_heat_removed / (4.184 * config.delta_t_ct_c)
    power_pump_ct = _PUMP_POWER(m_ct, config.pump_pressure_ct_pa, config.pump_e_ct, 1000.0)

    m_ct_air, water_evap, water_windage, water_drainoff = _cooling_tower(
        t_oa_c,
        rh_oa_pct,
        p_atm_pa,
        config.at_ct_c,
        ct_heat_removed,
        config.delta_t_ct_c,
        config.windage_p,
        config.cc,
        config.lg_ratio,
    )
    density_oa = 1.0 / HAPropsSI("Vha", "T", t_oa_c + 273.15, "RH", rh_oa_pct / 100.0, "P", p_atm_pa)
    power_fan_ct = _FAN_POWER(m_ct_air, density_oa, config.fan_pressure_ct_pa, config.fan_e_ct)

    power_comp = (
        power_it
        + (power_it / config.ups_e - power_it)
        + (power_it / (1 - config.pd_lr) - power_it)
        + power_it * config.l_percentage
        + power_fan_crac
        + power_pump_hd
        + power_chiller
        + power_pump_cw
        + power_pump_ct
        + power_fan_ct
    )
    pue = power_comp / power_it

    water_comp = (
        (hd_amount + hd_amount_ae)
        + np.maximum(water_evap, 0.0)
        + np.maximum(water_windage, 0.0)
        + np.maximum(water_drainoff, 0.0)
    )
    wue = water_comp * 3600.0 / power_it
    return pue, wue


def pue_wue_case2(
    t_oa_c: np.ndarray,
    rh_oa_pct: np.ndarray,
    p_atm_pa: np.ndarray | float,
    cop_gp,
    config: Case2Config = Case2Config(),
) -> tuple[np.ndarray, np.ndarray]:
    """PUE/WUE for Case 2 (large-scale, waterside economizer + water-cooled chiller).

    Note: the published source computes the cooling tower's air mass flow
    and its own water-use terms (evaporation/windage/blowdown) from the
    base IT+UPS+PD+lighting heat load ``Q`` rather than the fuller
    ``CT_heat_removed`` used for the water/pump mass-flow terms elsewhere
    in this same function (and used consistently in Cases 1 and 5) --
    ported as published rather than silently corrected, since (unlike the
    Case 1 indexing bug) it does not produce degenerate output and no
    paper table indicates it is wrong; flagged here as an observed
    asymmetry worth a Methods footnote if this case is reported in detail.
    """
    power_it = 1.0
    t_oa_c = np.asarray(t_oa_c, dtype=float)
    rh_oa_pct = np.asarray(rh_oa_pct, dtype=float)
    p_atm_pa = np.broadcast_to(np.asarray(p_atm_pa, dtype=float), t_oa_c.shape).copy()

    twb_oa = HAPropsSI("Twb", "T", t_oa_c + 273.15, "RH", rh_oa_pct / 100.0, "P", p_atm_pa) - 273.15
    q = (
        power_it
        + (power_it / config.ups_e - power_it)
        + (power_it / (1 - config.pd_lr) - power_it)
        + power_it * config.l_percentage
    )
    q_heat_latent = q / config.shr - q

    t_sa, d_sa, h_sa = _chiller_system(
        config.t_up_c, config.t_lw_c, config.dp_up_c, config.dp_lw_c, config.rh_up_pct, config.rh_lw_pct, p_atm_pa
    )
    h_sa = h_sa / 1000.0
    t_ra = t_sa + config.delta_t_air_c
    m_sa = np.full_like(t_oa_c, q / (1.01 * config.delta_t_air_c))
    density_sa = 1.0 / HAPropsSI("Vha", "T", t_sa + 273.15, "W", d_sa, "P", p_atm_pa)
    power_fan_crac = _FAN_POWER(m_sa, density_sa, config.fan_pressure_crac_pa, config.fan_e_crac)

    hd_amount = np.full_like(t_oa_c, max(q_heat_latent / 2266.0, 0.0))
    power_pump_hd = _PUMP_POWER(hd_amount, config.pump_pressure_hd_pa, config.pump_e_hd, 1000.0)

    cooling_required = q + q_heat_latent
    t_sfw = t_ra - (t_ra - t_sa) / config.hte
    t_rfw = t_sfw + config.delta_t_water_c
    we_use, we_heat_removed = _waterside_economizer(
        t_sfw, t_rfw, twb_oa, config.at_ct_c, config.at_he_c, np.full_like(t_oa_c, cooling_required)
    )
    m_sfw = we_heat_removed / (4.2 * config.delta_t_water_c)
    power_pump_we = _PUMP_POWER(m_sfw, config.pump_pressure_we_pa, config.pump_e_we, 1000.0)

    chiller_heat_removed = cooling_required - we_heat_removed
    m_sw = chiller_heat_removed / (4.2 * config.delta_t_water_c)
    power_pump_cw = _PUMP_POWER(m_sw, config.pump_pressure_cw_pa, config.pump_e_cw, 1000.0)

    cop_chiller = chiller_cop_water_cooled(twb_oa, config.at_ct_c, config.chiller_load, config.pcop, cop_gp)
    power_chiller = chiller_heat_removed / cop_chiller

    ct_heat_removed = cooling_required + power_chiller  # noqa: F841 (see docstring: source uses bare Q below, not this)
    m_ct = np.full_like(t_oa_c, q) / (4.184 * config.delta_t_ct_c)
    power_pump_ct = _PUMP_POWER(m_ct, config.pump_pressure_ct_pa, config.pump_e_ct, 1000.0)

    m_ct_air, water_evap, water_windage, water_drainoff = _cooling_tower(
        t_oa_c,
        rh_oa_pct,
        p_atm_pa,
        config.at_ct_c,
        np.full_like(t_oa_c, q),
        config.delta_t_ct_c,
        config.windage_p,
        config.cc,
        config.lg_ratio,
    )
    density_oa = 1.0 / HAPropsSI("Vha", "T", t_oa_c + 273.15, "RH", rh_oa_pct / 100.0, "P", p_atm_pa)
    power_fan_ct = _FAN_POWER(m_ct_air, density_oa, config.fan_pressure_ct_pa, config.fan_e_ct)

    power_comp = (
        power_it
        + (power_it / config.ups_e - power_it)
        + (power_it / (1 - config.pd_lr) - power_it)
        + power_it * config.l_percentage
        + power_fan_crac
        + power_pump_hd
        + power_pump_we
        + power_chiller
        + power_pump_cw
        + power_pump_ct
        + power_fan_ct
    )
    pue = power_comp / power_it

    water_comp = (
        hd_amount + np.maximum(water_evap, 0.0) + np.maximum(water_windage, 0.0) + np.maximum(water_drainoff, 0.0)
    )
    wue = water_comp * 3600.0 / power_it
    return pue, wue


def pue_wue_case5(
    t_oa_c: np.ndarray,
    rh_oa_pct: np.ndarray,
    p_atm_pa: np.ndarray | float,
    cop_gp,
    config: Case5Config = Case5Config(),
) -> tuple[np.ndarray, np.ndarray]:
    """PUE/WUE for Case 5 (midsize, water-cooled chiller, no economizer)."""
    power_it = 1.0
    t_oa_c = np.asarray(t_oa_c, dtype=float)
    rh_oa_pct = np.asarray(rh_oa_pct, dtype=float)
    p_atm_pa = np.broadcast_to(np.asarray(p_atm_pa, dtype=float), t_oa_c.shape).copy()

    q = (
        power_it
        + (power_it / config.ups_e - power_it)
        + (power_it / (1 - config.pd_lr) - power_it)
        + power_it * config.l_percentage
    )
    q_heat_latent = q / config.shr - q

    t_sa, d_sa, h_sa = _chiller_system(
        config.t_up_c, config.t_lw_c, config.dp_up_c, config.dp_lw_c, config.rh_up_pct, config.rh_lw_pct, p_atm_pa
    )
    h_sa = h_sa / 1000.0
    m_sa = np.full_like(t_oa_c, q / (1.01 * config.delta_t_air_c))
    density_sa = 1.0 / HAPropsSI("Vha", "T", t_sa + 273.15, "W", d_sa, "P", p_atm_pa)
    power_fan_crac = _FAN_POWER(m_sa, density_sa, config.fan_pressure_crac_pa, config.fan_e_crac)

    hd_amount = np.full_like(t_oa_c, max(q_heat_latent / 2266.0, 0.0))
    power_pump_hd = _PUMP_POWER(hd_amount, config.pump_pressure_hd_pa, config.pump_e_hd, 1000.0)

    chiller_heat_removed = np.full_like(t_oa_c, q_heat_latent + q)
    m_sw = chiller_heat_removed / (4.2 * config.delta_t_water_c)
    power_pump_cw = _PUMP_POWER(m_sw, config.pump_pressure_cw_pa, config.pump_e_cw, 1000.0)

    twb_oa = HAPropsSI("Twb", "T", t_oa_c + 273.15, "RH", rh_oa_pct / 100.0, "P", p_atm_pa) - 273.15
    cop_chiller = chiller_cop_water_cooled(twb_oa, config.at_ct_c, config.chiller_load, config.pcop, cop_gp)
    power_chiller = chiller_heat_removed / cop_chiller

    ct_heat_removed = chiller_heat_removed + power_chiller
    m_ct = ct_heat_removed / (4.184 * config.delta_t_ct_c)
    power_pump_ct = _PUMP_POWER(m_ct, config.pump_pressure_ct_pa, config.pump_e_ct, 1000.0)

    m_ct_air, water_evap, water_windage, water_drainoff = _cooling_tower(
        t_oa_c,
        rh_oa_pct,
        p_atm_pa,
        config.at_ct_c,
        ct_heat_removed,
        config.delta_t_ct_c,
        config.windage_p,
        config.cc,
        config.lg_ratio,
    )
    density_oa = 1.0 / HAPropsSI("Vha", "T", t_oa_c + 273.15, "RH", rh_oa_pct / 100.0, "P", p_atm_pa)
    power_fan_ct = _FAN_POWER(m_ct_air, density_oa, config.fan_pressure_ct_pa, config.fan_e_ct)

    power_comp = (
        power_it
        + (power_it / config.ups_e - power_it)
        + (power_it / (1 - config.pd_lr) - power_it)
        + power_it * config.l_percentage
        + power_fan_crac
        + power_pump_hd
        + power_chiller
        + power_pump_cw
        + power_pump_ct
        + power_fan_ct
    )
    pue = power_comp / power_it

    water_comp = (
        hd_amount + np.maximum(water_evap, 0.0) + np.maximum(water_windage, 0.0) + np.maximum(water_drainoff, 0.0)
    )
    wue = water_comp * 3600.0 / power_it
    return pue, wue


def check_pue_wue_plausible(
    pue: np.ndarray,
    wue: np.ndarray,
    pue_bounds: tuple[float, float] = (1.0, 3.0),
    wue_bounds: tuple[float, float] = (0.0, 6.0),
) -> None:
    """Sanity-check PUE/WUE ranges, mirroring ``pue.check_pue_plausible``.

    WUE upper bound of 6 L/kWh is deliberately generous relative to the
    paper's Fig. 6 case-1 range (up to ~4 L/kWh) to avoid flagging
    legitimate extreme-heat days while still catching a broken unit
    conversion or NaN propagation (e.g. the Case 1 indexing bug this module
    corrects, which produced ``WUE ~= 0`` unconditionally).
    """
    for name, arr, (lower, upper) in (("PUE", pue, pue_bounds), ("WUE", wue, wue_bounds)):
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            raise ValueError(f"No finite {name} values found for plausibility check.")
        if finite.min() < lower or finite.max() > upper:
            raise ValueError(
                f"{name} values outside plausible range [{lower}, {upper}]: "
                f"min={finite.min():.4f}, max={finite.max():.4f}."
            )
