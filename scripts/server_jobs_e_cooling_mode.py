"""Server job E (draft 4): Case 1 cooling-mode shares (supplemental chiller on vs. economizer only).

Where the decision is made: ``cooling_archetypes_fast.pue_wue_case1`` calls
``_air_side_economizer_fast``, which returns ``ae_use`` per pairing; the
chiller load is then ``chiller_heat_removed_pre = q - ae_use * m_cd * (h_ra - h_cd)``
(plus a latent term that is zero when that is zero). This script calls
``_air_side_economizer_fast`` with the Case 1 midpoint configuration
(``Case1Config()``, draft3_values.md section H) and repeats those few lines of
``pue_wue_case1`` to recover the chiller load. No existing code is changed.

Mode per pairing (hot = daily Tmax with RHmin, cool = Tmin with RHmax):
``chiller`` if ``chiller_heat_removed_pre > 1e-9 * q``, else ``economizer``
(free or adiabatic cooling only). Chiller pairings are split by cause:
``hot`` (T_oa >= t_up, branch C) or ``humid`` (branch A, outdoor humidity
ratio above the supply upper limit).

Shares per facility, model, period and calendar month (all days of that
month in the window); equal-weight ecoregion mean of facilities per model;
mean and p25/p75 across the 20 models. Two shares: of all pairings (two per
day), and of days with the hot pairing in chiller mode.

Validation: on one model's historical Case 1 inputs the full
``pue_wue_case1`` is rerun at one facility per grid cell and compared with
the stored ``pue_hot_pairing`` / ``pue_cool_pairing``, and pairings classed
``economizer`` must carry no chiller power.

Writes ``docs/results-figures/server_jobs_oct23/jobE_cooling_mode_shares.csv``
(all ecoregions, all months) and ``jobE_cooling_mode_JJAS_CP_WV.csv``.

Assumptions log
----------------
- Facilities in one LOCA2 cell share inputs and pressure, so the mode is
  computed once per cell (representative facility) and applied to every
  facility in the cell; ecoregion values are facility-weighted as in draft 3.
- The tolerance 1e-9 * q only absorbs floating-point residue: with
  ``ae_use = 1`` the economizer heat equals ``q`` algebraically
  (``m_sa = q / (h_ra - h_sa)``), so the chiller is either off or carries
  the full load. The script reports the count of partial-load pairings.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from climate_risk_dc.climate.cooling_archetypes import Case1Config
from climate_risk_dc.climate.cooling_archetypes_fast import _air_side_economizer_fast
from server_jobs_common import MONTHS, N_MODELS, ORDER, OUT_DIR, PERIOD_LABEL, PERIODS, daily_path, eco, facilities, \
    models

CFG = Case1Config()
COLS = ["facility_id", "date", "tasmax_c", "tasmin_c", "hursmax_pct", "hursmin_pct", "pue_hot_pairing", "pue_cool_pairing"]


def chiller_load(t_oa_c: np.ndarray, rh_oa_pct: np.ndarray, p_pa: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """(chiller_heat_removed_pre, ae_use, q) exactly as computed inside ``pue_wue_case1``."""
    c = CFG
    q = 1.0 + (1.0 / c.ups_e - 1.0) + (1.0 / (1 - c.pd_lr) - 1.0) + c.l_percentage
    ae = _air_side_economizer_fast(c.t_up_c, c.t_lw_c, c.dp_up_c, c.dp_lw_c, c.rh_up_pct, c.rh_lw_pct,
                                   t_oa_c, rh_oa_pct, p_pa, c.delta_t_air_c)
    h_sa_kj, h_cd_kj = ae["h_sa"] / 1000.0, ae["h_cd"] / 1000.0
    h_ra_kj = (ae["h_sa"] + 1010.0 * c.delta_t_air_c) / 1000.0
    m_sa = q / (h_ra_kj - h_sa_kj)
    m_cd = m_sa * (h_sa_kj - h_ra_kj) / (h_cd_kj - h_ra_kj)
    heat = ae["ae_use"] * m_cd * (h_ra_kj - h_cd_kj)
    return q - heat, ae["ae_use"], q


def _one(args) -> pd.DataFrame:
    path, reps, pressure, gcm, period = args
    df = pd.read_csv(path, usecols=COLS[:6])
    df = df[df["facility_id"].isin(reps)].reset_index(drop=True)
    assert df["hursmin_pct"].between(0, 100.5).all() and df["hursmax_pct"].between(0, 100.5).all()
    p = df["facility_id"].map(pressure).to_numpy()
    out = df[["facility_id"]].copy()
    out["month"] = df["date"].str.slice(5, 7).astype(int)
    stats = {}
    for pair, tcol, rcol in (("hot", "tasmax_c", "hursmin_pct"), ("cool", "tasmin_c", "hursmax_pct")):
        tc = df[tcol].to_numpy()
        load, ae_use, q = chiller_load(tc, df[rcol].to_numpy(), p)
        on = load > 1e-9 * q
        stats[f"partial_{pair}"] = int(((load > 1e-9 * q) & (load < q * (1 - 1e-9))).sum())
        stats[f"neg_{pair}"] = int((load < -1e-9 * q).sum())
        assert (on == (ae_use == 0)).all(), "chiller-on does not coincide with ae_use == 0"
        out[f"chiller_{pair}"] = on
        out[f"chiller_{pair}_hot"] = on & (tc >= CFG.t_up_c)
        out[f"chiller_{pair}_humid"] = on & (tc < CFG.t_up_c)
    out["chiller_pairings"] = (out["chiller_hot"].astype(int) + out["chiller_cool"].astype(int)) / 2.0
    out["chiller_pairings_hot_T"] = (out["chiller_hot_hot"].astype(int) + out["chiller_cool_hot"].astype(int)) / 2.0
    out["chiller_pairings_humid"] = (out["chiller_hot_humid"].astype(int) + out["chiller_cool_humid"].astype(int)) / 2.0
    agg = out.groupby(["facility_id", "month"])[
        ["chiller_pairings", "chiller_pairings_hot_T", "chiller_pairings_humid", "chiller_hot", "chiller_hot_hot",
         "chiller_hot_humid", "chiller_cool"]].mean().reset_index()
    return agg.assign(gcm=gcm, period=period, **stats)


def validate(root: Path, fac: pd.DataFrame, reps: pd.Series, gcm: str) -> None:
    """Rerun pue_wue_case1 at the 22 representative facilities; compare with stored pairings."""
    from climate_risk_dc.climate.cooling_archetypes import load_cop_water_cooled_gp
    from climate_risk_dc.climate.cooling_archetypes_fast import pue_wue_case1

    cop_gp = load_cop_water_cooled_gp(root / "data/lei_masanet_2022/COP_2.pkl")
    df = pd.read_csv(daily_path(root, gcm, "ae-chiller", "historical"), usecols=COLS)
    df = df[df["facility_id"].isin(set(reps))].reset_index(drop=True)
    p = df["facility_id"].map(dict(zip(fac["facility_id"], fac["pressure_pa"]))).to_numpy()
    for pair, tcol, rcol in (("hot", "tasmax_c", "hursmin_pct"), ("cool", "tasmin_c", "hursmax_pct")):
        pue, _ = pue_wue_case1(df[tcol].to_numpy(), df[rcol].to_numpy(), p, cop_gp, CFG)
        diff = np.abs(pue - df[f"pue_{pair}_pairing"].to_numpy())
        load, _, q = chiller_load(df[tcol].to_numpy(), df[rcol].to_numpy(), p)
        on = load > 1e-9 * q
        print(f"[validate] E {gcm} historical {pair} pairing: rerun vs stored PUE max |diff| = {diff.max():.2e}; "
              f"PUE economizer-only days {pue[~on].min():.4f}-{pue[~on].max():.4f}, "
              f"chiller days {pue[on].min():.4f}-{pue[on].max():.4f} (n = {on.sum()})" if on.any() else "")
        assert diff.max() < 1e-6, "reconstructed inputs do not reproduce stored PUE"


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    fac = facilities(root, with_pressure=True)
    reps = fac.sort_values("facility_id").groupby("grid_cell")["facility_id"].first()
    gcms = models(root)
    validate(root, fac, reps, gcms[0])

    pressure = dict(zip(fac["facility_id"], fac["pressure_pa"]))
    jobs = [(daily_path(root, g, "ae-chiller", p), set(reps), pressure, g, p) for g in gcms for p in PERIODS]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        parts = list(pool.map(_one, jobs))
    d = pd.concat(parts, ignore_index=True)
    print(f"[check] partial-load pairings: hot {d.drop_duplicates(['gcm', 'period'])['partial_hot'].sum()}, "
          f"cool {d.drop_duplicates(['gcm', 'period'])['partial_cool'].sum()}; negative loads: "
          f"{d.drop_duplicates(['gcm', 'period'])[['neg_hot', 'neg_cool']].sum().sum()}")
    assert d["gcm"].nunique() == N_MODELS

    cell = dict(zip(reps.values, reps.index))
    d["grid_cell"] = d["facility_id"].map(cell)
    d = d.drop(columns="facility_id")
    pm = fac[["facility_id", "abbr", "grid_cell"]].merge(d, on="grid_cell")

    shares = {"chiller_pairings": "share of pairings, chiller on",
              "chiller_pairings_hot_T": "share of pairings, chiller on (T >= t_up)",
              "chiller_pairings_humid": "share of pairings, chiller on (humidity)",
              "chiller_hot": "share of days, hot pairing chiller on",
              "chiller_hot_hot": "share of days, hot pairing chiller on (T >= t_up)",
              "chiller_hot_humid": "share of days, hot pairing chiller on (humidity)",
              "chiller_cool": "share of days, cool pairing chiller on"}
    rows = []
    for (period, month), g in pm.groupby(["period", "month"]):
        for col, label in shares.items():
            e = eco(g, col)
            for a in ORDER:
                rows.append({"period": PERIOD_LABEL[period], "month": month, "month_name": MONTHS[month - 1],
                             "ecoregion": a, "quantity": col, "description": label,
                             "mean": round(e.loc[a, "mean"], 4), "p25": round(e.loc[a, "p25"], 4),
                             "p75": round(e.loc[a, "p75"], 4)})
    out = pd.DataFrame(rows)
    out.to_csv(root / OUT_DIR / "jobE_cooling_mode_shares.csv", index=False)
    sel = out[out["ecoregion"].isin(["CP", "WV"]) & out["month"].between(6, 9)]
    sel.to_csv(root / OUT_DIR / "jobE_cooling_mode_JJAS_CP_WV.csv", index=False)

    # Models (of 20) in which July's share exceeds June's and August's
    agree = []
    for q in ("chiller_pairings", "chiller_hot"):
        m = pm.groupby(["period", "abbr", "gcm", "month"])[q].mean().unstack("month")
        for (period, abbr), g in m.groupby(level=["period", "abbr"]):
            if abbr in ("CP", "WV"):
                agree.append({"quantity": q, "ecoregion": abbr, "period": PERIOD_LABEL[period],
                              "n_models_jul_gt_jun": int((g[7] > g[6]).sum()),
                              "n_models_jul_gt_aug": int((g[7] > g[8]).sum()), "of": N_MODELS})
    agree = pd.DataFrame(agree)
    agree.to_csv(root / OUT_DIR / "jobE_july_agreement.csv", index=False)
    print(agree.to_string(index=False))

    with pd.option_context("display.width", 250, "display.max_rows", 500):
        for q in ("chiller_pairings", "chiller_hot", "chiller_hot_hot", "chiller_hot_humid"):
            s = sel[sel["quantity"] == q]
            s = s.assign(v=s.apply(lambda r: f"{100 * r['mean']:.1f} [{100 * r['p25']:.1f}, {100 * r['p75']:.1f}]", axis=1))
            print(f"\n{shares[q]} (%):")
            print(s.pivot_table(index=["ecoregion", "period"], columns="month_name", values="v", aggfunc="first")[
                ["Jun", "Jul", "Aug", "Sep"]])


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Server job E: Case 1 cooling-mode shares.")
    p.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    p.add_argument("--workers", type=int, default=16)
    return p.parse_args()


if __name__ == "__main__":
    run(parse_args())
