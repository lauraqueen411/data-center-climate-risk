"""Server job D (draft 4): wet-bulb method check.

Historical (1985-2014) afternoon wet-bulb (daily Tmax with RHmin) at the
facility's elevation-adjusted pressure, one representative facility per
LOCA2 grid cell (22 cells), all 20 models, by three methods:

1. ``cooling_archetypes_fast._wetbulb_from_rh`` (pipeline: Arden Buck
   saturation vapour pressure, ASHRAE enthalpy, 8 Newton iterations).
2. The same enthalpy-matching Newton solve with Hyland-Wexler saturation
   vapour pressure over liquid water (ASHRAE Fundamentals, ch. 1, eq. 6).
3. CoolProp ``HAPropsSI('B', ...)`` as the reference.

Validation first: (1) mapped back to every facility in its cell must
reproduce ``hist_wetbulb_C`` in ``draft3/table1.csv``.

Part 2 compares the three ACCESS-CM2 ``*_historical_coolprop.csv`` files with
their ``_fast`` counterparts (pue_daily, wue_daily).

Writes to ``docs/results-figures/server_jobs_oct23/``: ``jobD_wetbulb_summary.csv``
(method differences overall and by bin), ``jobD_wetbulb_ecoregion.csv``,
``jobD_coolprop_vs_fast.csv``.

Assumptions log
----------------
- **Representative facility per cell:** facilities in one cell share the
  daily inputs, the pressure and the ecoregion (checked), so a cell's value
  applies exactly to every facility in it. Ecoregion values are therefore
  computed as draft 3 does (equal weight per facility), not per cell.
- **Bins** are on the CoolProp wet-bulb (reference).
- **RH inputs** are used as stored. For CoolProp, RH is clipped to
  [1e-4, 1] because HAPropsSI rejects RH = 0 or > 1; the pipeline methods
  take RH unclipped, as the pipeline does. LOCA2 ``hurs`` is relative to
  liquid water above 0 °C; methods (1) and (2) treat it as over liquid water
  at all temperatures, while CoolProp uses ice saturation below 0 °C, so
  part of any sub-zero difference is a definition difference, not error.
- **Method (2) uses the same 8 Newton iterations** and finite-difference
  derivative as (1), so only the saturation formula differs.
- **CoolProp files are older than the current Case 1 and Case 5 runs.** The
  ``_coolprop`` files were written on 2026-09-11, before commit ac70bb8
  corrected the Case 1 ``rh_up_pct`` (60 -> 77.5) and the Case 5
  ``delta_t_water_c`` and ``fan_pressure_ct_pa`` midpoints. The like-for-like
  fast counterpart for Cases 1 and 5 is therefore the superseded (pre-
  correction) ``_fast`` file in ``superseded_rh_up60_2026-10-08/``; Case 2 was
  not rerun, so its current ``_fast`` file is like-for-like. Both are
  reported.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from climate_risk_dc.climate.cooling_archetypes_fast import _h_from_t_w, _wetbulb_from_rh
from server_jobs_common import DRAFT3, FAC_DIR, NAME, ORDER, OUT_DIR, compare, daily_path, eco, facilities, models, t

BINS = [-np.inf, 0.0, 5.0, 10.0, np.inf]
BIN_LABELS = ["< 0", "0–5", "5–10", "> 10"]


def _p_ws_hyland_wexler(t_c: np.ndarray) -> np.ndarray:
    """Saturation vapour pressure over liquid water (Pa), Hyland-Wexler, ASHRAE Fundamentals ch. 1 eq. 6."""
    tk = np.asarray(t_c, dtype=float) + 273.15
    return np.exp(-5.8002206e3 / tk + 1.3914993 - 4.8640239e-2 * tk + 4.1764768e-5 * tk**2
                  - 1.4452093e-8 * tk**3 + 6.5459673 * np.log(tk))


def _w_hw(t_c, rh_frac, p_pa):
    p_v = rh_frac * _p_ws_hyland_wexler(t_c)
    return 0.62198 * p_v / (p_pa - p_v)


def wetbulb_hyland_wexler(t_c: np.ndarray, rh_frac: np.ndarray, p_pa: np.ndarray, n_iter: int = 8) -> np.ndarray:
    """Same enthalpy-matching Newton solve as ``_wetbulb_from_rh``, with Hyland-Wexler p_ws."""
    h_oa = _h_from_t_w(t_c, _w_hw(t_c, rh_frac, p_pa))
    twb = np.asarray(t_c, dtype=float) * 0.7
    for _ in range(n_iter):
        h_sat = _h_from_t_w(twb, _w_hw(twb, 1.0, p_pa))
        eps = 1e-4
        df = (_h_from_t_w(twb + eps, _w_hw(twb + eps, 1.0, p_pa)) - h_sat) / eps
        twb = twb - (h_sat - h_oa) / df
    return twb


def _one_model(args) -> pd.DataFrame:
    from CoolProp.HumidAirProp import HAPropsSI

    path, reps, pressure = args
    df = pd.read_csv(path, usecols=["facility_id", "date", "tasmax_c", "hursmin_pct"])
    df = df[df["facility_id"].isin(reps)].reset_index(drop=True)
    assert df["hursmin_pct"].between(0, 100.5).all(), f"{path.name}: hursmin outside 0-100%"
    tc, rh, p = df["tasmax_c"].to_numpy(), df["hursmin_pct"].to_numpy() / 100.0, df["facility_id"].map(pressure).to_numpy()
    df["wb_pipeline"] = _wetbulb_from_rh(tc, rh, p)
    df["wb_hw"] = wetbulb_hyland_wexler(tc, rh, p)
    df["wb_coolprop"] = HAPropsSI("B", "T", tc + 273.15, "P", p, "R", np.clip(rh, 1e-4, 1.0)) - 273.15
    df["year"] = df["date"].str.slice(0, 4).astype(int)
    return df.drop(columns=["date"])


def wetbulb_part(root: Path, fac: pd.DataFrame, gcms: list[str], workers: int) -> tuple[list[str], list[str]]:
    reps = fac.sort_values("facility_id").groupby("grid_cell")["facility_id"].first()
    pressure = dict(zip(fac["facility_id"], fac["pressure_pa"]))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        parts = list(pool.map(_one_model, [(daily_path(root, g, "ae-chiller", "historical"), set(reps), pressure) for g in gcms]))
    d = pd.concat([p.assign(gcm=g) for g, p in zip(gcms, parts)], ignore_index=True)
    assert np.isfinite(d[["wb_pipeline", "wb_hw", "wb_coolprop"]]).all().all(), "non-finite wet-bulb"
    assert (d.groupby("gcm")["facility_id"].nunique() == 22).all() and d.groupby("gcm")["year"].nunique().eq(30).all()

    # --- point-level differences vs CoolProp
    d["bin"] = pd.cut(d["wb_coolprop"], BINS, labels=BIN_LABELS, right=False)
    rows = []
    for m, col in (("(1) pipeline Buck", "wb_pipeline"), ("(2) Hyland-Wexler", "wb_hw")):
        diff = d[col] - d["wb_coolprop"]
        for b, sel in [("all", slice(None))] + [(lab, d["bin"] == lab) for lab in BIN_LABELS]:
            s = diff[sel]
            rows.append({"method": m, "coolprop_wetbulb_bin_C": b, "n_points": len(s), "share_of_points": round(len(s) / len(d), 4),
                         "mean_diff_C": round(s.mean(), 3), "mean_abs_diff_C": round(s.abs().mean(), 3),
                         "max_abs_diff_C": round(s.abs().max(), 3), "p99_abs_diff_C": round(s.abs().quantile(0.99), 3)})
    summ = pd.DataFrame(rows)
    summ.to_csv(root / OUT_DIR / "jobD_wetbulb_summary.csv", index=False)
    print(summ.to_string(index=False))

    # --- ecoregion annual means: annual mean per cell and model, period mean, map to facilities
    cell = dict(zip(reps.values, reps.index))
    d["grid_cell"] = d["facility_id"].map(cell)
    per = d.groupby(["grid_cell", "gcm", "year"])[["wb_pipeline", "wb_hw", "wb_coolprop"]].mean().groupby(
        level=["grid_cell", "gcm"]).mean().reset_index()
    pm = fac[["facility_id", "abbr", "grid_cell"]].merge(per, on="grid_cell")
    erows, lines = [], []
    for m, col in (("(1) pipeline Buck", "wb_pipeline"), ("(2) Hyland-Wexler", "wb_hw"), ("(3) CoolProp", "wb_coolprop")):
        e = eco(pm, col)
        order = list(e["mean"].sort_values().index)
        for a in ORDER:
            erows.append({"method": m, "ecoregion": a, "mean": t(e.loc[a, "mean"], "temp"), "p25": t(e.loc[a, "p25"], "temp"),
                          "p75": t(e.loc[a, "p75"], "temp"), "mean_raw": round(e.loc[a, "mean"], 4)})
        gap = e.loc["WV", "mean"] - e.loc["KM", "mean"]
        lines.append(f"{m}: order {' < '.join(order)} (expected ECSF < BM < CP < KM < WV: "
                     f"{'holds' if order == ORDER else 'DOES NOT HOLD'}); WV - KM = {gap:+.3f} °C")
    ecodf = pd.DataFrame(erows)
    ecodf.to_csv(root / OUT_DIR / "jobD_wetbulb_ecoregion.csv", index=False)
    print(ecodf.pivot(index="ecoregion", columns="method", values="mean_raw").loc[ORDER])
    print("\n".join(lines))

    # --- validation: method (1) vs draft 3 Table 1 and the Figure 2 cache
    ref = pd.read_csv(root / DRAFT3 / "table1.csv", dtype=str)
    ref["abbr"] = ref["ecoregion"].map({v: k for k, v in NAME.items()})
    ours = ecodf[ecodf["method"] == "(1) pipeline Buck"].set_index("ecoregion")["mean"]
    mism = compare("D hist_wetbulb_C", ours, ref.drop_duplicates("abbr").set_index("abbr")["hist_wetbulb_C"])
    cache = pd.read_csv(root / FAC_DIR / "analysis/fig_show_everything/model_values_cache.csv")
    c = cache[cache["metric"] == "wetbulb_annual"].set_index(["facility_id", "gcm"])["value"]
    o = pm.set_index(["facility_id", "gcm"])["wb_pipeline"]
    print(f"[validate] D method (1) vs cached wetbulb_annual, all 109 facilities x 20 models: max |diff| = {(o - c).abs().max():.2e}")
    return mism, lines


def coolprop_part(root: Path, fac: pd.DataFrame) -> pd.DataFrame:
    sup = root / FAC_DIR / "superseded_rh_up60_2026-10-08"
    rows = []
    for arch, case in (("ae-chiller", "Case 1"), ("we-chiller", "Case 2"), ("chiller-only", "Case 5")):
        cp = pd.read_csv(daily_path(root, "ACCESS-CM2", arch, "historical", "coolprop"),
                         usecols=["facility_id", "date", "pue_daily", "wue_daily"])
        refs = {"current _fast": daily_path(root, "ACCESS-CM2", arch, "historical")}
        if arch != "we-chiller":
            refs["superseded _fast (same parameters as coolprop)"] = sup / daily_path(root, "ACCESS-CM2", arch, "historical").name
        for label, p in refs.items():
            fa = pd.read_csv(p, usecols=["facility_id", "date", "pue_daily", "wue_daily"])
            m = cp.merge(fa, on=["facility_id", "date"], suffixes=("_cp", "_fast"), validate="one_to_one")
            assert len(m) == len(cp) == len(fa)
            m["year"] = m["date"].str.slice(0, 4).astype(int)
            r = {"case": case, "fast_file": label, "n_days": len(m)}
            for v in ("pue", "wue"):
                dd = m[f"{v}_daily_fast"] - m[f"{v}_daily_cp"]
                r[f"{v}_mean_diff_fast_minus_cp"] = dd.mean()
                r[f"{v}_mean_abs_diff"] = dd.abs().mean()
                r[f"{v}_max_abs_diff"] = dd.abs().max()
            ann = m.groupby(["facility_id", "year"])[[f"{v}_daily_{s}" for v in ("pue", "wue") for s in ("cp", "fast")]].mean()
            per = ann.groupby(level="facility_id").mean().reset_index().merge(fac[["facility_id", "abbr"]], on="facility_id")
            e = per.groupby("abbr").mean(numeric_only=True).loc[ORDER]
            for a in ORDER:
                for v in ("pue", "wue"):
                    r[f"{a}_{v}_cp"] = e.loc[a, f"{v}_daily_cp"]
                    r[f"{a}_{v}_fast_minus_cp"] = e.loc[a, f"{v}_daily_fast"] - e.loc[a, f"{v}_daily_cp"]
            rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(root / OUT_DIR / "jobD_coolprop_vs_fast.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 50, "display.float_format", "{:.5f}".format):
        print(out.T)
    return out


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    fac = facilities(root, with_pressure=True)
    assert (fac.groupby("grid_cell")[["pressure_pa", "abbr"]].nunique() == 1).all().all()
    mism, _ = wetbulb_part(root, fac, models(root), args.workers)
    coolprop_part(root, fac)
    print("\nMISMATCHES:\n" + "\n".join(mism) if mism else "\nAll validations match draft 3.")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Server job D: wet-bulb method check.")
    p.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    p.add_argument("--workers", type=int, default=16)
    return p.parse_args()


if __name__ == "__main__":
    run(parse_args())
