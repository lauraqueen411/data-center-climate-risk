"""Server jobs A and B (draft 4): historical Tmax for Table 1 and absolute monthly PUE/WUE for Figure 4.

Job A: 1985-2014 annual mean daily Tmax by ecoregion (°C), mean and p25/p75
across the 20 models, from ``analysis/facility_annual_tmax.csv``. Validation:
the mid- and end-century changes must reproduce ``draft3/warming.csv``.

Job B: absolute monthly mean PUE and WUE, Cases 1, 2, 5, all three periods,
Columbia Plateau and Willamette Valley, from
``analysis/facility_month_summary.csv``. Validation: mid-century minus
historical reproduces ``draft3/monthly.csv``.

Writes ``docs/results-figures/server_jobs_oct23/jobA_hist_tmax.csv`` and
``jobB_monthly_pue_wue.csv``.

Assumptions log
----------------
- Conventions are draft 3's (``build_draft3_values.py``): annual means of
  daily Tmax; equal-weight ecoregion mean of facilities per model; mean and
  p25/p75 across models. Monthly PUE/WUE are the stored ``pue_mean`` and
  ``wue_mean`` (all days of that month in the window), the same values
  draft 3's Figure 4 differences were built from.
- Decimals: temperature 1, PUE 3, WUE 2 (absolute values). Raw unrounded
  means are kept in ``mean_raw`` so small month-to-month PUE differences
  can still be separated.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from server_jobs_common import DRAFT3, FAC_DIR, MONTHS, N_MODELS, OUT_DIR, PERIOD_LABEL, compare, eco, facilities, t

CASES = {"Case 1": "ae-chiller", "Case 2": "we-chiller", "Case 5": "chiller-only"}


def job_a(root: Path, fac: pd.DataFrame) -> list[str]:
    a = root / FAC_DIR / "analysis"
    tx = pd.concat([pd.read_csv(a / "facility_annual_tmax.csv"), pd.read_csv(a / "facility_annual_tmax_endcentury.csv")])
    assert tx["gcm"].nunique() == N_MODELS and tx["tasmax_c"].between(-30, 50).all()
    yrs = tx.groupby("period")["year"].agg(["min", "max"])
    print("[A] year ranges:\n", yrs)
    w = tx.groupby(["facility_id", "gcm", "period"])["tasmax_c"].mean().unstack("period").reset_index()
    for p in ("midcentury", "endcentury"):
        w[f"d_{p}"] = w[p] - w["historical"]
    w = w.merge(fac[["facility_id", "abbr"]], on="facility_id")

    rows, mism = [], []
    for col, label, diff in (("historical", "1985–2014", False), ("d_midcentury", "2045–2074", True),
                             ("d_endcentury", "2075–2100", True)):
        e = eco(w, col)
        for abbr in e.index:
            rows.append({"ecoregion": abbr, "quantity": "hist_tmax_C" if col == "historical" else "d_tmax_C",
                         "period": label, "mean": t(e.loc[abbr, "mean"], "temp"), "p25": t(e.loc[abbr, "p25"], "temp"),
                         "p75": t(e.loc[abbr, "p75"], "temp"), "mean_raw": e.loc[abbr, "mean"]})
    out = pd.DataFrame(rows)
    out.to_csv(root / OUT_DIR / "jobA_hist_tmax.csv", index=False)

    ref = pd.read_csv(root / DRAFT3 / "warming.csv", dtype=str).set_index(["ecoregion", "period"])
    ours = out[out["quantity"] == "d_tmax_C"].set_index(["ecoregion", "period"])
    for c in ("mean", "p25", "p75"):
        mism += compare(f"A warming {c}", ours[c], ref[c])
    print(out[out["quantity"] == "hist_tmax_C"].to_string(index=False))
    return mism


def job_b(root: Path, fac: pd.DataFrame) -> list[str]:
    mo = pd.read_csv(root / FAC_DIR / "analysis/facility_month_summary.csv")
    assert mo["gcm"].nunique() == N_MODELS
    assert mo["pue_mean"].between(1.0, 3.0).all() and mo["wue_mean"].between(0.0, 20.0).all()
    mo = mo.merge(fac[["facility_id", "abbr"]], on="facility_id")
    mo = mo[mo["abbr"].isin(["CP", "WV"])]
    mo["case"] = mo["archetype"].map({v: k for k, v in CASES.items()})

    # model-level ecoregion means: (case, abbr, period, gcm, month)
    em = mo.groupby(["case", "abbr", "period", "gcm", "month"])[["pue_mean", "wue_mean"]].mean()
    rows = []
    for (case, abbr, period), g in em.groupby(level=["case", "abbr", "period"]):
        for var, col, kind in (("PUE", "pue_mean", "pue"), ("WUE", "wue_mean", "wue")):
            s = g[col].droplevel(["case", "abbr", "period"])
            q = s.groupby(level="month")
            assert (q.size() == N_MODELS).all()
            for m, mean, p25, p75 in zip(q.mean().index, q.mean(), q.quantile(0.25), q.quantile(0.75)):
                rows.append({"case": case, "ecoregion": abbr, "period": PERIOD_LABEL[period], "month": m,
                             "month_name": MONTHS[m - 1], "variable": var, "mean": t(mean, kind), "p25": t(p25, kind),
                             "p75": t(p75, kind), "mean_raw": mean})
    out = pd.DataFrame(rows).sort_values(["case", "ecoregion", "variable", "period", "month"])
    out.to_csv(root / OUT_DIR / "jobB_monthly_pue_wue.csv", index=False)

    # Validation: mid-century minus historical, per facility and model, as draft 3 did.
    ref = pd.read_csv(root / DRAFT3 / "monthly.csv", dtype=str)
    ref = ref[ref["quantity"].isin(["d_pue_midcentury", "d_wue_midcentury", "hist_pue"])]
    ref = ref.set_index(["region", "quantity", "case", "month"])
    w = mo.pivot_table(index=["case", "abbr", "facility_id", "gcm", "month"], columns="period",
                       values=["pue_mean", "wue_mean"])
    chk = []
    for var, kind in (("pue", "pue"), ("wue", "wue")):
        d = (w[(f"{var}_mean", "midcentury")] - w[(f"{var}_mean", "historical")]).groupby(
            level=["case", "abbr", "gcm", "month"]).mean()
        for (case, abbr, month), s in d.groupby(level=["case", "abbr", "month"]):
            chk.append({"region": abbr, "quantity": f"d_{var}_midcentury", "case": case, "month": str(month),
                        "mean": t(s.mean(), kind, True), "p25": t(s.quantile(0.25), kind, True),
                        "p75": t(s.quantile(0.75), kind, True)})
    chk = pd.DataFrame(chk).set_index(["region", "quantity", "case", "month"])
    chk = chk[chk.index.isin(ref.index)]
    mism = []
    for c in ("mean", "p25", "p75"):
        mism += compare(f"B monthly diff {c}", chk[c], ref[c])
    # hist_pue (Case 1) is also in draft 3: absolute values
    h = out[(out["case"] == "Case 1") & (out["variable"] == "PUE") & (out["period"] == "1985–2014")]
    h = h.assign(region=h["ecoregion"], quantity="hist_pue", month=h["month"].astype(str)).set_index(
        ["region", "quantity", "case", "month"])
    for c in ("mean", "p25", "p75"):
        mism += compare(f"B hist_pue {c}", h[c], ref[c])
    return mism


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    (root / OUT_DIR).mkdir(parents=True, exist_ok=True)
    fac = facilities(root)
    mism = job_a(root, fac) + job_b(root, fac)
    print("\nMISMATCHES:" if mism else "\nAll validations match draft 3.")
    print("\n".join(mism))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Server jobs A and B (draft 4).")
    p.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return p.parse_args()


if __name__ == "__main__":
    run(parse_args())
