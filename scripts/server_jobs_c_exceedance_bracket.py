"""Server job C (draft 4): exceedance-day threshold bracket at 32, 35 and 40 °C.

Reads only ``facility_id``, ``date`` and ``tasmax_c`` from the Case 1
(``ae-chiller``) daily ``_fast`` files; tasmax is identical across cases.
Annual count of days with Tmax above the threshold, per facility and year,
then the mean over the years in the window; change = future minus
historical per facility and model; equal-weight ecoregion mean per model;
mean and p25/p75 across the 20 models.

Validation: the 35 °C values must reproduce ``draft3/table1.csv``
(``hist_exceedance_days``, ``d_exceedance`` and its p25/p75).

Writes ``docs/results-figures/server_jobs_oct23/jobC_exceedance_bracket.csv``
(threshold, ecoregion, quantity, period, mean, p25, p75) and
``jobC_contrast.csv`` (CP minus WV and ecoregion ranks per threshold).

Assumptions log
----------------
- **Strict exceedance** (Tmax > threshold), as in
  ``heat_indices.annual_tmax_threshold_days``, which produced the draft 3
  values (via the stored heat-ensemble outputs).
- Draft 3's exceedance days came from ``outputs/month2_heat_ensemble``
  (nearest LOCA2 cell from per-model GeoTIFFs); here they are recounted from
  the daily ``tasmax_c`` the PUE/WUE files store. Any 35 °C mismatch with
  draft 3 is reported, not adjusted.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from server_jobs_common import DRAFT3, N_MODELS, NAME, ORDER, OUT_DIR, PERIOD_LABEL, PERIODS, compare, daily_path, \
    eco, facilities, models, t

THRESHOLDS = (32.0, 35.0, 40.0)


def _read(path: Path) -> pd.DataFrame:
    """Period mean of annual exceedance counts per facility, one column per threshold."""
    df = pd.read_csv(path, usecols=["facility_id", "date", "tasmax_c"])
    assert df["tasmax_c"].between(-40, 60).all(), f"{path.name}: tasmax outside -40..60 °C"
    df["year"] = df["date"].str.slice(0, 4).astype(int)
    for th in THRESHOLDS:
        df[th] = df["tasmax_c"] > th
    annual = df.groupby(["facility_id", "year"])[list(THRESHOLDS)].sum()
    nyr = annual.index.get_level_values("year").nunique()
    return annual.groupby(level="facility_id").mean().assign(n_years=nyr)


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    fac = facilities(root)
    gcms = models(root)
    jobs = [(g, p) for g in gcms for p in PERIODS]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        parts = list(pool.map(_read, [daily_path(root, g, "ae-chiller", p) for g, p in jobs]))
    rows = []
    for (g, p), df in zip(jobs, parts):
        s, e = PERIODS[p]
        assert df["n_years"].iloc[0] == e - s + 1, (g, p)
        rows.append(df.drop(columns="n_years").reset_index().melt("facility_id", var_name="threshold").assign(gcm=g, period=p))
    long = pd.concat(rows)
    w = long.pivot_table(index=["threshold", "facility_id", "gcm"], columns="period", values="value").reset_index()
    for p in ("midcentury", "endcentury"):
        w[f"d_{p}"] = w[p] - w["historical"]
    w = w.merge(fac[["facility_id", "abbr"]], on="facility_id")

    out, contrast, raw = [], [], {}
    for th in THRESHOLDS:
        sub = w[w["threshold"] == th]
        for col, qty, period, diff in (("historical", "hist_days", "historical", False),
                                       ("d_midcentury", "d_days", "midcentury", True),
                                       ("d_endcentury", "d_days", "endcentury", True)):
            e = eco(sub, col)
            raw[(th, col)] = e["mean"]
            for a in ORDER:
                out.append({"threshold_C": th, "ecoregion": a, "quantity": qty, "period": PERIOD_LABEL[period],
                            "mean": t(e.loc[a, "mean"], "days", diff), "p25": t(e.loc[a, "p25"], "days", diff),
                            "p75": t(e.loc[a, "p75"], "days", diff), "mean_raw": e.loc[a, "mean"]})
            per_model = sub.groupby(["gcm", "abbr"])[col].mean().unstack("abbr")
            rank = e["mean"].rank(ascending=False).astype(int)
            contrast.append({"threshold_C": th, "quantity": col, "CP_minus_WV": round(e.loc["CP", "mean"] - e.loc["WV", "mean"], 1),
                             "n_models_CP_gt_WV": int((per_model["CP"] > per_model["WV"]).sum()), "of": N_MODELS,
                             "ranking_high_to_low": " > ".join(e["mean"].sort_values(ascending=False).index),
                             **{f"rank_{a}": rank[a] for a in ORDER}})
    out = pd.DataFrame(out)
    contrast = pd.DataFrame(contrast)
    out.to_csv(root / OUT_DIR / "jobC_exceedance_bracket.csv", index=False)
    contrast.to_csv(root / OUT_DIR / "jobC_contrast.csv", index=False)

    # Validate 35 °C against draft 3 Table 1
    ref = pd.read_csv(root / DRAFT3 / "table1.csv", dtype=str)
    ref["abbr"] = ref["ecoregion"].map({v: k for k, v in NAME.items()})
    o35 = out[out["threshold_C"] == 35.0]
    mism = []
    h = o35[o35["quantity"] == "hist_days"].set_index("ecoregion")["mean"]
    mism += compare("C hist 35", h, ref.drop_duplicates("abbr").set_index("abbr")["hist_exceedance_days"])
    d = o35[o35["quantity"] == "d_days"].set_index(["ecoregion", "period"])
    r = ref.set_index(["abbr", "period"])
    for c, rc in (("mean", "d_exceedance"), ("p25", "d_exceedance_p25"), ("p75", "d_exceedance_p75")):
        mism += compare(f"C d35 {c}", d[c], r[rc])
    # Cross-check against the daily recount draft 3 stored in the Figure 2 cache
    cache = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/fig_show_everything/model_values_cache.csv")
    chk = cache[cache["metric"] == "exceedance_days_check"].set_index(["facility_id", "gcm", "period"])["value"]
    ours = long[long["threshold"] == 35.0].set_index(["facility_id", "gcm", "period"])["value"]
    print(f"[validate] C 35 °C recount vs cached exceedance_days_check: max |diff| = {(ours - chk).abs().max():.2e}")
    stored = cache[cache["metric"] == "exceedance_days"].set_index(["facility_id", "gcm", "period"])["value"]
    print(f"[info] C 35 °C recount vs stored heat-ensemble exceedance_days: max |diff| = {(ours - stored).abs().max():.3f}")

    pd.set_option("display.width", 200)
    print(out.pivot_table(index=["threshold_C", "quantity", "period"], columns="ecoregion", values="mean", aggfunc="first")[ORDER])
    print(contrast.to_string(index=False))
    print("\nMISMATCHES:\n" + "\n".join(mism) if mism else "\nAll validations match draft 3.")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Server job C: exceedance-day threshold bracket.")
    p.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    p.add_argument("--workers", type=int, default=12)
    return p.parse_args()


if __name__ == "__main__":
    run(parse_args())
