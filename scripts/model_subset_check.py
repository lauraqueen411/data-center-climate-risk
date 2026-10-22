"""Compare the 20-model ensemble with the 9-model subset shared with Pierce & Cayan (2025).

Re-aggregates the saved per-model outputs for three model sets: ALL20 (every
``included_in_ensemble == yes`` model in ``outputs/ensemble_pue_wue/gcm_manifest.csv``),
SUB9 (the nine models below) and OTHER11 (ALL20 minus SUB9). Nothing is re-run
and no daily file is read.

Inputs and aggregation are those of ``build_draft3_values.py`` (its ``Data``
class, ecoregion assignment, rounding helper ``t`` and constants are imported):
annual means; change = future minus historical per facility and model;
ecoregion value = equal-weight mean of its facilities within each model; then
mean, p25 and p75 across the models of the set (pandas linear quantiles).

Step 1 validates this code path: the ALL20 values, rounded with draft 3's
``t``, must equal ``docs/results-figures/draft3/{table1,table_cases25,warming,
monthly,model_agreement,rh_sensitivity}.csv`` cell for cell, or the script
stops before writing anything else.

Writes only to ``docs/results-figures/model_subset_check/``:
- ``validation.csv``: every compared cell, draft 3 value, recomputed value.
- ``ecoregion_values_long.csv``: set, case, variable, period, ecoregion, mean, p25, p75 (unrounded).
- ``ecoregion_table_wide.csv``: draft3/table1.csv layout, one block per set and case, plus tmax columns.
- ``rank_order.csv``: ecoregion rank order of each change, per set, period, case.
- ``wue_ratio.csv``: future/historical Case 1-5 WUE ratio per ecoregion and set.
- ``model_agreement_by_set.csv``: draft3/model_agreement.csv extended to all cases,
  exceedance days, both setpoints and the three sets.
- ``pct_diff_vs_all20.csv``: every SUB9 and OTHER11 ensemble mean against ALL20, with flags.
- ``per_model.csv`` and ``per_model_wv_minus_cp_pue.png``.
- ``rh_sensitivity_by_set.csv``: the Figure 5 numbers per set.
- ``run_manifest.csv``.

Assumptions log
----------------
- **ALL20** is the manifest's ``included_in_ensemble == yes`` rows (20 of 27);
  the other seven have no SSP3-7.0 humidity or no SSP3-7.0 downscaling.
- **Percent difference** is (set mean - ALL20 mean) / |ALL20 mean| x 100, on
  unrounded values. Where the ALL20 value is tiny (Case 2 PUE change ~0.0001)
  a large percentage can be a negligible absolute shift; ``abs_diff`` is
  reported beside it so such flags can be judged.
- **Rank** is the ecoregion's position (1 = largest) among the five
  ecoregions for the same set, case, variable and period, ``method="min"``
  (ties share the better rank). ECSF and KM have one facility each.
- **"CP leads acute"** is read as CP having the largest exceedance-day change
  of the five ecoregions; **"WV leads chronic"** as WV having the largest PUE
  (or WUE) change. The WV-vs-CP comparison is reported separately.
- **WUE ratio** is (historical mean + change mean) / historical mean of the
  ensemble means, the reading of "quadruples" and "doubles" from Table 1.
- **Mean tmax change across facilities** (per-model CSV) is the equal-weight
  mean over all 109 facilities, not over ecoregions.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker  # noqa: F401
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_draft3_values import ABBR, CASES, MONTHS, NAME, ORDER, PERIOD_LABEL, Data, t  # noqa: E402

from climate_risk_dc import figstyle  # noqa: E402

SUB9 = ["ACCESS-CM2", "CESM2-LENS", "CNRM-ESM2-1", "EC-Earth3-Veg", "FGOALS-g3",
        "GFDL-ESM4", "MPI-ESM1-2-HR", "MRI-ESM2-0", "TaiESM1"]
SET_ORDER = ["ALL20", "SUB9", "OTHER11"]
FUTURE = list(PERIOD_LABEL)  # midcentury, endcentury
KIND = {"pue": "pue", "wue": "wue", "exceedance_days": "days", "tmax": "temp", "wetbulb": "temp"}
SET_COLOR = {"SUB9": figstyle.BLUE, "OTHER11": figstyle.ORANGE}
SET_MARKER = {"SUB9": "o", "OTHER11": "^"}


# ---------------------------------------------------------------------------
# Aggregation (one path for every set)
# ---------------------------------------------------------------------------


def eco_stats(pm: pd.DataFrame, col: str, models: list[str], keys: tuple[str, ...] = ()) -> pd.DataFrame:
    """Ecoregion mean of facilities within each model, then mean and IQR across ``models``."""
    pm = pm[pm["gcm"].isin(models)]
    m = pm.groupby(list(keys) + ["abbr", "gcm"])[col].mean()
    g = m.groupby(level=list(keys) + ["abbr"])
    out = pd.DataFrame({"mean": g.mean(), "p25": g.quantile(0.25), "p75": g.quantile(0.75), "n_models": g.size()})
    assert (out["n_models"] == len(models)).all(), f"{col}: every ecoregion must use all {len(models)} models"
    return out


def by_model(pm: pd.DataFrame, col: str, models: list[str]) -> pd.DataFrame:
    """gcm x ecoregion table of the ecoregion mean."""
    return pm[pm["gcm"].isin(models)].groupby(["gcm", "abbr"])[col].mean().unstack("abbr")


def tmax_per_model(d: Data) -> pd.DataFrame:
    tx = d.tmax.groupby(["facility_id", "gcm", "period"])["tasmax_c"].mean().unstack("period").reset_index()
    for p in FUTURE:
        tx[f"d_{p}"] = tx[p] - tx["historical"]
    return tx.merge(d.fac[["facility_id", "abbr"]], on="facility_id")


def sources(d: Data) -> list[tuple[str, str, pd.DataFrame]]:
    """(case, variable, per-model frame) for every quantity in the ecoregion tables."""
    out = [("all", "tmax", tmax_per_model(d)),
           ("all", "exceedance_days", d.per_model("exceedance_days", "all")),
           ("all", "wetbulb", d.per_model("wetbulb_annual", "ae-chiller"))]
    for case, arch in CASES.items():
        out += [(case, "pue", d.per_model("pue", arch)), (case, "wue", d.per_model("wue", arch))]
    return out


def long_table(d: Data, sets: dict[str, list[str]]) -> pd.DataFrame:
    n_fac = d.fac["abbr"].value_counts()
    rows = []
    for set_name, models in sets.items():
        for case, var, pm in sources(d):
            cols = ["historical"] + ([f"d_{p}" for p in FUTURE] if var != "wetbulb" else [])
            for col in cols:
                s = eco_stats(pm, col, models)
                for abbr in ORDER:
                    rows.append({"set": set_name, "case": case, "variable": var,
                                 "quantity": "historical" if col == "historical" else "change",
                                 "period": col.removeprefix("d_"), "ecoregion": abbr,
                                 "n_facilities": int(n_fac[abbr]), "n_models": int(s.loc[abbr, "n_models"]),
                                 "mean": s.loc[abbr, "mean"], "p25": s.loc[abbr, "p25"], "p75": s.loc[abbr, "p75"]})
    return pd.DataFrame(rows)


def lookup(long: pd.DataFrame) -> dict:
    return {(r.set, r.case, r.variable, r.period, r.ecoregion): r for r in long.itertuples()}


# ---------------------------------------------------------------------------
# Wide table in the draft3/table1.csv layout
# ---------------------------------------------------------------------------


def wide_table(long: pd.DataFrame, hist_iqr: bool = False, sets=SET_ORDER, cases=tuple(CASES)) -> pd.DataFrame:
    L = lookup(long)
    rows = []
    for s in sets:
        for case in cases:
            for abbr in ORDER:
                for p, label in PERIOD_LABEL.items():
                    g = lambda case_, var, per: L[(s, case_, var, per, abbr)]  # noqa: E731
                    r = {"set": s, "case": case, "ecoregion": NAME[abbr], "n_facilities": g("all", "tmax", "historical").n_facilities,
                         "hist_exceedance_days": t(g("all", "exceedance_days", "historical").mean, "days"),
                         "hist_wetbulb_C": t(g("all", "wetbulb", "historical").mean, "temp"),
                         "hist_pue": t(g(case, "pue", "historical").mean, "pue"),
                         "hist_wue": t(g(case, "wue", "historical").mean, "wue")}
                    if hist_iqr:
                        for k in ("pue", "wue"):
                            h = g(case, k, "historical")
                            r[f"hist_{k}_p25"], r[f"hist_{k}_p75"] = t(h.p25, k), t(h.p75, k)
                    r["period"] = label
                    for col, var, cs in (("d_exceedance", "exceedance_days", "all"), ("d_pue", "pue", case),
                                         ("d_wue", "wue", case), ("d_tmax_C", "tmax", "all")):
                        v = g(cs, var, p)
                        kind = KIND[var]
                        diff = kind != "temp"
                        r[col], r[f"{col}_p25"], r[f"{col}_p75"] = (t(x, kind, diff) for x in (v.mean, v.p25, v.p75))
                    r["hist_tmax_C"] = t(g("all", "tmax", "historical").mean, "temp")
                    r["n_models"] = g(case, "pue", p).n_models
                    rows.append(r)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Step 1: validation against draft 3
# ---------------------------------------------------------------------------


def warming_table(d: Data, models: list[str]) -> pd.DataFrame:
    tx = tmax_per_model(d)
    tx = tx[tx["gcm"].isin(models)]
    fac = tx.groupby(["facility_id", "abbr"])[[f"d_{p}" for p in FUTURE]].mean().reset_index()
    rows = []
    for p, label in PERIOD_LABEL.items():
        e = eco_stats(tx, f"d_{p}", models)
        for abbr in ORDER:
            s = fac[fac["abbr"] == abbr][f"d_{p}"]
            rows.append({"ecoregion": abbr, "period": label, "mean": t(e.loc[abbr, "mean"], "temp"),
                         "p25": t(e.loc[abbr, "p25"], "temp"), "p75": t(e.loc[abbr, "p75"], "temp"),
                         "facility_min": t(s.min(), "temp"), "facility_max": t(s.max(), "temp")})
    return pd.DataFrame(rows)


def monthly_table(d: Data, models: list[str]) -> pd.DataFrame:
    mo = d.month[d.month["gcm"].isin(models)].merge(d.fac[["facility_id", "abbr"]], on="facility_id")
    mo["case"] = mo["archetype"].map({v: k for k, v in CASES.items()})
    rows = []

    def add(series: pd.Series, region: str, quantity: str, case: str, kind: str, diff: bool) -> None:
        for m in range(1, 13):
            v = series.xs(m, level="month")
            assert len(v) == len(models)
            rows.append({"region": region, "quantity": quantity, "case": case, "month": m, "month_name": MONTHS[m - 1],
                         "mean": t(v.mean(), kind, diff), "p25": t(v.quantile(0.25), kind, diff),
                         "p75": t(v.quantile(0.75), kind, diff)})

    for region in ("CP", "WV"):
        r = mo[mo["abbr"] == region]
        for case in CASES:
            w = r[r["case"] == case].pivot_table(index=["facility_id", "gcm", "month"], columns="period",
                                                  values=["pue_mean", "wue_mean"])
            per_model = lambda s: s.groupby(level=["gcm", "month"]).mean()  # noqa: E731
            add(per_model(w[("pue_mean", "midcentury")] - w[("pue_mean", "historical")]), region, "d_pue_midcentury", case, "pue", True)
            if case == "Case 1":
                add(per_model(w[("wue_mean", "midcentury")] - w[("wue_mean", "historical")]), region, "d_wue_midcentury", case, "wue", True)
                add(per_model(w[("pue_mean", "historical")]), region, "hist_pue", case, "pue", False)
        mt = d.mtmax[d.mtmax["gcm"].isin(models)].merge(d.fac[["facility_id", "abbr"]], on="facility_id")
        mt = mt[(mt["abbr"] == region) & (mt["period"] == "historical")].groupby(["gcm", "month"])["tasmax_c"].mean()
        add(mt, region, "hist_tmax_C", "all", "temp", False)
    return pd.DataFrame(rows)


PAIRS = {"pue": (("WV", "CP"), ("WV", "KM"), ("CP", "KM")), "wue": (("WV", "CP"), ("WV", "KM"), ("CP", "KM")),
         "exceedance_days": (("CP", "WV"), ("WV", "CP"))}


def agreement_table(d: Data, sets: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for set_name, models in sets.items():
        for run in ("77.5%", "60%"):
            blocks = [("Case 1", m, d.per_model(m, CASES["Case 1"], run)) for m in ("pue", "wue")]
            if run == "77.5%":
                blocks += [(c, m, d.per_model(m, CASES[c])) for c in ("Case 2", "Case 5") for m in ("pue", "wue")]
                blocks.append(("all", "exceedance_days", d.per_model("exceedance_days", "all")))
            for case, metric, pm in blocks:
                for p, label in PERIOD_LABEL.items():
                    m = by_model(pm, f"d_{p}", models)
                    for a, b in PAIRS[metric]:
                        rows.append({"set": set_name, "setpoint": run if metric != "exceedance_days" else "n/a",
                                     "case": case, "metric": metric.upper() if metric != "exceedance_days" else "exceedance_days",
                                     "period": label, "comparison": f"{a} > {b}",
                                     "n_models": int((m[a] > m[b]).sum()), "of": len(models)})
    return pd.DataFrame(rows)


def rh_table(d: Data, models: list[str], set_name: str | None = None) -> pd.DataFrame:
    raw = {}
    for run in ("60%", "77.5%"):
        for metric in ("pue", "wue"):
            pm = d.per_model(metric, "ae-chiller", run)
            for col in ("historical", "d_midcentury", "d_endcentury"):
                raw[(run, metric, col)] = eco_stats(pm, col, models)
    rows = []
    for abbr in ORDER:
        for run in ("60%", "77.5%", "difference (60% − 77.5%)"):
            r = ({"set": set_name} if set_name else {}) | {"ecoregion": abbr, "setpoint": run}
            for metric in ("pue", "wue"):
                for col in ("historical", "d_midcentury", "d_endcentury"):
                    if run.startswith("difference"):
                        v = raw[("60%", metric, col)].loc[abbr, "mean"] - raw[("77.5%", metric, col)].loc[abbr, "mean"]
                        r[f"{metric}_{col}"] = t(v, metric, True)
                    else:
                        s = raw[(run, metric, col)].loc[abbr]
                        r[f"{metric}_{col}"] = t(s["mean"], metric, col != "historical")
                        if col != "historical":
                            r[f"{metric}_{col}_p25"], r[f"{metric}_{col}_p75"] = t(s["p25"], metric, True), t(s["p75"], metric, True)
            rows.append(r)
    return pd.DataFrame(rows), raw


def compare(name: str, ours: pd.DataFrame, theirs: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    ours = ours.astype(object).where(ours.notna(), "").astype(str)
    assert list(ours.columns) == list(theirs.columns), f"{name}: columns differ: {list(ours.columns)} vs {list(theirs.columns)}"
    assert len(ours) == len(theirs), f"{name}: {len(ours)} rows vs {len(theirs)}"
    rows = []
    for (_, a), (_, b) in zip(ours.iterrows(), theirs.iterrows()):
        for c in ours.columns:
            if c in keys:
                continue
            rows.append({"file": name, "row": " | ".join(str(b[k]) for k in keys), "column": c,
                         "draft3": b[c], "recomputed": a[c], "match": a[c] == b[c]})
    return pd.DataFrame(rows)


def validate(d: Data, long: pd.DataFrame, all20: list[str], draft3: Path) -> pd.DataFrame:
    read = lambda f: pd.read_csv(draft3 / f, dtype=str, keep_default_na=False)  # noqa: E731
    cols1 = list(read("table1.csv").columns)
    t1 = wide_table(long, sets=["ALL20"], cases=["Case 1"])[cols1]
    cols25 = list(read("table_cases25.csv").columns)
    t25 = wide_table(long, hist_iqr=True, sets=["ALL20"], cases=["Case 2", "Case 5"])[cols25]
    ag = agreement_table(d, {"ALL20": all20})
    ag = ag[(ag["case"] == "Case 1") & (ag["setpoint"] == "77.5%")][["setpoint", "metric", "period", "comparison", "n_models", "of"]]
    rh, _ = rh_table(d, all20)
    return pd.concat([
        compare("table1.csv", t1, read("table1.csv"), ["ecoregion", "period"]),
        compare("table_cases25.csv", t25, read("table_cases25.csv"), ["case", "ecoregion", "period"]),
        compare("warming.csv", warming_table(d, all20), read("warming.csv"), ["ecoregion", "period"]),
        compare("monthly.csv", monthly_table(d, all20), read("monthly.csv"), ["region", "quantity", "case", "month", "month_name"]),
        compare("model_agreement.csv", ag.reset_index(drop=True), read("model_agreement.csv"), ["setpoint", "metric", "period", "comparison"]),
        compare("rh_sensitivity.csv", rh, read("rh_sensitivity.csv"), ["ecoregion", "setpoint"]),
    ], ignore_index=True)


# ---------------------------------------------------------------------------
# Step 3: conclusion checks
# ---------------------------------------------------------------------------


def with_ranks(long: pd.DataFrame) -> pd.DataFrame:
    keys = ["set", "case", "variable", "quantity", "period"]
    out = long.copy()
    out["rank"] = out.groupby(keys)["mean"].rank(ascending=False, method="min").astype(int)
    return out


def rank_order(long: pd.DataFrame) -> pd.DataFrame:
    rk = with_ranks(long)
    rk = rk[(rk["quantity"] == "change") & rk["variable"].isin(["exceedance_days", "pue", "wue", "tmax"])]
    rows = []
    for (s, case, var, p), g in rk.groupby(["set", "case", "variable", "period"], sort=False):
        g = g.sort_values(["rank", "ecoregion"])
        v = g.set_index("ecoregion")["mean"]
        rows.append({"set": s, "period": PERIOD_LABEL[p], "case": case, "variable": var,
                     "rank_order": " > ".join(g["ecoregion"]), "leader": g["ecoregion"].iloc[0],
                     "CP_rank": int(g.set_index("ecoregion").loc["CP", "rank"]),
                     "WV_rank": int(g.set_index("ecoregion").loc["WV", "rank"]),
                     "larger_of_CP_WV": "CP" if v["CP"] > v["WV"] else "WV"})
    return pd.DataFrame(rows)


def wue_ratio(long: pd.DataFrame) -> pd.DataFrame:
    L = lookup(long)
    rows = []
    for s in SET_ORDER:
        for case in CASES:
            for abbr in ORDER:
                h = L[(s, case, "wue", "historical", abbr)].mean
                r = {"set": s, "case": case, "ecoregion": abbr, "hist_wue": h}
                for p in FUTURE:
                    r[f"ratio_{p}"] = (h + L[(s, case, "wue", p, abbr)].mean) / h
                rows.append(r)
    return pd.DataFrame(rows)


def pct_diff(long: pd.DataFrame) -> pd.DataFrame:
    rk = with_ranks(long)
    keys = ["case", "variable", "quantity", "period", "ecoregion"]
    base = rk[rk["set"] == "ALL20"][keys + ["mean", "rank"]].rename(columns={"mean": "all20_mean", "rank": "all20_rank"})
    other = rk[rk["set"] != "ALL20"][["set"] + keys + ["mean", "rank"]].rename(columns={"mean": "set_mean", "rank": "set_rank"})
    out = other.merge(base, on=keys)
    out["abs_diff"] = out["set_mean"] - out["all20_mean"]
    out["pct_diff"] = 100 * out["abs_diff"] / out["all20_mean"].abs()
    out["flag_gt10pct"] = out["pct_diff"].abs() > 10
    out["flag_rank_change"] = out["set_rank"] != out["all20_rank"]
    return out[["set"] + keys + ["all20_mean", "set_mean", "abs_diff", "pct_diff", "all20_rank", "set_rank",
                                 "flag_gt10pct", "flag_rank_change"]]


# ---------------------------------------------------------------------------
# Step 4: per-model values
# ---------------------------------------------------------------------------


def per_model_table(d: Data, all20: list[str]) -> pd.DataFrame:
    out = pd.DataFrame(index=pd.Index(sorted(all20), name="model"))
    member = d.manifest.set_index("model")["member"]
    out["member"] = member.reindex(out.index)
    out["in_SUB9"] = np.where(out.index.isin(SUB9), "yes", "no")
    pue_pm, ex_pm, tx = d.per_model("pue", CASES["Case 1"]), d.per_model("exceedance_days", "all"), tmax_per_model(d)
    for p in FUTURE:
        pue = by_model(pue_pm, f"d_{p}", all20)
        ex = by_model(ex_pm, f"d_{p}", all20)
        out[f"wv_minus_cp_d_pue_case1_{p}"] = pue["WV"] - pue["CP"]
        out[f"cp_minus_wv_d_exceedance_days_{p}"] = ex["CP"] - ex["WV"]
        out[f"mean_d_tmax_C_all_facilities_{p}"] = tx[tx["gcm"].isin(all20)].groupby("gcm")[f"d_{p}"].mean()
    return out.reset_index()


def plot_per_model(pm: pd.DataFrame, path: Path) -> None:
    figstyle.apply()
    pm = pm.sort_values("wv_minus_cp_d_pue_case1_midcentury").reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(figstyle.FIG_WIDTH_IN, 5.2), sharey=True)
    y = np.arange(len(pm))
    for ax, p in zip(axes, FUTURE):
        col = f"wv_minus_cp_d_pue_case1_{p}"
        ax.axvline(0, color=figstyle.ZERO_LINE_COLOR, lw=figstyle.ZERO_LINE_WIDTH)
        for s, flag in (("SUB9", "yes"), ("OTHER11", "no")):
            m = pm["in_SUB9"] == flag
            ax.scatter(pm.loc[m, col], y[m], s=36, marker=SET_MARKER[s], color=SET_COLOR[s],
                       edgecolor="white", linewidth=0.8, zorder=3, label=f"{s} model")
            ax.axvline(pm.loc[m, col].mean(), color=SET_COLOR[s], lw=1, ls="--", zorder=2, label=f"{s} mean")
        ax.set_title(f"{PERIOD_LABEL[p]} minus 1985–2014")
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    axes[0].set_yticks(y, pm["model"])
    axes[0].set_ylim(-0.7, len(pm) - 0.3)
    fig.supxlabel("Willamette Valley minus Columbia Plateau change in Case 1 PUE", fontsize=figstyle.LABEL_PT)
    handles, labels = axes[0].get_legend_handles_labels()
    order = [labels.index(k) for k in ("SUB9 model", "SUB9 mean", "OTHER11 model", "OTHER11 mean")]
    fig.legend([handles[i] for i in order], [labels[i] for i in order], loc="upper center", ncol=4,
               frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------------


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out = root / "docs/results-figures/model_subset_check"
    draft3 = root / "docs/results-figures/draft3"
    d = Data(root)
    man = d.manifest
    all20 = sorted(man.loc[man["included_in_ensemble"] == "yes", "model"])
    assert len(all20) == 20 and set(SUB9) <= set(all20), "SUB9 must be a subset of the 20 ensemble models"
    sets = {"ALL20": all20, "SUB9": sorted(SUB9), "OTHER11": sorted(set(all20) - set(SUB9))}
    assert len(sets["OTHER11"]) == 11

    # Step 1
    long = long_table(d, sets)
    val = validate(d, long, all20, draft3)
    out.mkdir(parents=True, exist_ok=True)
    val.to_csv(out / "validation.csv", index=False)
    bad = val[~val["match"]]
    print(f"Validation: {len(val)} cells compared, {len(bad)} mismatches")
    if len(bad):
        print(bad.to_string(index=False))
        sys.exit("Stopping: ALL20 does not reproduce draft 3 (see validation.csv).")

    # Steps 2-5
    long.to_csv(out / "ecoregion_values_long.csv", index=False)
    wide_table(long).to_csv(out / "ecoregion_table_wide.csv", index=False)
    ro = rank_order(long)
    ro.to_csv(out / "rank_order.csv", index=False)
    wr = wue_ratio(long)
    wr.to_csv(out / "wue_ratio.csv", index=False)
    ag = agreement_table(d, sets)
    ag.to_csv(out / "model_agreement_by_set.csv", index=False)
    pdiff = pct_diff(long)
    pdiff.to_csv(out / "pct_diff_vs_all20.csv", index=False)
    pm = per_model_table(d, all20)
    pm.to_csv(out / "per_model.csv", index=False)
    plot_per_model(pm, out / "per_model_wv_minus_cp_pue.png")
    rh = pd.concat([rh_table(d, m, s)[0] for s, m in sets.items()], ignore_index=True)
    rh.to_csv(out / "rh_sensitivity_by_set.csv", index=False)
    pd.DataFrame([{"script": "scripts/model_subset_check.py", "sets": "; ".join(f"{k}={','.join(v)}" for k, v in sets.items()),
                   "inputs": "build_draft3_values.Data (outputs/ensemble_pue_wue_facilities/analysis, superseded_rh_up60 cache, gcm_manifest)",
                   "validated_against": "draft3/table1, table_cases25, warming, monthly, model_agreement, rh_sensitivity",
                   "validation_cells": len(val), "validation_mismatches": len(bad)}]).to_csv(out / "run_manifest.csv", index=False)
    print(f"Wrote {out}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare the 20-model ensemble with the 9-model Pierce & Cayan subset.")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
