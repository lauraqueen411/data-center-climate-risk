"""Supporting tables for the heat chapter, built from saved outputs only.

Writes to ``outputs/chapter_tables/``: one tidy CSV and one markdown file per
table (Tables 1, 2, 3, 4, 6 and 7), ``ensemble_facts.md``, ``README.md`` and
``run_manifest.csv``. Nothing is re-run; no daily files are read.

Conventions (``docs/results-figures/methodology.md``): annual means; change =
future minus historical per facility and model; ecoregion mean of facilities
(equal weight) within each model; ensemble mean and interquartile range (25th
to 75th percentile) across the 20 models. Facility ranges use each facility's
ensemble-mean value.

Tidy columns: ``ecoregion, case, metric, period, statistic, value,
n_facilities, n_models``. Table 6 adds ``month`` and Table 7 adds ``setpoint``.

Inputs: the same local files as ``build_draft3_values.py`` (its ``Data``
class is reused), all under ``outputs/ensemble_pue_wue_facilities/``.

Assumptions log
----------------
- **Skipped, not approximated:** Table 5 (humidity and wet-bulb change), Table 8
  (cooling-mode shares) and daily minimum temperature. None is in the saved
  outputs; all need the daily files on the server.
- **Decimals:** PUE 3 and WUE 2 as the task specifies, but 4 and 3 for
  differences smaller than 0.01 (PUE) and 0.1 (WUE), the rule
  ``build_draft3_values.py`` uses. Monthly changes and the setpoint difference
  would otherwise round to zero. Exceedance days, temperature 1; shares 3.
- **August share of overhead** is day-weighted with fixed month lengths
  (February 28.25). The saved monthly means do not carry day counts, and the
  models' calendars differ, so this is nominal weighting. It is computed on
  the ecoregion-mean monthly PUE of each model, then averaged across models.
- **"All facilities"** rows are the equal-weight mean of the 109 facilities
  within each model.
- **Table 4 strict inequalities:** a model counts as "positive" only if the
  ecoregion-mean change is above zero, and as "WV > CP" only if strictly greater.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from build_draft3_values import CASES, MONTHS, N_MODELS, NAME, ORDER, Data

PERIODS = {"historical": "1985-2014", "midcentury": "2045-2074", "endcentury": "2075-2100"}
FUTURE = ("midcentury", "endcentury")
ALL = "All facilities"
DAYS_IN_MONTH = pd.Series([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], index=range(1, 13))
TIDY = ["ecoregion", "case", "metric", "period", "statistic", "value", "n_facilities", "n_models"]


@dataclass(frozen=True)
class QuotedValue:
    """A number quoted in the chapter text, checked against the computed one."""

    label: str
    quoted: str
    decimals: int


QUOTED = {
    "wv_aug": QuotedValue("Willamette Valley August Case 1 PUE change, mid-century", "+0.057", 3),
    "cp_sep": QuotedValue("Columbia Plateau September Case 1 PUE change, mid-century", "+0.040", 3),
    "sp_775": QuotedValue("Willamette Valley mid-century Case 1 PUE increase, lower quoted value", "0.015", 3),
    "sp_60": QuotedValue("Willamette Valley mid-century Case 1 PUE increase, higher quoted value", "0.021", 3),
    "sp_hist": QuotedValue("Historical Case 1 PUE difference between setpoints", "+0.0065", 4),
    "aug_share": QuotedValue("August share of historical annual overhead", "more than one third", 3),
    "warm_lo": QuotedValue("Mid-century tasmax warming, lowest ecoregion (°C)", "2.4", 1),
    "warm_hi": QuotedValue("Mid-century tasmax warming, highest ecoregion (°C)", "2.8", 1),
}


# ---------------------------------------------------------------------------
# Formatting and aggregation
# ---------------------------------------------------------------------------


def fmt(x: float, kind: str, diff: bool = False) -> str:
    """Fixed-decimal string, so trailing zeros survive in CSV and markdown."""
    if pd.isna(x):
        return ""
    x = float(x)
    d = {"pue": 3, "wue": 2, "days": 1, "temp": 1, "share": 3}[kind]
    if diff and kind == "pue" and abs(x) < 0.01:
        d = 4
    if diff and kind == "wue" and abs(x) < 0.1:
        d = 3
    return f"{x:+.{d}f}" if diff else f"{x:.{d}f}"


def by_model(df: pd.DataFrame, col: str, extra: tuple[str, ...] = ()) -> pd.Series:
    """Ecoregion mean of facilities within each model, plus an all-facility mean."""
    keys = list(extra) + ["abbr", "gcm"]
    eco = df.groupby(keys)[col].mean()
    everything = df.assign(abbr=ALL).groupby(keys)[col].mean()
    return pd.concat([eco, everything])


def stats(series_by_model: pd.Series) -> pd.DataFrame:
    """Mean and interquartile range across models, for every index level except ``gcm``."""
    levels = [n for n in series_by_model.index.names if n != "gcm"]
    g = series_by_model.groupby(level=levels)
    out = pd.DataFrame({"mean": g.mean(), "p25": g.quantile(0.25), "p75": g.quantile(0.75), "n_models": g.size()})
    assert (out["n_models"] == N_MODELS).all(), "every statistic must use all 20 models"
    return out


def facility_range(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """Min and max across facilities of each facility's ensemble-mean value."""
    fac = df.groupby(["facility_id", "abbr"])[col].mean().reset_index()
    return fac.groupby("abbr")[col].agg(["min", "max"])


class Rows:
    def __init__(self, n_fac: pd.Series) -> None:
        self.rows: list[dict] = []
        self.n_fac = n_fac

    def add(self, abbr: str, case: str, metric: str, period: str, statistic: str, value: str,
            n_models: int | str = N_MODELS, **extra) -> None:
        name = NAME.get(abbr, abbr)
        n = self.n_fac.get(abbr, "")
        self.rows.append({"ecoregion": name, "case": case, "metric": metric, "period": PERIODS.get(period, period),
                          "statistic": statistic, "value": value, "n_facilities": n, "n_models": n_models, **extra})

    def add_stats(self, s: pd.DataFrame, case: str, metric: str, period: str, kind: str, diff: bool,
                  regions: list[str], **extra) -> None:
        for abbr in regions:
            for stat in ("mean", "p25", "p75"):
                self.add(abbr, case, metric, period, stat, fmt(s.loc[abbr, stat], kind, diff), **extra)

    def frame(self, extra_cols: tuple[str, ...] = ()) -> pd.DataFrame:
        return pd.DataFrame(self.rows)[TIDY + list(extra_cols)]


def wide(df: pd.DataFrame, index: list[str]) -> pd.DataFrame:
    """Human-readable layout: one row per ``index``, one column per period and statistic."""
    d = df.copy()
    d["col"] = d["period"] + " " + d["statistic"]
    order = list(dict.fromkeys(d["col"]))
    w = d.pivot_table(index=index, columns="col", values="value", aggfunc="first", sort=False)
    return w[order].reset_index().fillna("")


def md_table(df: pd.DataFrame) -> str:
    cols = [str(c) for c in df.columns]
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------


def table1(d: Data, n_fac: pd.Series, found: dict) -> pd.DataFrame:
    """Change in annual mean daily maximum temperature."""
    tx = d.tmax.groupby(["facility_id", "gcm", "period"])["tasmax_c"].mean().unstack("period").reset_index()
    tx = tx.merge(d.fac[["facility_id", "abbr"]], on="facility_id")
    r = Rows(n_fac)
    for p in FUTURE:
        tx[f"d_{p}"] = tx[p] - tx["historical"]
        s = stats(by_model(tx, f"d_{p}"))
        r.add_stats(s, "all", "tasmax_change_C", p, "temp", False, ORDER + [ALL])
        lo, hi = s.loc[ORDER, "mean"].min(), s.loc[ORDER, "mean"].max()
        r.add("Range across the 5 ecoregions", "all", "tasmax_change_C", p, "ecoregion_min", fmt(lo, "temp"))
        r.add("Range across the 5 ecoregions", "all", "tasmax_change_C", p, "ecoregion_max", fmt(hi, "temp"))
        if p == "midcentury":
            found["warm_lo"], found["warm_hi"] = lo, hi
    return r.frame()


def table2(d: Data, n_fac: pd.Series) -> pd.DataFrame:
    """Historical baselines, with the range across facilities."""
    r = Rows(n_fac)

    def block(pm: pd.DataFrame, case: str, metric: str, kind: str, offset: float = 0.0) -> None:
        pm = pm.assign(historical=pm["historical"] - offset)
        s, rng = stats(by_model(pm, "historical")), facility_range(pm, "historical")
        for abbr in ORDER:
            r.add(abbr, case, metric, "historical", "mean", fmt(s.loc[abbr, "mean"], kind))
            r.add(abbr, case, metric, "historical", "facility_min", fmt(rng.loc[abbr, "min"], kind))
            r.add(abbr, case, metric, "historical", "facility_max", fmt(rng.loc[abbr, "max"], kind))

    block(d.per_model("exceedance_days", "all"), "all", "exceedance_days", "days")
    block(d.per_model("wetbulb_annual", "ae-chiller"), "all", "wetbulb_afternoon_C", "temp")
    for case, arch in CASES.items():
        block(d.per_model("pue", arch), case, "pue", "pue")
        block(d.per_model("pue", arch), case, "pue_overhead", "pue", offset=1.0)
        block(d.per_model("wue", arch), case, "wue_L_per_kWh", "wue")
    return r.frame()


def table3(d: Data, n_fac: pd.Series) -> pd.DataFrame:
    """Cases 2 and 5 in the layout of the Case 1 results table."""
    r = Rows(n_fac)
    for case in ("Case 2", "Case 5"):
        for metric, kind, label in (("pue", "pue", "pue"), ("wue", "wue", "wue_L_per_kWh")):
            pm = d.per_model(metric, CASES[case])
            r.add_stats(stats(by_model(pm, "historical")), case, label, "historical", kind, False, ORDER)
            if metric == "pue":
                over = stats(by_model(pm.assign(historical=pm["historical"] - 1.0), "historical"))
                r.add_stats(over, case, "pue_overhead", "historical", kind, False, ORDER)
            for p in FUTURE:
                r.add_stats(stats(by_model(pm, f"d_{p}")), case, f"{label}_change", p, kind, True, ORDER)
    return r.frame()


def table4(d: Data, n_fac: pd.Series) -> pd.DataFrame:
    """Number of models agreeing on the sign and on the WV/CP ordering of the change."""
    r = Rows(n_fac)
    blocks = [(case, metric, d.per_model(metric, arch)) for case, arch in CASES.items() for metric in ("pue", "wue")]
    blocks.append(("all", "exceedance_days", d.per_model("exceedance_days", "all")))
    for case, metric, pm in blocks:
        for p in FUTURE:
            m = by_model(pm, f"d_{p}").unstack("abbr")
            pairs = (("WV", "CP"), ("CP", "WV")) if metric == "exceedance_days" else (("WV", "CP"),)
            for a, b in pairs:
                r.add(f"{a} vs {b}", case, f"{metric}_change", p, f"n_models_{a}_gt_{b}", str(int((m[a] > m[b]).sum())))
            for abbr in ORDER:
                r.add(abbr, case, f"{metric}_change", p, "n_models_positive", str(int((m[abbr] > 0).sum())))
    return r.frame()


def table6(d: Data, n_fac: pd.Series, found: dict) -> tuple[pd.DataFrame, list[str]]:
    """Case 1 monthly PUE: historical value, change, and August's share of annual overhead."""
    mo = d.month[d.month["archetype"] == CASES["Case 1"]]
    w = mo.pivot_table(index=["facility_id", "gcm", "month"], columns="period", values="pue_mean").reset_index()
    w = w.merge(d.fac[["facility_id", "abbr"]], on="facility_id")
    w["overhead"] = w["historical"] - 1.0
    r = Rows(n_fac)
    means = {}
    for col, metric, period, diff in (("historical", "pue", "historical", False), ("overhead", "pue_overhead", "historical", False)):
        s = stats(by_model(w, col, ("month",)))
        for m in range(1, 13):
            r.add_stats(s.xs(m, level="month"), "Case 1", metric, period, "pue", diff, ORDER, month=m)
    for p in FUTURE:
        w[f"d_{p}"] = w[p] - w["historical"]
        s = stats(by_model(w, f"d_{p}", ("month",)))
        means[p] = s["mean"]
        for m in range(1, 13):
            r.add_stats(s.xs(m, level="month"), "Case 1", "pue_change", p, "pue", True, ORDER, month=m)

    # August share: per model, on the ecoregion-mean monthly values, day-weighted.
    def share(col: str) -> pd.DataFrame:
        e = by_model(w, col, ("month",)).unstack("month") * DAYS_IN_MONTH
        return stats(e[8] / e.sum(axis=1))

    shares = {"hist_overhead": share("overhead"), "d_midcentury": share("d_midcentury"), "d_endcentury": share("d_endcentury")}
    for key, metric, period in (("hist_overhead", "august_share_of_annual_overhead", "historical"),
                                ("d_midcentury", "august_share_of_annual_overhead_change", "midcentury"),
                                ("d_endcentury", "august_share_of_annual_overhead_change", "endcentury")):
        r.add_stats(shares[key], "Case 1", metric, period, "share", False, ORDER + [ALL], month=8)

    mid = means["midcentury"]
    found["wv_aug"], found["cp_sep"] = mid[(8, "WV")], mid[(9, "CP")]
    found["aug_share"] = shares["hist_overhead"]["mean"]
    cp = {m: mid[(m, "CP")] for m in (6, 7, 8, 9)}
    cpe = {m: means["endcentury"][(m, "CP")] for m in (6, 7, 8, 9)}
    dip = cp[7] < cp[6] and cp[7] < cp[8]
    notes = [
        f"- Willamette Valley, August, mid-century change: {fmt(found['wv_aug'], 'pue', True)}.",
        f"- Columbia Plateau, September, mid-century change: {fmt(found['cp_sep'], 'pue', True)}.",
        "- Columbia Plateau, mid-century change by month: " + ", ".join(f"{MONTHS[m - 1]} {fmt(v, 'pue', True)}" for m, v in cp.items())
        + f". July is {'below' if dip else 'not below'} both June and August.",
        "- Columbia Plateau, end-of-century change by month: " + ", ".join(f"{MONTHS[m - 1]} {fmt(v, 'pue', True)}" for m, v in cpe.items()) + ".",
        "- August share of historical annual overhead (PUE − 1): "
        + ", ".join(f"{a} {fmt(found['aug_share'][a], 'share')}" for a in ORDER + [ALL]) + ".",
        "- August share of the annual mid-century change in overhead: "
        + ", ".join(f"{a} {fmt(shares['d_midcentury']['mean'][a], 'share')}" for a in ORDER + [ALL]) + ".",
    ]
    found["aug_share_change"] = shares["d_midcentury"]["mean"]
    return r.frame(("month",)), notes


def table7(d: Data, n_fac: pd.Series, found: dict) -> pd.DataFrame:
    """Case 1 under the 60% and 77.5% upper supply-air RH setpoints."""
    r = Rows(n_fac)
    regions = ORDER + [ALL]
    for metric, kind, label in (("pue", "pue", "pue"), ("wue", "wue", "wue_L_per_kWh")):
        pm = {run: d.per_model(metric, CASES["Case 1"], run) for run in ("60%", "77.5%")}
        cols = ["historical", "d_midcentury", "d_endcentury"]
        both = pm["60%"].merge(pm["77.5%"], on=["facility_id", "gcm", "abbr", "grid_cell"], suffixes=("_60", "_775"))
        assert len(both) == len(pm["60%"]) == len(pm["77.5%"]), "the two runs must cover the same facilities and models"
        for col in cols:
            period = col.removeprefix("d_")
            name = label if col == "historical" else f"{label}_change"
            for run in ("60%", "77.5%"):
                s = stats(by_model(pm[run], col))
                r.add_stats(s, "Case 1", name, period, kind, col != "historical", regions, setpoint=run)
                if metric == "pue" and col == "d_midcentury":
                    found["sp_60" if run == "60%" else "sp_775"] = s.loc["WV", "mean"]
            both["diff"] = both[f"{col}_60"] - both[f"{col}_775"]
            s = stats(by_model(both, "diff"))
            r.add_stats(s, "Case 1", name, period, kind, True, regions, setpoint="60% minus 77.5%")
            if metric == "pue" and col == "historical":
                found["sp_hist"] = s["mean"]
    return r.frame(("setpoint",))


def ensemble_facts(d: Data) -> str:
    m = d.manifest[d.manifest["included_in_ensemble"] == "yes"]
    assert len(m) == N_MODELS and set(m["model"]) == set(d.cache["77.5%"]["gcm"].unique())
    excluded = d.manifest[d.manifest["included_in_ensemble"] != "yes"]
    eco = d.fac.groupby("abbr").agg(n_facilities=("facility_id", "size"), n_grid_cells=("grid_cell", "nunique")).loc[ORDER]
    eco.insert(0, "ecoregion", [NAME[a] for a in ORDER])
    cells = (d.fac.groupby(["grid_cell", "grid_lat", "grid_lon", "abbr"]).size().rename("n_facilities").reset_index()
             .sort_values(["n_facilities", "grid_cell"], ascending=[False, True]))
    cells["ecoregion"] = cells.pop("abbr").map(NAME)
    return "\n".join([
        "# Ensemble facts",
        "",
        f"Scenario SSP3-7.0. Periods: {', '.join(PERIODS.values())}.",
        "",
        f"## Climate models ({len(m)}, one member each)",
        "",
        md_table(m[["model", "member"]]),
        "",
        f"Candidates not used ({len(excluded)}): "
        + "; ".join(f"{x.model} ({x.exclusion_reason})" for x in excluded.itertuples()) + ".",
        "",
        f"## Facilities ({len(d.fac)}) and LOCA2 grid cells ({d.fac['grid_cell'].nunique()})",
        "",
        md_table(eco),
        "",
        "## Facilities per grid cell",
        "",
        "Grid cell IDs are the row and column indices used in `facility_grid_cells.csv`; latitude and longitude are the cell centre.",
        "",
        md_table(cells[["grid_cell", "grid_lat", "grid_lon", "ecoregion", "n_facilities"]]),
        "",
    ])


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------


def quoted_check(found: dict) -> pd.DataFrame:
    def row(key: str, computed: str, agrees: bool, note: str = "") -> dict:
        q = QUOTED[key]
        return {"quantity": q.label, "quoted": q.quoted, "computed": computed, "agrees": "yes" if agrees else "NO", "note": note}

    def same(key: str, value: float) -> bool:
        q = QUOTED[key]
        return round(float(value), q.decimals) == round(float(q.quoted), q.decimals)

    lo, hi = sorted([found["sp_775"], found["sp_60"]])
    hist = found["sp_hist"]
    share = found["aug_share"]
    rows = [
        row("wv_aug", fmt(found["wv_aug"], "pue", True), same("wv_aug", found["wv_aug"])),
        row("cp_sep", fmt(found["cp_sep"], "pue", True), same("cp_sep", found["cp_sep"])),
        row("sp_775", f"{lo:.3f}", same("sp_775", lo),
            f"This is the {'77.5%' if lo == found['sp_775'] else '60%'} setpoint."),
        row("sp_60", f"{hi:.3f}", same("sp_60", hi),
            f"This is the {'77.5%' if hi == found['sp_775'] else '60%'} setpoint."),
        row("sp_hist", "; ".join(f"{a} {fmt(hist[a], 'pue', True)}" for a in ORDER + [ALL]),
            any(same("sp_hist", abs(hist[a])) for a in ORDER + [ALL]),
            "60% minus 77.5%. Agreement is judged on the magnitude, for any ecoregion or all facilities."),
        row("aug_share", "; ".join(f"{a} {fmt(share[a], 'share')}" for a in ORDER + [ALL]),
            bool((share > 1 / 3).any()),
            "Day-weighted August (PUE − 1) over the annual sum. August's share of the mid-century change in overhead is "
            + "; ".join(f"{a} {fmt(found['aug_share_change'][a], 'share')}" for a in ORDER + [ALL]) + "."),
        row("warm_lo", f"{found['warm_lo']:.1f}", same("warm_lo", found["warm_lo"])),
        row("warm_hi", f"{found['warm_hi']:.1f}", same("warm_hi", found["warm_hi"])),
    ]
    return pd.DataFrame(rows)


def readme(check: pd.DataFrame, notes6: list[str]) -> str:
    return "\n".join([
        "# Chapter supporting tables",
        "",
        "Built by `scripts/build_chapter_tables.py` from saved outputs on the local machine. No climate or PUE/WUE "
        "model was re-run and no daily file was read. 109 facilities, 20 climate models, SSP3-7.0; historical "
        "1985–2014, mid-century 2045–2074, end-of-century 2075–2100.",
        "",
        "**Aggregation order, all tables.** Annual means per facility and model; change is the future value minus "
        "that facility's own 1985–2014 value in the same model; the ecoregion value is the equal-weight mean of its "
        "facilities within each model; the reported `mean`, `p25` and `p75` are the mean and the 25th and 75th "
        "percentiles of those ecoregion values across the 20 models. This is the order in `methodology.md`. The "
        "ensemble mean is used, not the median. `facility_min` and `facility_max` are taken over each facility's "
        "ensemble-mean value.",
        "",
        "**Decimals.** PUE 3, WUE 2, exceedance days and temperature 1, shares 3. Differences smaller than 0.01 (PUE) "
        "or 0.1 (WUE) carry one more decimal so that monthly changes and the setpoint difference do not round to zero.",
        "",
        "**Files.** Each table has a tidy `.csv` (`ecoregion, case, metric, period, statistic, value, n_facilities, "
        "n_models`; Table 6 adds `month`, Table 7 adds `setpoint`) and a readable `.md`.",
        "",
        "## Tables",
        "",
        "**Table 1, `table1_warming`.** Change in annual mean daily maximum temperature by ecoregion, mid- and "
        "end-of-century, with the lowest and highest ecoregion mean. Daily minimum temperature is skipped (see below).",
        "",
        "**Table 2, `table2_historical_baselines`.** Historical ensemble mean and the range across facilities of "
        "annual exceedance days (Tmax > 35 °C), annual mean afternoon wet-bulb temperature, and PUE, PUE overhead "
        "(PUE − 1) and WUE for Cases 1, 2 and 5. ECSF and KM have one facility each, so their range is a single value.",
        "",
        "**Table 3, `table3_cases_2_and_5`.** Historical PUE, PUE overhead and WUE with the mid- and end-of-century "
        "change and interquartile range, for Cases 2 and 5. Recomputed from the same cache as "
        "`docs/results-figures/draft3/table_cases25.csv` and in tidy form.",
        "",
        "**Table 4, `table4_model_agreement`.** For Cases 1, 2 and 5, PUE and WUE, and for exceedance days: the "
        "number of the 20 models in which the Willamette Valley change exceeds the Columbia Plateau change (and the "
        "reverse for exceedance days), and the number in which each ecoregion's change is positive.",
        "",
        "**Table 6, `table6_monthly_pue_case1`.** Case 1 monthly historical PUE and overhead, and monthly mid- and "
        "end-of-century change, for all five ecoregions. Also August's share of annual overhead, day-weighted with "
        "fixed month lengths and computed on each model's ecoregion-mean monthly values.",
        "",
        *notes6,
        "",
        "**Table 7, `table7_setpoint_sensitivity`.** Case 1 PUE and WUE under the 60% and 77.5% upper supply-air "
        "relative humidity setpoints: historical value, mid- and end-of-century change, and the difference "
        "(60% minus 77.5%) computed per facility and model before aggregation.",
        "",
        "**`ensemble_facts.md`.** The 20 models and members, facilities and grid cells per ecoregion, and facilities "
        "per grid cell.",
        "",
        "## Skipped",
        "",
        "- **Daily minimum temperature (Table 1).** Only daily maximum temperature was saved per facility and model.",
        "- **Table 5, humidity and wet-bulb change.** Saved humidity is historical only (afternoon relative humidity "
        "and wet-bulb). Future values, and daily mean or maximum relative humidity, need the daily files on the server.",
        "- **Table 8, cooling-mode shares.** The pipeline does not store the economizer branch or chiller mode of "
        "each day. This needs the daily weather on the server, and possibly a rerun.",
        "",
        "## Quoted values against computed values",
        "",
        md_table(check),
        "",
    ])


# ---------------------------------------------------------------------------


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out = root / "outputs/chapter_tables"
    out.mkdir(parents=True, exist_ok=True)
    d = Data(root)
    n_fac = d.fac["abbr"].value_counts()
    n_fac[ALL] = len(d.fac)
    found: dict = {}

    t1 = table1(d, n_fac, found)
    t2 = table2(d, n_fac)
    t3 = table3(d, n_fac)
    t4 = table4(d, n_fac)
    t6, notes6 = table6(d, n_fac, found)
    t7 = table7(d, n_fac, found)

    tables = {
        "table1_warming": (t1, ["ecoregion", "n_facilities"], "Table 1. Change in annual mean daily maximum temperature (°C)"),
        "table2_historical_baselines": (t2, ["ecoregion", "n_facilities", "case", "metric"], "Table 2. Historical baselines (1985–2014)"),
        "table3_cases_2_and_5": (t3, ["case", "ecoregion", "n_facilities", "metric"], "Table 3. Cases 2 and 5 by ecoregion"),
        "table4_model_agreement": (t4, ["case", "metric", "ecoregion"], "Table 4. Model agreement (number of 20 models)"),
        "table6_monthly_pue_case1": (t6, ["ecoregion", "month", "metric"], "Table 6. Case 1 monthly PUE"),
        "table7_setpoint_sensitivity": (t7, ["ecoregion", "metric", "setpoint"], "Table 7. Case 1 supply-air RH setpoint sensitivity"),
    }
    for name, (df, index, title) in tables.items():
        df.to_csv(out / f"{name}.csv", index=False)
        (out / f"{name}.md").write_text(f"# {title}\n\n{md_table(wide(df, index))}\n")

    (out / "ensemble_facts.md").write_text(ensemble_facts(d))
    check = quoted_check(found)
    text = readme(check, notes6)
    (out / "README.md").write_text(text)
    pd.DataFrame([{
        "script": "scripts/build_chapter_tables.py", "n_facilities": len(d.fac), "n_models": N_MODELS, "scenario": "ssp370",
        "periods": "; ".join(f"{k}={v}" for k, v in PERIODS.items()), "cases": "; ".join(f"{k}={v}" for k, v in CASES.items()),
        "ensemble_statistic": "mean, p25, p75 across models", "tables_built": "1 (tasmax only); 2; 3; 4; 6; 7",
        "tables_skipped": "5; 8; tasmin", "inputs": "outputs/ensemble_pue_wue_facilities (saved summaries only)",
    }]).to_csv(out / "run_manifest.csv", index=False)

    print(f"Wrote {out}")
    print("\n".join(notes6))
    print(check.to_string(index=False))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the heat chapter's supporting tables from saved outputs.")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
