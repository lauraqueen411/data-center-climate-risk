"""Export the values needed to write revised draft 3 (docs/tasks/draft3-data-tables.md).

Writes everything to ``docs/results-figures/draft3/``: ``draft3_values.md``,
``table1.csv``, ``table_cases25.csv``, ``monthly.csv``, ``rh_sensitivity.csv``,
``warming.csv``, ``facility_facts.csv``, ``model_agreement.csv``, plus copies
of the final Figure 2-5 PNGs and ``methodology.md``.

Conventions (``docs/results-figures/methodology.md``): annual means; change =
future minus historical per facility and model; ecoregion mean of facilities
(equal weight per record) within each model; ensemble mean and interquartile
range across the 20 models. Facility-level facts use each facility's
ensemble-mean value.

Decimals: exceedance days 1; PUE 3 (4 for differences under 0.01); WUE 2 (3 for
differences under 0.1); temperature 1.

Inputs (all local, nothing re-run):
- ``outputs/ensemble_pue_wue_facilities/analysis/``: ``fig_show_everything/model_values_cache.csv``
  and ``facility_grid_cells.csv``, ``facility_month_summary.csv``,
  ``facility_annual_tmax.csv`` (historical, mid-century),
  ``facility_annual_tmax_endcentury.csv`` and ``facility_monthly_tmax.csv``
  (both read on the server from the Case 1 daily ``tasmax_c`` column on 2026-10-08).
- ``.../superseded_rh_up60_2026-10-08/analysis/fig_show_everything/model_values_cache.csv``: the 60% run.
- ``outputs/ensemble_pue_wue/gcm_manifest.csv``: models and members.

Assumptions log
----------------
- **Cooling-mode shares** (task F) are not computed: the PUE/WUE pipeline does
  not store which economizer branch or chiller mode each pairing uses.
- **Ties** in facility minima and maxima are common (facilities sharing a
  LOCA2 cell have identical values); every facility at the extreme is counted.
"""

from __future__ import annotations

import argparse
import dataclasses
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from climate_risk_dc.climate.cooling_archetypes import Case1Config, Case2Config, Case5Config

N_MODELS = 20
ORDER = ["ECSF", "BM", "CP", "KM", "WV"]
ABBR = {
    "Eastern Cascades Slopes and Foothills": "ECSF",
    "Blue Mountains": "BM",
    "Columbia Plateau": "CP",
    "Klamath Mountains/California High North Coast Range": "KM",
    "Willamette Valley": "WV",
}
NAME = {"ECSF": "Eastern Cascades Slopes and Foothills", "BM": "Blue Mountains", "CP": "Columbia Plateau",
        "KM": "Klamath Mountains", "WV": "Willamette Valley"}
CASES = {"Case 1": "ae-chiller", "Case 2": "we-chiller", "Case 5": "chiller-only"}
PERIOD_LABEL = {"midcentury": "2045–2074", "endcentury": "2075–2100"}
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


# ---------------------------------------------------------------------------
# Rounding
# ---------------------------------------------------------------------------


def rnd(x: float, kind: str, diff: bool = False) -> float:
    if pd.isna(x):
        return x
    d = {"days": 1, "temp": 1, "pue": 3, "wue": 2}[kind]
    if diff and kind == "pue" and abs(x) < 0.01:
        d = 4
    if diff and kind == "wue" and abs(x) < 0.1:
        d = 3
    return round(float(x), d)


def _decimals(x: float, kind: str, diff: bool) -> int:
    d = {"days": 1, "temp": 1, "pue": 3, "wue": 2}[kind]
    if diff and kind == "pue" and abs(x) < 0.01:
        d = 4
    if diff and kind == "wue" and abs(x) < 0.1:
        d = 3
    return d


def t(x: float, kind: str, diff: bool = False) -> str:
    """Table cell: fixed decimals, no sign, so trailing zeros survive in CSV and markdown."""
    if pd.isna(x):
        return ""
    x = float(x)
    return f"{round(x, _decimals(x, kind, diff)):.{_decimals(x, kind, diff)}f}"


def f(x: float, kind: str, diff: bool = False) -> str:
    x = float(x)
    v = rnd(x, kind, diff)
    d = {"days": 1, "temp": 1, "pue": 3, "wue": 2}[kind]
    if diff and kind == "pue" and abs(x) < 0.01:
        d = 4
    if diff and kind == "wue" and abs(x) < 0.1:
        d = 3
    return f"{v:+.{d}f}" if diff else f"{v:.{d}f}"


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


class Data:
    def __init__(self, root: Path) -> None:
        self.a = root / "outputs/ensemble_pue_wue_facilities/analysis"
        self.fac = pd.read_csv(self.a / "fig_show_everything/facility_grid_cells.csv")
        self.fac["abbr"] = self.fac["ecoregion"].map(ABBR)
        assert len(self.fac) == 109 and self.fac["abbr"].notna().all()
        self.cache = {"77.5%": pd.read_csv(self.a / "fig_show_everything/model_values_cache.csv"),
                      "60%": pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/superseded_rh_up60_2026-10-08/"
                                         "analysis/fig_show_everything/model_values_cache.csv")}
        self.month = pd.read_csv(self.a / "facility_month_summary.csv")
        t = pd.concat([pd.read_csv(self.a / "facility_annual_tmax.csv"),
                       pd.read_csv(self.a / "facility_annual_tmax_endcentury.csv")])
        self.tmax = t
        self.mtmax = pd.read_csv(self.a / "facility_monthly_tmax.csv")
        self.manifest = pd.read_csv(root / "outputs/ensemble_pue_wue/gcm_manifest.csv")
        for name, df in (("cache", self.cache["77.5%"]), ("cache60", self.cache["60%"]), ("month", self.month),
                         ("tmax", self.tmax), ("monthly tmax", self.mtmax)):
            assert df["gcm"].nunique() == N_MODELS, f"{name}: expected {N_MODELS} models"

    def per_model(self, metric: str, archetype: str, run: str = "77.5%") -> pd.DataFrame:
        """facility_id, abbr, gcm, historical, d_midcentury, d_endcentury."""
        c = self.cache[run]
        c = c[(c["metric"] == metric) & (c["archetype"] == archetype)]
        w = c.pivot_table(index=["facility_id", "gcm"], columns="period", values="value").reset_index()
        for p in PERIOD_LABEL:
            if p in w:
                w[f"d_{p}"] = w[p] - w["historical"]
        return w.merge(self.fac[["facility_id", "abbr", "grid_cell"]], on="facility_id")


def eco(per_model: pd.DataFrame, col: str) -> pd.DataFrame:
    """Ecoregion mean within each model, then mean and IQR across models."""
    m = per_model.groupby(["abbr", "gcm"])[col].mean()
    g = m.groupby(level="abbr")
    out = pd.DataFrame({"mean": g.mean(), "p25": g.quantile(0.25), "p75": g.quantile(0.75), "n_models": g.size()})
    assert (out["n_models"] == N_MODELS).all()
    return out.loc[ORDER]


def eco_by_model(per_model: pd.DataFrame, col: str) -> pd.DataFrame:
    return per_model.groupby(["gcm", "abbr"])[col].mean().unstack("abbr")


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------


def table_main(d: Data, cases: list[str], hist_iqr: bool) -> pd.DataFrame:
    n = d.fac["abbr"].value_counts()
    ex = d.per_model("exceedance_days", "all")
    wb = d.per_model("wetbulb_annual", "ae-chiller")
    rows = []
    for case in cases:
        pue, wue = d.per_model("pue", CASES[case]), d.per_model("wue", CASES[case])
        h = {"pue": eco(pue, "historical"), "wue": eco(wue, "historical"),
             "ex": eco(ex, "historical"), "wb": eco(wb, "historical")}
        for abbr in ORDER:
            for p, label in PERIOD_LABEL.items():
                r = {"case": case, "ecoregion": NAME[abbr], "n_facilities": int(n[abbr]),
                     "hist_exceedance_days": t(h["ex"].loc[abbr, "mean"], "days"),
                     "hist_wetbulb_C": t(h["wb"].loc[abbr, "mean"], "temp"),
                     f"hist_pue": t(h["pue"].loc[abbr, "mean"], "pue"), f"hist_wue": t(h["wue"].loc[abbr, "mean"], "wue")}
                if hist_iqr:
                    for k in ("pue", "wue"):
                        r[f"hist_{k}_p25"], r[f"hist_{k}_p75"] = t(h[k].loc[abbr, "p25"], k), t(h[k].loc[abbr, "p75"], k)
                r["period"] = label
                for col, src, kind in (("d_exceedance", ex, "days"), ("d_pue", pue, "pue"), ("d_wue", wue, "wue")):
                    s = eco(src, f"d_{p}").loc[abbr]
                    r[col], r[f"{col}_p25"], r[f"{col}_p75"] = (t(s[k], kind, True) for k in ("mean", "p25", "p75"))
                rows.append(r)
    return pd.DataFrame(rows)


def facility_facts(d: Data) -> tuple[pd.DataFrame, list[str]]:
    pue = d.per_model("pue", "ae-chiller")
    fac = pue.groupby(["facility_id", "abbr", "grid_cell"])[["historical", "d_midcentury"]].mean().reset_index()
    wb = d.per_model("wetbulb_annual", "ae-chiller").groupby("facility_id")["historical"].mean().rename("wetbulb")
    fac = fac.merge(wb, on="facility_id")
    lines, rows = [], []

    def extreme(col: str, which: str, sub: pd.DataFrame, scope: str, kind: str, diff: bool) -> None:
        v = sub[col].min() if which == "min" else sub[col].max()
        at = sub[np.isclose(sub[col], v, atol=1e-12)]
        rows.append({"quantity": col, "scope": scope, "which": which, "value": t(v, kind, diff),
                     "n_facilities_at_value": len(at), "ecoregions": "; ".join(sorted(at["abbr"].unique())),
                     "grid_cells": "; ".join(sorted(at["grid_cell"].unique()))})

    for which in ("min", "max"):
        extreme("historical", which, fac, "all", "pue", False)
        extreme("d_midcentury", which, fac, "all", "pue", True)
    for abbr in ORDER:
        for which in ("min", "max"):
            extreme("d_midcentury", which, fac[fac["abbr"] == abbr], abbr, "pue", True)
    tbl = pd.DataFrame(rows)
    s = fac.loc[fac["d_midcentury"].idxmin()]
    wv_low = fac[fac["abbr"] == s["abbr"]]["wetbulb"].min()
    lines.append(f"- **Smallest mid-century change:** facility {int(s['facility_id'])}, {NAME[s['abbr']]} ({s['abbr']}), "
                 f"grid cell {s['grid_cell']}, change {f(s['d_midcentury'], 'pue', True)}, historical wet-bulb "
                 f"{f(s['wetbulb'], 'temp')} °C (lowest in its ecoregion: "
                 f"{'yes' if np.isclose(s['wetbulb'], wv_low) else 'no'}, ecoregion minimum {f(wv_low, 'temp')} °C).")
    cells = d.fac.groupby("abbr")["grid_cell"].nunique().loc[ORDER]
    lines.append(f"- **Distinct LOCA2 grid cells:** {d.fac['grid_cell'].nunique()} in total; "
                 + ", ".join(f"{a} {cells[a]}" for a in ORDER) + ".")
    return tbl, lines


def agreement(d: Data, run: str = "77.5%", pairs=(("WV", "CP"), ("WV", "KM"), ("CP", "KM"))) -> pd.DataFrame:
    rows = []
    for metric in ("pue", "wue"):
        pm = d.per_model(metric, "ae-chiller", run)
        for p, label in PERIOD_LABEL.items():
            m = eco_by_model(pm, f"d_{p}")
            for a, b in pairs:
                rows.append({"setpoint": run, "metric": metric.upper(), "period": label, "comparison": f"{a} > {b}",
                             "n_models": int((m[a] > m[b]).sum()), "of": N_MODELS})
    return pd.DataFrame(rows)


def monthly(d: Data) -> tuple[pd.DataFrame, list[str]]:
    mo = d.month.merge(d.fac[["facility_id", "abbr"]], on="facility_id")
    mo = mo[mo["abbr"].isin(["CP", "WV"])]
    mo["case"] = mo["archetype"].map({v: k for k, v in CASES.items()})
    rows = []

    def add(series_by_model: pd.Series, region: str, quantity: str, case: str, kind: str, diff: bool) -> None:
        g = series_by_model.groupby(level="month")
        for m in range(1, 13):
            v = series_by_model.xs(m, level="month")
            rows.append({"region": region, "quantity": quantity, "case": case, "month": m, "month_name": MONTHS[m - 1],
                         "mean": t(v.mean(), kind, diff), "p25": t(v.quantile(0.25), kind, diff),
                         "p75": t(v.quantile(0.75), kind, diff), "_raw": v.mean()})

    for region in ("CP", "WV"):
        r = mo[mo["abbr"] == region]
        for case in CASES:
            c = r[r["case"] == case]
            w = c.pivot_table(index=["facility_id", "gcm", "month"], columns="period", values=["pue_mean", "wue_mean"])
            dp = (w[("pue_mean", "midcentury")] - w[("pue_mean", "historical")]).groupby(level=["gcm", "month"]).mean()
            add(dp, region, "d_pue_midcentury", case, "pue", True)
            if case == "Case 1":
                dw = (w[("wue_mean", "midcentury")] - w[("wue_mean", "historical")]).groupby(level=["gcm", "month"]).mean()
                add(dw, region, "d_wue_midcentury", case, "wue", True)
                hp = w[("pue_mean", "historical")].groupby(level=["gcm", "month"]).mean()
                add(hp, region, "hist_pue", case, "pue", False)
        mt = d.mtmax.merge(d.fac[["facility_id", "abbr"]], on="facility_id")
        mt = mt[(mt["abbr"] == region) & (mt["period"] == "historical")].groupby(["gcm", "month"])["tasmax_c"].mean()
        add(mt, region, "hist_tmax_C", "all", "temp", False)
    out = pd.DataFrame(rows)

    lines = []
    for region in ("CP", "WV"):
        for quantity, kind in (("d_pue_midcentury", "pue"), ("d_wue_midcentury", "wue")):
            for case in CASES:
                s = out[(out["region"] == region) & (out["quantity"] == quantity) & (out["case"] == case)]
                if s.empty:
                    continue
                pk = s.loc[s["_raw"].idxmax()]
                lines.append(f"- {region}, {case}, {quantity.split('_')[1].upper()} change: peak {pk['month_name']} "
                             f"{f(pk['_raw'], kind, True)} [{f(pk['p25'], kind, True)}, {f(pk['p75'], kind, True)}]")
    cp = out[(out["region"] == "CP") & (out["quantity"] == "d_pue_midcentury") & (out["case"] == "Case 1")].set_index("month")["_raw"]
    is_min = cp[7] < cp[6] and cp[7] < cp[8]
    lines.append(f"- CP's July Case 1 PUE change is {'a' if is_min else 'not a'} local minimum "
                 f"(Jun {f(cp[6], 'pue', True)}, Jul {f(cp[7], 'pue', True)}, Aug {f(cp[8], 'pue', True)}).")
    for region in ("CP", "WV"):
        tm = out[(out["region"] == region) & (out["quantity"] == "hist_tmax_C")].set_index("month")["_raw"]
        hp = out[(out["region"] == region) & (out["quantity"] == "hist_pue")].set_index("month")["_raw"]
        lines.append(f"- {region}: hottest month by historical mean Tmax is {MONTHS[tm.idxmax() - 1]} "
                     f"(Jul {f(tm[7], 'temp')} °C, Aug {f(tm[8], 'temp')} °C); highest historical Case 1 PUE is "
                     f"{MONTHS[hp.idxmax() - 1]} (Jul {hp[7]:.4f}, Aug {hp[8]:.4f}; four decimals to separate them).")
    return out.drop(columns="_raw"), lines


def rh_sensitivity(d: Data) -> tuple[pd.DataFrame, list[str]]:
    rows, raw = [], {}
    for run in ("60%", "77.5%"):
        for metric in ("pue", "wue"):
            pm = d.per_model(metric, "ae-chiller", run)
            for col in ("historical", "d_midcentury", "d_endcentury"):
                raw[(run, metric, col)] = eco(pm, col)
    for abbr in ORDER:
        for run in ("60%", "77.5%", "difference (60% − 77.5%)"):
            r = {"ecoregion": abbr, "setpoint": run}
            for metric in ("pue", "wue"):
                for col in ("historical", "d_midcentury", "d_endcentury"):
                    diff = col != "historical" or run.startswith("difference")
                    if run.startswith("difference"):
                        v = raw[("60%", metric, col)].loc[abbr, "mean"] - raw[("77.5%", metric, col)].loc[abbr, "mean"]
                        r[f"{metric}_{col}"] = t(v, metric, True)
                    else:
                        s = raw[(run, metric, col)].loc[abbr]
                        r[f"{metric}_{col}"] = t(s["mean"], metric, diff)
                        if col != "historical":
                            r[f"{metric}_{col}_p25"], r[f"{metric}_{col}_p75"] = t(s["p25"], metric, True), t(s["p75"], metric, True)
            rows.append(r)
    ag = agreement(d, "60%", pairs=(("WV", "CP"),))
    lines = ["- WV > CP model agreement at 60%: " + "; ".join(
        f"{r.metric} {r.period} {r.n_models}/20" for r in ag.itertuples())]
    return pd.DataFrame(rows), lines


def warming(d: Data) -> pd.DataFrame:
    tx = d.tmax.groupby(["facility_id", "gcm", "period"])["tasmax_c"].mean().unstack("period").reset_index()
    for p in PERIOD_LABEL:
        tx[f"d_{p}"] = tx[p] - tx["historical"]
    tx = tx.merge(d.fac[["facility_id", "abbr"]], on="facility_id")
    fac = tx.groupby(["facility_id", "abbr"])[["d_midcentury", "d_endcentury"]].mean().reset_index()
    rows = []
    for p, label in PERIOD_LABEL.items():
        e = eco(tx, f"d_{p}")
        for abbr in ORDER:
            s = fac[fac["abbr"] == abbr][f"d_{p}"]
            rows.append({"ecoregion": abbr, "period": label, "mean": t(e.loc[abbr, "mean"], "temp"),
                         "p25": t(e.loc[abbr, "p25"], "temp"), "p75": t(e.loc[abbr, "p75"], "temp"),
                         "facility_min": t(s.min(), "temp"), "facility_max": t(s.max(), "temp")})
    return pd.DataFrame(rows)


def methods_facts(d: Data) -> list[str]:
    lines = []
    for case, cfg in (("Case 1", Case1Config), ("Case 2", Case2Config), ("Case 5", Case5Config)):
        fields = [x for x in dataclasses.fields(cfg) if x.name != "p_atm_pa"]
        lines.append(f"- **{case}: {len(fields)} parameters** at Table B.1 midpoints: "
                     + ", ".join(f"`{x.name}` = {x.default:g}" for x in fields) + ".")
    lines.append("- `p_atm_pa` is not a Table B.1 parameter; it is replaced per facility by the elevation-adjusted "
                 "pressure (`elevation_adjusted_pressure_pa`, standard barometric formula on the LOCA elevation field).")
    lines.append("- **Wet-bulb:** `cooling_archetypes_fast._wetbulb_from_rh`, eight Newton iterations solving for the "
                 "temperature at which saturated-air enthalpy equals the outdoor-air enthalpy, with Arden Buck "
                 "saturation vapour pressure and the ASHRAE moist-air enthalpy formula, at the facility's "
                 "elevation-adjusted pressure. The model uses it for the chiller COP input; the figures use the same "
                 "function on the hot pairing (daily Tmax with RHmin) for historical afternoon wet-bulb. Within 0.21 °C "
                 "of CoolProp where wet-bulb is above 5 °C (audit B4).")
    m = d.manifest[d.manifest["included_in_ensemble"] == "yes"]
    lines.append(f"- **Models ({len(m)}), SSP3-7.0, one member each:** "
                 + ", ".join(f"{r.model} ({r.member})" for r in m.itertuples()) + ".")
    lines.append("- **Midpoints:** every equipment parameter in the current (77.5%) run is the Table B.1 midpoint "
                 "(audit B1). The 60% run differs from it in Case 1 only by `rh_up_pct` (60.0, the lower bound of "
                 "60–95%); its Case 5 also used `delta_t_water_c` 5.0 and `fan_pressure_ct_pa` 250, which do not "
                 "enter Case 1. Confirmed from the code history: commit ac70bb8 changes only these three defaults, "
                 "and the one other later commit to these modules (0875509) changes a docstring and removes an "
                 "unused function.")
    return lines


# ---------------------------------------------------------------------------


def md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(out)


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out = root / "docs/results-figures/draft3"
    out.mkdir(parents=True, exist_ok=True)
    d = Data(root)

    t1 = table_main(d, ["Case 1"], hist_iqr=False).drop(columns="case")
    t25 = table_main(d, ["Case 2", "Case 5"], hist_iqr=True).drop(
        columns=["hist_exceedance_days", "hist_wetbulb_C", "d_exceedance", "d_exceedance_p25", "d_exceedance_p75"])
    facts, fact_lines = facility_facts(d)
    agree = agreement(d)
    mon, mon_lines = monthly(d)
    rh, rh_lines = rh_sensitivity(d)
    warm = warming(d)
    meth = methods_facts(d)

    t1.to_csv(out / "table1.csv", index=False)
    t25.to_csv(out / "table_cases25.csv", index=False)
    facts.to_csv(out / "facility_facts.csv", index=False)
    agree.to_csv(out / "model_agreement.csv", index=False)
    mon.to_csv(out / "monthly.csv", index=False)
    rh.to_csv(out / "rh_sensitivity.csv", index=False)
    warm.to_csv(out / "warming.csv", index=False)

    pv = agree.pivot_table(index=["metric", "comparison"], columns="period", values="n_models", aggfunc="first").reset_index()
    pv[list(PERIOD_LABEL.values())] = pv[list(PERIOD_LABEL.values())].astype(int)
    md = [
        "# Values for revised draft 3",
        "",
        "Generated by `scripts/build_draft3_values.py` from the final, audited outputs. Conventions as in "
        "`methodology.md`: annual means; change = future minus historical per facility and model; ecoregion mean "
        "of facilities within each model; ensemble mean with the interquartile range across the 20 models in "
        "brackets. Decimals: exceedance days 1, PUE 3 (4 below 0.01 for differences), WUE 2 (3 below 0.1 for "
        "differences), temperature 1.",
        "",
        "## A. Table 1 (Case 1) — `table1.csv`",
        "",
        md_table(t1),
        "",
        "## B. Cases 2 and 5 — `table_cases25.csv`",
        "",
        md_table(t25),
        "",
        "## C. Facility-level facts (Case 1) — `facility_facts.csv`",
        "",
        md_table(facts),
        "",
        *fact_lines,
        "",
        "## D. Model agreement (Case 1, number of models of 20) — `model_agreement.csv`",
        "",
        md_table(pv),
        "",
        "## E. Monthly change (Figure 4) — `monthly.csv`",
        "",
        "Peak month of the mid-century change, mean [p25, p75]:",
        "",
        *mon_lines,
        "",
        "## F. Humidity setpoint sensitivity (Figure 5) — `rh_sensitivity.csv`",
        "",
        md_table(rh),
        "",
        *rh_lines,
        "- `rh_up` is the only Case 1 parameter that differs between the two runs (see H).",
        "- **Cooling-mode shares:** not available. The PUE/WUE pipeline does not store which economizer branch or "
        "chiller mode each pairing uses, so the June–September share in supplemental-chiller mode cannot be read "
        "from the outputs. It could be computed on the server from the daily weather (~30 min) while access lasts.",
        "",
        "## G. Warming (change in annual mean daily Tmax, °C) — `warming.csv`",
        "",
        md_table(warm),
        "",
        "## H. Methods facts",
        "",
        *meth,
        "",
    ]
    (out / "draft3_values.md").write_text("\n".join(md) + "\n")

    figs = root / "docs/results-figures"
    for name in ("fig2_show_everything.png", "fig3_maps.png", "fig4_monthly.png", "fig5_rh_setpoint.png", "methodology.md"):
        shutil.copy2(figs / name, out / name)
    print(f"Wrote {out}")
    print("\n".join(mon_lines + rh_lines + fact_lines))
    print(pv.to_string(index=False))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export the values for revised draft 3.")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
