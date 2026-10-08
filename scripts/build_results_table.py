"""Table 1 and the Results number check for the heat risk chapter
(docs/tasks/results-table.md).

Uses the same analysis outputs and conventions as Figures 2-5
(docs/results-figures/methodology.md): annual means; change = future minus
historical per facility and model; ecoregion mean of facilities within each
model; ensemble mean with the interquartile range across the 20 models.
No figure is changed.

Inputs (``outputs/ensemble_pue_wue_facilities/analysis/``):
- ``fig_show_everything/model_values_cache.csv``: period-mean annual PUE,
  WUE, exceedance days and historical wet-bulb per facility and model.
- ``fig_show_everything/facility_grid_cells.csv``: ecoregion, grid cell,
  lon/lat.
- ``facility_month_summary.csv``: monthly mean PUE/WUE per facility and model.
- ``facility_annual_tmax.csv``: annual mean daily Tmax per facility and model.

Outputs (``docs/results-figures/``):
- ``table1.csv``: Case 1, one row per ecoregion and future period.
- ``table1_cases.csv``: the same layout for Cases 2 and 5 (reference only).
- ``results_number_check.csv``: item, draft value, computed value, match.
- ``saturation_check.csv``: historical monthly Case 1 PUE for CP and WV.

Assumptions log
----------------
- **Rounding for the check.** A draft value matches when the computed value
  rounded to the draft's decimals equals it (half away from zero). Ranges
  match when both ends match. The computed column shows one more decimal
  than the draft so near misses are visible.
- **"Everywhere" and "every ecoregion"** claims are tested on ecoregion
  means (the plotted aggregate), both future periods unless the item names
  one. The facility-level maximum is given alongside for reference.
- **Ecoregion order** compares the ranking of ecoregion-mean changes.
- **Parameter counts** are the fields of ``Case1Config``/``Case2Config``/
  ``Case5Config`` except ``p_atm_pa``, which is overridden by the
  elevation-adjusted pressure and is not a Table B.1 parameter.
- **Cooling modes** are not stored by the PUE/WUE pipeline (the daily CSVs
  hold PUE/WUE and weather only), so the supplemental-chiller share of the
  saturation check is not computed here.
"""

from __future__ import annotations

import argparse
import dataclasses
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import numpy as np
import pandas as pd

from climate_risk_dc.climate.cooling_archetypes import Case1Config, Case2Config, Case5Config

N_MODELS = 20
ORDER = ["ECSF", "BM", "CP", "KM", "WV"]
ECOREGION_ABBR = {
    "Eastern Cascades Slopes and Foothills": "ECSF",
    "Blue Mountains": "BM",
    "Columbia Plateau": "CP",
    "Klamath Mountains/California High North Coast Range": "KM",
    "Willamette Valley": "WV",
}
TABLE_NAME = {
    "ECSF": "Eastern Cascades Slopes and Foothills", "BM": "Blue Mountains", "CP": "Columbia Plateau",
    "KM": "Klamath Mountains", "WV": "Willamette Valley",
}
CASES = {"ae-chiller": "Case 1", "we-chiller": "Case 2", "chiller-only": "Case 5"}
PERIOD_LABEL = {"midcentury": "2045–2074", "endcentury": "2075–2100"}
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def rnd(x: float, d: int) -> float:
    return float(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP))


def fmt(x: float, d: int) -> str:
    return f"{rnd(x, d):.{d}f}"


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


def load(root: Path) -> dict[str, pd.DataFrame]:
    a = root / "outputs/ensemble_pue_wue_facilities/analysis"
    fac = pd.read_csv(a / "fig_show_everything/facility_grid_cells.csv")
    fac["abbr"] = fac["ecoregion"].map(ECOREGION_ABBR)
    assert len(fac) == 109 and fac["abbr"].notna().all()
    cache = pd.read_csv(a / "fig_show_everything/model_values_cache.csv")
    cache["case"] = cache["archetype"].map(CASES).fillna("all")
    month = pd.read_csv(a / "facility_month_summary.csv")
    month["case"] = month["archetype"].map(CASES)
    tmax = pd.read_csv(a / "facility_annual_tmax.csv")
    for name, df in (("cache", cache), ("month", month), ("tmax", tmax)):
        assert df["gcm"].nunique() == N_MODELS, f"{name}: expected {N_MODELS} models"
    return {"fac": fac, "cache": cache, "month": month, "tmax": tmax}


def per_model_changes(cache: pd.DataFrame, fac: pd.DataFrame) -> pd.DataFrame:
    """facility_id, gcm, case, metric, period, change (future minus historical)."""
    v = cache[cache["metric"].isin(["pue", "wue", "exceedance_days"])]
    w = v.pivot_table(index=["facility_id", "gcm", "case", "metric"], columns="period", values="value")
    out = pd.concat([(w[p] - w["historical"]).rename("change").reset_index().assign(period=p) for p in PERIOD_LABEL],
                    ignore_index=True)
    return out.merge(fac[["facility_id", "abbr"]], on="facility_id")


def historical(cache: pd.DataFrame, fac: pd.DataFrame) -> pd.DataFrame:
    """Per-facility historical values (mean across models): exceedance days, wet-bulb, PUE/WUE per case."""
    h = cache[cache["period"] == "historical"]
    h = h.groupby(["facility_id", "metric", "case"])["value"].mean().reset_index()
    return h.merge(fac[["facility_id", "abbr"]], on="facility_id")


def eco_stats(per_model: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """Ecoregion mean within each model, then ensemble mean and IQR across models."""
    m = per_model.groupby(keys + ["abbr", "gcm"])["change"].mean()
    g = m.groupby(level=keys + ["abbr"])
    out = pd.DataFrame({"mean": g.mean(), "p25": g.quantile(0.25), "p75": g.quantile(0.75), "n_models": g.size()})
    assert (out["n_models"] == N_MODELS).all()
    return out.reset_index()


# ---------------------------------------------------------------------------
# Table 1
# ---------------------------------------------------------------------------


def build_table(changes: pd.DataFrame, hist: pd.DataFrame, fac: pd.DataFrame, cases: list[str]) -> pd.DataFrame:
    e = eco_stats(changes, ["case", "metric", "period"]).set_index(["case", "metric", "period", "abbr"])
    hx = hist[hist["metric"] == "exceedance_days"].groupby("abbr")["value"].mean()
    hw = hist[hist["metric"] == "wetbulb_annual"].groupby("abbr")["value"].mean()
    n = fac["abbr"].value_counts()
    rows = []
    for case in cases:
        for abbr in ORDER:
            for period, label in PERIOD_LABEL.items():
                row = {"case": case, "ecoregion": TABLE_NAME[abbr], "n_facilities": int(n[abbr]),
                       "hist_exceedance_days": rnd(hx[abbr], 1), "hist_wetbulb_C": rnd(hw[abbr], 1), "period": label}
                for col, (c, metric, d) in {"d_exceedance": ("all", "exceedance_days", 1), "d_pue": (case, "pue", 3),
                                            "d_wue": (case, "wue", 2)}.items():
                    s = e.loc[(c, metric, period, abbr)]
                    row[col], row[f"{col}_p25"], row[f"{col}_p75"] = rnd(s["mean"], d), rnd(s["p25"], d), rnd(s["p75"], d)
                rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Number check
# ---------------------------------------------------------------------------


class Check:
    def __init__(self) -> None:
        self.rows: list[dict] = []

    def value(self, section: str, item: str, draft: float, computed: float, d: int) -> None:
        self.rows.append({"section": section, "item": item, "draft_value": f"{draft:.{d}f}",
                          "computed_value": f"{computed:.{d + 1}f}", "match": "yes" if rnd(computed, d) == draft else "no"})

    def range(self, section: str, item: str, draft: tuple[float, float], lo: float, hi: float, d: int) -> None:
        ok = rnd(lo, d) == draft[0] and rnd(hi, d) == draft[1]
        self.rows.append({"section": section, "item": item, "draft_value": f"{draft[0]:.{d}f}–{draft[1]:.{d}f}",
                          "computed_value": f"{lo:.{d + 1}f}–{hi:.{d + 1}f}", "match": "yes" if ok else "no"})

    def claim(self, section: str, item: str, draft: str, computed: str, ok: bool) -> None:
        self.rows.append({"section": section, "item": item, "draft_value": draft, "computed_value": computed,
                          "match": "yes" if ok else "no"})


def number_check(data: dict[str, pd.DataFrame], changes: pd.DataFrame, hist: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    fac = data["fac"]
    c = Check()
    e = eco_stats(changes, ["case", "metric", "period"]).set_index(["case", "metric", "period", "abbr"])["mean"]
    hx = hist[hist["metric"] == "exceedance_days"].groupby("abbr")["value"].mean()

    # 3.1
    for abbr, v in {"CP": 19.2, "KM": 17.6, "BM": 7.8, "WV": 4.1, "ECSF": 1.3}.items():
        c.value("3.1", f"Historical exceedance days, {abbr}", v, hx[abbr], 1)
    for period, drafts in (("midcentury", {"CP": 30.4, "WV": 10.9}), ("endcentury", {"CP": 48.9, "WV": 22.6})):
        for abbr, v in drafts.items():
            c.value("3.1", f"Change in exceedance days, {PERIOD_LABEL[period]}, {abbr}", v, e[("all", "exceedance_days", period, abbr)], 1)
    mid_ex = {a: e[("all", "exceedance_days", "midcentury", a)] for a in ORDER}
    smallest = min(mid_ex, key=mid_ex.get)
    c.claim("3.1", "WV has the smallest mid-century increase in exceedance days", "WV", smallest, smallest == "WV")
    for abbr, v in {"WV": 0.021, "KM": 0.019, "CP": 0.017, "BM": 0.011, "ECSF": 0.011}.items():
        c.value("3.1", f"Case 1 change in PUE, 2045–2074, {abbr}", v, e[("Case 1", "pue", "midcentury", abbr)], 3)
    for abbr, v in {"WV": 0.27, "KM": 0.24, "CP": 0.21, "BM": 0.13, "ECSF": 0.13}.items():
        c.value("3.1", f"Case 1 change in WUE, 2045–2074, {abbr}", v, e[("Case 1", "wue", "midcentury", abbr)], 2)
    for abbr, v in {"WV": 0.037, "CP": 0.025}.items():
        c.value("3.1", f"Case 1 change in PUE, 2075–2100, {abbr}", v, e[("Case 1", "pue", "endcentury", abbr)], 3)

    # Model agreement: WV ecoregion-mean change > CP, per model.
    m = changes[changes["abbr"].isin(["CP", "WV"])].groupby(["case", "metric", "period", "abbr", "gcm"])["change"].mean().unstack("abbr")
    agree = (m["WV"] > m["CP"]).groupby(level=["case", "metric", "period"]).sum().rename("n_models_wv_gt_cp").reset_index()
    agree = agree[agree["case"] != "all"]

    wv1 = changes[(changes["abbr"] == "WV") & (changes["case"] == "Case 1") & (changes["metric"] == "pue") & (changes["period"] == "midcentury")]
    wv_fac = wv1.groupby("facility_id")["change"].mean()
    c.range("3.1", "Facility range of Case 1 change in PUE within WV, 2045–2074", (0.014, 0.023), wv_fac.min(), wv_fac.max(), 3)

    # Case 2
    c2_end = {a: e[("Case 2", "pue", "endcentury", a)] for a in ORDER}
    c.claim("3.1", "Case 2 change in PUE below 0.006 in every ecoregion, 2075–2100", "< 0.006",
            f"max {max(c2_end.values()):.4f} ({max(c2_end, key=c2_end.get)})", max(c2_end.values()) < 0.006)
    for period in PERIOD_LABEL:
        o1 = sorted(ORDER, key=lambda a: e[("Case 1", "pue", period, a)])
        o2 = sorted(ORDER, key=lambda a: e[("Case 2", "pue", period, a)])
        c.claim("3.1", f"Case 2 PUE change has the same ecoregion order as Case 1, {PERIOD_LABEL[period]}", "yes",
                f"Case 1 {'<'.join(o1)}; Case 2 {'<'.join(o2)}", o1 == o2)
    c2_wue = max(e[("Case 2", "wue", p, a)] for p in PERIOD_LABEL for a in ORDER)
    c2_wue_fac = changes[(changes["case"] == "Case 2") & (changes["metric"] == "wue")].groupby(["facility_id", "period"])["change"].mean().max()
    c.claim("3.1", "Case 2 change in WUE below 0.01 L/kWh everywhere", "< 0.01",
            f"ecoregion max {c2_wue:.4f}; facility max {c2_wue_fac:.4f}", c2_wue < 0.01)
    # Case 5
    for period, rng in (("midcentury", (0.006, 0.007)), ("endcentury", (0.010, 0.012))):
        vals = [e[("Case 5", "pue", period, a)] for a in ORDER]
        ok = all(rng[0] <= rnd(v, 3) <= rng[1] for v in vals)
        c.claim("3.1", f"Case 5 ecoregion change in PUE within {rng[0]:.3f}–{rng[1]:.3f}, {PERIOD_LABEL[period]}",
                f"{rng[0]:.3f}–{rng[1]:.3f}", f"{min(vals):.4f}–{max(vals):.4f}", ok)
    c5_wue = max(e[("Case 5", "wue", p, a)] for p in PERIOD_LABEL for a in ORDER)
    c5_wue_fac = changes[(changes["case"] == "Case 5") & (changes["metric"] == "wue")].groupby(["facility_id", "period"])["change"].mean().max()
    c.claim("3.1", "Case 5 change in WUE below 0.03 L/kWh everywhere", "< 0.03",
            f"ecoregion max {c5_wue:.4f}; facility max {c5_wue_fac:.4f}", c5_wue < 0.03)

    # 3.2
    hp = hist[(hist["metric"] == "pue") & (hist["case"] == "Case 1")]["value"]
    c.range("3.2", "Facility range of historical Case 1 PUE", (1.125, 1.139), hp.min(), hp.max(), 3)
    fp = changes[(changes["case"] == "Case 1") & (changes["metric"] == "pue") & (changes["period"] == "midcentury")].groupby("facility_id")["change"].mean()
    c.range("3.2", "Facility range of Case 1 change in PUE, 2045–2074", (0.011, 0.023), fp.min(), fp.max(), 3)

    # 3.3 monthly, mid-century, ecoregion mean within each model
    mo = data["month"].merge(fac[["facility_id", "abbr"]], on="facility_id")
    mo = mo[mo["period"].isin(["historical", "midcentury"])]
    mw = mo.pivot_table(index=["facility_id", "abbr", "gcm", "case", "month"], columns="period", values=["pue_mean", "wue_mean"])
    mch = pd.DataFrame({"pue": mw[("pue_mean", "midcentury")] - mw[("pue_mean", "historical")],
                        "wue": mw[("wue_mean", "midcentury")] - mw[("wue_mean", "historical")]}).reset_index()
    mch = mch.melt(id_vars=["facility_id", "abbr", "gcm", "case", "month"], var_name="metric", value_name="change")
    me = eco_stats(mch, ["case", "metric", "month"]).set_index(["case", "metric", "abbr", "month"])["mean"]
    for abbr, drafts in (("WV", {7: 0.060, 8: 0.066, 10: 0.013, 11: 0.015}), ("CP", {6: 0.043, 7: 0.032, 9: 0.049})):
        for mon, v in drafts.items():
            c.value("3.3", f"Case 1 monthly PUE change, {abbr} {MONTHS[mon - 1]}", v, me[("Case 1", "pue", abbr, mon)], 3)

    def peak(case: str, metric: str, abbr: str) -> tuple[int, float]:
        s = me.loc[(case, metric, abbr)]
        return int(s.idxmax()), float(s.max())

    for abbr, (mon, v) in {"WV": (8, 0.82), "CP": (9, 0.63)}.items():
        pm, pv = peak("Case 1", "wue", abbr)
        c.claim("3.3", f"Case 1 WUE peak, {abbr}", f"{MONTHS[mon - 1]} {v:.2f}", f"{MONTHS[pm - 1]} {pv:.3f}",
                pm == mon and rnd(pv, 2) == v)
    for abbr in ("CP", "WV"):
        pm, pv = peak("Case 5", "pue", abbr)
        c.claim("3.3", f"Case 5 PUE peak, {abbr}", "Jul–Aug, about 0.011", f"{MONTHS[pm - 1]} {pv:.4f}",
                pm in (7, 8) and rnd(pv, 3) == 0.011)
    for abbr, v in {"WV": 0.011, "CP": 0.006}.items():
        c.value("3.3", f"Case 2 monthly PUE change, {abbr} Aug", v, me[("Case 2", "pue", abbr, 8)], 3)

    # 3.4
    t = data["tmax"]
    tw = t.groupby(["facility_id", "gcm", "period"])["tasmax_c"].mean().unstack("period")
    tch = (tw["midcentury"] - tw["historical"]).groupby(level="facility_id").mean().rename("tmax").reset_index()
    tch = tch.merge(fac[["facility_id", "abbr"]], on="facility_id")
    for group, rng in ((["ECSF", "BM", "CP"], (2.7, 2.8)), (["KM", "WV"], (2.3, 2.5))):
        s = tch[tch["abbr"].isin(group)]["tmax"]
        detail = "; ".join(f"{a} {tch.loc[tch['abbr'] == a, 'tmax'].min():.2f}–{tch.loc[tch['abbr'] == a, 'tmax'].max():.2f}" for a in group)
        ok = rnd(s.min(), 1) == rng[0] and rnd(s.max(), 1) == rng[1]
        c.claim("3.4", f"Facility range of change in annual mean Tmax, {', '.join(group)} (°C)", f"{rng[0]:.1f}–{rng[1]:.1f}",
                f"{s.min():.2f}–{s.max():.2f} ({detail})", ok)
    wb = hist[hist["metric"] == "wetbulb_annual"][["facility_id", "abbr", "value"]].rename(columns={"value": "wb"})
    w_wv = wb[wb["abbr"] == "WV"]["wb"]
    c.range("3.4", "Historical annual afternoon wet-bulb, WV facilities (°C)", (10.7, 11.7), w_wv.min(), w_wv.max(), 1)
    for abbr, v in {"BM": 8.1, "ECSF": 6.8}.items():
        s = wb[wb["abbr"] == abbr]["wb"]
        c.claim("3.4", f"Historical annual afternoon wet-bulb, {abbr} (°C)", f"{v:.1f}",
                f"{s.min():.2f}–{s.max():.2f}", rnd(s.min(), 1) == v and rnd(s.max(), 1) == v)
    gap = wb[(wb["wb"].map(lambda x: rnd(x, 1)) > 8.1) & (wb["wb"].map(lambda x: rnd(x, 1)) < 9.8)]
    above = wb[wb["wb"] > 8.5]["wb"].min()
    below = wb[wb["wb"] < 9.0]["wb"].max()
    c.claim("3.4", "No facility between 8.1 and 9.8 °C wet-bulb", "none", f"{len(gap)} facilities; gap {below:.2f}–{above:.2f}", len(gap) == 0)
    sel = fp.reset_index().merge(fac, on="facility_id").merge(wb[["facility_id", "wb"]], on="facility_id")
    sel_wv = sel[sel["abbr"] == "WV"]
    lowest_pue = sel_wv.loc[sel_wv["change"].idxmin()]
    lowest_wb = sel_wv.loc[sel_wv["wb"].idxmin()]
    same_cell = sel_wv[sel_wv["grid_cell"] == lowest_pue["grid_cell"]]
    c.claim("3.4", "WV facility with smallest Case 1 PUE change has the lowest WV wet-bulb", "yes (PUE change 0.014)",
            f"id {lowest_pue['facility_id']} ({lowest_pue['facility_name'] if pd.notna(lowest_pue['facility_name']) else 'unnamed'}), "
            f"{lowest_pue['lat']:.4f} N {abs(lowest_pue['lon']):.4f} W, cell {lowest_pue['grid_cell']} "
            f"({len(same_cell)} facilities in cell); PUE change {lowest_pue['change']:.4f}; wet-bulb {lowest_pue['wb']:.2f} "
            f"(WV min {lowest_wb['wb']:.2f}, cell {lowest_wb['grid_cell']})",
            lowest_pue["grid_cell"] == lowest_wb["grid_cell"])

    # Methods
    n_cells = fac["grid_cell"].nunique()
    c.claim("Methods", "Distinct LOCA2 grid cells", "22", str(n_cells), n_cells == 22)
    counts = {name: len([f for f in dataclasses.fields(cfg) if f.name != "p_atm_pa"])
              for name, cfg in (("Case 1", Case1Config), ("Case 2", Case2Config), ("Case 5", Case5Config))}
    c.claim("Methods", "Equipment parameters at Table B.1 midpoints", "–",
            "; ".join(f"{k} {v}" for k, v in counts.items()), True)
    return pd.DataFrame(c.rows), agree


def saturation(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    mo = data["month"].merge(data["fac"][["facility_id", "abbr"]], on="facility_id")
    mo = mo[(mo["case"] == "Case 1") & mo["abbr"].isin(["CP", "WV"])]
    eco = mo.groupby(["abbr", "period", "gcm", "month"])["pue_mean"].mean().groupby(level=["abbr", "period", "month"]).mean()
    out = eco.unstack("period")[["historical", "midcentury"]].reset_index()
    out["change"] = out["midcentury"] - out["historical"]
    out["month_name"] = out["month"].map(lambda m: MONTHS[m - 1])
    return out


# ---------------------------------------------------------------------------


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out = root / "docs/results-figures"
    data = load(root)
    changes = per_model_changes(data["cache"], data["fac"])
    hist = historical(data["cache"], data["fac"])

    t1 = build_table(changes, hist, data["fac"], ["Case 1"]).drop(columns="case")
    t1.to_csv(out / "table1.csv", index=False)
    build_table(changes, hist, data["fac"], ["Case 2", "Case 5"]).to_csv(out / "table1_cases.csv", index=False)
    check, agree = number_check(data, changes, hist)
    check.to_csv(out / "results_number_check.csv", index=False)
    sat = saturation(data)
    sat.to_csv(out / "saturation_check.csv", index=False)

    pd.set_option("display.width", 250, "display.max_colwidth", 200)
    print("=== Table 1 (Case 1) ===")
    print(t1.to_string(index=False))
    print("\n=== Mismatches ===")
    mism = check[check["match"] == "no"]
    print(mism.to_string(index=False) if len(mism) else "none")
    print("\n=== Model agreement: models (of 20) with WV change > CP change ===")
    print(agree.pivot_table(index=["metric", "case"], columns="period", values="n_models_wv_gt_cp")[list(PERIOD_LABEL)].to_string())
    print("\n=== Parameters at Table B.1 midpoints ===")
    print(check.loc[check["item"].str.startswith("Equipment"), "computed_value"].iloc[0])
    print("\n=== Saturation: Case 1 monthly PUE, ecoregion mean, ensemble mean ===")
    for abbr, g in sat.groupby("abbr"):
        print(f"{abbr}: " + "  ".join(f"{r.month_name} {r.historical:.3f}/{r.change:+.3f}" for r in g.itertuples())
              + f"   highest historical: {g.loc[g['historical'].idxmax(), 'month_name']}")
    print(f"\nWrote table1.csv, table1_cases.csv, results_number_check.csv, saturation_check.csv to {out}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Table 1 and the Results number check.")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
