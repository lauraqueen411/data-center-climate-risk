"""Revised Figure 5 of the heat risk chapter: sensitivity of Case 1 PUE to
the supply-air relative humidity upper setpoint (the mechanism test).

Two panels, ecoregions on the x-axis in wet-bulb order:
- (a) historical (1985-2014) annual mean Case 1 PUE
- (b) mid-century change (2045-2074 minus 1985-2014) in annual mean Case 1 PUE

At each ecoregion, one point per setpoint (60%, the run superseded on
2026-10-08, and 77.5%, the Table B.1 midpoint used for all other results):
the ecoregion mean of facilities, averaged over the 20 models, with a thin
vertical bar for the interquartile range across models of the ecoregion
mean.

Inputs (nothing re-run):
- 77.5%: ``outputs/ensemble_pue_wue_facilities/analysis/fig_show_everything/model_values_cache.csv``
- 60%: ``outputs/ensemble_pue_wue_facilities/superseded_rh_up60_2026-10-08/analysis/fig_show_everything/model_values_cache.csv``
- ``.../analysis/fig_show_everything/facility_grid_cells.csv``: ecoregion.

Outputs:
- ``docs/results-figures/fig5_rh_setpoint.{png,pdf}`` and ``fig5_rh_setpoint_caption.txt``.
- ``outputs/ensemble_pue_wue_facilities/analysis/fig5_rh_setpoint/``:
  ``ecoregion_values.csv``, ``run_manifest.csv``.

Assumptions log
----------------
- **The two runs differ only in Case 1 ``rh_up_pct``** (60 vs 77.5). The
  2026-10-08 correction also changed two Case 5 parameters, which do not
  affect Case 1. Both runs used the same weather, models and code otherwise.
- **Values** follow the other figures: per facility and model, the period
  mean of annual means; change is mid-century minus historical per facility
  and model; ecoregion mean of facilities within each model (equal weight
  per record); then the mean and interquartile range across the 20 models.
- **Ecoregion order** is ascending historical annual afternoon wet-bulb, the
  ecoregion mean of facilities, as in Figure 2. Wet-bulb does not depend on
  the setpoint.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from climate_risk_dc import figstyle

N_MODELS = 20
PERIODS = {"historical": (1985, 2014), "midcentury": (2045, 2074)}
ECOREGION_ABBR = {
    "Eastern Cascades Slopes and Foothills": "ECSF",
    "Blue Mountains": "BM",
    "Columbia Plateau": "CP",
    "Klamath Mountains/California High North Coast Range": "KM",
    "Willamette Valley": "WV",
}
RUNS = {  # label -> (cache folder relative to analysis/, color)
    "60%": ("../superseded_rh_up60_2026-10-08/analysis/fig_show_everything", "#9e9e9e"),
    "77.5%": ("fig_show_everything", "#252525"),
}
OFFSET = 0.13
MARKER_SIZE = 4.0
BAR_WIDTH = 1.0
FIG_HEIGHT_IN = 3.2
X_LABEL = "Increasing wet-bulb temperature →"
LEGEND_TITLE = "Supply-air relative humidity upper setpoint:"
BOUNDS = (1.0, 3.0)


def _analysis_dir(root: Path) -> Path:
    return root / "outputs/ensemble_pue_wue_facilities/analysis"


def _per_model(root: Path, folder: str) -> pd.DataFrame:
    """Case 1 PUE per facility and model: historical value and mid-century change."""
    c = pd.read_csv(_analysis_dir(root) / folder / "model_values_cache.csv")
    c = c[(c["metric"] == "pue") & (c["archetype"] == "ae-chiller") & c["period"].isin(PERIODS)]
    if c["value"].isna().any() or not c["value"].between(*BOUNDS).all():
        raise ValueError(f"{folder}: Case 1 PUE outside {BOUNDS} or NaN")
    w = c.pivot_table(index=["facility_id", "gcm"], columns="period", values="value")
    out = pd.concat([w["historical"].rename("value").reset_index().assign(kind="historical"),
                     (w["midcentury"] - w["historical"]).rename("value").reset_index().assign(kind="change")])
    assert (out.groupby(["facility_id", "kind"])["gcm"].nunique() == N_MODELS).all(), f"{folder}: missing models"
    return out


def _ecoregion_order(root: Path, fac: pd.DataFrame) -> list[str]:
    c = pd.read_csv(_analysis_dir(root) / "fig_show_everything/model_values_cache.csv")
    wb = c[(c["metric"] == "wetbulb_annual") & (c["period"] == "historical")].groupby("facility_id")["value"].mean()
    eco = fac.set_index("facility_id")["abbr"]
    return list(wb.groupby(eco).mean().sort_values().index)


def ecoregion_values(root: Path) -> tuple[pd.DataFrame, list[str]]:
    fac = pd.read_csv(_analysis_dir(root) / "fig_show_everything/facility_grid_cells.csv", usecols=["facility_id", "ecoregion"])
    fac["abbr"] = fac["ecoregion"].map(ECOREGION_ABBR)
    parts = []
    for label, (folder, _) in RUNS.items():
        pm = _per_model(root, folder).merge(fac[["facility_id", "abbr"]], on="facility_id")
        m = pm.groupby(["kind", "abbr", "gcm"])["value"].mean()
        g = m.groupby(level=["kind", "abbr"])
        parts.append(pd.DataFrame({"mean": g.mean(), "p25": g.quantile(0.25), "p75": g.quantile(0.75),
                                   "n_models": g.size()}).reset_index().assign(rh_up=label))
    out = pd.concat(parts, ignore_index=True)
    assert (out["n_models"] == N_MODELS).all()
    n = fac["abbr"].value_counts()
    out["n_facilities"] = out["abbr"].map(n)
    return out, _ecoregion_order(root, fac)


def plot(vals: pd.DataFrame, order: list[str], path: Path) -> None:
    figstyle.apply()
    fig, axes = plt.subplots(1, 2, figsize=(figstyle.FIG_WIDTH_IN, FIG_HEIGHT_IN),
                             gridspec_kw={"left": 0.115, "right": 0.985, "top": 0.91, "bottom": 0.255, "wspace": 0.33})
    titles = {"historical": "(a) Historical (1985–2014)", "change": "(b) Mid-century change"}
    ylabels = {"historical": "Annual mean PUE", "change": "Change in annual mean PUE"}
    x = {abbr: i for i, abbr in enumerate(order)}
    for ax, kind in zip(axes, ("historical", "change")):
        for k, (label, (_, color)) in enumerate(RUNS.items()):
            d = vals[(vals["kind"] == kind) & (vals["rh_up"] == label)].set_index("abbr").loc[order]
            xs = [x[a] + (k - 0.5) * 2 * OFFSET for a in order]
            ax.vlines(xs, d["p25"], d["p75"], color=color, linewidth=BAR_WIDTH, zorder=2)
            ax.plot(xs, d["mean"], linestyle="none", marker="o", markersize=MARKER_SIZE, markeredgewidth=0,
                    color=color, zorder=3)
        if kind == "change":
            ax.axhline(0, color=figstyle.ZERO_LINE_COLOR, linewidth=figstyle.ZERO_LINE_WIDTH, zorder=1)
        ax.set_xticks(range(len(order)), order)
        ax.set_xlim(-0.5, len(order) - 0.5)
        ax.set_xlabel(X_LABEL)
        ax.set_ylabel(ylabels[kind])
        ax.set_title(titles[kind], loc="left")
    handles = [Line2D([], [], linestyle="none")] + [
        Line2D([], [], linestyle="none", marker="o", markersize=MARKER_SIZE + 1, markeredgewidth=0, color=c)
        for _, c in RUNS.values()]
    leg = fig.legend(handles, [LEGEND_TITLE] + list(RUNS), loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3,
                     handlelength=0.8, handletextpad=0.3, columnspacing=1.2, borderaxespad=0.3)
    leg.legend_handles[0].set_visible(False)
    for ext in ("png", "pdf"):
        fig.savefig(path.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)


def caption(vals: pd.DataFrame) -> str:
    n = vals.drop_duplicates("abbr").set_index("abbr")["n_facilities"]
    return (
        "Case 1 annual mean power usage effectiveness (PUE) with the supply-air relative humidity upper setpoint at "
        "60% (grey) and 77.5% (black), the lower bound and midpoint of the Lei and Masanet (2022) range. (a) Historical "
        "(1985–2014) value. (b) Mid-century change (2045–2074 minus 1985–2014). Ecoregions are ordered from left to "
        "right by historical annual mean afternoon wet-bulb temperature: Eastern Cascades Slopes and Foothills "
        f"(ECSF; n = {n['ECSF']}), Blue Mountains (BM; n = {n['BM']}), Columbia Plateau (CP; n = {n['CP']}), Klamath "
        f"Mountains (KM; n = {n['KM']}) and Willamette Valley (WV; n = {n['WV']}). Points are the ecoregion mean of "
        "facilities averaged over 20 climate models (SSP3-7.0); bars are the interquartile range (25th–75th "
        "percentile) of the ecoregion mean across models. All other results use the 77.5% setpoint."
    )


def print_summary(vals: pd.DataFrame, order: list[str]) -> None:
    w = vals.pivot_table(index=["kind", "abbr"], columns="rh_up", values="mean")
    w["difference (60 - 77.5)"] = w["60%"] - w["77.5%"]
    for kind in ("historical", "change"):
        print(f"\n=== {kind}: ecoregion mean Case 1 PUE by setpoint ===")
        print(w.loc[kind].loc[order].round(4).to_string())
    for label in RUNS:
        c = w.loc["change", label]
        print(f"   {label}: WV minus CP change = {c['WV'] - c['CP']:+.4f}")


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out_dir = _analysis_dir(root) / "fig5_rh_setpoint"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = root / "docs" / "results-figures"
    vals, order = ecoregion_values(root)
    vals.sort_values(["kind", "abbr", "rh_up"]).to_csv(out_dir / "ecoregion_values.csv", index=False)
    pd.DataFrame([{"n_models": N_MODELS, "scenario": "ssp370", "case": "Case 1 = ae-chiller",
                   "rh_up_runs": "60% (superseded_rh_up60_2026-10-08); 77.5% (current)",
                   "periods": "; ".join(f"{k}={a}-{b}" for k, (a, b) in PERIODS.items()),
                   "band": "iqr", "ecoregion_order": "; ".join(order)}]).to_csv(out_dir / "run_manifest.csv", index=False)
    plot(vals, order, fig_dir / "fig5_rh_setpoint.png")
    (fig_dir / "fig5_rh_setpoint_caption.txt").write_text(caption(vals) + "\n")
    print_summary(vals, order)
    print(f"\nWrote {fig_dir / 'fig5_rh_setpoint.png'} (+ .pdf) and CSVs under {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the revised Figure 5 (RH setpoint sensitivity).")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
