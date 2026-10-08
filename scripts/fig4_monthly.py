"""Figure 4 of the heat risk chapter: when in the year the change happens
(docs/tasks/result-figures.md).

PUE only, mid-century only. Two panels on a shared y-axis: (a) Columbia
Plateau and (b) Willamette Valley, each showing the monthly change in the
mean of that ecoregion's facilities for Cases 1 (solid), 2 (dashed) and 5
(dotted). Line color is the region (Okabe-Ito bluish green for CP, reddish
purple for WV, as in Figure 5; blue and vermillion are kept for the periods
in the other figures).

Each line is the mean across the 20 models; its band is the interquartile
range across models of the ecoregion mean.

Inputs (nothing re-run):
- ``outputs/ensemble_pue_wue_facilities/analysis/facility_month_summary.csv``
  (``build_pue_wue_ensemble_summary.py``): monthly mean of daily PUE/WUE per
  facility, model, case and period.
- ``.../analysis/fig_show_everything/facility_grid_cells.csv``: facility
  Level III ecoregion.

Outputs:
- ``docs/results-figures/fig4_monthly.{png,pdf}`` and ``fig4_monthly_caption.txt``.
- ``outputs/ensemble_pue_wue_facilities/analysis/fig4_monthly/``:
  ``monthly_changes.csv`` (plotted values and bands), ``peaks.csv``,
  ``run_manifest.csv``.

Assumptions log
----------------
- **Monthly mean** is the mean of all daily values in that calendar month
  over the 30 (historical, mid-century) or 26 (end-of-century) years, as
  stored in ``facility_month_summary.csv``. Every year contributes the same
  number of days to a month except February in leap years, so this equals
  the mean of yearly monthly means to within that one day.
- **Change** is future minus historical per facility, model and month. The
  statewide mean weights the 109 facilities equally; so do the CP and WV
  means.
- **Peak month** is the month with the largest ensemble-mean change.
- **Plotted content.** The task specified statewide PUE and WUE panels for
  both future periods plus Case 1 CP/WV panels. At Laura's request
  (2026-10-07) the figure shows PUE, mid-century (``PLOT_PERIOD``), CP and WV
  only, all three cases. Statewide, WUE and end-of-century values are still
  computed and written to the CSVs and the printed summary.
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
N_FACILITIES = 109
PERIODS = {"historical": (1985, 2014), "midcentury": (2045, 2074), "endcentury": (2075, 2100)}
FUTURE_PERIODS = ["midcentury", "endcentury"]
PERIOD_LABEL = {"midcentury": "Mid-century (2045–2074)", "endcentury": "End-of-century (2075–2100)"}
PLOT_PERIOD = "midcentury"
CASES = {"ae-chiller": "Case 1", "we-chiller": "Case 2", "chiller-only": "Case 5"}
CASE_STYLE = {"Case 1": "-", "Case 2": "--", "Case 5": ":"}
REGIONS = {"Columbia Plateau": "CP", "Willamette Valley": "WV"}
REGION_COLOR = {"CP": figstyle.BLUISH_GREEN, "WV": figstyle.REDDISH_PURPLE}
REGION_TITLE = {"CP": "Columbia Plateau", "WV": "Willamette Valley"}
PLOT_METRIC = "pue"
MONTH_TICKS = list("JFMAMJJASOND")
BAND_ALPHA = 0.3
LINE_WIDTH = 1.3
FIG_HEIGHT_IN = 3.3
BOUNDS = {"pue_mean": (1.0, 3.0), "wue_mean": (0.0, 20.0)}


def _analysis_dir(root: Path) -> Path:
    return root / "outputs/ensemble_pue_wue_facilities/analysis"


def _load(root: Path) -> pd.DataFrame:
    month = pd.read_csv(_analysis_dir(root) / "facility_month_summary.csv")
    for col, (lo, hi) in BOUNDS.items():
        v = month[col]
        if v.isna().any() or (v < lo).any() or (v > hi).any():
            raise ValueError(f"{col}: values outside [{lo}, {hi}] or NaN (range {v.min()}-{v.max()})")
    eco = pd.read_csv(_analysis_dir(root) / "fig_show_everything" / "facility_grid_cells.csv", usecols=["facility_id", "ecoregion"])
    month = month.merge(eco, on="facility_id", validate="many_to_one")
    month["case"] = month["archetype"].map(CASES)
    counts = month.groupby(["case", "period", "month"]).agg(fac=("facility_id", "nunique"), gcm=("gcm", "nunique"))
    assert (counts["fac"] == N_FACILITIES).all() and (counts["gcm"] == N_MODELS).all(), "Missing facilities or models"
    return month


def _changes(month: pd.DataFrame) -> pd.DataFrame:
    """Per facility, model, case, month, metric and future period: future minus historical."""
    long = month.melt(id_vars=["facility_id", "ecoregion", "gcm", "case", "period", "month"],
                      value_vars=["pue_mean", "wue_mean"], var_name="metric", value_name="value")
    long["metric"] = long["metric"].str.replace("_mean", "")
    wide = long.pivot_table(index=["facility_id", "ecoregion", "gcm", "case", "month", "metric"], columns="period", values="value")
    out = pd.concat([(wide[p] - wide["historical"]).rename("change").reset_index().assign(period=p) for p in FUTURE_PERIODS],
                    ignore_index=True)
    assert out["change"].notna().all()
    return out


def _aggregate(changes: pd.DataFrame) -> pd.DataFrame:
    """Statewide and CP/WV means (all cases): ensemble mean and IQR across models of the per-model aggregate."""
    keys = ["metric", "case", "period", "month"]
    state = changes.groupby(keys + ["gcm"])["change"].mean().reset_index().assign(region="Statewide")
    regional = changes[changes["ecoregion"].isin(REGIONS)].groupby(keys + ["ecoregion", "gcm"])["change"].mean().reset_index()
    regional = regional.assign(region=regional["ecoregion"].map(REGIONS)).drop(columns="ecoregion")
    per_model = pd.concat([state, regional], ignore_index=True)
    g = per_model.groupby(["region"] + keys)["change"]
    out = g.mean().rename("mean").to_frame()
    out["p25"], out["p75"] = g.quantile(0.25), g.quantile(0.75)
    out["n_models"] = g.size()
    return out.reset_index()


def _peaks(agg: pd.DataFrame) -> pd.DataFrame:
    idx = agg.groupby(["region", "metric", "case", "period"])["mean"].idxmax()
    return agg.loc[idx, ["region", "metric", "case", "period", "month", "mean", "p25", "p75"]].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------


def _line(ax, d: pd.DataFrame, color: str, style: str) -> None:
    d = d.sort_values("month")
    ax.fill_between(d["month"], d["p25"], d["p75"], color=color, alpha=BAND_ALPHA, linewidth=0)
    ax.plot(d["month"], d["mean"], color=color, linestyle=style, linewidth=LINE_WIDTH)


def _style_axis(ax: plt.Axes) -> None:
    ax.axhline(0, color=figstyle.ZERO_LINE_COLOR, linewidth=figstyle.ZERO_LINE_WIDTH, zorder=0.5)
    ax.set_xlim(0.6, 12.4)
    ax.set_xticks(range(1, 13), MONTH_TICKS)


def plot(agg: pd.DataFrame, n_facilities: dict[str, int], path: Path) -> None:
    figstyle.apply()
    fig, axes = plt.subplots(1, 2, figsize=(figstyle.FIG_WIDTH_IN, FIG_HEIGHT_IN), sharex=True, sharey=True,
                             gridspec_kw={"left": 0.115, "right": 0.985, "top": 0.91, "bottom": 0.27, "wspace": 0.08})
    d = agg[(agg["period"] == PLOT_PERIOD) & (agg["metric"] == PLOT_METRIC)]
    for ax, letter, region in zip(axes, "ab", REGION_COLOR):
        for case, style in CASE_STYLE.items():
            _line(ax, d[(d["region"] == region) & (d["case"] == case)], REGION_COLOR[region], style)
        ax.set_title(f"({letter}) {REGION_TITLE[region]} (n = {n_facilities[region]})", loc="left")
        _style_axis(ax)
        ax.set_xlabel("Month")
    axes[0].set_ylabel("Change in monthly mean PUE")

    handles = [Line2D([], [], color=figstyle.INK, linestyle=s, linewidth=LINE_WIDTH) for s in CASE_STYLE.values()]
    fig.legend(handles, list(CASE_STYLE), loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, handlelength=2.6,
               columnspacing=2.0, borderaxespad=0.3)
    for ext in ("png", "pdf"):
        fig.savefig(path.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)


def caption(n_facilities: dict[str, int]) -> str:
    return (
        "Monthly change from 1985–2014 to mid-century (2045–2074) in mean power usage effectiveness (PUE) for "
        f"facilities in (a) the Columbia Plateau (n = {n_facilities['CP']}) and (b) the Willamette Valley "
        f"(n = {n_facilities['WV']}), for cooling-system Cases 1 (solid), 2 (dashed) and 5 (dotted). Lines are the mean "
        "of the ecoregion's facilities, averaged over 20 climate models (SSP3-7.0); shading is the interquartile range "
        "(25th–75th percentile) across models. Both panels share the y-axis."
    )


def print_summary(peaks: pd.DataFrame) -> None:
    print("\n=== Peak month (largest ensemble-mean change) [IQR across models] ===")
    month_name = dict(enumerate(["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"], 1))
    for _, r in peaks.sort_values(["region", "metric", "case", "period"], ascending=[False, True, True, False]).iterrows():
        fmt = "{:.4f}" if r["metric"] == "pue" else "{:.3f}"
        print(f"   {r['region']:9s} {r['metric'].upper()} {r['case']} {r['period']:10s} {month_name[r['month']]}  "
              f"+{fmt.format(r['mean'])} [{fmt.format(r['p25'])}, {fmt.format(r['p75'])}]")


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out_dir = _analysis_dir(root) / "fig4_monthly"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = root / "docs" / "results-figures"

    month = _load(root)
    n_facilities = {abbr: month.loc[month["ecoregion"] == name, "facility_id"].nunique() for name, abbr in REGIONS.items()}
    agg = _aggregate(_changes(month))
    assert (agg["n_models"] == N_MODELS).all()
    peaks = _peaks(agg)
    agg.sort_values(["region", "metric", "case", "period", "month"]).to_csv(out_dir / "monthly_changes.csv", index=False)
    peaks.to_csv(out_dir / "peaks.csv", index=False)
    pd.DataFrame([{
        "n_models": N_MODELS, "n_facilities": N_FACILITIES, "scenario": "ssp370",
        "periods": "; ".join(f"{k}={a}-{b}" for k, (a, b) in PERIODS.items()),
        "cases": "; ".join(f"{v}={k}" for k, v in CASES.items()), "band": "iqr", "plotted": "PUE, midcentury, CP and WV, Cases 1/2/5",
    }]).to_csv(out_dir / "run_manifest.csv", index=False)

    plot(agg, n_facilities, fig_dir / "fig4_monthly.png")
    (fig_dir / "fig4_monthly_caption.txt").write_text(caption(n_facilities) + "\n")
    print_summary(peaks)
    print(f"\nWrote {fig_dir / 'fig4_monthly.png'} (+ .pdf) and CSVs under {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Figure 4 (monthly change) of the heat chapter.")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
