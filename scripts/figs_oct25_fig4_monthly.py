"""Figure 4 redraw (reviewer comment, October 2026): one row per cooling
case, one column per ecoregion, seasonal cycles for all three periods.

Layout: 3 rows (Cases 1, 2, 5) x 2 columns (Columbia Plateau, Willamette
Valley); 6.5 in wide (``figstyle.FIG_WIDTH_IN``), 300 dpi, figstyle fonts,
as ``scripts/fig4_monthly.py``. Each panel shows absolute monthly mean PUE
for 1985-2014, 2045-2074 and 2075-2100: line = mean across the 20 models of
the ecoregion facility mean; band = interquartile range across models (the
aggregation of the current Figure 4 and of ``draft3/monthly.csv``). The
y-axis is shared within a row and differs between rows.

Two versions:
- ``fig4_monthly_3x2.{png,pdf}``: absolute values, with an inset per panel
  showing the mid-century change (2045-2074 minus 1985-2014) and its IQR;
  inset y-axes are shared within a row.
- ``fig4_monthly_3x2_change.{png,pdf}`` (alternative): mid- and end-century
  change only, same layout.

Inputs: ``analysis/facility_month_summary.csv`` (the changes' IQR needs
per-model differences, which ``jobB_monthly_pue_wue.csv`` does not hold) and
``facility_grid_cells.csv``. Checks, written to ``fig4_check.csv``:
absolute values reproduce ``server_jobs_oct23/jobB_monthly_pue_wue.csv``;
mid-century changes reproduce ``draft3/monthly.csv``. Plotted values go to
``fig4_plotted_values.csv``.

Assumptions log
----------------
- **Period colors** follow Figure 2 (Okabe-Ito blue for 2045-2074,
  vermillion for 2075-2100), with black for 1985-2014.
- **Inset placement** is the upper left of each panel, where the absolute
  curves are lowest (January-May) in every case and ecoregion; each
  row's main y-limit is raised so no band (interpolated between months)
  reaches the inset or its month labels.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

from climate_risk_dc import figstyle
from figs_oct25_common import OUT_DIR, apply_figstyle, side_by_side

N_MODELS = 20
CASES = {"ae-chiller": "Case 1", "we-chiller": "Case 2", "chiller-only": "Case 5"}
REGIONS = {"Columbia Plateau": "CP", "Willamette Valley": "WV"}
PERIOD_LABEL = {"historical": "1985–2014", "midcentury": "2045–2074", "endcentury": "2075–2100"}
PERIOD_COLOR = {"historical": figstyle.INK, "midcentury": figstyle.BLUE, "endcentury": figstyle.VERMILLION}
MONTH_TICKS = list("JFMAMJJASOND")
BAND_ALPHA = 0.22
LINE_WIDTH = 1.2
FIG_HEIGHT_IN = 7.6
INSET_RECT = [0.035, 0.49, 0.40, 0.40]  # axes fraction: x, y, w, h
INSET_PT = 7.0
INSET_TICK_CLEAR = 0.08  # axes fraction below the inset frame taken by its month labels
CHECK_PT = {"CP": 7, "WV": 8}  # months listed in the check
IQR_LABEL = "IQR across 20 models"


def _aggregate(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    a = root / "outputs/ensemble_pue_wue_facilities/analysis"
    mo = pd.read_csv(a / "facility_month_summary.csv")
    assert mo["pue_mean"].between(1.0, 3.0).all()
    eco = pd.read_csv(a / "fig_show_everything/facility_grid_cells.csv", usecols=["facility_id", "ecoregion"])
    mo = mo.merge(eco, on="facility_id", validate="many_to_one")
    mo = mo[mo["ecoregion"].isin(REGIONS)].assign(region=lambda d: d["ecoregion"].map(REGIONS),
                                                    case=lambda d: d["archetype"].map(CASES))
    w = mo.pivot_table(index=["case", "region", "facility_id", "gcm", "month"], columns="period", values="pue_mean")
    parts = {"absolute": {p: w[p] for p in PERIOD_LABEL},
             "change": {p: w[p] - w["historical"] for p in ("midcentury", "endcentury")}}
    rows = []
    for kind, d in parts.items():
        for period, s in d.items():
            pm = s.groupby(level=["case", "region", "gcm", "month"]).mean()
            g = pm.groupby(level=["case", "region", "month"])
            assert (g.size() == N_MODELS).all()
            rows.append(pd.DataFrame({"mean": g.mean(), "p25": g.quantile(0.25), "p75": g.quantile(0.75)})
                        .reset_index().assign(kind=kind, period=period))
    agg = pd.concat(rows, ignore_index=True)
    n = mo.groupby("region")["facility_id"].nunique().to_dict()
    assert n == {"CP": 61, "WV": 33}, n
    return agg, pd.Series(n)


def _check(root: Path, agg: pd.DataFrame) -> pd.DataFrame:
    rows = []
    # absolute vs jobB (raw mean, and formatted strings for mean/p25/p75)
    b = pd.read_csv(root / "docs/results-figures/server_jobs_oct23/jobB_monthly_pue_wue.csv")
    b = b[b["variable"] == "PUE"].assign(period=lambda d: d["period"].map({v: k for k, v in PERIOD_LABEL.items()}))
    m = agg[agg["kind"] == "absolute"].merge(b, left_on=["case", "region", "period", "month"],
                                             right_on=["case", "ecoregion", "period", "month"], suffixes=("", "_b"))
    assert len(m) == 3 * 2 * 3 * 12
    rows.append({"check": "absolute mean vs jobB mean_raw", "n": len(m), "max_abs_diff": (m["mean"] - m["mean_raw"]).abs().max()})
    for c in ("mean", "p25", "p75"):
        rows.append({"check": f"absolute {c} vs jobB {c} (3 dp)", "n": len(m),
                     "n_match": int((m[c].map(lambda x: f"{x:.3f}") == m[f"{c}_b"].map(lambda x: f"{x:.3f}")).sum()),
                     "max_abs_diff": (m[c] - m[f"{c}_b"]).abs().max()})
    # mid-century change vs draft3/monthly.csv
    d3 = pd.read_csv(root / "docs/results-figures/draft3/monthly.csv")
    d3 = d3[d3["quantity"] == "d_pue_midcentury"]
    c = agg[(agg["kind"] == "change") & (agg["period"] == "midcentury")].merge(d3, on=["region", "case", "month"], suffixes=("", "_d3"))
    assert len(c) == 3 * 2 * 12
    for col in ("mean", "p25", "p75"):
        dec = c[f"{col}_d3"].map(lambda v: 4 if abs(v) < 0.01 else 3)
        same = [round(x, int(k)) == round(y, int(k)) for x, y, k in zip(c[col], c[f"{col}_d3"], dec)]
        rows.append({"check": f"mid-century change {col} vs draft3/monthly.csv (draft 3 decimals)", "n": len(c),
                     "n_match": int(sum(same)), "max_abs_diff": (c[col] - c[f"{col}_d3"]).abs().max()})
    return pd.DataFrame(rows)


def _band(ax, d: pd.DataFrame, color: str, lw: float = LINE_WIDTH) -> None:
    d = d.sort_values("month")
    ax.fill_between(d["month"], d["p25"], d["p75"], color=color, alpha=BAND_ALPHA, linewidth=0)
    ax.plot(d["month"], d["mean"], color=color, linewidth=lw)


def _months(ax, labels: bool, size: float | None = None) -> None:
    ax.set_xlim(0.6, 12.4)
    ax.set_xticks(range(1, 13), MONTH_TICKS if labels else [""] * 12)
    if size:
        ax.tick_params(labelsize=size, length=2, pad=1)


def _grid(n: pd.Series):
    apply_figstyle()
    fig, axes = plt.subplots(3, 2, figsize=(figstyle.FIG_WIDTH_IN, FIG_HEIGHT_IN), sharex=True, sharey="row",
                             gridspec_kw={"left": 0.155, "right": 0.985, "top": 0.955, "bottom": 0.105,
                                          "wspace": 0.07, "hspace": 0.2})
    for j, (region, title) in enumerate((("CP", "Columbia Plateau"), ("WV", "Willamette Valley"))):
        axes[0, j].set_title(f"{title} (n = {n[region]})", loc="center", pad=14, fontweight="bold")
    for i, case in enumerate(CASES.values()):
        fig.text(0.012, np.mean([axes[i, 0].get_position().y0, axes[i, 0].get_position().y1]), case,
                 rotation=90, ha="left", va="center", fontsize=figstyle.LABEL_PT, fontweight="bold")
    for k, ax in enumerate(axes.flat):
        ax.set_title(f"({'abcdef'[k]})", loc="left", pad=3)
        _months(ax, labels=True)
    for ax in axes[-1]:
        ax.set_xlabel("Month")
    return fig, axes


def plot_absolute(agg: pd.DataFrame, n: pd.Series, path: Path) -> None:
    fig, axes = _grid(n)
    all_insets = []
    ab, ch = agg[agg["kind"] == "absolute"], agg[(agg["kind"] == "change") & (agg["period"] == "midcentury")]
    for i, case in enumerate(CASES.values()):
        insets = []
        for j, region in enumerate(REGIONS.values()):
            ax = axes[i, j]
            for period, color in PERIOD_COLOR.items():
                _band(ax, ab[(ab["case"] == case) & (ab["region"] == region) & (ab["period"] == period)], color)
            ins = ax.inset_axes(INSET_RECT)
            ins.axhline(0, color=figstyle.ZERO_LINE_COLOR, linewidth=figstyle.ZERO_LINE_WIDTH, zorder=0.5)
            _band(ins, ch[(ch["case"] == case) & (ch["region"] == region)], PERIOD_COLOR["midcentury"], lw=1.0)
            _months(ins, labels=True, size=INSET_PT - 1)
            ins.set_xticks(range(1, 13, 2), MONTH_TICKS[::2])
            # y tick labels inside the inset, so they do not collide with the main axis labels
            ins.tick_params(axis="y", labelsize=INSET_PT - 1, direction="in", pad=-2, length=2)
            ins.yaxis.set_major_locator(MaxNLocator(nbins=3, min_n_ticks=2))
            for s in ins.spines.values():
                s.set_linewidth(0.6)
            ins.set_title("Change, 2045–2074", fontsize=INSET_PT, pad=1.5)
            ins.patch.set_alpha(0.9)
            insets.append(ins)
        # shared inset y-range within the row
        lo = min(x.get_ylim()[0] for x in insets)
        hi = max(x.get_ylim()[1] for x in insets)
        for x in insets:
            x.set_ylim(lo, hi)
        all_insets += insets
        # main y-range: data below the inset's bottom edge
        sub = ab[ab["case"] == case]
        dmin, dmax = sub["p25"].min(), sub["p75"].max()
        pad = 0.04 * (dmax - dmin)
        y0 = dmin - pad
        # highest band value (interpolated between months) under the inset, which spans x <= x_right
        x_right = 0.6 + (INSET_RECT[0] + INSET_RECT[2]) * 11.8
        xs = np.linspace(1, x_right, 200)
        left = max(np.interp(xs, g["month"], g["p75"]).max()
                   for _, g in sub.sort_values("month").groupby(["region", "period"]))
        span = max((dmax + pad - y0), (left - y0) / (INSET_RECT[1] - INSET_TICK_CLEAR))
        axes[i, 0].set_ylim(y0, y0 + span)
        axes[i, 0].yaxis.set_major_locator(MaxNLocator(nbins=5))
        axes[i, 0].set_ylabel("Monthly mean PUE")
    handles = [Line2D([], [], color=c, linewidth=LINE_WIDTH) for c in PERIOD_COLOR.values()]
    handles.append(Patch(facecolor="grey", alpha=0.35, edgecolor="none"))
    fig.canvas.draw()  # tick labels exist only after limits are final
    for ins in all_insets:
        for lab in ins.get_yticklabels():
            lab.set_horizontalalignment("left")
            lab.set_bbox({"facecolor": "white", "edgecolor": "none", "pad": 0.3, "alpha": 0.8})
    fig.legend(handles, list(PERIOD_LABEL.values()) + [IQR_LABEL], loc="lower center",
               bbox_to_anchor=(0.55, 0.0), ncol=4, handlelength=1.8, columnspacing=1.2, fontsize=figstyle.TICK_PT)
    for ext in ("png", "pdf"):
        fig.savefig(path.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)


def plot_change(agg: pd.DataFrame, n: pd.Series, path: Path) -> None:
    fig, axes = _grid(n)
    ch = agg[agg["kind"] == "change"]
    for i, case in enumerate(CASES.values()):
        for j, region in enumerate(REGIONS.values()):
            ax = axes[i, j]
            ax.axhline(0, color=figstyle.ZERO_LINE_COLOR, linewidth=figstyle.ZERO_LINE_WIDTH, zorder=0.5)
            for period in ("midcentury", "endcentury"):
                _band(ax, ch[(ch["case"] == case) & (ch["region"] == region) & (ch["period"] == period)], PERIOD_COLOR[period])
        axes[i, 0].yaxis.set_major_locator(MaxNLocator(nbins=5))
        axes[i, 0].set_ylabel("Change in monthly\nmean PUE")
    handles = [Line2D([], [], color=PERIOD_COLOR[p], linewidth=LINE_WIDTH) for p in ("midcentury", "endcentury")]
    handles.append(Patch(facecolor="grey", alpha=0.35, edgecolor="none"))
    fig.legend(handles, [f"{PERIOD_LABEL[p]} minus 1985–2014" for p in ("midcentury", "endcentury")]
               + [IQR_LABEL], loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3,
               handlelength=1.5, columnspacing=0.9, handletextpad=0.5, fontsize=figstyle.TICK_PT)
    for ext in ("png", "pdf"):
        fig.savefig(path.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out = root / OUT_DIR
    out.mkdir(parents=True, exist_ok=True)
    agg, n = _aggregate(root)
    agg.sort_values(["kind", "case", "region", "period", "month"]).to_csv(out / "fig4_plotted_values.csv", index=False)
    chk = _check(root, agg)
    chk.to_csv(out / "fig4_check.csv", index=False)
    plot_absolute(agg, n, out / "fig4_monthly_3x2.png")
    plot_change(agg, n, out / "fig4_monthly_3x2_change.png")
    old = root / "docs/results-figures/fig4_monthly.png"
    side_by_side(old, out / "fig4_monthly_3x2.png", out / "preview_fig4_old_vs_new.png")
    side_by_side(old, out / "fig4_monthly_3x2_change.png", out / "preview_fig4_old_vs_change.png")

    with pd.option_context("display.width", 200, "display.max_colwidth", 80):
        print(chk.to_string(index=False))
        for region, month in CHECK_PT.items():
            s = agg[(agg["region"] == region) & (agg["month"] == month)]
            print(f"\n{region}, month {month}:")
            print(s.pivot_table(index="case", columns=["kind", "period"], values="mean").round(4).to_string())


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Figure 4 redraw: 3x2 seasonal cycles by case and ecoregion.")
    p.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return p.parse_args()


if __name__ == "__main__":
    run(parse_args())
