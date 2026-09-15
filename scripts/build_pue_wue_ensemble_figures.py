"""Phase 2/3 of docs/tasks/pue_wue_analysis.md: figures and tables from the PUE/WUE ensemble summary.

Reads the summary tables ``scripts/build_pue_wue_ensemble_summary.py`` wrote
under ``outputs/ensemble_pue_wue_facilities/analysis/`` and produces
Figures 1-6 and 8 (Figure 7 needs the heat-index ensemble, built separately
by ``scripts/aggregate_heat_ensemble.py``, and paired in a follow-up step;
Figure 9 is contingent on a sub-region decision not yet made, per the task
doc -- not built speculatively) plus Tables 1-4, per
``docs/data-center-chapter-outline.md``'s candidate figure/table list.

Archetype labels use the corrected Case numbers (Phase 0 found "Case 3" in
the task doc and outline was a typo for Case 2 -- confirmed with the user):
Case 1 (``ae-chiller``), Case 2 (``we-chiller``), Case 5 (``chiller-only``).

Outputs under ``outputs/ensemble_pue_wue_facilities/``:
- ``figures/figure{1..6,8}_*.png``
- ``tables/table{1..4}_*.csv``
"""

from __future__ import annotations

import argparse
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

# --------------------------------------------------------------------------
# Palette (see facility_PUE_explore.ipynb / dataviz skill references/palette.md)
# --------------------------------------------------------------------------
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
AXIS_LINE = "#c3c2b7"
SERIES_BLUE = "#2a78d6"
SERIES_ORANGE = "#eb6834"
SERIES_AQUA = "#1baf7a"

PERIOD_COLORS = {"historical": SERIES_BLUE, "midcentury": SERIES_ORANGE, "endcentury": SERIES_AQUA}
PERIOD_LABELS = {"historical": "Historical\n(1985-2014)", "midcentury": "Mid-century\n(2045-2074)", "endcentury": "End-of-century\n(2075-2100)"}
PERIOD_ORDER = ["historical", "midcentury", "endcentury"]

ARCHETYPE_LABELS = {
    "ae-chiller": "Case 1: AE + adiabatic",
    "we-chiller": "Case 2: waterside economizer",
    "chiller-only": "Case 5: chiller only",
}
ARCHETYPE_ORDER = ["ae-chiller", "we-chiller", "chiller-only"]
ARCHETYPE_MARKERS = {"ae-chiller": "o", "we-chiller": "^", "chiller-only": "s"}

SEQUENTIAL_BLUE_STEPS = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
    "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b",
]
SEQUENTIAL_BLUE_CMAP = LinearSegmentedColormap.from_list("sequential_blue", SEQUENTIAL_BLUE_STEPS)

SPOTCHECK_IDS = {19988712: "Meta Prineville", 231769626: "Google The Dalles"}


def _apply_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 140,
        "axes.edgecolor": AXIS_LINE,
        "axes.labelcolor": INK_SECONDARY,
        "axes.titlecolor": INK_PRIMARY,
        "text.color": INK_PRIMARY,
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "grid.color": GRIDLINE,
        "axes.grid": True,
        "grid.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 10,
    })


def _oregon_axes(fig, pos=111):
    ax = fig.add_subplot(pos, projection=ccrs.PlateCarree())
    ax.set_extent([-125, -116, 42, 46.5], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale("50m"), edgecolor=AXIS_LINE, linewidth=0.7, zorder=1)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), edgecolor=INK_SECONDARY, linewidth=0.8, zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), edgecolor=AXIS_LINE, linewidth=0.8, zorder=1)
    ax.grid(False)
    return ax


def figure1_site_map(by_site: pd.DataFrame, fig_dir: Path) -> None:
    sub = by_site[(by_site["archetype"] == "ae-chiller") & (by_site["period"] == "historical") & (by_site["metric"] == "pue_mean")]
    fig = plt.figure(figsize=(8, 4.8))
    ax = _oregon_axes(fig)
    sc = ax.scatter(sub["lon"], sub["lat"], c=sub["mean"], cmap=SEQUENTIAL_BLUE_CMAP, s=45, edgecolor="white", linewidth=0.4, zorder=2, transform=ccrs.PlateCarree())
    cbar = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("Historical baseline PUE (Case 1, ensemble mean)", color=INK_SECONDARY)
    cbar.ax.tick_params(color=INK_MUTED, labelcolor=INK_MUTED)
    ax.set_title("Figure 1. Oregon data center locations, historical baseline PUE", color=INK_PRIMARY, fontsize=11)
    fig.tight_layout()
    fig.savefig(fig_dir / "figure1_site_map.png", bbox_inches="tight")
    plt.close(fig)


def _grouped_bar(ax, statewide: pd.DataFrame, metric: str, ylabel: str) -> None:
    sub = statewide[statewide["metric"] == metric]
    x = np.arange(len(ARCHETYPE_ORDER))
    width = 0.25
    for i, period in enumerate(PERIOD_ORDER):
        vals = sub[sub["period"] == period].set_index("archetype").reindex(ARCHETYPE_ORDER)
        lo = vals["mean"] - vals["p05"]
        hi = vals["p95"] - vals["mean"]
        ax.bar(x + (i - 1) * width, vals["mean"], width=width, color=PERIOD_COLORS[period], label=PERIOD_LABELS[period].split("\n")[0], edgecolor=AXIS_LINE, linewidth=0.6)
        ax.errorbar(x + (i - 1) * width, vals["mean"], yerr=[lo, hi], fmt="none", ecolor=INK_PRIMARY, elinewidth=1.1, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels([ARCHETYPE_LABELS[a] for a in ARCHETYPE_ORDER], fontsize=8.5)
    ax.set_ylabel(ylabel)


def figure2_pue_by_period_archetype(statewide: pd.DataFrame, fig_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    _grouped_bar(ax, statewide, "pue_mean", "Fleet-mean PUE")
    ax.legend(frameon=False, fontsize=8.5)
    ax.set_title("Figure 2. Fleet-mean PUE by period and archetype\n(ensemble mean, 5th-95th pct. whiskers, n=20 GCMs)", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(fig_dir / "figure2_pue_by_period_archetype.png", bbox_inches="tight")
    plt.close(fig)


def figure3_wue_by_period_archetype(statewide: pd.DataFrame, fig_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    _grouped_bar(ax, statewide, "wue_mean", "Fleet-mean WUE (L/kWh)")
    ax.legend(frameon=False, fontsize=8.5)
    hist = statewide[(statewide["metric"] == "wue_mean") & (statewide["period"] == "historical")].set_index("archetype")
    holds = hist.loc["ae-chiller", "mean"] < hist.loc["we-chiller", "mean"]
    note = "Case 1 < Case 2 WUE (paper's finding) holds here" if holds else "Case 1 < Case 2 WUE does NOT hold here"
    ax.set_title(f"Figure 3. Fleet-mean WUE by period and archetype\n{note}", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(fig_dir / "figure3_wue_by_period_archetype.png", bbox_inches="tight")
    plt.close(fig)
    print(f"  Figure 3 check: {note} (Case 1={hist.loc['ae-chiller', 'mean']:.4f}, Case 2={hist.loc['we-chiller', 'mean']:.4f})")


def figure4_pue_vs_wue_scatter(by_site: pd.DataFrame, fig_dir: Path) -> None:
    pue = by_site[by_site["metric"] == "pue_mean"][["facility_id", "archetype", "period", "mean"]].rename(columns={"mean": "pue"})
    wue = by_site[by_site["metric"] == "wue_mean"][["facility_id", "archetype", "period", "mean"]].rename(columns={"mean": "wue"})
    merged = pue.merge(wue, on=["facility_id", "archetype", "period"])

    fig, ax = plt.subplots(figsize=(7.5, 6))
    for archetype in ARCHETYPE_ORDER:
        for period in PERIOD_ORDER:
            sub = merged[(merged["archetype"] == archetype) & (merged["period"] == period)]
            ax.scatter(sub["pue"], sub["wue"], color=PERIOD_COLORS[period], marker=ARCHETYPE_MARKERS[archetype], s=22, alpha=0.75, linewidth=0)
    period_handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=PERIOD_COLORS[p], label=PERIOD_LABELS[p].split("\n")[0], markersize=7) for p in PERIOD_ORDER]
    archetype_handles = [plt.Line2D([0], [0], marker=ARCHETYPE_MARKERS[a], color=INK_MUTED, label=ARCHETYPE_LABELS[a], markersize=7, linestyle="") for a in ARCHETYPE_ORDER]
    leg1 = ax.legend(handles=period_handles, title="Period", loc="upper left", frameon=False, fontsize=8)
    ax.add_artist(leg1)
    ax.legend(handles=archetype_handles, title="Archetype", loc="lower right", frameon=False, fontsize=8)
    ax.set_xlabel("Ensemble-mean PUE (per facility)")
    ax.set_ylabel("Ensemble-mean WUE (L/kWh, per facility)")
    ax.set_title("Figure 4. PUE vs. WUE across sites, archetypes, and periods", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(fig_dir / "figure4_pue_vs_wue_scatter.png", bbox_inches="tight")
    plt.close(fig)


def figure5_monthly_climatology(month_summary: pd.DataFrame, fig_dir: Path) -> None:
    fleet_month = month_summary.groupby(["archetype", "period", "month"], as_index=False)["pue_mean"].mean()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)
    for ax, archetype in zip(axes, ARCHETYPE_ORDER):
        sub = fleet_month[fleet_month["archetype"] == archetype]
        for period in PERIOD_ORDER:
            s = sub[sub["period"] == period].sort_values("month")
            ax.plot(s["month"], s["pue_mean"], color=PERIOD_COLORS[period], linewidth=2, marker="o", markersize=3.5, label=PERIOD_LABELS[period].split("\n")[0])
        ax.set_xticks(range(1, 13))
        ax.set_xticklabels(["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"], fontsize=7.5)
        ax.set_title(ARCHETYPE_LABELS[archetype], fontsize=9.5)
    axes[0].set_ylabel("Fleet-mean daily PUE")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Figure 5. Monthly PUE climatology (fleet mean across ensemble)", y=1.03, color=INK_PRIMARY)
    fig.tight_layout()
    fig.savefig(fig_dir / "figure5_monthly_climatology.png", bbox_inches="tight")
    plt.close(fig)

    peak_month = fleet_month.loc[fleet_month.groupby(["archetype", "period"])["pue_mean"].idxmax()]
    print("  Figure 5 peak-month check (July=7, September=9):")
    print(peak_month.sort_values(["archetype", "period"])[["archetype", "period", "month", "pue_mean"]].to_string(index=False))


def figure6_delta_map(delta: pd.DataFrame, by_site: pd.DataFrame, fig_dir: Path) -> None:
    ensemble_mean_delta = delta.groupby(["facility_id", "archetype"], as_index=False)["pue_mean_delta_endcentury"].mean()
    coords = by_site[["facility_id", "lon", "lat"]].drop_duplicates("facility_id")
    ensemble_mean_delta = ensemble_mean_delta.merge(coords, on="facility_id")

    vmin, vmax = ensemble_mean_delta["pue_mean_delta_endcentury"].min(), ensemble_mean_delta["pue_mean_delta_endcentury"].max()
    fig = plt.figure(figsize=(15, 6))
    for i, archetype in enumerate(ARCHETYPE_ORDER):
        ax = _oregon_axes(fig, pos=131 + i)
        sub = ensemble_mean_delta[ensemble_mean_delta["archetype"] == archetype]
        sc = ax.scatter(sub["lon"], sub["lat"], c=sub["pue_mean_delta_endcentury"], cmap=SEQUENTIAL_BLUE_CMAP, vmin=vmin, vmax=vmax, s=40, edgecolor="white", linewidth=0.4, zorder=2, transform=ccrs.PlateCarree())
        ax.set_title(ARCHETYPE_LABELS[archetype], fontsize=9.5)
    cbar = fig.colorbar(sc, ax=fig.axes, shrink=0.7, pad=0.02)
    cbar.set_label("Ensemble-mean PUE delta (end-of-century - historical)", color=INK_SECONDARY)
    fig.suptitle("Figure 6. Projected PUE change by facility", y=1.0, color=INK_PRIMARY)
    fig.savefig(fig_dir / "figure6_delta_map.png", bbox_inches="tight")
    plt.close(fig)


def figure7_heat_pue_trend(statewide: pd.DataFrame, heat_summary_path: Path, fig_dir: Path) -> None:
    """Extreme heat index + PUE trend together, both on SSP3-7.0 / the same 20-GCM ensemble.

    Both metrics come from ``scripts/aggregate_heat_ensemble.py`` /
    ``scripts/build_pue_wue_ensemble_summary.py`` -- the same GCM list read
    from the same ``gcm_manifest.csv`` -- so unlike the original 1.1
    (SSP5-8.5, single-GCM) vs. PUE/WUE (SSP3-7.0, 20-GCM ensemble) mismatch
    the task doc anticipated, no scenario-mismatch caveat is needed here.
    """
    heat = pd.read_csv(heat_summary_path)
    heat_days = heat[heat["metric"] == "tmax_threshold_days"].groupby("period", as_index=False)["mean"].mean()
    heat_days = heat_days.set_index("period").reindex(PERIOD_ORDER)["mean"]

    pue = statewide[(statewide["metric"] == "pue_mean") & (statewide["archetype"] == "ae-chiller")]
    pue = pue.set_index("period").reindex(PERIOD_ORDER)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    x = np.arange(len(PERIOD_ORDER))
    axes[0].bar(x, heat_days, color=SERIES_ORANGE, edgecolor=AXIS_LINE, linewidth=0.6)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([PERIOD_LABELS[p].split("\n")[0] for p in PERIOD_ORDER], fontsize=8.5)
    axes[0].set_ylabel("Fleet-mean threshold-exceedance days/yr")
    axes[0].set_title("Extreme heat (ASHRAE A2 threshold)", fontsize=10)

    lo = pue["mean"] - pue["p05"]
    hi = pue["p95"] - pue["mean"]
    axes[1].bar(x, pue["mean"], color=SERIES_BLUE, edgecolor=AXIS_LINE, linewidth=0.6)
    axes[1].errorbar(x, pue["mean"], yerr=[lo, hi], fmt="none", ecolor=INK_PRIMARY, elinewidth=1.1, capsize=3)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([PERIOD_LABELS[p].split("\n")[0] for p in PERIOD_ORDER], fontsize=8.5)
    axes[1].set_ylabel("Fleet-mean PUE (Case 1)")
    axes[1].set_title("Cooling efficiency (PUE, ensemble mean + 5th-95th pct.)", fontsize=10)

    fig.suptitle("Figure 7. Extreme heat and PUE trend, same SSP3-7.0 / 20-GCM ensemble", y=1.03, color=INK_PRIMARY)
    fig.tight_layout()
    fig.savefig(fig_dir / "figure7_heat_pue_trend.png", bbox_inches="tight")
    plt.close(fig)

    out = pd.DataFrame({"period": PERIOD_ORDER, "fleet_mean_threshold_days": heat_days.values, "fleet_mean_pue_case1": pue["mean"].values})
    return out


def figure8_ensemble_spread(period_summary: pd.DataFrame, fig_dir: Path) -> None:
    fleet_gcm = period_summary.groupby(["gcm", "archetype", "period"], as_index=False)["pue_mean"].mean()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)
    x = np.arange(len(PERIOD_ORDER))
    for ax, archetype in zip(axes, ARCHETYPE_ORDER):
        sub = fleet_gcm[fleet_gcm["archetype"] == archetype]
        for gcm in sorted(sub["gcm"].unique()):
            s = sub[sub["gcm"] == gcm].set_index("period").reindex(PERIOD_ORDER)["pue_mean"]
            ax.plot(x, s, color=INK_MUTED, linewidth=0.8, alpha=0.55, zorder=1)
        ens_mean = sub.groupby("period")["pue_mean"].mean().reindex(PERIOD_ORDER)
        ax.plot(x, ens_mean, color=SERIES_BLUE, linewidth=2.5, marker="o", zorder=2, label="Ensemble mean")
        ax.set_xticks(x)
        ax.set_xticklabels([PERIOD_LABELS[p].split("\n")[0] for p in PERIOD_ORDER], fontsize=8)
        ax.set_title(ARCHETYPE_LABELS[archetype], fontsize=9.5)
    axes[0].set_ylabel("Fleet-mean PUE")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Figure 8. Ensemble spread (individual GCMs, thin lines) vs. ensemble mean", y=1.03, color=INK_PRIMARY)
    fig.tight_layout()
    fig.savefig(fig_dir / "figure8_ensemble_spread.png", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# Tables
# --------------------------------------------------------------------------

def table1_data_sources(tab_dir: Path) -> None:
    rows = [
        {"Data": "Temperature (tasmax/tasmin)", "Source": "LOCA2 CMIP6, 20 of 27 GCMs with SSP3-7.0 downscaled (see gcm_manifest.csv)", "Use": "PUE/WUE and extreme-heat drivers"},
        {"Data": "Humidity (hursmax/hursmin)", "Source": "LOCA CMIP6 humidity companion (Pierce & Cayan, Scripps/UCSD), not MACA -- confirmed directly from loca2_io.py", "Use": "Cooling-system performance (PUE/WUE)"},
        {"Data": "Facility inventory", "Source": "IM3/PNNL open-source data center atlas, 109 Oregon facilities", "Use": "Site locations for all facility-level analysis"},
        {"Data": "Elevation", "Source": "LOCA elevation field (elevation.LOCA_2016-04-02.nc)", "Use": "Elevation-adjusted atmospheric pressure"},
        {"Data": "PUE/WUE model", "Source": "Lei & Masanet (2022), Cases 1/2/5 (Table B.1 midpoint parameters)", "Use": "Cooling-system physics"},
    ]
    pd.DataFrame(rows).to_csv(tab_dir / "table1_data_sources.csv", index=False)


def table2_archetype_definitions(tab_dir: Path) -> None:
    rows = [
        {"Case": 1, "Key": "ae-chiller", "Size class": "Large-scale", "Configuration": "Airside economizer + adiabatic cooling, water-cooled chiller supplemental", "Direct evaporation": "yes", "Space humidification": "adiabatic", "Cooling tower": "yes"},
        {"Case": 2, "Key": "we-chiller", "Size class": "Large-scale", "Configuration": "Waterside economizer + water-cooled chiller", "Direct evaporation": "n/a", "Space humidification": "n/a", "Cooling tower": "yes"},
        {"Case": 5, "Key": "chiller-only", "Size class": "Midsize", "Configuration": "Water-cooled chiller only, no economizer", "Direct evaporation": "n/a", "Space humidification": "adiabatic", "Cooling tower": "yes"},
    ]
    pd.DataFrame(rows).to_csv(tab_dir / "table2_archetype_definitions.csv", index=False)


def table3_historical_validation(period_summary: pd.DataFrame, by_site: pd.DataFrame, tab_dir: Path) -> None:
    rows = []
    for facility_id, name in SPOTCHECK_IDS.items():
        ens = by_site[(by_site["facility_id"] == facility_id) & (by_site["archetype"] == "ae-chiller") & (by_site["period"] == "historical") & (by_site["metric"] == "pue_mean")]
        single_gcm = period_summary[(period_summary["facility_id"] == facility_id) & (period_summary["archetype"] == "ae-chiller") & (period_summary["period"] == "historical") & (period_summary["gcm"] == "ACCESS-CM2")]
        rows.append({
            "facility": name,
            "archetype": "Case 1 (ae-chiller)",
            "ensemble_mean_pue": float(ens["mean"].iloc[0]) if len(ens) else np.nan,
            "ensemble_p05_pue": float(ens["p05"].iloc[0]) if len(ens) else np.nan,
            "ensemble_p95_pue": float(ens["p95"].iloc[0]) if len(ens) else np.nan,
            "single_gcm_ACCESS-CM2_pue": float(single_gcm["pue_mean"].iloc[0]) if len(single_gcm) else np.nan,
            "paper_plausibility_range": "1.0-1.2 (Fig. 6 AE-archetype contour)",
        })
    pd.DataFrame(rows).to_csv(tab_dir / "table3_historical_validation.csv", index=False)


def table4_summary_statistics(statewide: pd.DataFrame, tab_dir: Path) -> None:
    out = statewide.copy()
    out["mean_p05_p95"] = out.apply(lambda r: f"{r['mean']:.4f} ({r['p05']:.4f}-{r['p95']:.4f})", axis=1)
    out = out[["metric", "archetype", "period", "mean_p05_p95", "mean", "median", "p05", "p95", "n"]]
    out.to_csv(tab_dir / "table4_summary_statistics.csv", index=False)


def run(args: argparse.Namespace) -> None:
    _apply_style()
    root = Path(args.repo_root).resolve()
    analysis_dir = root / "outputs" / "ensemble_pue_wue_facilities" / "analysis"
    fig_dir = root / "outputs" / "ensemble_pue_wue_facilities" / "figures"
    tab_dir = root / "outputs" / "ensemble_pue_wue_facilities" / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    period_summary = pd.read_csv(analysis_dir / "facility_period_summary.csv")
    month_summary = pd.read_csv(analysis_dir / "facility_month_summary.csv")
    statewide = pd.read_csv(analysis_dir / "ensemble_summary_statewide.csv")
    by_site = pd.read_csv(analysis_dir / "ensemble_summary_by_site.csv")
    delta = pd.read_csv(analysis_dir / "ensemble_delta.csv")

    print("Figure 1 (site map)...")
    figure1_site_map(by_site, fig_dir)
    print("Figure 2 (PUE by period x archetype)...")
    figure2_pue_by_period_archetype(statewide, fig_dir)
    print("Figure 3 (WUE by period x archetype)...")
    figure3_wue_by_period_archetype(statewide, fig_dir)
    print("Figure 4 (PUE vs WUE scatter)...")
    figure4_pue_vs_wue_scatter(by_site, fig_dir)
    print("Figure 5 (monthly climatology)...")
    figure5_monthly_climatology(month_summary, fig_dir)
    print("Figure 6 (delta map)...")
    figure6_delta_map(delta, by_site, fig_dir)
    heat_summary_path = root / "outputs" / "month2_heat_ensemble" / "facility_heat_indices_ensemble_summary.csv"
    if heat_summary_path.exists():
        print("Figure 7 (heat index + PUE trend)...")
        fig7_table = figure7_heat_pue_trend(statewide, heat_summary_path, fig_dir)
        fig7_table.to_csv(tab_dir / "figure7_heat_pue_trend_data.csv", index=False)
    else:
        print(f"Figure 7 skipped -- {heat_summary_path} not found (run aggregate_heat_ensemble.py first)")
    print("Figure 8 (ensemble spread)...")
    figure8_ensemble_spread(period_summary, fig_dir)

    print("Table 1 (data sources)...")
    table1_data_sources(tab_dir)
    print("Table 2 (archetype definitions)...")
    table2_archetype_definitions(tab_dir)
    print("Table 3 (historical validation)...")
    table3_historical_validation(period_summary, by_site, tab_dir)
    print("Table 4 (summary statistics)...")
    table4_summary_statistics(statewide, tab_dir)

    print(f"Wrote figures to {fig_dir}, tables to {tab_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Phase 2/3 figures and tables for the PUE/WUE ensemble analysis.")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
