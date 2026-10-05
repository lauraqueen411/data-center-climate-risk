"""Figures for the restructured Results chapter (docs/results-draft.md,
2026-09-15 redesign). Replaces the original raw-value figure set
(docs/data-center-chapter-outline.md's Figures 1-9) with four delta-first
figures, arrived at by prototyping each one against a single question --
"does this figure deliver its one-sentence claim at a glance" -- rather than
by choosing an encoding that could hold every available dimension. See the
conversation this script comes out of for the two designs that failed that
test before landing here (a 3-row % change dot-whisker chart and a %-change
facility scatter that turned out to hide a correlation-of-ratios artifact).

Four figures, four distinct jobs -- deliberately non-overlapping so no two
figures are answering the same question in different clothes:

- **Figure 1 (orientation map)**: where are Oregon's 109 facilities, and how
  much hotter is Columbia Plateau already, today, than Willamette Valley --
  sets up Figure 3's "the milder region grows fastest" finding as a genuine
  twist rather than assumed background.
- **Figure 2 (seasonal shape)**: the response concentrates in a few summer
  months, and archetype (economizer design) changes the *shape* of that
  concentration, not just its size -- placed before Figure 3 so the
  statewide change-factor chart's later switch to a May-Oct-only PUE/WUE
  baseline (rather than an annual mean) already has its justification on
  the page, instead of asserting it cold.
- **Figure 3 (statewide change-factor chart)**: heat exposure is the driver,
  PUE and WUE are two responses to it that move by very different amounts --
  shown as three lines expressed as a change factor relative to their own
  historical baseline (historical = 1.0), not as three peers ranked on a
  shared % axis (that version crushed PUE and required reading three
  separate axis labels to see the divergence; a common baseline fixes both).
  "Change factor" (future value / historical value) is this project's own
  existing term for this calculation -- see methods-document.md Section
  2.1's flood-depth scaling, which computes the same ratio -- so this figure
  reuses that vocabulary rather than importing an economics term ("index")
  for an identical calculation.
- **Figure 4 (facility-level scatter)**: today's hottest facilities are not
  tomorrow's most-degraded ones -- absolute PUE/WUE delta vs. historical CDD,
  each facility a point. Deliberately **absolute** delta, not % change: %
  change on this specific comparison is a correlation-of-ratios trap,
  because historical WUE is itself extremely strongly correlated with CDD
  (r=0.97) -- dividing by it inflated the apparent WUE-CDD correlation to
  r=-0.90 versus the true r=-0.36 (matching the original draft) that
  absolute delta recovers.

Ecoregion aggregation conventions (kept consistent with what's already
established elsewhere in this project, not invented here):
- PUE/WUE at ecoregion granularity (none needed here beyond Figure 1's
  per-facility view) would use the single representative facility per
  ecoregion (``ecoregion_representative_facilities.csv``), matching the
  original Figure 9's convention -- Jensen's-inequality risk from averaging
  across facilities first, given the model's thresholded economizer/chiller
  switching (see ``build_pue_wue_ensemble_summary.py``'s Assumptions log).
- Heat indices use the full within-ecoregion facility mean instead (Figure
  1's cross-hatch calculation) -- direct climate-index calculations, no
  nonlinear model in the loop, so the averaging-bias concern doesn't apply
  and using every facility is simply more robust.
- Statewide aggregates (Figure 3) average facilities within each GCM first,
  divide each GCM's future value by its own historical value, then take the
  mean and 5th-95th percentile of those change factors across the 20 GCMs,
  so bands show inter-model (climate) uncertainty alone and collapse to
  exactly 1.0 at the historical baseline. This differs from
  ``ensemble_summary_statewide.csv``, which pools every facility x GCM value.

Assumptions log
----------------
- **Figure 1**: color = historical ensemble-mean threshold-exceedance days
  (orange sequential, matching Figure 3's heat=orange convention). Facilities
  in the top quartile of relative ensemble spread ((p95-p05)/mean, cutoff
  ~0.37) are cross-hatched rather than shown at reduced alpha -- alpha was
  tried first and rejected because it confounds with the color-value
  encoding (a low-value point is already pale; "pale because uncertain" and
  "pale because mild" became indistinguishable). The cross-hatched points
  cluster heavily in Willamette Valley (mean relative spread 0.40 vs.
  Columbia Plateau's 0.24) -- plausibly because its historical exceedance-day
  count sits close to the ASHRAE threshold itself, where small model-to-model
  differences swing the day count by a large relative amount. Worth a
  sentence in the chapter text, not just a figure footnote.
- **Figure 2**: end-of-century minus historical only, not both future
  periods -- a 6-line-per-panel version (3 archetypes x 2 periods) was
  judged likely to re-introduce the overload this redesign has been cutting;
  flag if mid-century turns out to tell a materially different shape. (Now
  mid-century, matching Figure 1, with a 5th-95th percentile band across the
  20 GCMs' fleet-mean deltas.)
- **Figure 3**: statewide, Case 1 only, no ecoregion facet -- a by-ecoregion
  small-multiple version was prototyped and rejected because a shared y-axis
  across all five ecoregions let the two single-facility regions' extreme
  values (a baseline-near-zero artifact, see below) dominate the scale and
  flatten the three real multi-facility regions into an undifferentiated
  cluster near the bottom.
- **Figure 4**: end-of-century only, not both periods -- a static scatter
  has no clean way to show a per-facility trajectory without either 218
  overlapping marks or 109 little arrows. No per-point uncertainty whiskers
  (Figure 3 already carries the ensemble-uncertainty message at the
  statewide level). Point alpha (not jitter) reveals overlapping facilities:
  many share the same or an adjacent LOCA2 grid cell (~6km resolution) --
  only ~9 of Columbia Plateau's 61 facilities and ~8 of Willamette Valley's
  33 have visually distinct historical CDD values.
- Both n=1 ecoregions (Eastern Cascades Slopes and Foothills, Klamath
  Mountains/California High North Coast Range) are shown, not dropped, but
  visually de-emphasized (hollow markers) everywhere they appear -- their
  historical-near-zero baselines produce artificially extreme % changes
  elsewhere in this chapter (results-draft.md's own Section 1 flags this),
  and de-emphasis keeps them visible without letting them set the terms for
  the rest of the figure.

Output: docs/results-figures/figure{1..4}_*.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator

from climate_risk_dc.config import load_dataset_paths

INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
AXIS_LINE = "#c3c2b7"
GRIDLINE = "#e1e0d9"
SERIES_BLUE = "#2a78d6"
SERIES_ORANGE = "#eb6834"
SERIES_AQUA = "#1baf7a"
SERIES_PURPLE = "#8858b0"

SEQUENTIAL_ORANGE_STEPS = [
    "#fde3d3", "#fbd0b3", "#f9bc93", "#f6a874", "#f39457", "#ef803d",
    "#eb6834", "#d85a2b", "#c14d24", "#a8401d", "#8f3417", "#762911", "#5e1f0c",
]
SEQUENTIAL_ORANGE_CMAP = LinearSegmentedColormap.from_list("sequential_orange", SEQUENTIAL_ORANGE_STEPS)

# Same blue ramp as build_pue_wue_ensemble_figures.py's SEQUENTIAL_BLUE_CMAP,
# reused here rather than inventing a second "PUE blue" so the two figure
# sets share one color vocabulary.
SEQUENTIAL_BLUE_STEPS = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
    "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b",
]
SEQUENTIAL_BLUE_CMAP = LinearSegmentedColormap.from_list("sequential_blue", SEQUENTIAL_BLUE_STEPS)

SEQUENTIAL_GREEN_STEPS = [
    "#d7f5e9", "#c2eedd", "#ade7d2", "#98e0c6", "#83d9bb", "#6ed2af",
    "#59cba3", "#1baf7a", "#189c6c", "#15895e", "#127650", "#0f6342", "#0c5034",
]
SEQUENTIAL_GREEN_CMAP = LinearSegmentedColormap.from_list("sequential_green", SEQUENTIAL_GREEN_STEPS)

METRIC_COLORS = {"heat": SERIES_ORANGE, "pue": SERIES_BLUE, "wue": SERIES_PURPLE}
METRIC_LABELS = {"heat": "Threshold-exceedance days", "pue": "PUE (Case 1)", "wue": "WUE (Case 1)"}

ARCHETYPE_LABELS = {"ae-chiller": "Case 1", "we-chiller": "Case 2", "chiller-only": "Case 5"}
ARCHETYPE_ORDER = ["ae-chiller", "we-chiller", "chiller-only"]
# chiller-only is purple, not the aqua/green used elsewhere in this file --
# Figure 2 plots it on the same axes as we-chiller's orange line, the same
# red/green-adjacent colorblind problem already fixed in Figure 3 (Erica
# Fleishman review, comment on that chart's orange/green pairing).
ARCHETYPE_COLORS = {"ae-chiller": SERIES_BLUE, "we-chiller": SERIES_ORANGE, "chiller-only": SERIES_PURPLE}

PERIOD_ORDER = ["historical", "midcentury", "endcentury"]
PERIOD_LABELS = ["Historical\n(1985-2014)", "Mid-century\n(2045-2074)", "End-of-century\n(2075-2100)"]
MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
HEAT_METRIC = "tmax_threshold_days"

ECOREGION_ORDER = [
    "Columbia Plateau", "Willamette Valley", "Blue Mountains",
    "Eastern Cascades Slopes and Foothills", "Klamath Mountains/California High North Coast Range",
]
DEEMPHASIZED = {"Eastern Cascades Slopes and Foothills", "Klamath Mountains/California High North Coast Range"}
# Same hues as the study-area map (Figure M0), saturated for use on points.
# Each ecoregion also gets its own marker shape, so the scatter doesn't rely
# on hue alone (red/green-safe, per Erica Fleishman's colorblind comment).
ECOREGION_COLORS = {
    "Columbia Plateau": "#3a9d5d",
    "Willamette Valley": "#d0559b",
    "Blue Mountains": "#3b6fd0",
    "Eastern Cascades Slopes and Foothills": "#e8931f",
    "Klamath Mountains/California High North Coast Range": "#1a9e96",
}
ECOREGION_MARKERS = {
    "Columbia Plateau": "o",
    "Willamette Valley": "s",
    "Blue Mountains": "^",
    "Eastern Cascades Slopes and Foothills": "D",
    "Klamath Mountains/California High North Coast Range": "v",
}
ECOREGION_SHORT_NAMES = {"Klamath Mountains/California High North Coast Range": "Klamath Mountains"}


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


# --------------------------------------------------------------------------
# Figure 1: orientation map (historical) + three change-factor maps
# --------------------------------------------------------------------------

# Same choice as build_pue_wue_ensemble_figures.py's OREGON_ALBERS: plain
# PlateCarree renders Oregon ~2:1 wide:tall against its true ~1.4:1 shape,
# since it plots degrees of longitude and latitude at equal visual scale.
OREGON_ALBERS = ccrs.AlbersEqualArea(central_longitude=-120.5, central_latitude=43.0, standard_parallels=(43.0, 45.5))


def _oregon_axes(fig, pos) -> plt.Axes:
    ax = fig.add_subplot(pos, projection=OREGON_ALBERS)
    ax.set_extent([-125, -116, 42, 46.5], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale("50m"), edgecolor=AXIS_LINE, linewidth=0.7, zorder=1)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), edgecolor=INK_SECONDARY, linewidth=0.8, zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), edgecolor=AXIS_LINE, linewidth=0.8, zorder=1)
    ax.grid(False)
    return ax


def _facility_metric_table(root: Path) -> pd.DataFrame:
    """One row per facility: ecoregion, historical value, and mid-century
    absolute change (future - historical, native units) for heat, PUE, and
    WUE (Case 1).

    Mid-century, not end-of-century, per feedback on Figure 1 specifically --
    methods-document.md calls mid-century "the primary analysis period"
    (matched to data-center mechanical/cooling equipment's service life), so
    this orientation figure leads with that horizon; end-of-century remains
    the headline comparison in Figures 3/4's text.

    Absolute change, not change factor (future / historical) -- a ratio
    blows up for the two single-facility ecoregions (Eastern Cascades Slopes
    and Foothills, Klamath Mountains/California High North Coast Range),
    whose historical baselines sit close to zero: a small absolute increase
    divides into a huge ratio. Figure 4 made the same call for the same
    reason (see its own docstring). Unlike Figure 3's shared-axis line chart,
    each panel here has its own independent color scale, so there's no
    cross-metric comparability lost by dropping the common
    ratio-relative-to-1.0 baseline."""
    by_site = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/ensemble_summary_by_site.csv")
    coords = by_site[["facility_id", "lon", "lat", "US_L3NAME"]].drop_duplicates("facility_id").rename(columns={"US_L3NAME": "ecoregion"})

    heat = pd.read_csv(root / "outputs/month2_heat_ensemble/facility_heat_indices_ensemble_summary.csv")
    heat = heat[heat["metric"] == HEAT_METRIC]
    heat_hist = heat[heat["period"] == "historical"].set_index("facility_id")
    heat_mid = heat[heat["period"] == "midcentury"].set_index("facility_id")
    heat_delta = pd.DataFrame({
        "heat_historical": heat_hist["mean"],
        "heat_rel_spread": (heat_hist["p95"] - heat_hist["p05"]) / heat_hist["mean"],
        "heat_delta": heat_mid["mean"] - heat_hist["mean"],
    })

    delta = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/ensemble_delta.csv")
    delta = delta[delta["archetype"] == "ae-chiller"].groupby("facility_id")[
        ["pue_mean_historical", "pue_mean_delta_midcentury", "wue_mean_historical", "wue_mean_delta_midcentury"]
    ].mean()
    delta = delta.rename(columns={
        "pue_mean_historical": "pue_historical", "wue_mean_historical": "wue_historical",
        "pue_mean_delta_midcentury": "pue_delta", "wue_mean_delta_midcentury": "wue_delta",
    })

    table = coords.set_index("facility_id").join(heat_delta).join(delta[["pue_historical", "wue_historical", "pue_delta", "wue_delta"]]).reset_index()
    return table


def _draw_panel(fig, pos, table: pd.DataFrame, value_col: str, cmap, colorbar_label: str, uncertainty: pd.Series | None = None, robust: bool = False, cbar_ax=None) -> plt.Axes:
    """One map panel: ecoregion polygons filled by the mean of that
    ecoregion's own facility-level outputs (opaque), all facilities plotted
    on top as points on the same color scale. Sharing one norm between the
    polygon fill and the points keeps "same color = same value" true across
    both layers.

    Mean of facility outputs, not one representative facility (this figure's
    original convention, and ``build_pue_wue_ensemble_summary.py``'s for the
    model-input stage): Jensen's-inequality bias only applies when averaging
    *inputs* to a nonlinear function (climate-input averaging, upstream, was
    rightly avoided for exactly that reason) -- averaging already-computed
    per-facility *outputs* is unbiased, uses every facility instead of one
    point estimate, and doesn't carry the "is this representative facility
    actually typical?" caveat Blue Mountains raised under the old convention.

    ``robust=True`` computes vmin/vmax from the multi-facility ecoregions
    only (excluding the two single-facility ones in DEEMPHASIZED), then
    clips the color norm so the excluded regions still render, just pinned
    to the extreme color rather than stretching the scale to fit them. This
    matters specifically for the change-factor panels: Eastern Cascades'
    near-zero heat baseline produces a change factor of 25.6, which would
    otherwise compress Columbia Plateau (3.3-3.6) and Willamette Valley
    (5.1-6.6)'s real, multi-facility difference into an indistinguishable
    pale band at the bottom of the scale. The historical panel doesn't
    have this problem (its de-emphasized regions aren't scale-dominating
    outliers) and doesn't use this option."""
    ax = _oregon_axes(fig, pos)
    scale_source = table[~table["ecoregion"].isin(DEEMPHASIZED)] if robust else table
    vmin, vmax = scale_source[value_col].min(), scale_source[value_col].max()
    norm = plt.Normalize(vmin=vmin, vmax=vmax, clip=True)

    ecoregions = gpd.read_file(load_dataset_paths().ecoregion_l3_shp).to_crs("EPSG:4326")
    eco_mean = table.groupby("ecoregion")[value_col].mean()
    for row in ecoregions.itertuples():
        eco_name = row.US_L3NAME
        if eco_name in eco_mean.index:
            face = cmap(norm(eco_mean[eco_name]))
        else:
            face = "none"
        ax.add_geometries([row.geometry], crs=ccrs.PlateCarree(), facecolor=face, edgecolor=AXIS_LINE, linewidth=0.8, zorder=0)
        if robust and eco_name in DEEMPHASIZED:
            ax.add_geometries([row.geometry], crs=ccrs.PlateCarree(), facecolor="none", edgecolor=INK_PRIMARY, linewidth=0.5, hatch="....", zorder=1)

    sc = ax.scatter(table["lon"], table["lat"], c=table[value_col], cmap=cmap, norm=norm, s=32, edgecolor=INK_PRIMARY, linewidth=0.5, zorder=2, transform=ccrs.PlateCarree())
    if uncertainty is not None:
        cutoff = uncertainty.quantile(0.75)
        flagged = table[uncertainty >= cutoff]
        ax.scatter(flagged["lon"], flagged["lat"], facecolor="none", edgecolor=INK_PRIMARY, linewidth=0.5, hatch="xxxx", s=32, zorder=3, transform=ccrs.PlateCarree())
    if cbar_ax is not None:
        cbar = fig.colorbar(sc, cax=cbar_ax, orientation="horizontal", extend="max" if robust else "neither")
    else:
        cbar = fig.colorbar(sc, ax=ax, shrink=0.75, pad=0.02, extend="max" if robust else "neither")
    cbar.locator = MaxNLocator(nbins=5)
    cbar.update_ticks()
    cbar.set_label(colorbar_label, color=INK_PRIMARY, fontsize=16)
    cbar.ax.tick_params(color=INK_PRIMARY, labelcolor=INK_PRIMARY, labelsize=14)
    return ax


def figure1_orientation_map(root: Path, fig_dir: Path) -> None:
    """First pass: bare 2x3 grid (historical top row, change-factor bottom
    row; heat/PUE/WUE columns), all six panels the same size and each with
    its own horizontal colorbar directly beneath it. Region labels,
    cross-hatching, and the ensemble-spread flag are deliberately left out
    of this pass -- added back once the base layout is confirmed."""
    table = _facility_metric_table(root)

    fig = plt.figure(figsize=(13, 9.5), constrained_layout=True)
    gs = fig.add_gridspec(4, 3, height_ratios=[1, 0.05, 1, 0.05], wspace=0.02)

    cmaps = [SEQUENTIAL_ORANGE_CMAP, SEQUENTIAL_BLUE_CMAP, SEQUENTIAL_GREEN_CMAP]
    var_names = ["Exceedance Days", "PUE", "WUE"]
    units = ["days/yr", "PUE", "L/kWh"]
    rows = [
        (0, 1, "Historical (1985-2014)", ["heat_historical", "pue_historical", "wue_historical"]),
        (2, 3, "Mid-Century Delta (Case 1)", ["heat_delta", "pue_delta", "wue_delta"]),
    ]
    letters = iter("abcdef")

    for map_row, cbar_row, row_title, value_cols in rows:
        for col, (value_col, cmap, var_name, unit) in enumerate(zip(value_cols, cmaps, var_names, units)):
            cax = fig.add_subplot(gs[cbar_row, col])
            ax = _draw_panel(fig, gs[map_row, col], table, value_col, cmap, unit, cbar_ax=cax)
            ax.text(0.03, 0.97, f"({next(letters)}) {var_name}", transform=ax.transAxes, ha="left", va="top", fontsize=16, color=INK_PRIMARY, zorder=4)
            if col == 1:
                ax.set_title(row_title, fontsize=18, color=INK_PRIMARY)

    fig.savefig(fig_dir / "figure1_orientation_map.png", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 3: statewide change-factor chart (helpers)
# --------------------------------------------------------------------------

def _change_factor(per_gcm: pd.DataFrame) -> pd.DataFrame:
    """Per-model future / historical ratio, then mean and 5th-95th percentile
    across models. ``per_gcm``: one row per GCM (fleet-mean values), one
    column per period in PERIOD_ORDER.

    Each model is divided by its *own* historical value, so model baseline
    bias cancels: historical = 1.0 for every model (zero-width band at the
    baseline) and the future bands show only inter-model disagreement about
    the change itself. Same calculation as methods-document.md Section 2.1's
    flood-depth change factor, applied to heat/PUE/WUE."""
    ratio = per_gcm[PERIOD_ORDER].div(per_gcm["historical"], axis=0)
    return pd.DataFrame({"cf_mean": ratio.mean(), "cf_p05": ratio.quantile(0.05), "cf_p95": ratio.quantile(0.95)})


def _statewide_heat_change_factor(root: Path) -> pd.DataFrame:
    facilities_dir = root / "outputs/month2_heat_ensemble/facilities"
    gcms = sorted({p.name.split("_")[3] for p in facilities_dir.glob("facility_heat_indices_*_ssp370_historical.csv")})
    # Fleet mean per GCM first (facilities averaged), so the band is
    # inter-model (climate) uncertainty only, not between-facility spread.
    per_gcm = pd.DataFrame(
        {period: [pd.read_csv(facilities_dir / f"facility_heat_indices_{gcm}_ssp370_{period}.csv")[HEAT_METRIC].mean() for gcm in gcms] for period in PERIOD_ORDER},
        index=gcms,
    )
    return _change_factor(per_gcm)


# --------------------------------------------------------------------------
# Figure 2: seasonal shape by archetype
# --------------------------------------------------------------------------

def figure2_seasonal_delta(root: Path, fig_dir: Path) -> None:
    """Monthly mid-century minus historical PUE/WUE, fleet mean, by archetype.

    Mid-century (not end-of-century) to match Figure 1's horizon. Panel
    letters replace subplot titles; the mid-century basis lives in the
    y-axis labels, since in-figure titles are dropped as redundant with the
    manuscript caption.

    Shaded band = interquartile range (25th-75th percentile) across the
    20-GCM ensemble of each GCM's fleet-mean delta (facilities averaged
    first, then percentile across GCMs) -- pure inter-model spread. IQR
    rather than 5th-95th, matching OCA7's own convention for this figure
    type (e.g. Figures 8, 17: ensemble mean line + model dots/envelope),
    because with only 20 models the 5th/95th tails are each defined by a
    single model, so which model briefly leads at a given month (e.g.
    KACE-1-0-G vs. CanESM5 for Case 1 in June-August) produces a
    month-to-month pinch in the band that isn't present in the mean; IQR
    draws on the middle 10 of 20 models and is far less sensitive to any
    one model's rank swapping. This differs from Figure 3's band, which
    uses 5th-95th and pools every facility x GCM value into one
    distribution, so also carries between-facility spread (roughly 1.5-2x
    wider for Case 1 in August).
    """
    month_summary = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/facility_month_summary.csv")
    fleet_gcm = month_summary.groupby(["archetype", "gcm", "month", "period"], as_index=False)[["pue_mean", "wue_mean"]].mean()

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    panels = [("pue_mean", "Mid-century ΔPUE", "(a) PUE"), ("wue_mean", "Mid-century ΔWUE (L/kWh)", "(b) WUE")]
    for ax, (metric, ylabel, letter) in zip(axes, panels):
        wide = fleet_gcm.pivot_table(index=["archetype", "gcm", "month"], columns="period", values=metric)
        delta = (wide["midcentury"] - wide["historical"]).rename("delta").reset_index()
        for archetype in ARCHETYPE_ORDER:
            d = delta[delta["archetype"] == archetype].groupby("month")["delta"]
            mean, lo, hi = d.mean(), d.quantile(0.25), d.quantile(0.75)
            color = ARCHETYPE_COLORS[archetype]
            ax.fill_between(range(1, 13), lo.reindex(range(1, 13)), hi.reindex(range(1, 13)), color=color, alpha=0.15, linewidth=0, zorder=1)
            ax.plot(range(1, 13), mean.reindex(range(1, 13)), color=color, linewidth=2.2, marker="o", markersize=4, label=ARCHETYPE_LABELS[archetype], zorder=2)
        ax.axhline(0, color=INK_MUTED, linewidth=0.8, linestyle=":")
        ax.set_xticks(range(1, 13))
        ax.set_xticklabels(MONTH_LABELS, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=14)
        # Erica Fleishman review: black axis text in place of this style's
        # default muted gray, same fix already applied to Figure 3.
        ax.tick_params(axis="both", colors=INK_PRIMARY, labelsize=14)
        ax.tick_params(axis="x", labelsize=12)
        ax.yaxis.label.set_color(INK_PRIMARY)
        ax.text(0.03, 0.97, letter, transform=ax.transAxes, ha="left", va="top", fontsize=16, color=INK_PRIMARY)

    axes[0].legend(frameon=False, fontsize=14, loc="upper left", bbox_to_anchor=(0.0, 0.90))
    fig.tight_layout()
    fig.savefig(fig_dir / "figure2_seasonal_delta.png", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 3: statewide change-factor chart
# --------------------------------------------------------------------------

WARM_SEASON_MONTHS = range(5, 11)  # May-Oct (6 months), per Figure 2's own seasonal-concentration finding


def _warm_season_pue_wue_by_gcm(root: Path) -> pd.DataFrame:
    """Fleet-mean PUE/WUE over the warm-season months only, one row per
    (gcm, archetype, period): facilities are averaged first, so downstream
    spread across GCMs is climate uncertainty alone.

    Warm season = May-Oct. Annual PUE/WUE dilutes the real signal with ~6
    near-flat winter months (Figure 2 shows this directly) and, for WUE
    specifically, an annual historical baseline close to zero inflates the
    annual change factor via a small-denominator effect (2.38 annual vs. 2.29
    May-Oct). Heat exposure (threshold-exceedance days) stays on its native
    annual definition: it is already effectively warm-season-only (Oregon
    essentially never exceeds 35C outside summer), and no monthly breakdown
    of it is computed.
    """
    month = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/facility_month_summary.csv")
    warm = month[month["month"].isin(WARM_SEASON_MONTHS)]
    return warm.groupby(["gcm", "archetype", "period"], as_index=False)[["pue_mean", "wue_mean"]].mean()


def figure3_statewide_change_factor(root: Path, fig_dir: Path) -> None:
    warm_by_gcm = _warm_season_pue_wue_by_gcm(root)

    def pue_wue_change_factor(metric: str) -> pd.DataFrame:
        sub = warm_by_gcm[warm_by_gcm["archetype"] == "ae-chiller"].pivot(index="gcm", columns="period", values=metric)
        return _change_factor(sub)

    series = {"heat": _statewide_heat_change_factor(root), "pue": pue_wue_change_factor("pue_mean"), "wue": pue_wue_change_factor("wue_mean")}

    fig, ax = plt.subplots(figsize=(7.5, 6))
    x = range(len(PERIOD_ORDER))
    for key, df in series.items():
        color = METRIC_COLORS[key]
        ax.plot(x, df["cf_mean"], color=color, linewidth=2.2, marker="o", markersize=5, label=METRIC_LABELS[key], zorder=2)
        ax.fill_between(x, df["cf_p05"], df["cf_p95"], color=color, alpha=0.15, zorder=1, linewidth=0)
    ax.axhline(1.0, color=INK_MUTED, linewidth=0.8, linestyle=":", zorder=0)
    ax.set_xticks(list(x))
    ax.set_xticklabels(PERIOD_LABELS, fontsize=14)
    # The annual-vs-warm-season basis distinction (heat exposure: annual;
    # PUE/WUE: May-Oct mean, per _warm_season_pue_wue_by_gcm's docstring)
    # lives in the manuscript caption, not here -- kept off the in-figure
    # label itself per feedback to keep this concise.
    ax.set_ylabel("Change factor rel. to 1985-2014 baseline", fontsize=14)
    # Erica Fleishman review: red/green (orange heat-exposure line vs. green
    # WUE line, with overlapping shaded uncertainty bands) is a colorblind
    # accessibility problem here specifically, since the two bands visually
    # overlap on a shared axis -- WUE gets a purple that isn't confusable
    # with the orange heat-exposure line, rather than recoloring WUE
    # site-wide (Figure 1's WUE map panel is a separate, non-overlapping
    # context and keeps its own green). She also asked for black axis
    # labels/ticks in place of this style's default muted gray.
    ax.tick_params(axis="both", colors=INK_PRIMARY, labelsize=14)
    ax.yaxis.label.set_color(INK_PRIMARY)
    ax.legend(frameon=False, fontsize=14, loc="upper left")
    fig.tight_layout()
    fig.savefig(fig_dir / "figure3_statewide_change_factor.png", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 4: facility-level scatter
# --------------------------------------------------------------------------

def _facility_historical_humidity(root: Path) -> pd.DataFrame:
    """Per-facility historical (1985-2014) warm-season (May-Oct) mean afternoon
    relative humidity and wet-bulb temperature, plus mean temperature,
    averaged across the 20 GCMs.

    "Afternoon" = the model's hot pairing: daily Tmax with the day's minimum
    RH (hursmin), the same inputs the PUE/WUE model sees. Wet-bulb comes from
    the model's own ``_wetbulb_from_rh`` at each facility's elevation-adjusted
    pressure. Read from the ae-chiller historical daily files (the climate
    inputs are identical across archetypes). Cached to
    ``outputs/ensemble_pue_wue_facilities/analysis/facility_historical_humidity.csv``.
    """
    cache = root / "outputs/ensemble_pue_wue_facilities/analysis/facility_historical_humidity.csv"
    if cache.exists():
        return pd.read_csv(cache)
    import xarray as xr
    from climate_risk_dc.climate.cooling_archetypes import elevation_adjusted_pressure_pa
    from climate_risk_dc.climate.cooling_archetypes_fast import _wetbulb_from_rh
    from climate_risk_dc.geo import load_oregon_data_centers

    fac = load_oregon_data_centers(root / "data/data_centers_im3_pnnl/im3_open_source_data_center_atlas_v2026.02.09.csv")
    elev = xr.open_dataset(root / "data/elevation.LOCA_2016-04-02.nc")["Elevation"]
    lon, lat = fac["lon"].to_numpy(), fac["lat"].to_numpy()
    lon_vals, lat_vals = np.asarray(elev["lon"].values), np.asarray(elev["lat"].values)
    lon_q = np.where(lon < 0, lon + 360.0, lon) if lon_vals.max() > 180.0 else lon
    lat_idx = np.abs(lat_vals[:, None] - lat[None, :]).argmin(axis=0)
    lon_idx = np.abs(lon_vals[:, None] - lon_q[None, :]).argmin(axis=0)
    pressure = pd.Series(elevation_adjusted_pressure_pa(elev.values[lat_idx, lon_idx]), index=fac["id"].to_numpy())

    facilities_dir = root / "outputs/ensemble_pue_wue_facilities"
    parts = []
    for f in sorted(facilities_dir.glob("facility_pue_wue_*_ssp370_ae-chiller_historical_fast.csv")):
        d = pd.read_csv(f, usecols=["facility_id", "date", "tasmax_c", "tasmin_c", "hursmin_pct"])
        d = d[pd.to_datetime(d["date"]).dt.month.isin(WARM_SEASON_MONTHS)]
        d["tmean_c"] = (d["tasmax_c"] + d["tasmin_c"]) / 2
        d["wetbulb_c"] = _wetbulb_from_rh(d["tasmax_c"].to_numpy(), d["hursmin_pct"].to_numpy() / 100.0, d["facility_id"].map(pressure).to_numpy())
        parts.append(d.groupby("facility_id")[["hursmin_pct", "wetbulb_c", "tmean_c"]].mean())
    out = pd.concat(parts).groupby(level=0).mean().reset_index().rename(columns={"hursmin_pct": "rh_afternoon_pct"})
    out.to_csv(cache, index=False)
    return out


def figure4_facility_scatter(root: Path, fig_dir: Path, horizon: str = "endcentury", x_metric: str = "cdd") -> None:
    """Facility-level absolute PUE/WUE change vs. historical CDD, one point per
    facility, colored/shaped by ecoregion.

    No regression line or correlation coefficient: the points form distinct
    ecoregion clusters (a between-cluster pattern, not a within-cluster linear
    trend), so a single fitted line and Pearson r overstate a general
    relationship. The finding shown is descriptive -- historically cooler
    regions (Willamette Valley) see the largest absolute changes in cooling
    efficiency. ``horizon`` is "endcentury" or "midcentury"; mid-century is
    saved with a ``_midcentury`` filename suffix. ``x_metric`` is "cdd"
    (cooling degree days) or "exceedance" (annual threshold-exceedance days,
    matching Figure 1's measure of baseline heat); the latter gets an
    ``_exceedance`` suffix; "rh" / "wetbulb" use warm-season historical
    afternoon relative humidity / wet-bulb temperature (see
    ``_facility_historical_humidity``).
    """
    label = {"endcentury": "End-of-century", "midcentury": "Mid-century"}[horizon]
    by_site = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/ensemble_summary_by_site.csv")
    eco_map = by_site[["facility_id", "US_L3NAME"]].drop_duplicates().rename(columns={"US_L3NAME": "ecoregion"})
    xlabel = {"cdd": "Historical CDD (\u00b0C-days)", "exceedance": "Historical exceedance days (days/yr)",
              "rh": "Historical May\u2013Oct afternoon RH (%)", "wetbulb": "Historical May\u2013Oct afternoon wet-bulb (\u00b0C)"}[x_metric]
    if x_metric in ("rh", "wetbulb"):
        hum = _facility_historical_humidity(root)
        cdd = hum[["facility_id", {"rh": "rh_afternoon_pct", "wetbulb": "wetbulb_c"}[x_metric]]].set_axis(["facility_id", "cdd_historical"], axis=1)
    else:
        heat_metric = {"cdd": "cdd", "exceedance": HEAT_METRIC}[x_metric]
        cdd = pd.read_csv(root / "outputs/month2_heat_ensemble/facility_heat_indices_ensemble_summary.csv")
        cdd = cdd[(cdd["metric"] == heat_metric) & (cdd["period"] == "historical")][["facility_id", "mean"]].rename(columns={"mean": "cdd_historical"})
    delta = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/ensemble_delta.csv")
    delta = delta[delta["archetype"] == "ae-chiller"]
    per_facility = delta.groupby("facility_id")[[f"pue_mean_delta_{horizon}", f"wue_mean_delta_{horizon}"]].mean().reset_index()
    per_facility = per_facility.rename(columns={f"pue_mean_delta_{horizon}": "pue_delta", f"wue_mean_delta_{horizon}": "wue_delta"})
    df = per_facility.merge(cdd, on="facility_id").merge(eco_map, on="facility_id")

    def panel(ax, ycol: str, ylabel: str, letter: str) -> None:
        # All markers filled, including the two single-facility ecoregions
        # (Erica Fleishman review: hollow markers unnecessary since every
        # facility is plotted as its own point).
        for eco, color in ECOREGION_COLORS.items():
            sub = df[df["ecoregion"] == eco]
            if sub.empty:
                continue
            ax.scatter(sub["cdd_historical"], sub[ycol], s=70, marker=ECOREGION_MARKERS[eco], facecolor=color, edgecolor=INK_PRIMARY, linewidth=0.6, alpha=0.85, label=ECOREGION_SHORT_NAMES.get(eco, eco), zorder=2)
        # Headroom above the data so the panel letter never lands on a point.
        y0, y1 = ax.get_ylim()
        ax.set_ylim(y0, y1 + 0.10 * (y1 - y0))
        if x_metric == "exceedance":
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.text(0.03, 0.97, letter, transform=ax.transAxes, ha="left", va="top", fontsize=16, color=INK_PRIMARY)
        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.tick_params(axis="both", colors=INK_PRIMARY, labelsize=14)
        ax.xaxis.label.set_color(INK_PRIMARY)
        ax.yaxis.label.set_color(INK_PRIMARY)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
    panel(axes[0], "pue_delta", f"{label} \u0394PUE", "(a) PUE")
    panel(axes[1], "wue_delta", f"{label} \u0394WUE (L/kWh)", "(b) WUE")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=13, bbox_to_anchor=(0.5, -0.12), markerscale=1.2)
    fig.tight_layout()
    suffix = (f"_{x_metric}" if x_metric != "cdd" else "") + ("_midcentury" if horizon == "midcentury" else "")
    fig.savefig(fig_dir / f"figure4_facility_scatter{suffix}.png", bbox_inches="tight")
    plt.close(fig)


def _facility_warming(root: Path, horizon: str = "midcentury") -> pd.DataFrame:
    """Per-facility warm-season (May-Oct) mean-temperature change, ``horizon``
    minus historical, averaged across the 20 GCMs (each GCM differenced
    first). Cached to ``.../analysis/facility_warming_{horizon}.csv``."""
    cache = root / f"outputs/ensemble_pue_wue_facilities/analysis/facility_warming_{horizon}.csv"
    if cache.exists():
        return pd.read_csv(cache)
    facilities_dir = root / "outputs/ensemble_pue_wue_facilities"

    def warm_tmean(period: str, gcm_file: Path) -> pd.Series:
        d = pd.read_csv(gcm_file, usecols=["facility_id", "date", "tasmax_c", "tasmin_c"])
        d = d[pd.to_datetime(d["date"]).dt.month.isin(WARM_SEASON_MONTHS)]
        return ((d["tasmax_c"] + d["tasmin_c"]) / 2).groupby(d["facility_id"]).mean()

    diffs = []
    for hist_file in sorted(facilities_dir.glob("facility_pue_wue_*_ssp370_ae-chiller_historical_fast.csv")):
        fut_file = Path(str(hist_file).replace("_historical_fast", f"_{horizon}_fast"))
        diffs.append(warm_tmean(horizon, fut_file) - warm_tmean("historical", hist_file))
    out = pd.concat(diffs, axis=1).mean(axis=1).rename("warming_c").reset_index()
    out.to_csv(cache, index=False)
    return out


def figure4_hazard_vs_response(root: Path, fig_dir: Path, horizon: str = "midcentury") -> None:
    """Draft Figure 4: (a) warm-season warming and (b) PUE change (Case 1), each
    against historical May-Oct afternoon wet-bulb temperature.

    Message: the climate hazard (warming) is nearly uniform across Oregon's
    facilities, but the cooling-efficiency response differs about 2x, and
    tracks baseline wet-bulb -- vulnerability, not hazard, varies. Facilities
    form ecoregion clusters (effective n of about five), so no fit line.
    """
    label = {"endcentury": "End-of-century", "midcentury": "Mid-century"}[horizon]
    by_site = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/ensemble_summary_by_site.csv")
    eco = by_site[["facility_id", "US_L3NAME"]].drop_duplicates().rename(columns={"US_L3NAME": "ecoregion"})
    delta = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/ensemble_delta.csv")
    delta = delta[delta["archetype"] == "ae-chiller"].groupby("facility_id")[f"pue_mean_delta_{horizon}"].mean().rename("pue_delta").reset_index()
    df = delta.merge(_facility_historical_humidity(root), on="facility_id").merge(_facility_warming(root, horizon), on="facility_id").merge(eco, on="facility_id")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), sharex=True)
    specs = [("warming_c", f"{label} warming (\u00b0C)", "(a) Warming"), ("pue_delta", f"{label} \u0394PUE", "(b) PUE change")]
    for ax, (ycol, ylabel, letter) in zip(axes, specs):
        for name, color in ECOREGION_COLORS.items():
            sub = df[df["ecoregion"] == name]
            ax.scatter(sub["wetbulb_c"], sub[ycol], s=70, marker=ECOREGION_MARKERS[name], facecolor=color, edgecolor=INK_PRIMARY, linewidth=0.6, alpha=0.85, label=ECOREGION_SHORT_NAMES.get(name, name), zorder=2)
        if ycol == "warming_c":
            ax.set_ylim(2.0, 4.0)
        else:
            y0, y1 = ax.get_ylim()
            ax.set_ylim(y0, y1 + 0.10 * (y1 - y0))
        ax.text(0.03, 0.97, letter, transform=ax.transAxes, ha="left", va="top", fontsize=16, color=INK_PRIMARY)
        ax.set_xlabel("Historical May\u2013Oct afternoon wet-bulb (\u00b0C)", fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.tick_params(axis="both", colors=INK_PRIMARY, labelsize=14)
        ax.xaxis.label.set_color(INK_PRIMARY)
        ax.yaxis.label.set_color(INK_PRIMARY)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=13, bbox_to_anchor=(0.5, -0.12), markerscale=1.2)
    fig.tight_layout()
    fig.savefig(fig_dir / f"figure4_hazard_vs_response{'_' + horizon if horizon != 'endcentury' else ''}.png", bbox_inches="tight")
    plt.close(fig)


def figure4_temperature_vs_humidity(root: Path, fig_dir: Path, horizon: str = "midcentury") -> None:
    """Draft Figure 4: ΔPUE (top row) and ΔWUE (bottom row), Case 1, against
    historical May-Oct mean temperature (left column) and historical May-Oct
    afternoon wet-bulb temperature (right column).

    Message: temperature alone gives an unstable, uninformative picture of
    which facilities see the largest cooling-efficiency change; humidity-aware
    wet-bulb temperature orders them consistently. Rank correlations across
    the 109 facilities (Spearman): temperature measures -0.1 to -0.3, wet-bulb
    ~+0.9. No fit line -- facilities form ecoregion clusters (effective n
    of about five).
    """
    label = {"endcentury": "End-of-century", "midcentury": "Mid-century"}[horizon]
    by_site = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/ensemble_summary_by_site.csv")
    eco = by_site[["facility_id", "US_L3NAME"]].drop_duplicates().rename(columns={"US_L3NAME": "ecoregion"})
    delta = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/ensemble_delta.csv")
    delta = delta[delta["archetype"] == "ae-chiller"].groupby("facility_id")[[f"pue_mean_delta_{horizon}", f"wue_mean_delta_{horizon}"]].mean()
    delta.columns = ["pue_delta", "wue_delta"]
    df = delta.reset_index().merge(_facility_historical_humidity(root), on="facility_id").merge(eco, on="facility_id")

    fig, axes = plt.subplots(2, 2, figsize=(13, 10.5), sharex="col", sharey="row")
    xcols = [("tmean_c", "Historical May\u2013Oct mean temperature (\u00b0C)"), ("wetbulb_c", "Historical May\u2013Oct afternoon wet-bulb (\u00b0C)")]
    ycols = [("pue_delta", f"{label} \u0394PUE", "PUE"), ("wue_delta", f"{label} \u0394WUE (L/kWh)", "WUE")]
    letters = iter("abcd")
    for row, (ycol, ylabel, var) in enumerate(ycols):
        for col, (xcol, xlabel) in enumerate(xcols):
            ax = axes[row, col]
            for name, color in ECOREGION_COLORS.items():
                sub = df[df["ecoregion"] == name]
                ax.scatter(sub[xcol], sub[ycol], s=70, marker=ECOREGION_MARKERS[name], facecolor=color, edgecolor=INK_PRIMARY, linewidth=0.6, alpha=0.85, label=ECOREGION_SHORT_NAMES.get(name, name), zorder=2)
            y0, y1 = ax.get_ylim()
            ax.set_ylim(y0, y1 + 0.10 * (y1 - y0))
            ax.text(0.03, 0.97, f"({next(letters)}) {var}", transform=ax.transAxes, ha="left", va="top", fontsize=16, color=INK_PRIMARY)
            ax.tick_params(axis="both", colors=INK_PRIMARY, labelsize=14)
            if row == 1:
                ax.set_xlabel(xlabel, fontsize=14)
                ax.xaxis.label.set_color(INK_PRIMARY)
            if col == 0:
                ax.set_ylabel(ylabel, fontsize=14)
                ax.yaxis.label.set_color(INK_PRIMARY)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=13, bbox_to_anchor=(0.5, -0.07), markerscale=1.2)
    fig.tight_layout()
    fig.savefig(fig_dir / f"figure4_temperature_vs_humidity{'_' + horizon if horizon != 'endcentury' else ''}.png", bbox_inches="tight")
    plt.close(fig)


def run(args: argparse.Namespace) -> None:
    _apply_style()
    root = Path(args.repo_root).resolve()
    fig_dir = root / "docs" / "results-figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print("Figure 1 (orientation map)...")
    figure1_orientation_map(root, fig_dir)
    print("Figure 2 (seasonal shape by archetype)...")
    figure2_seasonal_delta(root, fig_dir)
    print("Figure 3 (statewide change-factor chart)...")
    figure3_statewide_change_factor(root, fig_dir)
    print("Figure 4 (facility-level scatter)...")
    figure4_facility_scatter(root, fig_dir)
    figure4_facility_scatter(root, fig_dir, horizon="midcentury")
    print(f"Wrote figures to {fig_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the four delta-first Results figures (2026-09-15 redesign).")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
