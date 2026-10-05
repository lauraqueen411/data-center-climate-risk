"""The "show everything" results figure for the heat risk chapter
(docs/tasks/task-show-everything-figure.md).

One figure showing every cooling case (1, 2, 5), both future windows and
every ecoregion that contains a facility, laid out after Lei and Masanet
(2022), Fig. 5. Replaces the statewide relative-change figure.

- Panel (a): change in annual threshold-exceedance days (Tmax > 35 C).
- Panels (b-d): change in annual mean PUE, Cases 1, 2, 5.
- Panels (e-g): change in annual mean WUE (L/kWh), Cases 1, 2, 5.

Per ecoregion and future window, one strip: facility points (mean of the
20 models' change), the ecoregion mean of those points (equal weight per
facility), and a shaded bar giving the spread across the 20 models of the
per-model ecoregion mean change (``BAND``).

Inputs (all already on disk, nothing re-run):
- Daily PUE/WUE per facility, per model:
  ``outputs/ensemble_pue_wue_facilities/facility_pue_wue_{gcm}_ssp370_{archetype}_{period}_fast.csv``
  (Case 1 = ``ae-chiller``, Case 2 = ``we-chiller``, Case 5 = ``chiller-only``).
- Period-mean annual exceedance-day counts per facility, per model:
  ``outputs/month2_heat_ensemble/facilities/facility_heat_indices_{gcm}_ssp370_{period}.csv``
  (``tmax_threshold_days``).
- Historical afternoon wet-bulb: recomputed here from the daily Tmax and
  RHmin in the Case 1 historical files, at each facility's elevation-adjusted
  pressure, with the same functions the PUE/WUE model uses.

Outputs:
- ``docs/results-figures/fig_show_everything.{png,pdf}`` and
  ``fig_show_everything_caption.txt``.
- ``outputs/ensemble_pue_wue_facilities/analysis/fig_show_everything/``:
  ``facility_changes.csv``, ``ecoregion_changes.csv``,
  ``model_agreement.csv``, ``model_values_cache.csv`` (per facility, model,
  case, period values read from the 180 daily files) and ``run_manifest.csv``.

Assumptions log
----------------
- **Annual means, then period means.** PUE/WUE: the mean of daily values
  within each calendar year, then the mean over the years in the window.
  This differs slightly from the pooled all-days mean in
  ``analysis/facility_period_summary.csv`` because leap years carry one
  extra day; the task fixes the annual-then-period convention.
- **Exceedance days come from the stored heat-ensemble outputs** (period
  mean of annual counts, Tmax > 35 C, nearest LOCA2 cell sampled from the
  per-model GeoTIFFs). Year-by-year counts are not stored, but the period
  mean of annual counts is exactly the plotted metric. As a check, the same
  count is recomputed from the daily ``tasmax_c`` in the PUE/WUE files and
  the maximum difference is printed.
- **Wet-bulb window** is set by ``WETBULB_WINDOW``. ``"annual"`` (the
  default) averages all days; ``"may_oct"`` reproduces the draft's
  warm-season values. Either way: per-year mean, then mean over 1985-2014,
  then mean over the 20 models. Both are computed so the ecoregion order
  under each can be compared.
- **LOCA2 grid cell** is not stored by the PUE/WUE pipeline; it is
  recovered with the same nearest-cell rule
  (``loca2_io.sample_datacube_at_points``) on the Oregon LOCA2 grid in the
  humidity cache. Verified: 22 cells, matching 22 distinct daily input
  series across the 109 facilities.
- **Band percentiles** use ``np.percentile`` with default linear
  interpolation over 20 values. All three band choices are written to the
  ecoregion CSV whatever ``BAND`` is.
- **Ecoregion order** is computed (ascending ecoregion-mean historical
  wet-bulb, equal weight per facility), not hard-coded.
"""

from __future__ import annotations

import argparse
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------------
BAND = "iqr"  # "iqr" | "p05_p95" | "minmax"
WETBULB_WINDOW = "annual"  # "annual" | "may_oct"
SHARE_Y_ROWS = True

BAND_COLUMNS = {"iqr": ("p25", "p75"), "p05_p95": ("p05", "p95"), "minmax": ("min", "max")}
BAND_LEGEND = {
    "iqr": "Interquartile range",
    "p05_p95": "5th–95th percentile",
    "minmax": "Full range",
}
BAND_CAPTION = {
    "iqr": "the interquartile range (25th–75th percentile)",
    "p05_p95": "the 5th–95th percentile range",
    "minmax": "the full range (minimum to maximum)",
}

# ---------------------------------------------------------------------------
# Fixed conventions
# ---------------------------------------------------------------------------
N_MODELS = 20
THRESHOLD_C = 35.0
PERIODS = {"historical": (1985, 2014), "midcentury": (2045, 2074), "endcentury": (2075, 2100)}
FUTURE_PERIODS = ["midcentury", "endcentury"]
PERIOD_LEGEND = {"midcentury": "Mid-century (2045–2074)", "endcentury": "End-of-century (2075–2100)"}
CASES = {"ae-chiller": "Case 1", "we-chiller": "Case 2", "chiller-only": "Case 5"}
WARM_SEASON_MONTHS = range(5, 11)

ECOREGION_SHORT = {
    "Eastern Cascades Slopes and Foothills": "Cascades",
    "Blue Mountains": "Blue Mtns",
    "Columbia Plateau": "Columbia",
    "Klamath Mountains/California High North Coast Range": "Klamath",
    "Willamette Valley": "Willamette",
}

# Styling (task spec)
STRIP_OFFSET = 0.19
STRIP_COLORS = {"midcentury": "#6da7ec", "endcentury": "#184f95"}
BAR_WIDTH = 0.26
BAR_ALPHA = 0.22
JITTER = 0.075
JITTER_SEED = 20261003
MARKER_AREA = 12
MEAN_HALF_WIDTH = 0.12
MEAN_COLOR = "#1f1f1d"
AXIS_COLOR = "#6b6a66"
GRID_COLOR = "#e4e3df"
FONT_SIZE = 7.5
TICK_FONT_SIZE = 7
# Condensed sans: five three-line ecoregion labels must fit under each
# ~1.6 in panel at >= 7 pt; DejaVu Sans (matplotlib default) collides.
FONT_FAMILY = ["Liberation Sans Narrow", "Nimbus Sans Narrow", "DejaVu Sans Condensed"]

FNAME_RE = re.compile(
    r"facility_pue_wue_(?P<gcm>.+)_ssp370_(?P<archetype>ae-chiller|we-chiller|chiller-only)_"
    r"(?P<period>historical|midcentury|endcentury)_fast\.csv"
)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def _facility_table(root: Path) -> pd.DataFrame:
    """Facility id, name, lon/lat, ecoregion, LOCA2 grid cell and pressure."""
    import geopandas as gpd

    from climate_risk_dc.climate.cooling_archetypes import elevation_adjusted_pressure_pa
    from climate_risk_dc.config import load_dataset_paths
    from climate_risk_dc.geo import assign_us_l3_ecoregion, load_oregon_data_centers

    fac = load_oregon_data_centers(root / "data/data_centers_im3_pnnl/im3_open_source_data_center_atlas_v2026.02.09.csv")
    labeled = assign_us_l3_ecoregion(fac, gpd.read_file(load_dataset_paths().ecoregion_l3_shp))
    assert labeled["US_L3NAME"].notna().all(), "Facilities outside every Level III ecoregion polygon"
    out = pd.DataFrame(labeled.drop(columns="geometry"))[["id", "name", "lon", "lat", "US_L3NAME"]]
    out = out.rename(columns={"id": "facility_id", "name": "facility_name", "US_L3NAME": "ecoregion"})

    # LOCA2 grid cell: same nearest-cell rule as loca2_io.sample_datacube_at_points,
    # on the Oregon LOCA2 grid the PUE/WUE run sampled from.
    grid_file = next((root / "outputs/ensemble_pue_wue/humidity_cache").glob("hursmax_hursmin_*_historical_historical.nc"))
    with xr.open_dataset(grid_file) as grid:
        lat_vals, lon_vals = grid["lat"].values.astype(float), grid["lon"].values.astype(float)
    lon_q = np.where(out["lon"] < 0, out["lon"] + 360.0, out["lon"]) if lon_vals.max() > 180 else out["lon"].to_numpy()
    ilat = np.abs(lat_vals[:, None] - out["lat"].to_numpy()[None, :]).argmin(axis=0)
    ilon = np.abs(lon_vals[:, None] - lon_q[None, :]).argmin(axis=0)
    out["grid_cell"] = [f"{a}_{b}" for a, b in zip(ilat, ilon)]
    out["grid_lat"] = lat_vals[ilat]
    out["grid_lon"] = lon_vals[ilon]

    # Elevation-adjusted pressure, as in run_pue_wue_facilities.sample_elevation_at_points
    elev = xr.open_dataset(root / "data/elevation.LOCA_2016-04-02.nc")["Elevation"]
    e_lat, e_lon = np.asarray(elev["lat"].values), np.asarray(elev["lon"].values)
    lon_e = np.where(out["lon"] < 0, out["lon"] + 360.0, out["lon"]) if e_lon.max() > 180 else out["lon"].to_numpy()
    eilat = np.abs(e_lat[:, None] - out["lat"].to_numpy()[None, :]).argmin(axis=0)
    eilon = np.abs(e_lon[:, None] - lon_e[None, :]).argmin(axis=0)
    out["pressure_pa"] = elevation_adjusted_pressure_pa(elev.values[eilat, eilon])
    return out.reset_index(drop=True)


def _read_one_file(path: Path, pressure: dict[int, float]) -> pd.DataFrame:
    """Per-facility values from one daily file, long format (facility_id, metric, value).

    Every file: annual mean PUE/WUE, averaged over the years in the window.
    Case 1 files additionally: annual Tmax > 35 C counts (cross-check), and,
    for the historical window, annual and May-Oct mean afternoon wet-bulb.
    """
    from climate_risk_dc.climate.cooling_archetypes_fast import _wetbulb_from_rh

    meta = FNAME_RE.match(path.name).groupdict()
    is_case1 = meta["archetype"] == "ae-chiller"
    cols = ["facility_id", "date", "pue_daily", "wue_daily"] + (["tasmax_c", "hursmin_pct"] if is_case1 else [])
    df = pd.read_csv(path, usecols=cols)
    df["year"] = df["date"].str.slice(0, 4).astype(int)
    df["month"] = df["date"].str.slice(5, 7).astype(int)
    start, end = PERIODS[meta["period"]]
    if df["year"].min() != start or df["year"].max() != end or df["year"].nunique() != end - start + 1:
        raise ValueError(f"{path.name}: years {df['year'].min()}-{df['year'].max()} do not match {start}-{end}")

    agg = {"pue": ("pue_daily", "mean"), "wue": ("wue_daily", "mean")}
    if is_case1:
        df["hot"] = df["tasmax_c"] > THRESHOLD_C
        agg["exceedance_days_check"] = ("hot", "sum")
        if meta["period"] == "historical":
            df["wetbulb"] = _wetbulb_from_rh(
                df["tasmax_c"].to_numpy(), df["hursmin_pct"].to_numpy() / 100.0, df["facility_id"].map(pressure).to_numpy()
            )
            agg["wetbulb_annual"] = ("wetbulb", "mean")
    annual = df.groupby(["facility_id", "year"]).agg(**agg)
    if is_case1 and meta["period"] == "historical":
        warm = df[df["month"].isin(WARM_SEASON_MONTHS)].groupby(["facility_id", "year"])["wetbulb"].mean()
        annual["wetbulb_may_oct"] = warm
    period_mean = annual.groupby(level="facility_id").mean()
    long = period_mean.reset_index().melt(id_vars="facility_id", var_name="metric", value_name="value")
    long["gcm"], long["archetype"], long["period"] = meta["gcm"], meta["archetype"], meta["period"]
    return long


def _model_values(root: Path, facilities: pd.DataFrame, workers: int, refresh: bool) -> pd.DataFrame:
    """Per (facility, gcm, archetype, period, metric) values, cached after the first read."""
    cache = _out_dir(root) / "model_values_cache.csv"
    if cache.exists() and not refresh:
        print(f"Reading cached per-model values: {cache}")
        return pd.read_csv(cache)

    files = sorted((root / "outputs/ensemble_pue_wue_facilities").glob("facility_pue_wue_*_fast.csv"))
    if len(files) != N_MODELS * 9:
        raise FileNotFoundError(f"Expected {N_MODELS * 9} daily PUE/WUE files, found {len(files)}")
    pressure = dict(zip(facilities["facility_id"], facilities["pressure_pa"]))
    print(f"Reading {len(files)} daily PUE/WUE files with {workers} workers...")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        parts = list(pool.map(_read_one_file, files, [pressure] * len(files)))
    pue_wue = pd.concat(parts, ignore_index=True)

    heat_dir = root / "outputs/month2_heat_ensemble/facilities"
    heat_parts = []
    for f in sorted(heat_dir.glob("facility_heat_indices_*_ssp370_*.csv")):
        gcm, period = re.match(r"facility_heat_indices_(.+)_ssp370_(\w+)\.csv", f.name).groups()
        h = pd.read_csv(f, usecols=["facility_id", "tmax_threshold_days"]).rename(columns={"tmax_threshold_days": "value"})
        h["metric"], h["gcm"], h["archetype"], h["period"] = "exceedance_days", gcm, "all", period
        heat_parts.append(h)
    heat = pd.concat(heat_parts, ignore_index=True)
    if heat["gcm"].nunique() != N_MODELS:
        raise ValueError(f"Expected {N_MODELS} models in heat outputs, found {heat['gcm'].nunique()}")

    out = pd.concat([pue_wue, heat], ignore_index=True)
    out.to_csv(cache, index=False)
    print(f"Wrote {cache} ({len(out)} rows)")
    return out


def _check_units(values: pd.DataFrame) -> None:
    """Plausible-range checks on load, following the repo's _log_unit_sanity pattern."""
    bounds = {"pue": (1.0, 3.0), "wue": (0.0, 20.0), "exceedance_days": (0.0, 366.0),
              "exceedance_days_check": (0.0, 366.0), "wetbulb_annual": (-20.0, 35.0), "wetbulb_may_oct": (-20.0, 35.0)}
    for metric, (lo, hi) in bounds.items():
        v = values.loc[values["metric"] == metric, "value"]
        if v.isna().any() or (v < lo).any() or (v > hi).any():
            raise ValueError(f"{metric}: values outside [{lo}, {hi}] or NaN (range {v.min()}-{v.max()})")


# ---------------------------------------------------------------------------
# Computation
# ---------------------------------------------------------------------------


def _out_dir(root: Path) -> Path:
    d = root / "outputs/ensemble_pue_wue_facilities/analysis/fig_show_everything"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _per_model_changes(values: pd.DataFrame) -> pd.DataFrame:
    """Future minus historical, per facility and model. Columns: facility_id, gcm, case, metric, period, change."""
    v = values[values["metric"].isin(["pue", "wue", "exceedance_days"])].copy()
    v["case"] = v["archetype"].map(CASES).fillna("all")
    wide = v.pivot_table(index=["facility_id", "gcm", "case", "metric"], columns="period", values="value")
    changes = pd.concat(
        [(wide[p] - wide["historical"]).rename("change").reset_index().assign(period=p) for p in FUTURE_PERIODS],
        ignore_index=True,
    )
    assert changes["change"].notna().all()
    counts = changes.groupby(["facility_id", "case", "metric", "period"])["gcm"].nunique()
    assert (counts == N_MODELS).all(), "Not every facility/case/metric/period has all 20 models"
    return changes


def _wetbulb(values: pd.DataFrame, facilities: pd.DataFrame) -> pd.DataFrame:
    """Per-facility historical afternoon wet-bulb (mean across models), both windows."""
    wb = values[values["metric"].isin(["wetbulb_annual", "wetbulb_may_oct"])]
    wb = wb.pivot_table(index=["facility_id", "gcm"], columns="metric", values="value").groupby(level="facility_id").mean()
    return facilities.merge(wb.reset_index(), on="facility_id")


def _ecoregion_order(fac_wb: pd.DataFrame, column: str) -> pd.Series:
    return fac_wb.groupby("ecoregion")[column].mean().sort_values()


def _summaries(changes: pd.DataFrame, fac_wb: pd.DataFrame, wb_col: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    keys = ["case", "metric", "period"]
    ch = changes.merge(fac_wb[["facility_id", "ecoregion", "grid_cell"]], on="facility_id")

    facility = ch.groupby(["facility_id", "ecoregion", "grid_cell"] + keys, as_index=False)["change"].mean()
    facility = facility.merge(fac_wb[["facility_id", "facility_name"]], on="facility_id")

    per_model = ch.groupby(["ecoregion", "gcm"] + keys, as_index=False)["change"].mean()

    def spread(s: pd.Series) -> pd.Series:
        a = s.to_numpy()
        assert a.size == N_MODELS
        return pd.Series({
            "p25": np.percentile(a, 25), "p75": np.percentile(a, 75),
            "p05": np.percentile(a, 5), "p95": np.percentile(a, 95),
            "min": a.min(), "max": a.max(),
        })

    band = per_model.groupby(["ecoregion"] + keys)["change"].apply(spread).unstack().reset_index()
    mean = facility.groupby(["ecoregion"] + keys, as_index=False)["change"].mean().rename(columns={"change": "mean"})
    eco_meta = fac_wb.groupby("ecoregion").agg(
        n_facilities=("facility_id", "size"), n_grid_cells=("grid_cell", "nunique"), wetbulb_c=(wb_col, "mean")
    ).reset_index()
    ecoregion = mean.merge(band, on=["ecoregion"] + keys).merge(eco_meta, on="ecoregion")
    ecoregion = ecoregion[["ecoregion"] + keys + ["mean", "p25", "p75", "p05", "p95", "min", "max", "n_facilities", "n_grid_cells", "wetbulb_c"]]
    return facility, ecoregion, per_model


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

PANELS = [
    ("a", "Exceedance days (all designs)", "all", "exceedance_days"),
    ("b", "PUE, Case 1", "Case 1", "pue"),
    ("c", "PUE, Case 2", "Case 2", "pue"),
    ("d", "PUE, Case 5", "Case 5", "pue"),
    ("e", "WUE, Case 1", "Case 1", "wue"),
    ("f", "WUE, Case 2", "Case 2", "wue"),
    ("g", "WUE, Case 5", "Case 5", "wue"),
]
Y_LABELS = {
    "exceedance_days": "Change in annual threshold-exceedance days",
    "pue": "Change in annual mean PUE",
    "wue": "Change in annual mean WUE (L/kWh)",
}


def _style_axis(ax: plt.Axes) -> None:
    ax.axhline(0, color=AXIS_COLOR, linewidth=0.8, zorder=0.5)
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.6, zorder=0)
    ax.grid(axis="x", visible=False)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS_COLOR)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=AXIS_COLOR, labelcolor="#1f1f1d", labelsize=TICK_FONT_SIZE, width=0.8, length=3)


def _draw_panel(ax, facility, ecoregion, case, metric, order, rng) -> None:
    for xi, eco in enumerate(order):
        for period, sign in zip(FUTURE_PERIODS, (-1, 1)):
            xc = xi + sign * STRIP_OFFSET
            color = STRIP_COLORS[period]
            row = ecoregion[(ecoregion["ecoregion"] == eco) & (ecoregion["case"] == case)
                            & (ecoregion["metric"] == metric) & (ecoregion["period"] == period)].iloc[0]
            lo, hi = row[BAND_COLUMNS[BAND][0]], row[BAND_COLUMNS[BAND][1]]
            ax.bar(xc, hi - lo, bottom=lo, width=BAR_WIDTH, color=color, alpha=BAR_ALPHA, linewidth=0, zorder=1)
            pts = facility[(facility["ecoregion"] == eco) & (facility["case"] == case)
                           & (facility["metric"] == metric) & (facility["period"] == period)]["change"].to_numpy()
            jitter = rng.uniform(-JITTER, JITTER, size=pts.size)
            ax.scatter(xc + jitter, pts, s=MARKER_AREA, color=color, edgecolors="white", linewidths=0.6, zorder=3)
            ax.plot([xc - MEAN_HALF_WIDTH, xc + MEAN_HALF_WIDTH], [row["mean"]] * 2, color=MEAN_COLOR,
                    linewidth=1.6, solid_capstyle="round", zorder=4)


def plot(facility: pd.DataFrame, ecoregion: pd.DataFrame, order: list[str], fig_dir: Path) -> None:
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": FONT_FAMILY, "font.size": FONT_SIZE, "axes.titlesize": FONT_SIZE, "axes.labelsize": FONT_SIZE,
                         "legend.fontsize": FONT_SIZE, "pdf.fonttype": 42, "axes.labelcolor": "#1f1f1d"})
    rng = np.random.default_rng(JITTER_SEED)
    meta = ecoregion.drop_duplicates("ecoregion").set_index("ecoregion")
    tick_labels = [f"{ECOREGION_SHORT[e]}\n{meta.loc[e, 'wetbulb_c']:.1f} °C\nn = {meta.loc[e, 'n_facilities']}" for e in order]

    fig = plt.figure(figsize=(7.5, 5.5))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.0, 3.17], wspace=0.24, left=0.058, right=0.997, top=0.955, bottom=0.165)
    ax_a = fig.add_subplot(outer[0, 0])
    inner = outer[0, 1].subgridspec(2, 3, wspace=0.05, hspace=0.42)
    axes = {"a": ax_a}
    for r, letters in enumerate(("bcd", "efg")):
        first = None
        for c, letter in enumerate(letters):
            ax = fig.add_subplot(inner[r, c], sharey=first if (SHARE_Y_ROWS and first is not None) else None)
            axes[letter] = ax
            first = first or ax

    for letter, title, case, metric in PANELS:
        ax = axes[letter]
        _draw_panel(ax, facility, ecoregion, case, metric, order, rng)
        _style_axis(ax)
        ax.set_title(f"({letter}) {title}", loc="left", pad=4)
        ax.set_xlim(-0.5, len(order) - 0.5)
        ax.set_xticks(range(len(order)))
        if letter in "bcd":
            ax.tick_params(axis="x", labelbottom=False)
        else:
            ax.set_xticklabels(tick_labels, fontsize=TICK_FONT_SIZE, linespacing=1.15)
        if letter in "abe":
            ax.set_ylabel(Y_LABELS[metric])
        elif SHARE_Y_ROWS:
            ax.tick_params(axis="y", labelleft=False)

    handles = [
        Line2D([], [], linestyle="none", marker="o", markersize=4.5, markerfacecolor=STRIP_COLORS["midcentury"],
               markeredgecolor="white", markeredgewidth=0.6),
        Line2D([], [], linestyle="none", marker="o", markersize=4.5, markerfacecolor=STRIP_COLORS["endcentury"],
               markeredgecolor="white", markeredgewidth=0.6),
        Line2D([], [], color=MEAN_COLOR, linewidth=1.6, solid_capstyle="round"),
        Patch(facecolor=AXIS_COLOR, alpha=BAR_ALPHA + 0.08, linewidth=0),
    ]
    labels = [PERIOD_LEGEND["midcentury"], PERIOD_LEGEND["endcentury"], "Ecoregion mean of facilities",
              f"{BAND_LEGEND[BAND]} of ecoregion mean across {N_MODELS} models"]
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0),
               handlelength=1.6, columnspacing=1.4, handletextpad=0.5)

    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_show_everything.{ext}", dpi=300)
    plt.close(fig)


def caption() -> str:
    return (
        "Projected change in (a) annual number of threshold-exceedance days, (b–d) power usage effectiveness "
        "(PUE), and (e–g) water usage effectiveness (WUE, L/kWh) for cooling-system Cases 1, 2, and 5. Changes "
        "are between the historical period (1985–2014) and mid-century (2045–2074; light blue) or "
        "end-of-century (2075–2100; dark blue). PUE and WUE are annual means of daily values. Ecoregions are "
        "ordered by historical "
        + ("annual mean" if WETBULB_WINDOW == "annual" else "May–October mean")
        + " afternoon wet-bulb temperature, shown below each name with the number of data centers (n). Points are "
        f"facilities (mean of {N_MODELS} climate models), jittered horizontally; facilities within the same LOCA2 "
        "grid cell have identical values. Black lines are ecoregion means, and shaded bars are "
        f"{BAND_CAPTION[BAND]} of the ecoregion mean across the {N_MODELS} models."
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

DRAFT_TABLE1 = {  # draft Table 1 (docs/results-draft.md), Case 1 mid-century
    ("Columbia Plateau", "pue"): 0.017, ("Columbia Plateau", "wue"): 0.21,
    ("Willamette Valley", "pue"): 0.021, ("Willamette Valley", "wue"): 0.27,
}


def report(values, changes, facility, ecoregion, per_model, fac_wb, order, out_dir) -> None:
    pd.set_option("display.width", 200)
    wb_col = f"wetbulb_{WETBULB_WINDOW}"

    print("\n=== Ecoregion order (ascending historical afternoon wet-bulb) ===")
    meta = fac_wb.groupby("ecoregion").agg(
        wetbulb_annual=("wetbulb_annual", "mean"), wetbulb_may_oct=("wetbulb_may_oct", "mean"),
        n_facilities=("facility_id", "size"), n_grid_cells=("grid_cell", "nunique"))
    print(meta.loc[order].round(2).to_string())
    alt = "may_oct" if WETBULB_WINDOW == "annual" else "annual"
    alt_order = list(_ecoregion_order(fac_wb, f"wetbulb_{alt}").index)
    print(f"Order using {WETBULB_WINDOW}: {[ECOREGION_SHORT[e] for e in order]}")
    print(f"Order using {alt}: {[ECOREGION_SHORT[e] for e in alt_order]}  -> {'same' if alt_order == order else 'DIFFERENT'}")

    print(f"\n=== Ecoregion mean and plotted band ({BAND}) per panel and period ===")
    lo, hi = BAND_COLUMNS[BAND]
    for letter, title, case, metric in PANELS:
        print(f"({letter}) {title}")
        sub = ecoregion[(ecoregion["case"] == case) & (ecoregion["metric"] == metric)].set_index(["ecoregion", "period"])
        for eco in order:
            line = "   ".join(
                f"{p[:3]}: {sub.loc[(eco, p), 'mean']:+.4g} [{sub.loc[(eco, p), lo]:+.4g}, {sub.loc[(eco, p), hi]:+.4g}]"
                for p in FUTURE_PERIODS)
            print(f"   {ECOREGION_SHORT[eco]:<11} {line}")

    print("\n=== Exceedance-day cross-check: stored heat-ensemble vs recomputed from daily tasmax ===")
    ex = values[values["metric"].isin(["exceedance_days", "exceedance_days_check"])]
    ex = ex.pivot_table(index=["facility_id", "gcm", "period"], columns="metric", values="value")
    diff = (ex["exceedance_days"] - ex["exceedance_days_check"]).abs()
    print(f"max |difference| = {diff.max():.2e} days/yr over {len(diff)} facility-model-period values "
          f"({(diff > 1e-3).sum()} differ by more than 0.001 day)")

    print("\n=== Reconciliation with draft Table 1 (Case 1, mid-century) ===")
    pooled = pd.read_csv(out_dir.parent / "ensemble_delta.csv")
    pooled = pooled[pooled["archetype"] == "ae-chiller"].merge(fac_wb[["facility_id", "ecoregion"]], on="facility_id")
    for (eco, metric), draft in DRAFT_TABLE1.items():
        new = ecoregion[(ecoregion["ecoregion"] == eco) & (ecoregion["case"] == "Case 1") & (ecoregion["metric"] == metric)
                        & (ecoregion["period"] == "midcentury")]["mean"].iloc[0]
        pool = pooled.groupby(["ecoregion", "facility_id"])[f"{metric}_mean_delta_midcentury"].mean().loc[eco].mean()
        decimals = 3 if metric == "pue" else 2
        status = "match" if round(new, decimals) == draft else "MISMATCH"
        print(f"   {eco:<18} {metric.upper()}: draft {draft:+.{decimals}f}  new annual {new:+.4f} "
              f"(rounds to {new:+.{decimals}f}, {status}); pooled all-days mean {pool:+.4f}")

    print("\n=== Model agreement: models (of 20) where ecoregion A's change exceeds ecoregion B's ===")
    pm = per_model.pivot_table(index=["case", "metric", "period", "gcm"], columns="ecoregion", values="change")
    rows = []
    for a, b in [("Willamette Valley", "Columbia Plateau"), ("Columbia Plateau", "Blue Mountains")]:
        cnt = (pm[a] > pm[b]).groupby(level=["case", "metric", "period"]).sum()
        for (case, metric, period), n in cnt.items():
            rows.append({"comparison": f"{ECOREGION_SHORT[a]} > {ECOREGION_SHORT[b]}", "case": case, "metric": metric,
                         "period": period, "n_models": int(n), "of": N_MODELS})
    agree = pd.DataFrame(rows)
    agree.to_csv(out_dir / "model_agreement.csv", index=False)
    print(agree.pivot_table(index=["comparison", "metric", "case"], columns="period", values="n_models")
          [FUTURE_PERIODS].to_string())


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out_dir = _out_dir(root)
    fig_dir = root / "docs" / "results-figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    facilities = _facility_table(root)
    values = _model_values(root, facilities, args.workers, args.refresh_cache)
    _check_units(values)

    wb_col = f"wetbulb_{WETBULB_WINDOW}"
    fac_wb = _wetbulb(values, facilities)
    order = list(_ecoregion_order(fac_wb, wb_col).index)
    changes = _per_model_changes(values)
    facility, ecoregion, per_model = _summaries(changes, fac_wb, wb_col)

    facility[["facility_id", "facility_name", "ecoregion", "case", "metric", "period", "change", "grid_cell"]].sort_values(
        ["metric", "case", "period", "ecoregion", "facility_id"]).to_csv(out_dir / "facility_changes.csv", index=False)
    eco_out = ecoregion.assign(order=ecoregion["ecoregion"].map({e: i for i, e in enumerate(order)}))
    eco_out.sort_values(["metric", "case", "period", "order"]).drop(columns="order").to_csv(
        out_dir / "ecoregion_changes.csv", index=False)
    facilities.drop(columns="pressure_pa").to_csv(out_dir / "facility_grid_cells.csv", index=False)

    plot(facility, ecoregion, order, fig_dir)
    (fig_dir / "fig_show_everything_caption.txt").write_text(caption() + "\n")

    pd.DataFrame([{
        "band": BAND, "wetbulb_window": WETBULB_WINDOW, "share_y_rows": SHARE_Y_ROWS, "n_models": N_MODELS,
        "threshold_c": THRESHOLD_C, "periods": "; ".join(f"{k}={a}-{b}" for k, (a, b) in PERIODS.items()),
        "cases": "; ".join(f"{v}={k}" for k, v in CASES.items()), "jitter_seed": JITTER_SEED,
        "ecoregion_order": "; ".join(order), "scenario": "ssp370",
    }]).to_csv(out_dir / "run_manifest.csv", index=False)

    report(values, changes, facility, ecoregion, per_model, fac_wb, order, out_dir)
    print(f"\nCaption:\n{caption()}")
    print(f"\nWrote {fig_dir / 'fig_show_everything.png'} (+ .pdf) and CSVs under {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the 'show everything' heat-chapter results figure.")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    parser.add_argument("--workers", type=int, default=6, help="Processes for reading the 180 daily files (default 6)")
    parser.add_argument("--refresh-cache", action="store_true", help="Re-read the daily files even if the cache exists")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
