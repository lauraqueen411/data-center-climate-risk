"""Figure 3 of the heat risk chapter: maps of acute and chronic heat risk
(docs/tasks/result-figures.md).

Message: acute risk (exceedance days) grows where it already is, while
chronic risk (annual mean PUE, Case 1) shifts towards the Willamette Valley.

- (a) historical annual exceedance days (Tmax > 35 C), 1985-2014
- (b) historical annual mean PUE, Case 1
- (c) change in exceedance days, 2045-2074 minus 1985-2014
- (d) change in annual mean PUE, Case 1, 2045-2074 minus 1985-2014

Each ecoregion that contains a facility is filled with its ecoregion mean;
the other Level III ecoregions are grey with outlines. Facilities are dots
colored by their own value on the same scale as the fill.

Inputs (written by ``scripts/fig_show_everything.py``, nothing re-run):
- ``outputs/ensemble_pue_wue_facilities/analysis/fig_show_everything/model_values_cache.csv``:
  period-mean values per facility, model, case and period.
- ``.../fig_show_everything/facility_grid_cells.csv``: facility lon/lat and
  Level III ecoregion.
- Level III ecoregion polygons (``load_dataset_paths().ecoregion_l3_shp``)
  and Natural Earth 50m state lines and coastline (via cartopy).

Outputs:
- ``docs/results-figures/fig3_maps.{png,pdf}`` and ``fig3_maps_caption.txt``.
- ``outputs/ensemble_pue_wue_facilities/analysis/fig3_maps/``:
  ``facility_values.csv``, ``ecoregion_values.csv``, ``run_manifest.csv``.

Assumptions log
----------------
- **Values.** Per facility and model, the period mean of annual values (as
  in Figure 2). Change is future minus historical per facility and model,
  then the mean across the 20 models.
- **Ecoregion fill** is the mean of that ecoregion's facility values, equal
  weight per facility (the Figure 2 ecoregion mean). ECSF and KM contain one
  facility each, so their fill equals that facility's value.
- **Color scales** run from the minimum to the maximum facility value in
  each panel (ecoregion means lie inside that range), so every dot is on
  scale and no colorbar needs an extension.
- **Ecoregion labels** sit at each polygon's representative point, moved by
  ``LABEL_OFFSET`` (degrees) where the point falls on facility dots.
- **Colormaps** drop their palest 15% (``CMAP_START``) so the lowest-valued
  ecoregion is not mistaken for the grey no-facility fill.
- **Ecoregion IQR** across models (of the per-model ecoregion mean) is
  written to ``ecoregion_values.csv`` for the text; maps show no spread.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from climate_risk_dc import figstyle
from climate_risk_dc.config import load_dataset_paths

N_MODELS = 20
PERIODS = {"historical": (1985, 2014), "midcentury": (2045, 2074)}
ECOREGION_ABBR = {
    "Eastern Cascades Slopes and Foothills": "ECSF",
    "Blue Mountains": "BM",
    "Columbia Plateau": "CP",
    "Klamath Mountains/California High North Coast Range": "KM",
    "Willamette Valley": "WV",
}
# (metric, archetype) in model_values_cache.csv
METRICS = {"exceedance_days": ("exceedance_days", "all"), "pue_case1": ("pue", "ae-chiller")}
PANELS = [
    # letter, metric, kind, title, colorbar label, cmap
    ("a", "exceedance_days", "historical", "Exceedance days", "Days per year with Tmax > 35 °C", "Oranges"),
    ("b", "pue_case1", "historical", "PUE, Case 1", "Annual mean PUE", "Purples"),
    ("c", "exceedance_days", "change", "Exceedance days", "Change in days per year", "Oranges"),
    ("d", "pue_case1", "change", "PUE, Case 1", "Change in annual mean PUE", "Purples"),
]
ROW_LABELS = {"historical": "Historical, 1985–2014", "change": "Change, 2045–2074 minus 1985–2014"}
LABEL_OFFSET = {"ECSF": (0.0, -0.55), "WV": (0.0, -0.4), "CP": (0.3, 0.05)}
CBAR_NBINS = {"exceedance_days": 5, "pue_case1": 4}
CMAP_START = 0.15  # fraction of the colormap skipped at the pale end
CBAR_FRACTION = 0.85  # colorbar width as a fraction of the map width
# Layout, inches
L = {"top": 0.05, "side": 0.1, "col_gap": 0.3, "row_label": 0.26, "title": 0.24,
     "cbar_gap": 0.08, "cbar_h": 0.11, "cbar_text": 0.52, "row_gap": 0.12}

OTHER_FILL = "#e6e6e6"
ECO_EDGE = "#555555"
ECO_EDGE_WIDTH = 0.4
DOT_AREA = 10
DOT_EDGE_WIDTH = 0.3

BOUNDS = {"exceedance_days": (0.0, 366.0), "pue": (1.0, 3.0)}


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


def _analysis_dir(root: Path) -> Path:
    return root / "outputs/ensemble_pue_wue_facilities/analysis"


def _load(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    src = _analysis_dir(root) / "fig_show_everything"
    values = pd.read_csv(src / "model_values_cache.csv")
    facilities = pd.read_csv(src / "facility_grid_cells.csv")
    for metric, (lo, hi) in BOUNDS.items():
        v = values.loc[values["metric"] == metric, "value"]
        if v.isna().any() or (v < lo).any() or (v > hi).any():
            raise ValueError(f"{metric}: values outside [{lo}, {hi}] or NaN (range {v.min()}-{v.max()})")
    return values, facilities


def _per_model(values: pd.DataFrame) -> pd.DataFrame:
    """Per facility and model: historical value and mid-century change. Columns: facility_id, gcm, metric, kind, value."""
    parts = []
    for name, (metric, archetype) in METRICS.items():
        v = values[(values["metric"] == metric) & (values["archetype"] == archetype) & values["period"].isin(PERIODS)]
        wide = v.pivot_table(index=["facility_id", "gcm"], columns="period", values="value")
        hist = wide["historical"].rename("value").reset_index().assign(metric=name, kind="historical")
        change = (wide["midcentury"] - wide["historical"]).rename("value").reset_index().assign(metric=name, kind="change")
        parts += [hist, change]
    out = pd.concat(parts, ignore_index=True)
    assert out["value"].notna().all()
    assert (out.groupby(["facility_id", "metric", "kind"])["gcm"].nunique() == N_MODELS).all(), "Missing models"
    return out


def _summaries(per_model: pd.DataFrame, facilities: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Facility values (ensemble mean) and ecoregion values (mean of facilities, plus IQR across models)."""
    pm = per_model.merge(facilities[["facility_id", "ecoregion"]], on="facility_id")
    facility = pm.groupby(["facility_id", "metric", "kind"], as_index=False)["value"].mean()
    facility = facilities[["facility_id", "facility_name", "lon", "lat", "ecoregion", "grid_cell"]].merge(facility, on="facility_id")

    eco_mean = facility.groupby(["ecoregion", "metric", "kind"])["value"].agg(["mean", "size"]).rename(columns={"size": "n_facilities"})
    per_model_eco = pm.groupby(["ecoregion", "metric", "kind", "gcm"])["value"].mean()
    q = per_model_eco.groupby(level=["ecoregion", "metric", "kind"]).quantile([0.25, 0.75]).unstack()
    q.columns = ["p25", "p75"]
    ecoregion = eco_mean.join(q).reset_index()
    ecoregion["abbr"] = ecoregion["ecoregion"].map(ECOREGION_ABBR)
    return facility, ecoregion


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------


def _panel(fig, rect, cax_rect, polygons, facility, ecoregion, letter, metric, kind, title, cbar_label, cmap_name):
    ax = figstyle.oregon_map_axes(fig, rect)
    cmap = _truncated(cmap_name)
    fac = facility[(facility["metric"] == metric) & (facility["kind"] == kind)]
    eco = ecoregion[(ecoregion["metric"] == metric) & (ecoregion["kind"] == kind)].set_index("ecoregion")["mean"]
    norm = plt.Normalize(fac["value"].min(), fac["value"].max())

    # Project with geopandas rather than letting cartopy reproject each polygon at draw time.
    projected = polygons.to_crs(ax.projection.proj4_init)
    for name, geom in zip(projected["US_L3NAME"], projected.geometry):
        face = cmap(norm(eco[name])) if name in eco.index else OTHER_FILL
        ax.add_geometries([geom], crs=ax.projection, facecolor=face, edgecolor=ECO_EDGE, linewidth=ECO_EDGE_WIDTH, zorder=1)

    # Draw the highest values last so the extremes stay visible where dots overlap.
    fac = fac.sort_values("value")
    sc = ax.scatter(fac["lon"], fac["lat"], c=fac["value"], cmap=cmap, norm=norm, s=DOT_AREA,
                    edgecolors=figstyle.INK, linewidths=DOT_EDGE_WIDTH, transform=ccrs.PlateCarree(), zorder=3)

    for name, abbr in ECOREGION_ABBR.items():
        pt = polygons.loc[polygons["US_L3NAME"] == name].geometry.iloc[0].representative_point()
        dx, dy = LABEL_OFFSET.get(abbr, (0.0, 0.0))
        ax.text(pt.x + dx, pt.y + dy, abbr, transform=ccrs.PlateCarree(), ha="center", va="center",
                fontsize=figstyle.TICK_PT, fontweight="bold", color=figstyle.INK, zorder=4,
                path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])

    ax.set_title(f"({letter}) {title}", loc="left", pad=3)
    cax = fig.add_axes(cax_rect)
    cbar = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cbar.set_label(cbar_label, labelpad=2)
    cbar.outline.set_linewidth(figstyle.SPINE_WIDTH)
    cbar.ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=CBAR_NBINS[metric], min_n_ticks=3))
    return ax


def _truncated(cmap_name: str) -> matplotlib.colors.Colormap:
    """The named colormap without its palest end, so the lowest value does not read as an empty polygon."""
    base = plt.get_cmap(cmap_name)
    return matplotlib.colors.ListedColormap(base(np.linspace(CMAP_START, 1.0, 256)), name=f"{cmap_name}_trunc")


def _map_aspect() -> float:
    """Height / width of the map frame (projected extent of figstyle.oregon_map_axes)."""
    fig = plt.figure()
    x0, x1, y0, y1 = figstyle.oregon_map_axes(fig, 111).get_extent()
    plt.close(fig)
    return (y1 - y0) / (x1 - x0)


def plot(facility: pd.DataFrame, ecoregion: pd.DataFrame, polygons: gpd.GeoDataFrame, path: Path) -> None:
    """Explicit layout in inches: the map frame has a fixed aspect, so gridspec cells would leave gaps."""
    figstyle.apply()
    W = figstyle.FIG_WIDTH_IN
    map_w = (W - 2 * L["side"] - L["col_gap"]) / 2
    map_h = map_w * _map_aspect()
    row_h = L["row_label"] + L["title"] + map_h + L["cbar_gap"] + L["cbar_h"] + L["cbar_text"]
    H = L["top"] + 2 * row_h + L["row_gap"]
    fig = plt.figure(figsize=(W, H))

    def rect(x, y, w, h):  # inches from the bottom-left to figure fraction
        return [x / W, y / H, w / W, h / H]

    for i, (letter, metric, kind, title, cbar_label, cmap) in enumerate(PANELS):
        r, c = divmod(i, 2)
        row_top = H - L["top"] - r * (row_h + L["row_gap"])
        x = L["side"] + c * (map_w + L["col_gap"])
        map_y = row_top - L["row_label"] - L["title"] - map_h
        cbar_y = map_y - L["cbar_gap"] - L["cbar_h"]
        inset = map_w * (1 - CBAR_FRACTION) / 2
        _panel(fig, rect(x, map_y, map_w, map_h), rect(x + inset, cbar_y, map_w - 2 * inset, L["cbar_h"]),
               polygons, facility, ecoregion, letter, metric, kind, title, cbar_label, cmap)
        if c == 0:
            fig.text(0.5, (row_top - L["row_label"] / 2) / H, ROW_LABELS[kind], ha="center", va="center",
                     fontsize=figstyle.LABEL_PT, fontweight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(path.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)


def caption(ecoregion: pd.DataFrame) -> str:
    n = ecoregion.drop_duplicates("ecoregion").set_index("abbr")["n_facilities"]
    return (
        "Historical (1985–2014) values and mid-century change (2045–2074 minus 1985–2014) in (a, c) annual number "
        "of threshold-exceedance days (daily maximum temperature above 35 °C) and (b, d) annual mean power usage "
        "effectiveness (PUE) for cooling-system Case 1. Values are the mean of 20 climate models (SSP3-7.0). "
        "Ecoregions containing facilities are filled with the mean of their facilities: Eastern Cascades Slopes and "
        f"Foothills (ECSF; n = {n['ECSF']}), Blue Mountains (BM; n = {n['BM']}), Columbia Plateau (CP; n = {n['CP']}), "
        f"Klamath Mountains (KM; n = {n['KM']}) and Willamette Valley (WV; n = {n['WV']}). Dots are facilities, "
        "colored on the same scale; other Level III ecoregions are shown in grey."
    )


def print_summary(ecoregion: pd.DataFrame, facility: pd.DataFrame) -> None:
    print("\n=== Ecoregion historical value and mid-century change (ensemble mean [IQR across models]) ===")
    for metric, fmt in (("exceedance_days", "{:.1f}"), ("pue_case1", "{:.4f}")):
        print(f"{metric}:")
        e = ecoregion[ecoregion["metric"] == metric].set_index(["abbr", "kind"])
        for abbr in ECOREGION_ABBR.values():
            h, c = e.loc[(abbr, "historical")], e.loc[(abbr, "change")]
            print(f"   {abbr:5s} n={int(h['n_facilities']):3d}  historical {fmt.format(h['mean'])} "
                  f"[{fmt.format(h['p25'])}, {fmt.format(h['p75'])}]   change +{fmt.format(c['mean'])} "
                  f"[{fmt.format(c['p25'])}, {fmt.format(c['p75'])}]")
    print("\n=== Facility range per panel ===")
    for (metric, kind), g in facility.groupby(["metric", "kind"]):
        print(f"   {metric:16s} {kind:10s} {g['value'].min():.4f} – {g['value'].max():.4f}")


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out_dir = _analysis_dir(root) / "fig3_maps"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = root / "docs" / "results-figures"

    values, facilities = _load(root)
    facility, ecoregion = _summaries(_per_model(values), facilities)
    polygons = gpd.read_file(load_dataset_paths().ecoregion_l3_shp).to_crs("EPSG:4326")
    missing = set(ECOREGION_ABBR) - set(polygons["US_L3NAME"])
    if missing:
        raise ValueError(f"Ecoregions missing from the shapefile: {missing}")

    facility.sort_values(["metric", "kind", "ecoregion", "facility_id"]).to_csv(out_dir / "facility_values.csv", index=False)
    ecoregion.sort_values(["metric", "kind", "ecoregion"]).to_csv(out_dir / "ecoregion_values.csv", index=False)
    pd.DataFrame([{
        "n_models": N_MODELS, "scenario": "ssp370",
        "periods": "; ".join(f"{k}={a}-{b}" for k, (a, b) in PERIODS.items()),
        "metrics": "exceedance_days (Tmax > 35 C); pue Case 1 = ae-chiller",
        "fill": "ecoregion mean of facility ensemble means", "color_scale": "facility min to max per panel",
    }]).to_csv(out_dir / "run_manifest.csv", index=False)

    plot(facility, ecoregion, polygons, fig_dir / "fig3_maps.png")
    (fig_dir / "fig3_maps_caption.txt").write_text(caption(ecoregion) + "\n")
    print_summary(ecoregion, facility)
    print(f"\nWrote {fig_dir / 'fig3_maps.png'} (+ .pdf) and CSVs under {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Figure 3 (heat-risk maps) of the heat chapter.")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
