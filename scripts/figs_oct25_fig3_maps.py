"""Figure 3 redraw (reviewer comment, October 2026): discrete colorbars and
facility dots that stay visible on every fill.

Same data, panels, labels, layout, figure size (6.5 in wide) and dpi (300)
as ``scripts/fig3_maps.py``, whose data functions and layout constants are
imported (``_load``, ``_per_model``, ``_summaries``, ``L``, ``PANELS``, ...),
so the plotted values are computed by the same code as the current figure.

Changes:
- **Binned colorbars** (``BoundaryNorm``) with round bin edges per panel
  (``EDGES``); no extensions, since every facility value lies inside the edges
  (asserted).
- **Colormaps** are sequential, colorblind-safe and monotonic in lightness,
  and stop at a mid-tone so no bin is near-black.
- **Dots**, two variants:
  - A: ColorBrewer YlOrBr (days) and YlGnBu (PUE); facilities drawn as white
    dots with a black edge (location only; the fill carries the value).
  - B: magma (days) and viridis (PUE), light end = low value; facilities
    colored by their own value on the same bins, with a black edge and a
    white halo.

Inputs as ``fig3_maps.py``. Outputs (``docs/results-figures/figures_oct25/``):
``fig3_maps_A.{png,pdf}``, ``fig3_maps_B.{png,pdf}``,
``fig3_bin_contrast.csv`` (per bin: color, luminance, contrast against the
dot's fill/edge/halo), ``fig3_check.csv`` and previews.

Assumptions log
----------------
- **Bin edges** are round numbers spanning the facility range of each panel
  (printed and in ``fig3_bin_contrast.csv``). A value equal to an edge goes
  in the upper bin (``BoundaryNorm`` convention).
- **Contrast test** uses the WCAG 2 contrast ratio between relative
  luminances. A dot counts as visible on a bin if its outermost ring (A:
  white fill or black edge; B: white halo or black edge) reaches a ratio of at least 3
  against the bin color.
- **Variant A drops the per-facility color**: dots mark locations only, so
  within-ecoregion facility differences are not shown (they are small
  next to the between-ecoregion differences). Variant B keeps them.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from climate_risk_dc import figstyle
from climate_risk_dc.config import load_dataset_paths
from fig3_maps import (
    CBAR_FRACTION, ECO_EDGE, ECO_EDGE_WIDTH, ECOREGION_ABBR, L, LABEL_OFFSET, OTHER_FILL, PANELS, ROW_LABELS,
    _load, _map_aspect, _per_model, _summaries,
)
from figs_oct25_common import OUT_DIR, apply_figstyle, side_by_side

EDGES = {
    ("exceedance_days", "historical"): np.arange(0, 24.01, 3),
    ("exceedance_days", "change"): np.arange(5, 35.01, 5),
    ("pue_case1", "historical"): np.round(np.arange(1.120, 1.13401, 0.002), 3),
    ("pue_case1", "change"): np.round(np.arange(0.008, 0.01701, 0.001), 3),
}
# label every n-th edge (all edges get a tick)
LABEL_EVERY = {("exceedance_days", "historical"): 1, ("exceedance_days", "change"): 1,
               ("pue_case1", "historical"): 2, ("pue_case1", "change"): 2}
TICK_FMT = {"exceedance_days": "{:.0f}", "pue_case1": "{:.3f}"}

VARIANTS = {
    "A": {"cmap": {"exceedance_days": ("YlOrBr", 0.08, 0.72), "pue_case1": ("YlGnBu", 0.08, 0.72)},
          "dot": "white"},
    "B": {"cmap": {"exceedance_days": ("magma_r", 0.02, 0.62), "pue_case1": ("viridis_r", 0.0, 0.62)},
          "dot": "value"},
}
DOT_AREA = 15
DOT_EDGE = 0.6
HALO_AREA = 38
MIN_CONTRAST = 3.0


def _luminance(rgb) -> float:
    c = np.asarray(mcolors.to_rgb(rgb))
    c = np.where(c <= 0.03928, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    return float(0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2])


def _contrast(a, b) -> float:
    la, lb = sorted((_luminance(a), _luminance(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


def _bin_colors(metric: str, kind: str, variant: str) -> tuple[np.ndarray, list]:
    name, lo, hi = VARIANTS[variant]["cmap"][metric]
    edges = EDGES[(metric, kind)]
    n = len(edges) - 1
    base = plt.get_cmap(name)
    return edges, [base(x) for x in np.linspace(lo, hi, n)]


def _panel(fig, rect, cax_rect, polygons, facility, ecoregion, letter, metric, kind, title, cbar_label, variant):
    ax = figstyle.oregon_map_axes(fig, rect)
    edges, colors = _bin_colors(metric, kind, variant)
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(edges, cmap.N)
    fac = facility[(facility["metric"] == metric) & (facility["kind"] == kind)].sort_values("value")
    eco = ecoregion[(ecoregion["metric"] == metric) & (ecoregion["kind"] == kind)].set_index("ecoregion")["mean"]
    assert fac["value"].between(edges[0], edges[-1], inclusive="left").all(), f"{metric} {kind}: value outside bins"

    projected = polygons.to_crs(ax.projection.proj4_init)
    for name, geom in zip(projected["US_L3NAME"], projected.geometry):
        face = cmap(norm(eco[name])) if name in eco.index else OTHER_FILL
        ax.add_geometries([geom], crs=ax.projection, facecolor=face, edgecolor=ECO_EDGE, linewidth=ECO_EDGE_WIDTH, zorder=1)

    pc = ccrs.PlateCarree()
    if VARIANTS[variant]["dot"] == "white":
        ax.scatter(fac["lon"], fac["lat"], s=DOT_AREA, facecolor="white", edgecolors=figstyle.INK, linewidths=DOT_EDGE,
                   transform=pc, zorder=3)
    else:
        ax.scatter(fac["lon"], fac["lat"], s=HALO_AREA, facecolor="white", edgecolors="none", transform=pc, zorder=2.9)
        ax.scatter(fac["lon"], fac["lat"], s=DOT_AREA, c=fac["value"], cmap=cmap, norm=norm, edgecolors=figstyle.INK,
                   linewidths=DOT_EDGE, transform=pc, zorder=3)

    for name, abbr in ECOREGION_ABBR.items():
        pt = polygons.loc[polygons["US_L3NAME"] == name].geometry.iloc[0].representative_point()
        dx, dy = LABEL_OFFSET.get(abbr, (0.0, 0.0))
        ax.text(pt.x + dx, pt.y + dy, abbr, transform=pc, ha="center", va="center", fontsize=figstyle.TICK_PT,
                fontweight="bold", color=figstyle.INK, zorder=4, path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])
    ax.set_title(f"({letter}) {title}", loc="left", pad=3)

    cax = fig.add_axes(cax_rect)
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal", spacing="uniform", ticks=edges)
    every = LABEL_EVERY[(metric, kind)]
    cbar.ax.set_xticklabels([TICK_FMT[metric].format(e) if i % every == 0 else "" for i, e in enumerate(edges)])
    cbar.set_label(cbar_label, labelpad=2)
    cbar.outline.set_linewidth(figstyle.SPINE_WIDTH)
    cbar.dividers = None
    return ax


def plot(facility, ecoregion, polygons, path: Path, variant: str) -> None:
    """Layout copied from fig3_maps.plot (explicit inches)."""
    apply_figstyle()
    W = figstyle.FIG_WIDTH_IN
    map_w = (W - 2 * L["side"] - L["col_gap"]) / 2
    map_h = map_w * _map_aspect()
    row_h = L["row_label"] + L["title"] + map_h + L["cbar_gap"] + L["cbar_h"] + L["cbar_text"]
    H = L["top"] + 2 * row_h + L["row_gap"]
    fig = plt.figure(figsize=(W, H))

    def rect(x, y, w, h):
        return [x / W, y / H, w / W, h / H]

    for i, (letter, metric, kind, title, cbar_label, _cmap) in enumerate(PANELS):
        r, c = divmod(i, 2)
        row_top = H - L["top"] - r * (row_h + L["row_gap"])
        x = L["side"] + c * (map_w + L["col_gap"])
        map_y = row_top - L["row_label"] - L["title"] - map_h
        cbar_y = map_y - L["cbar_gap"] - L["cbar_h"]
        inset = map_w * (1 - CBAR_FRACTION) / 2
        _panel(fig, rect(x, map_y, map_w, map_h), rect(x + inset, cbar_y, map_w - 2 * inset, L["cbar_h"]),
               polygons, facility, ecoregion, letter, metric, kind, title, cbar_label, variant)
        if c == 0:
            fig.text(0.5, (row_top - L["row_label"] / 2) / H, ROW_LABELS[kind], ha="center", va="center",
                     fontsize=figstyle.LABEL_PT, fontweight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(path.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)


def contrast_table() -> pd.DataFrame:
    rows = []
    for variant in VARIANTS:
        for (metric, kind), edges in EDGES.items():
            _, colors = _bin_colors(metric, kind, variant)
            for i, col in enumerate(colors):
                r = {"variant": variant, "metric": metric, "kind": kind, "bin": f"{edges[i]:g}–{edges[i + 1]:g}",
                     "color": mcolors.to_hex(col), "luminance": round(_luminance(col), 3),
                     "contrast_white": round(_contrast(col, "white"), 2), "contrast_black": round(_contrast(col, "black"), 2)}
                # A: white fill + black edge; B: white halo + black edge (the dot fill may equal the bin color)
                r["dot_visible"] = max(r["contrast_white"], r["contrast_black"]) >= MIN_CONTRAST
                rows.append(r)
    t = pd.DataFrame(rows)
    mono = t.groupby(["variant", "metric", "kind"])["luminance"].apply(lambda s: bool((np.diff(s.to_numpy()) < 0).all()))
    t = t.merge(mono.rename("luminance_monotonic").reset_index(), on=["variant", "metric", "kind"])
    return t


def check(root: Path, facility: pd.DataFrame, ecoregion: pd.DataFrame) -> pd.DataFrame:
    """Plotted values vs an independent recomputation from the cache and vs draft 3 Table 1 / facility facts."""
    cache = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/fig_show_everything/model_values_cache.csv")
    rows = []
    for name, (m, arch) in {"exceedance_days": ("exceedance_days", "all"), "pue_case1": ("pue", "ae-chiller")}.items():
        w = cache[(cache["metric"] == m) & (cache["archetype"] == arch)].pivot_table(
            index=["facility_id", "gcm"], columns="period", values="value")
        ind = {"historical": w["historical"].groupby(level="facility_id").mean(),
               "change": (w["midcentury"] - w["historical"]).groupby(level="facility_id").mean()}
        for kind, s in ind.items():
            p = facility[(facility["metric"] == name) & (facility["kind"] == kind)].set_index("facility_id")["value"]
            rows.append({"check": f"facility {name} {kind} vs cache recomputation", "max_abs_diff": float((p - s).abs().max())})
    t1 = pd.read_csv(root / "docs/results-figures/draft3/table1.csv")
    t1 = t1[t1["period"] == "2045–2074"].assign(abbr=lambda d: d["ecoregion"].str.split("/").str[0].map(
        {"Eastern Cascades Slopes and Foothills": "ECSF", "Blue Mountains": "BM", "Columbia Plateau": "CP",
         "Klamath Mountains": "KM", "Willamette Valley": "WV"})).set_index("abbr")
    e = ecoregion.set_index(["abbr", "metric", "kind"])["mean"]
    for col, key, dec in (("hist_exceedance_days", ("exceedance_days", "historical"), 1),
                          ("d_exceedance", ("exceedance_days", "change"), 1),
                          ("hist_pue", ("pue_case1", "historical"), 3), ("d_pue", ("pue_case1", "change"), 4)):
        diffs = [abs(e[(a, *key)] - t1.loc[a, col]) for a in t1.index]
        rows.append({"check": f"ecoregion {key[0]} {key[1]} vs draft3/table1.csv {col} (rounded to {dec} dp there)",
                     "max_abs_diff": float(max(diffs))})
    ff = pd.read_csv(root / "docs/results-figures/draft3/facility_facts.csv")
    ff = ff[ff["scope"] == "all"]
    for _, r in ff.iterrows():
        kind = "historical" if r["quantity"] == "historical" else "change"
        p = facility[(facility["metric"] == "pue_case1") & (facility["kind"] == kind)]["value"]
        v = p.min() if r["which"] == "min" else p.max()
        rows.append({"check": f"facility PUE {kind} {r['which']} vs draft3/facility_facts.csv ({r['value']})",
                     "max_abs_diff": float(abs(v - r["value"]))})
    return pd.DataFrame(rows)


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out = root / OUT_DIR
    out.mkdir(parents=True, exist_ok=True)
    values, facilities = _load(root)
    facility, ecoregion = _summaries(_per_model(values), facilities)
    polygons = gpd.read_file(load_dataset_paths().ecoregion_l3_shp).to_crs("EPSG:4326")

    ct = contrast_table()
    ct.to_csv(out / "fig3_bin_contrast.csv", index=False)
    chk = check(root, facility, ecoregion)
    chk.to_csv(out / "fig3_check.csv", index=False)
    for v in VARIANTS:
        plot(facility, ecoregion, polygons, out / f"fig3_maps_{v}.png", v)
        side_by_side(root / "docs/results-figures/fig3_maps.png", out / f"fig3_maps_{v}.png", out / f"preview_fig3_{v}_old_vs_new.png")

    for (m, k), e in EDGES.items():
        print(f"{m:16s} {k:10s} edges {', '.join(f'{x:g}' for x in e)}")
    print(ct.groupby(["variant", "metric", "kind"]).agg(all_visible=("dot_visible", "all"), monotonic=("luminance_monotonic", "first"),
                                                         min_white=("contrast_white", "min"), darkest_lum=("luminance", "min")))
    with pd.option_context("display.max_colwidth", 120, "display.width", 200):
        print(chk.to_string(index=False))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Figure 3 redraw: binned colorbars, high-contrast facility dots.")
    p.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return p.parse_args()


if __name__ == "__main__":
    run(parse_args())
