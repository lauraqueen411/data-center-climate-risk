"""Draft replacement for Figure 3: acute vs. chronic change by ecoregion, one dot per climate model.

Panel (a): change in threshold-exceedance days per year (acute).
Panel (b): change in annual-mean PUE, Case 1 (chronic; annual basis to match Table 1 / Figures 1 and 4).

Per model, facilities are averaged within each multi-facility ecoregion first, so
the spread across dots is inter-model (climate) uncertainty only. Ecoregions are
ordered by historical exceedance days (descending), which makes the reversal
between panels visible. Single-facility ecoregions are excluded.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_results_delta_figures import (  # noqa: E402
    AXIS_LINE,
    ECOREGION_COLORS,
    ECOREGION_MARKERS,
    HEAT_METRIC,
    INK_MUTED,
    INK_PRIMARY,
    INK_SECONDARY,
    _apply_style,
)

ECOREGIONS = ["Columbia Plateau", "Blue Mountains", "Willamette Valley"]
HORIZONS = [("midcentury", "Mid-century"), ("endcentury", "End-of-century")]


def _heat_deltas(root: Path, site: pd.DataFrame) -> pd.DataFrame:
    fac_dir = root / "outputs/month2_heat_ensemble/facilities"
    rows = []
    for f in sorted(fac_dir.glob("facility_heat_indices_*_ssp370_historical.csv")):
        gcm = f.name.split("_")[3]
        hist = pd.read_csv(f).set_index("facility_id")[HEAT_METRIC]
        for key, _ in HORIZONS:
            fut = pd.read_csv(fac_dir / f"facility_heat_indices_{gcm}_ssp370_{key}.csv").set_index("facility_id")[HEAT_METRIC]
            for fid, v in (fut - hist).items():
                rows.append({"facility_id": fid, "gcm": gcm, "horizon": key, "value": v})
    df = pd.DataFrame(rows).merge(site, on="facility_id")
    return df.groupby(["ecoregion", "gcm", "horizon"], as_index=False)["value"].mean()


def _pue_deltas(root: Path, site: pd.DataFrame) -> pd.DataFrame:
    d = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/ensemble_delta.csv")
    d = d[d["archetype"] == "ae-chiller"]
    rows = []
    for key, _ in HORIZONS:
        sub = d[["facility_id", "gcm", f"pue_mean_delta_{key}"]].rename(columns={f"pue_mean_delta_{key}": "value"})
        sub["horizon"] = key
        rows.append(sub)
    df = pd.concat(rows).merge(site, on="facility_id")
    return df.groupby(["ecoregion", "gcm", "horizon"], as_index=False)["value"].mean()


def build(root: Path, out: Path) -> None:
    by_site = pd.read_csv(root / "outputs/ensemble_pue_wue_facilities/analysis/ensemble_summary_by_site.csv")
    site = by_site[["facility_id", "US_L3NAME"]].drop_duplicates("facility_id").rename(columns={"US_L3NAME": "ecoregion"})

    panels = [
        ("(a)", _heat_deltas(root, site), "Change in exceedance days per year"),
        ("(b)", _pue_deltas(root, site), "Change in annual-mean PUE"),
    ]

    _apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
    rng = np.random.default_rng(7)
    group_gap, within = 3.0, 1.0
    for ax, (letter, df, ylabel) in zip(axes, panels):
        centers = []
        for i, eco in enumerate(ECOREGIONS):
            color, marker = ECOREGION_COLORS[eco], ECOREGION_MARKERS[eco]
            for j, (key, label) in enumerate(HORIZONS):
                x0 = i * group_gap + j * within
                vals = df[(df["ecoregion"] == eco) & (df["horizon"] == key)]["value"].to_numpy()
                jitter = rng.uniform(-0.22, 0.22, size=len(vals))
                face = color if key == "endcentury" else "white"
                ax.scatter(x0 + jitter, vals, s=34, marker=marker, facecolor=face, edgecolor=color, linewidth=1.3, alpha=0.9, zorder=2)
                ax.hlines(vals.mean(), x0 - 0.36, x0 + 0.36, color=INK_PRIMARY, linewidth=2.2, zorder=3)
                ax.text(x0, -0.02, "Mid" if key == "midcentury" else "End", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=11, color=INK_SECONDARY)
            centers.append(i * group_gap + within / 2)
        ax.set_xticks(centers)
        ax.set_xticklabels(ECOREGIONS, fontsize=13)
        ax.tick_params(axis="x", length=0, pad=34, labelcolor=INK_PRIMARY)
        ax.tick_params(axis="y", labelsize=12, colors=INK_PRIMARY)
        ax.set_ylabel(ylabel, fontsize=13, color=INK_PRIMARY)
        ax.set_xlim(-0.9, 2 * group_gap + within + 0.9)
        ax.grid(axis="x", visible=False)
        ax.spines["bottom"].set_color(AXIS_LINE)
        ax.text(0.02, 0.97, letter, transform=ax.transAxes, ha="left", va="top", fontsize=15, color=INK_PRIMARY)
    axes[1].set_ylim(bottom=0)
    axes[0].set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", default=".")
    args = p.parse_args()
    root = Path(args.repo_root).resolve()
    build(root, root / "docs/results-figures/figure3_ecoregion_model_spread_draft.png")
