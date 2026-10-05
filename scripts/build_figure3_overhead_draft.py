"""Draft Figure 3 variant: statewide change factor with cooling overhead (PUE - 1) in place of PUE.

PUE - 1 removes the fixed IT load from the denominator, so its change factor
reflects the change in cooling and other overhead energy. Everything else
matches ``figure3_statewide_change_factor``: Case 1, exceedance days annual,
PUE/WUE May-October mean, per-model ratio (each model divided by its own
historical value), 5th-95th percentile band across the 20 models.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_results_delta_figures import (  # noqa: E402
    INK_MUTED,
    INK_PRIMARY,
    METRIC_COLORS,
    PERIOD_LABELS,
    PERIOD_ORDER,
    SERIES_BLUE,
    _apply_style,
    _change_factor,
    _statewide_heat_change_factor,
    _warm_season_pue_wue_by_gcm,
)


def build(root: Path, out: Path) -> None:
    warm = _warm_season_pue_wue_by_gcm(root)
    case1 = warm[warm["archetype"] == "ae-chiller"]
    pue = case1.pivot(index="gcm", columns="period", values="pue_mean") - 1.0
    wue = case1.pivot(index="gcm", columns="period", values="wue_mean")
    series = [
        ("Threshold-exceedance days", _statewide_heat_change_factor(root), METRIC_COLORS["heat"]),
        ("Cooling energy overhead (PUE − 1)", _change_factor(pue), SERIES_BLUE),
        ("WUE", _change_factor(wue), METRIC_COLORS["wue"]),
    ]

    _apply_style()
    fig, ax = plt.subplots(figsize=(7.5, 6))
    x = range(len(PERIOD_ORDER))
    for label, df, color in series:
        ax.plot(x, df["cf_mean"], color=color, linewidth=2.2, marker="o", markersize=5, label=label, zorder=2)
        ax.fill_between(x, df["cf_p05"], df["cf_p95"], color=color, alpha=0.15, zorder=1, linewidth=0)
    ax.axhline(1.0, color=INK_MUTED, linewidth=0.8, linestyle=":", zorder=0)
    ax.set_xticks(list(x))
    ax.set_xticklabels(PERIOD_LABELS, fontsize=14)
    ax.set_ylabel("Change factor rel. to 1985-2014 baseline", fontsize=14)
    ax.tick_params(axis="both", colors=INK_PRIMARY, labelsize=14)
    ax.yaxis.label.set_color(INK_PRIMARY)
    ax.legend(frameon=False, fontsize=14, loc="upper left")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", default=".")
    args = p.parse_args()
    root = Path(args.repo_root).resolve()
    build(root, root / "docs/results-figures/figure3_statewide_change_factor_overhead_draft.png")
