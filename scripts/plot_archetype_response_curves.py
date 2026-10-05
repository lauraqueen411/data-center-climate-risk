"""Methods-document figure: how each cooling archetype's PUE/WUE model responds to temperature.

Sweeps outdoor dry-bulb temperature through each of the three cooling
archetypes' actual PUE/WUE model (``climate_risk_dc.climate.cooling_archetypes_fast``,
the same vectorized implementation the full 109-facility ensemble run uses),
at a few fixed relative-humidity levels, and plots the resulting PUE/WUE
"response curve" per archetype. Intended as a Methods-document visual aid
showing how the underlying physics-based model translates a climate driver
into the two headline metrics -- not a new analysis result.

Assumptions log
----------------
- **Temperature is not the model's only input** -- PUE/WUE here are a
  function of both outdoor temperature and relative humidity jointly (the
  economizer/adiabatic-cooling decision depends on wet-bulb conditions, which
  mix the two), plus atmospheric pressure/elevation. A single T-only curve
  would misrepresent the model as 1-D. Instead this script sweeps
  temperature at **three fixed RH levels (20%, 50%, 80%)**, chosen to span
  Oregon's plausible summer range from dry (east-of-Cascades high desert) to
  humid (Willamette Valley/coastal), so the figure is honest about the
  model's actual two-input structure.
- **Elevation held at sea level** (``p_atm_pa=101325``, the model's own
  default) -- elevation varies per facility, not per climate driver, so
  holding it fixed keeps this a temperature/humidity response-curve figure
  rather than a third sweep dimension. Higher elevation shifts the curves
  slightly (lower pressure changes psychrometrics and fan/pump power) but
  does not change the qualitative shape discussed here.
- **Table B.1 midpoint equipment parameters** (``Case1Config()``,
  ``Case2Config()``, ``Case5Config()`` defaults) -- the same defaults used
  throughout the rest of this project's PUE/WUE pipeline, not re-tuned for
  this figure.
- **Temperature range -10 to 45 degC**, 220 points -- comfortably spans
  Oregon's historical-to-end-of-century LOCA2 tasmax/tasmin range with
  margin at both ends, so both archetypes' full economizer-to-chiller
  transition is visible, not just the observed climate's slice of it.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from climate_risk_dc.climate.cooling_archetypes import Case1Config, Case2Config, Case5Config
from climate_risk_dc.climate.cooling_archetypes_fast import (
    load_cop_water_cooled_gp,
    pue_wue_case1,
    pue_wue_case2,
    pue_wue_case5,
)

# Palette -- matches scripts/build_pue_wue_ensemble_figures.py
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
AXIS_LINE = "#c3c2b7"

RH_COLORS = {20.0: "#eb6834", 50.0: "#2a78d6", 80.0: "#1baf7a"}
RH_LEVELS = [20.0, 50.0, 80.0]

ARCHETYPES = {
    "ae-chiller": ("Case 1: AE + adiabatic", pue_wue_case1, Case1Config()),
    "we-chiller": ("Case 2: waterside economizer", pue_wue_case2, Case2Config()),
    "chiller-only": ("Case 5: chiller only", pue_wue_case5, Case5Config()),
}
ARCHETYPE_ORDER = ["ae-chiller", "we-chiller", "chiller-only"]

T_RANGE_C = np.linspace(-10.0, 45.0, 220)
P_ATM_PA = 101325.0


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


def _compute_curves(cop_gp) -> dict:
    """Returns {archetype: {rh: (pue_array, wue_array)}}."""
    out = {}
    for key, (_, func, config) in ARCHETYPES.items():
        out[key] = {}
        for rh in RH_LEVELS:
            rh_arr = np.full_like(T_RANGE_C, rh)
            pue, wue = func(T_RANGE_C, rh_arr, P_ATM_PA, cop_gp, config)
            out[key][rh] = (pue, wue)
    return out


def _plot_metric(curves: dict, metric_idx: int, ylabel: str, title: str, fig_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharex=True)
    for ax, key in zip(axes, ARCHETYPE_ORDER, strict=True):
        label, _, _ = ARCHETYPES[key]
        for rh in RH_LEVELS:
            values = curves[key][rh][metric_idx]
            ax.plot(T_RANGE_C, values, color=RH_COLORS[rh], linewidth=2, label=f"{rh:.0f}% RH")
        ax.set_title(label, fontsize=9.5)
        ax.set_xlabel("Outdoor dry-bulb temperature (°C)")
    axes[0].set_ylabel(ylabel)
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle(title, y=1.03, color=INK_PRIMARY)
    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {fig_path}")


def run(args: argparse.Namespace) -> None:
    _apply_style()
    root = Path(args.repo_root).resolve()
    cop_path = root / "data" / "lei_masanet_2022" / "COP_2.pkl"
    if not cop_path.exists():
        raise FileNotFoundError(f"{cop_path} not present (gitignored data file, see CLAUDE.md)")
    cop_gp = load_cop_water_cooled_gp(cop_path)

    out_dir = root / "outputs" / "archetype_curves"
    out_dir.mkdir(parents=True, exist_ok=True)

    curves = _compute_curves(cop_gp)

    _plot_metric(
        curves, 0, "PUE", "PUE response to outdoor temperature, by archetype and relative humidity",
        out_dir / "archetype_pue_vs_temperature.png",
    )
    _plot_metric(
        curves, 1, "WUE (L/kWh)", "WUE response to outdoor temperature, by archetype and relative humidity",
        out_dir / "archetype_wue_vs_temperature.png",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot each cooling archetype's PUE/WUE response-curve vs. temperature (methods-document figure).")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
