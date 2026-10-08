"""Figure 5 of the heat risk chapter: why the change differs between
facilities (docs/tasks/result-figures.md).

Mid-century only, one point per facility, marker shape by ecoregion.

- (a) x: historical annual afternoon wet-bulb (C); y: change in annual mean
  daily Tmax (C). Expected: a nearly flat cloud (the hazard is uniform).
- (b) same x; y: change in annual mean PUE, Case 1.

Marker shape and color both give the ecoregion (``figstyle.ECOREGION_COLOR``).

Inputs (nothing re-run):
- ``outputs/ensemble_pue_wue_facilities/analysis/facility_annual_tmax.csv``:
  annual mean daily Tmax per facility, model and year (historical and
  mid-century), reduced on the server from ``tasmax_c`` in the Case 1 daily
  PUE/WUE files.
- ``.../analysis/fig_show_everything/model_values_cache.csv``: historical
  annual afternoon wet-bulb and period-mean annual PUE per facility and model.
- ``.../analysis/fig_show_everything/facility_grid_cells.csv``: ecoregion and
  LOCA2 grid cell.

Outputs:
- ``docs/results-figures/fig5_mechanism.{png,pdf}`` and ``fig5_mechanism_caption.txt``.
- ``outputs/ensemble_pue_wue_facilities/analysis/fig5_mechanism/``:
  ``facility_values.csv``, ``fits.csv``, ``run_manifest.csv``.

Assumptions log
----------------
- **Values.** Per facility and model: the mean over the window's years of
  annual means; change is mid-century minus historical; then the mean
  across the 20 models. Wet-bulb is the historical value, also averaged
  across models.
- **No fit lines** (Laura, 2026-10-07): the task asked for OLS lines with
  slope and R^2, but the 109 facilities fall in 22 LOCA2 grid cells (and
  about five ecoregion clusters), and facilities sharing a cell have
  identical values. OLS over the facilities and over the 22 cells is still
  written to ``fits.csv`` and the printed summary for reference only.
- **Ecoregion colors** are redundant with marker shape so the figure still
  reads in greyscale.
- **Case 1 only** (Laura, 2026-10-07): the task also asked for Case 5 as
  open markers in (b); it was dropped because Figure 4 already compares the
  cases.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

from climate_risk_dc import figstyle

N_MODELS = 20
PERIODS = {"historical": (1985, 2014), "midcentury": (2045, 2074)}
CASES = {"ae-chiller": "Case 1"}
ECOREGION_ABBR = {
    "Eastern Cascades Slopes and Foothills": "ECSF",
    "Blue Mountains": "BM",
    "Columbia Plateau": "CP",
    "Klamath Mountains/California High North Coast Range": "KM",
    "Willamette Valley": "WV",
}
MARKERS = {"ECSF": "D", "BM": "^", "CP": "o", "KM": "v", "WV": "s"}
MARKER_SIZE = 22
MARKER_EDGE = 0.6
FIG_HEIGHT_IN = 3.6
BOUNDS = {"tasmax_c": (-30.0, 50.0), "pue": (1.0, 3.0), "wetbulb_annual": (-20.0, 35.0)}


def _analysis_dir(root: Path) -> Path:
    return root / "outputs/ensemble_pue_wue_facilities/analysis"


def _check(name: str, v: pd.Series) -> None:
    lo, hi = BOUNDS[name]
    if v.isna().any() or (v < lo).any() or (v > hi).any():
        raise ValueError(f"{name}: values outside [{lo}, {hi}] or NaN (range {v.min()}-{v.max()})")


def _facility_values(root: Path) -> pd.DataFrame:
    """One row per facility: ecoregion, grid cell, wet-bulb, Tmax change and PUE change per case (ensemble means)."""
    a = _analysis_dir(root)
    tmax = pd.read_csv(a / "facility_annual_tmax.csv")
    _check("tasmax_c", tmax["tasmax_c"])
    for period, (start, end) in PERIODS.items():
        years = tmax.loc[tmax["period"] == period, "year"]
        if years.min() != start or years.max() != end or years.nunique() != end - start + 1:
            raise ValueError(f"facility_annual_tmax.csv {period}: years {years.min()}-{years.max()}, expected {start}-{end}")
    t = tmax.groupby(["facility_id", "gcm", "period"])["tasmax_c"].mean().unstack("period")
    assert t.groupby(level="facility_id").size().eq(N_MODELS).all(), "Missing models in Tmax"
    out = (t["midcentury"] - t["historical"]).groupby(level="facility_id").mean().rename("tmax_change_c").to_frame()

    cache = pd.read_csv(a / "fig_show_everything" / "model_values_cache.csv")
    wb = cache[(cache["metric"] == "wetbulb_annual") & (cache["period"] == "historical")]
    _check("wetbulb_annual", wb["value"])
    out["wetbulb_hist_c"] = wb.groupby("facility_id")["value"].mean()

    pue = cache[(cache["metric"] == "pue") & cache["archetype"].isin(CASES) & cache["period"].isin(PERIODS)]
    _check("pue", pue["value"])
    w = pue.pivot_table(index=["facility_id", "gcm", "archetype"], columns="period", values="value")
    ch = (w["midcentury"] - w["historical"]).groupby(level=["facility_id", "archetype"]).mean().unstack("archetype")
    for arch, case in CASES.items():
        out[f"pue_change_{case.replace(' ', '').lower()}"] = ch[arch]

    fac = pd.read_csv(a / "fig_show_everything" / "facility_grid_cells.csv",
                      usecols=["facility_id", "facility_name", "ecoregion", "grid_cell"])
    out = fac.merge(out.reset_index(), on="facility_id", validate="one_to_one")
    out["abbr"] = out["ecoregion"].map(ECOREGION_ABBR)
    assert out.drop(columns="facility_name").notna().all().all() and len(out) == 109  # 21 atlas names are blank
    return out


def _fits(fac: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for y in ["tmax_change_c", "pue_change_case1"]:
        for unit, d in (("facility", fac), ("grid_cell", fac.drop_duplicates("grid_cell"))):
            r = stats.linregress(d["wetbulb_hist_c"], d[y])
            rows.append({"y": y, "unit": unit, "n": len(d), "slope_per_c": r.slope, "intercept": r.intercept,
                         "r2": r.rvalue ** 2, "p_value": r.pvalue, "slope_stderr": r.stderr})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------


def _scatter(ax, fac: pd.DataFrame, y: str) -> None:
    for abbr, marker in MARKERS.items():
        d = fac[fac["abbr"] == abbr]
        ax.scatter(d["wetbulb_hist_c"], d[y], marker=marker, s=MARKER_SIZE, linewidths=MARKER_EDGE,
                   facecolors=figstyle.ECOREGION_COLOR[abbr], edgecolors=figstyle.INK, zorder=3)


def plot(fac: pd.DataFrame, path: Path) -> None:
    figstyle.apply()
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(figstyle.FIG_WIDTH_IN, FIG_HEIGHT_IN), sharex=True,
                                     gridspec_kw={"left": 0.1, "right": 0.985, "top": 0.92, "bottom": 0.31, "wspace": 0.36})
    _scatter(ax_a, fac, "tmax_change_c")
    ax_a.set_ylabel("Change in annual mean\ndaily Tmax (°C)")
    ax_a.set_title("(a) Warming", loc="left")
    # A zero-anchored y range so a flat cloud reads as flat, not as noise stretched to fill the panel.
    ax_a.set_ylim(0, fac["tmax_change_c"].max() * 1.15)

    _scatter(ax_b, fac, "pue_change_case1")
    ax_b.set_ylabel("Change in annual mean PUE")
    ax_b.set_title("(b) PUE response, Case 1", loc="left")
    ax_b.set_ylim(0, None)
    for ax in (ax_a, ax_b):
        ax.set_xlabel("Historical annual afternoon\nwet-bulb temperature (°C)")

    eco = [Line2D([], [], linestyle="none", marker=m, markersize=5.5, markerfacecolor=figstyle.ECOREGION_COLOR[k],
                  markeredgecolor=figstyle.INK, markeredgewidth=MARKER_EDGE) for k, m in MARKERS.items()]
    fig.legend(eco, list(MARKERS), loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=5, handlelength=1.0, handletextpad=0.4, columnspacing=1.2, borderaxespad=0.3)
    for ext in ("png", "pdf"):
        fig.savefig(path.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)


def caption(fac: pd.DataFrame) -> str:
    n = fac["abbr"].value_counts()
    return (
        "Mid-century (2045–2074 minus 1985–2014) change in (a) annual mean daily maximum temperature and (b) annual "
        "mean power usage effectiveness (PUE) for cooling-system Case 1, against historical annual mean afternoon wet-bulb temperature. Each marker is one "
        "facility (mean of 20 climate models, SSP3-7.0); marker shape and color give the ecoregion: Eastern Cascades Slopes and "
        f"Foothills (ECSF; n = {n['ECSF']}), Blue Mountains (BM; n = {n['BM']}), Columbia Plateau (CP; n = {n['CP']}), "
        f"Klamath Mountains (KM; n = {n['KM']}) and Willamette Valley (WV; n = {n['WV']}). Facilities in the same LOCA2 grid cell share identical values."
    )


def print_summary(fac: pd.DataFrame, fits: pd.DataFrame) -> None:
    t = fac["tmax_change_c"]
    print("\n=== Tmax change across facilities (mid-century, ensemble mean) ===")
    print(f"   range {t.min():.2f} – {t.max():.2f} °C (spread {t.max() - t.min():.2f}), mean {t.mean():.2f}")
    print(fac.groupby("abbr")["tmax_change_c"].agg(["mean", "min", "max", "size"]).round(2).to_string())
    print("\n=== OLS against historical annual wet-bulb (reference only; not plotted) ===")
    for _, r in fits.iterrows():
        print(f"   {r['y']:17s} {r['unit']:9s} n={r['n']:3d}  slope {r['slope_per_c']:+.5f}/°C "
              f"(SE {r['slope_stderr']:.5f})  R² {r['r2']:.3f}  p {r['p_value']:.2g}")


def run(args: argparse.Namespace) -> None:
    root = Path(args.repo_root).resolve()
    out_dir = _analysis_dir(root) / "fig5_mechanism"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = root / "docs" / "results-figures"

    fac = _facility_values(root)
    fits = _fits(fac)
    fac.sort_values(["ecoregion", "facility_id"]).to_csv(out_dir / "facility_values.csv", index=False)
    fits.to_csv(out_dir / "fits.csv", index=False)
    pd.DataFrame([{
        "n_models": N_MODELS, "scenario": "ssp370",
        "periods": "; ".join(f"{k}={a}-{b}" for k, (a, b) in PERIODS.items()),
        "cases": "; ".join(f"{v}={k}" for k, v in CASES.items()), "wetbulb_window": "annual",
        "plotted_fit": "none",
    }]).to_csv(out_dir / "run_manifest.csv", index=False)

    plot(fac, fig_dir / "fig5_mechanism.png")
    (fig_dir / "fig5_mechanism_caption.txt").write_text(caption(fac) + "\n")
    print_summary(fac, fits)
    print(f"\nWrote {fig_dir / 'fig5_mechanism.png'} (+ .pdf) and CSVs under {out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Figure 5 (mechanism) of the heat chapter.")
    parser.add_argument("--repo-root", default=".", help="Repository root path (default: current directory)")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
