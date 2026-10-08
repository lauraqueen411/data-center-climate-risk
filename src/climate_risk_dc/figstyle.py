"""Shared matplotlib style for the heat-chapter figures.

Sizes are in points at the width the figure is inserted in the Word
document (``FIG_WIDTH_IN`` = 6.5 in), so a figure built at that width and
saved without ``bbox_inches="tight"`` prints at these sizes exactly.

Used by ``scripts/fig_show_everything.py`` and ``scripts/fig3_maps.py``. The other chapter
figures (``build_results_delta_figures.py``, ``plot_archetype_response_curves.py``)
still set their own style; see "Assumptions log".

Assumptions log
----------------
- **Font sizes.** The other chapter figures are built ~13 in wide with
  14 pt ticks/labels and 16 pt panel letters, then inserted at 6.5 in, so
  they print at ~7 pt (ticks, axis labels) and ~8 pt (panel letters). This
  module sets a 9 pt floor for ticks and 10.5 pt for titles and axis
  labels, larger than those figures by design.
- **Font family.** The other figures set no family and render in
  matplotlib's default DejaVu Sans. This module uses Arial, then Helvetica,
  falling back to DejaVu Sans where neither is installed.
- **Ink.** All text, spines and ticks are black (the other figures use
  grey ink, grey gridlines and open top/right spines).
- **Map projection and extent** (``oregon_map_axes``) are copied from the
  orientation map in ``build_results_delta_figures.py`` (``OREGON_ALBERS``,
  extent 125-116 W, 42-46.5 N, Natural Earth 50m states, borders and
  coastline); only the ink is changed to black. cartopy is imported inside
  the function so the non-map figures do not depend on it.
"""

from __future__ import annotations

import matplotlib.pyplot as plt

FIG_WIDTH_IN = 6.5
TICK_PT = 9.0
LABEL_PT = 10.5
TITLE_PT = 10.5
LEGEND_PT = 10.0
INK = "black"
SPINE_WIDTH = 0.8
ZERO_LINE_WIDTH = 0.6
ZERO_LINE_COLOR = "#555555"
FONT_FAMILY = ["Arial", "Helvetica", "DejaVu Sans"]

# Okabe-Ito
BLUE = "#0072B2"
VERMILLION = "#D55E00"
BLUISH_GREEN = "#009E73"
REDDISH_PURPLE = "#CC79A7"
ORANGE = "#E69F00"
SKY_BLUE = "#56B4E9"

# Ecoregion colors (Okabe-Ito). CP and WV match Figure 4. Yellow is left out
# because open (Case 5) markers draw the color as a thin edge on white.
ECOREGION_COLOR = {"ECSF": ORANGE, "BM": SKY_BLUE, "CP": BLUISH_GREEN, "KM": INK, "WV": REDDISH_PURPLE}


def apply() -> None:
    """Set the shared rcParams. Call once before building a figure."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": FONT_FAMILY,
        "font.size": TICK_PT,
        "axes.titlesize": TITLE_PT,
        "axes.labelsize": LABEL_PT,
        "xtick.labelsize": TICK_PT,
        "ytick.labelsize": TICK_PT,
        "legend.fontsize": LEGEND_PT,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "axes.edgecolor": INK,
        "axes.linewidth": SPINE_WIDTH,
        "axes.grid": False,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": SPINE_WIDTH,
        "ytick.major.width": SPINE_WIDTH,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.top": False,
        "ytick.right": False,
        "lines.solid_capstyle": "butt",
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "savefig.dpi": 300,
    })


def oregon_map_axes(fig, subplot_spec):
    """Map axes over Oregon in the chapter's Albers projection, with state and coast lines."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    proj = ccrs.AlbersEqualArea(central_longitude=-120.5, central_latitude=43.0, standard_parallels=(43.0, 45.5))
    ax = fig.add_subplot(subplot_spec, projection=proj)
    ax.set_extent([-125, -116, 42, 46.5], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale("50m"), facecolor="none", edgecolor=INK, linewidth=0.6, zorder=2)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), edgecolor=INK, linewidth=0.6, zorder=2)
    ax.spines["geo"].set_edgecolor(INK)
    ax.spines["geo"].set_linewidth(SPINE_WIDTH)
    return ax
