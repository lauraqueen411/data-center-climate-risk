"""Shared helpers for the October 25 figure redraws (``figs_oct25_*.py``).

Output folder and an old-vs-new side-by-side preview. No data handling.
"""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

OUT_DIR = Path("docs/results-figures/figures_oct25")
PREVIEW_HEIGHT_PX = 700


def side_by_side(old: Path, new: Path, out: Path, height: int = PREVIEW_HEIGHT_PX) -> None:
    """Old (left) and new (right), scaled to a common height, labelled 'Old' and 'New'."""
    ims = []
    for p in (old, new):
        im = Image.open(p).convert("RGB")
        ims.append(im.resize((round(im.width * height / im.height), height), Image.LANCZOS))
    pad, head = 20, 36
    canvas = Image.new("RGB", (sum(i.width for i in ims) + 3 * pad, height + head + pad), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 24)
    except OSError:
        font = ImageFont.load_default()
    x = pad
    for label, im in zip(("Old", "New"), ims):
        draw.text((x, 6), label, fill="black", font=font)
        canvas.paste(im, (x, head))
        x += im.width + pad
    canvas.save(out)


def apply_figstyle() -> None:
    """``figstyle.apply()`` plus Liberation Sans (metrically identical to Arial) ahead of the DejaVu fallback.

    The chapter figures were built where Arial is installed; on the server it
    is not, and figstyle would fall back to DejaVu Sans, which is wider. Arial
    still wins wherever it exists.
    """
    import matplotlib.pyplot as plt

    from climate_risk_dc import figstyle

    figstyle.apply()
    fams = [f for f in figstyle.FONT_FAMILY if f != "DejaVu Sans"] + ["Liberation Sans", "DejaVu Sans"]
    plt.rcParams["font.sans-serif"] = fams
