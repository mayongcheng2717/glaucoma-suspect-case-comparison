"""Journal-grade matplotlib styling for the manuscript figures."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

# Nature-style typographic palette
PAL = {
    "blue": "#0F4D92",
    "red": "#B2182B",
    "orange": "#E08214",
    "teal": "#1B7837",
    "purple": "#6A3D9A",
    "grey": "#8C8C8C",
    "dark": "#252525",
    "sand": "#D9C9A3",
    "cyan": "#2C7FB8",
    "gold": "#C8A200",
}


def style(width_mm: float = 180, height_mm: float = 120, dpi: int = 400):
    """Return (fig, ax) with a single-column/双栏 sized canvas."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.5,
        "axes.titlesize": 8.5,
        "axes.labelsize": 7.5,
        "xtick.labelsize": 6.8,
        "ytick.labelsize": 6.8,
        "legend.fontsize": 6.8,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "lines.linewidth": 0.9,
        "patch.linewidth": 0.5,
        "savefig.dpi": dpi,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "legend.frameon": False,
    })
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4))
    return fig, ax


def panel_label(ax, letter: str, x: float = -0.13, y: float = 1.06):
    ax.text(x, y, letter.upper(), transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", va="top", ha="left")


def panel(ax, letter: str, x: float = -0.12, y: float = 1.08):
    """Journal-style panel letter: uppercase, bold, outside the axes, no title."""
    ax.text(x, y, letter.upper(), transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", va="bottom", ha="left")


def strip_titles(axes):
    """Remove in-panel titles: journal figures carry the description in the caption."""
    for ax in axes:
        ax.set_title("")


def save(fig, stem: str):
    """Save as PDF (vector) + 600-dpi PNG."""
    fig.savefig(stem + ".pdf")
    fig.savefig(stem + ".png", dpi=600)
    plt.close(fig)
