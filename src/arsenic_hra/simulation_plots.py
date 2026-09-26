"""Phase 7 figures: district HI distributions (design doc section 14.2).

HI spans many orders of magnitude, so each panel shows the share of iterations in
log-spaced bins on a logarithmic x-axis (a density per unit HI would be dominated by
the near-zero draws). The figure is a visual aid only;
every number reported comes from the iteration table, never from the plot.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e6e5e1"
SERIES = "#2a78d6"     # the single simulated series
P95_COLOR = "#eb6834"  # reference marker for the 95th percentile


HI_FLOOR = 1e-3  # HI below this is collected into the first bin


def _panel(ax, hi: np.ndarray, title: str, edges: np.ndarray) -> None:
    ax.set_facecolor(SURFACE)
    shown = np.clip(hi, edges[0], edges[-1])
    ax.hist(shown, bins=edges, weights=np.full(hi.size, 100.0 / hi.size), color=SERIES,
            edgecolor=SURFACE, linewidth=0.8)
    p95 = float(np.quantile(hi, 0.95))
    p_gt_1 = float(np.mean(hi > 1.0))
    below = float(np.mean(hi <= HI_FLOOR))
    ymax = ax.get_ylim()[1] * 1.18
    ax.set_ylim(0, ymax)
    box = {"facecolor": SURFACE, "edgecolor": "none", "pad": 1.5}
    ax.axvline(1.0, color=TEXT_SECONDARY, linestyle=(0, (4, 3)), linewidth=1.5)
    ax.axvline(p95, color=P95_COLOR, linewidth=2)
    ax.text(1.0, ymax * 0.97, "HI = 1 ", color=TEXT_SECONDARY, fontsize=9, va="top", ha="right", bbox=box)
    ax.text(p95, ymax * 0.97, f" P95 = {p95:.3g}", color=TEXT_PRIMARY, fontsize=9, va="top", ha="left", bbox=box)
    ax.set_xscale("log")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_title(title, loc="left", fontsize=10.5, color=TEXT_PRIMARY)
    ax.text(0.045, 0.97, f"P(HI > 1) = {100 * p_gt_1:.1f}%\n{100 * below:.1f}% of iterations at HI <= {HI_FLOOR:g}"
            f"\nN = {hi.size:,}", transform=ax.transAxes, ha="left", va="top", fontsize=8.5,
            color=TEXT_SECONDARY, bbox=box)
    ax.set_ylabel("share of iterations (%)", color=TEXT_SECONDARY, fontsize=9)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8.5)


def plot_district_hi(district: str, hi_adult: np.ndarray, hi_child: np.ndarray, out_dir: Path) -> Path:
    top = max(hi_adult.max(), hi_child.max())
    edges = np.logspace(np.log10(HI_FLOOR), np.log10(top), 50)
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True, facecolor=SURFACE)
    _panel(axes[0], hi_adult, f"(a) Adults - {district}", edges)
    _panel(axes[1], hi_child, f"(b) Children - {district}", edges)
    axes[1].set_xlabel("Hazard index, HI (dimensionless, log scale)", color=TEXT_SECONDARY, fontsize=9)
    fig.suptitle(f"Simulated HI distribution, {district} (primary run)", x=0.02, ha="left",
                 fontsize=12, color=TEXT_PRIMARY)
    fig.text(0.02, 0.005, f"First bin also holds every iteration with HI <= {HI_FLOOR:g} (arsenic below the "
             "detection limit, where the fitted distribution is extrapolated).", fontsize=7.5, color=TEXT_SECONDARY)
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    path = out_dir / f"simulation_hi_distribution_{district}.png"
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return path
