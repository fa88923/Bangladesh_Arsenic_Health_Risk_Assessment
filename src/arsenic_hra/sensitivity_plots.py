"""Phase 9 figures: tornado plots (design doc 16.3) and parameter-uncertainty intervals."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from arsenic_hra.convergence_plots import POP_COLORS  # noqa: E402
from arsenic_hra.simulation_plots import GRID, SURFACE, TEXT_PRIMARY, TEXT_SECONDARY  # noqa: E402

# Diverging pair for signed correlations (increases HI / decreases HI).
POSITIVE = "#e34948"
NEGATIVE = "#2a78d6"
INPUT_LABELS = {"C_mg_per_L": "Arsenic C", "IR_L_per_day": "Ingestion rate IR", "BW_kg": "Body weight BW",
                "EF_days_per_year": "Exposure frequency EF", "ET_h_per_day": "Exposure time ET",
                "SA_m2": "Skin area SA"}


def _style(ax) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8.5, length=0)


def plot_tornado(table: pd.DataFrame, district: str, out_dir: Path) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), facecolor=SURFACE, sharex=True)
    for ax, pop, tag in zip(axes, ("adult", "child"), ("(a) Adults", "(b) Children")):
        _style(ax)
        g = table[(table.district == district) & (table.population == pop) & (table.output == "HI")]
        g = g.sort_values("abs_rho")
        y = np.arange(len(g))
        colors = [POSITIVE if r > 0 else NEGATIVE for r in g.spearman_rho]
        ax.barh(y, g.spearman_rho, color=colors, height=0.62, edgecolor=SURFACE, linewidth=1)
        ax.set_yticks(y, [INPUT_LABELS[i] for i in g.input], color=TEXT_PRIMARY, fontsize=8.5)
        for yi, r in zip(y, g.spearman_rho):
            ax.text(r + (0.02 if r >= 0 else -0.02), yi, f"{r:+.2f}", va="center",
                    ha="left" if r >= 0 else "right", fontsize=8, color=TEXT_PRIMARY)
        ax.axvline(0, color=TEXT_SECONDARY, linewidth=0.8)
        ax.set_xlim(-1.05, 1.05)
        ax.set_title(f"{tag} - {district}", loc="left", fontsize=10, color=TEXT_PRIMARY)
        ax.set_xlabel("Spearman rank correlation with HI", fontsize=8.5, color=TEXT_SECONDARY)
    fig.suptitle(f"Sensitivity of HI to sampled inputs, {district} (primary run, N = 10,000)", x=0.01,
                 ha="left", fontsize=11, color=TEXT_PRIMARY)
    fig.text(0.01, 0.01, "Red: HI tends to rise with the input.  Blue: HI tends to fall.  "
             "Fixed ED, Kp, CF, RfD, CSF are excluded (zero variance).", fontsize=7.5, color=TEXT_SECONDARY)
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    path = out_dir / f"sensitivity_tornado_{district}.png"
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return path


def plot_parameter_uncertainty(summary: pd.DataFrame, statistic: str, label: str, districts: list[str],
                               out_dir: Path, log: bool = True) -> Path:
    fig, ax = plt.subplots(figsize=(8.5, 5.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8.5, length=0)
    y0 = np.arange(len(districts))[::-1]
    for offset, pop in ((0.17, "adult"), (-0.17, "child")):
        g = summary[(summary.statistic == statistic) & (summary.population == pop)].set_index("district")
        g = g.loc[districts]
        y = y0 + offset
        ax.hlines(y, g.param_unc_q025, g.param_unc_q975, color=POP_COLORS[pop], linewidth=2.5, alpha=0.55)
        ax.plot(g.primary_value, y, "o", color=POP_COLORS[pop], markersize=7, markeredgecolor=SURFACE,
                markeredgewidth=1.2, label=f"{pop}: primary estimate, 95% parameter-uncertainty interval")
    ax.set_yticks(y0, districts, color=TEXT_PRIMARY, fontsize=9)
    if log:
        ax.set_xscale("log")
    if statistic.startswith("P95_HI") or statistic == "mean_HI":
        ax.axvline(1.0, color=TEXT_SECONDARY, linewidth=1, linestyle=(0, (4, 3)))
    ax.set_xlabel(label, fontsize=9, color=TEXT_SECONDARY)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1, frameon=False, fontsize=8)
    ax.set_title(f"{label}: fitted-parameter uncertainty (bootstrap of wells and survey PSUs)", loc="left",
                 fontsize=10.5, color=TEXT_PRIMARY)
    fig.tight_layout()
    path = out_dir / f"sensitivity_parameter_uncertainty_{statistic}.png"
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return path
