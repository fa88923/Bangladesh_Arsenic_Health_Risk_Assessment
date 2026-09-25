"""Phase 8 figures: statistic vs N (across-seed mean with seed range), and MC error decay."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from arsenic_hra.simulation_plots import GRID, SURFACE, TEXT_PRIMARY, TEXT_SECONDARY  # noqa: E402

POP_COLORS = {"adult": "#2a78d6", "child": "#eb6834"}  # categorical slots 1 and 2, fixed order
LABELS = {"mean_HI": "Mean HI", "P50_HI": "Median (P50) HI", "P95_HI": "P95 HI",
          "P_HI_gt_1": "P(HI > 1)", "P_HI_gt_2": "P(HI > 2)", "P95_ELCR": "P95 ELCR (lifetime AT)"}


def _style(ax) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=7.5)


def plot_metric_vs_N(agg: pd.DataFrame, metric: str, districts: list[str], out_dir: Path) -> Path:
    fig, axes = plt.subplots(2, 5, figsize=(13, 5.6), facecolor=SURFACE)
    for ax, d in zip(axes.flat, districts):
        _style(ax)
        for pop, color in POP_COLORS.items():
            g = agg[(agg.district == d) & (agg.population == pop) & (agg.metric == metric)].sort_values("N")
            ax.fill_between(g.N, g.seed_min, g.seed_max, color=color, alpha=0.18, linewidth=0)
            ax.plot(g.N, g.across_seed_mean, color=color, linewidth=2, marker="o", markersize=4.5,
                    markeredgecolor=SURFACE, markeredgewidth=1, label=pop)
        ax.set_xscale("log")
        ax.set_xticks([1000, 5000, 10000, 20000], ["1k", "5k", "10k", "20k"])
        ax.minorticks_off()
        ax.axvline(10000, color=TEXT_SECONDARY, linewidth=0.8, linestyle=(0, (3, 3)))
        ax.set_title(d, loc="left", fontsize=9.5, color=TEXT_PRIMARY)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=2, frameon=False, fontsize=9)
    fig.suptitle(f"{LABELS[metric]} vs simulation size N: across-seed mean (line) and range of 5 seeds (band)",
                 x=0.01, ha="left", fontsize=11, color=TEXT_PRIMARY)
    fig.text(0.01, 0.01, "Dashed line: primary N = 10,000. Each panel has its own y-scale.",
             fontsize=7.5, color=TEXT_SECONDARY)
    fig.tight_layout(rect=(0, 0.03, 1, 0.94))
    path = out_dir / f"convergence_{metric}_vs_N.png"
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return path


def plot_error_decay(agg: pd.DataFrame, out_dir: Path) -> Path:
    """Between-seed spread vs N on log-log axes; Monte Carlo error should fall as N^-1/2."""
    metrics = list(LABELS)
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.4), facecolor=SURFACE)
    Ns = np.array([1000, 5000, 10000, 20000])
    for ax, metric in zip(axes.flat, metrics):
        _style(ax)
        g = agg[agg.metric == metric]
        ycol = "between_seed_sd" if metric.startswith("P_HI") else "between_seed_cv_pct"
        for (d, pop), gg in g.groupby(["district", "population"]):
            gg = gg.sort_values("N")
            ax.plot(gg.N, gg[ycol], color=POP_COLORS[pop], alpha=0.35, linewidth=1)
        med = g.groupby("N")[ycol].median()
        ref = med.loc[1000] * np.sqrt(1000 / Ns)
        ax.plot(Ns, ref, color=TEXT_PRIMARY, linewidth=1.5, linestyle=(0, (4, 3)), label="N^-1/2 from median at 1k")
        tol = 5.0 if ycol == "between_seed_cv_pct" else 0.01
        ax.axhline(tol, color=TEXT_SECONDARY, linewidth=1, linestyle=(0, (1, 2)))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(Ns, ["1k", "5k", "10k", "20k"])
        ax.minorticks_off()
        ax.set_title(LABELS[metric], loc="left", fontsize=9.5, color=TEXT_PRIMARY)
        ax.set_ylabel("between-seed SD" if ycol == "between_seed_sd" else "between-seed CV (%)",
                      fontsize=8, color=TEXT_SECONDARY)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=POP_COLORS["adult"], label="adult (one line per district)"),
               Line2D([], [], color=POP_COLORS["child"], label="child (one line per district)"),
               Line2D([], [], color=TEXT_PRIMARY, linestyle=(0, (4, 3)), label="N^-1/2 reference"),
               Line2D([], [], color=TEXT_SECONDARY, linestyle=(0, (1, 2)), label="tolerance (5% CV or 0.01)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8.5)
    fig.suptitle("Monte Carlo error vs N (5 independent seeds per N)", x=0.01, ha="left", fontsize=11,
                 color=TEXT_PRIMARY)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    path = out_dir / "convergence_mc_error_decay.png"
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return path
