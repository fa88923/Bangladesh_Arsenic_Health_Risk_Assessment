"""Diagnostic figures for the weighted body-weight fits."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .fitting import FAMILIES, PARAM_NAMES, frozen, weighted_ecdf_midpoints

# Categorical slots 1-4 of the reference palette, in fixed order; empirical data in neutral ink.
FAMILY_COLORS = {"Normal": "#2a78d6", "Lognormal": "#eb6834", "Gamma": "#1baf7a", "Triangular": "#eda100"}
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
EMPIRICAL_FILL = "#d6d5cf"
GRID = "#e6e5e0"

plt.rcParams.update({
    "axes.edgecolor": INK_SECONDARY,
    "axes.labelcolor": INK,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": INK_SECONDARY,
    "ytick.color": INK_SECONDARY,
    "legend.frameon": False,
    "legend.fontsize": 9,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
})


def _fits(fit_table: pd.DataFrame) -> dict:
    out = {}
    for row in fit_table.itertuples(index=False):
        row = row._asdict()
        params = {name: row[f"ls_{name}"] for name in PARAM_NAMES[row["distribution"]]}
        out[row["distribution"]] = (frozen(row["distribution"], params), bool(row["admissible"]))
    return out


def _label(family: str, admissible: bool, selected: str) -> str:
    if family == selected:
        return f"{family} (selected)"
    return family if admissible else f"{family} (inadmissible)"


def _style(family: str, admissible: bool, selected: str) -> dict:
    return {
        "color": FAMILY_COLORS[family],
        "linewidth": 2.2 if family == selected else 1.6,
        "linestyle": "-" if admissible else (0, (4, 2)),
    }


def plot_population(x, w, fit_table: pd.DataFrame, selected: str, population: str, out_dir: Path) -> list[Path]:
    """Write histogram/PDF, ECDF/CDF, per-family Q-Q, and tail figures; return their paths."""
    out_dir.mkdir(parents=True, exist_ok=True)
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    ecdf = weighted_ecdf_midpoints(x, w)
    fits = _fits(fit_table)
    title_pop = "Adults 18-69 (STEPS 2018)" if population == "adult" else "Children 0-59 months (MICS7)"
    grid = np.linspace(max(x.min() * 0.8, 1e-6), x.max() * 1.05, 800)
    paths = []

    # 1. Weighted histogram + fitted PDFs
    fig, ax = plt.subplots(figsize=(7, 4.3))
    width = 1.0 if population == "adult" else 0.5
    bins = np.arange(np.floor(x.min()) - width / 2, x.max() + width, width)
    ax.hist(x, bins=bins, weights=w, density=True, color=EMPIRICAL_FILL, edgecolor="white",
            linewidth=0.4, label="Survey-weighted data")
    for fam in FAMILIES:
        dist, ok = fits[fam]
        ax.plot(grid, dist.pdf(grid), label=_label(fam, ok, selected), **_style(fam, ok, selected))
    ax.set(xlabel="Body weight (kg)", ylabel="Density", title=f"{title_pop}: weighted histogram and LS-fitted PDFs")
    ax.set_xlim(np.quantile(x, 0.0005) * 0.9, np.quantile(x, 0.9995) * 1.1)
    ax.legend()
    paths.append(out_dir / "weighted_histogram_pdf.png")
    fig.savefig(paths[-1])
    plt.close(fig)

    # 2. Weighted ECDF + fitted CDFs
    fig, ax = plt.subplots(figsize=(7, 4.3))
    ax.step(ecdf.x, ecdf.F, where="post", color=INK_SECONDARY, linewidth=1.4, label="Weighted ECDF")
    for fam in FAMILIES:
        dist, ok = fits[fam]
        ax.plot(grid, dist.cdf(grid), label=_label(fam, ok, selected), **_style(fam, ok, selected))
    ax.set(xlabel="Body weight (kg)", ylabel="Cumulative probability", ylim=(0, 1.01),
           title=f"{title_pop}: weighted ECDF and LS-fitted CDFs")
    ax.set_xlim(np.quantile(x, 0.0005) * 0.9, np.quantile(x, 0.9995) * 1.1)
    ax.legend(loc="lower right")
    paths.append(out_dir / "weighted_ecdf_cdf.png")
    fig.savefig(paths[-1])
    plt.close(fig)

    # 3. Weighted Q-Q per family
    for fam in FAMILIES:
        dist, ok = fits[fam]
        theo = dist.ppf(ecdf.p)
        fig, ax = plt.subplots(figsize=(5, 5))
        lim = (min(theo.min(), ecdf.x.min()) * 0.95, max(theo.max(), ecdf.x.max()) * 1.02)
        ax.plot(lim, lim, color=INK_SECONDARY, linewidth=1, linestyle=(0, (3, 3)), label="y = x")
        ax.scatter(theo, ecdf.x, s=np.clip(ecdf.q * 4e3, 8, 60), color=FAMILY_COLORS[fam],
                   edgecolor="white", linewidth=0.6, label="distinct observed weights (size = survey mass)")
        ax.set(xlim=lim, ylim=lim, xlabel=f"{fam} theoretical quantile (kg)", ylabel="Weighted empirical quantile (kg)",
               title=f"{title_pop}\nWeighted Q-Q: {_label(fam, ok, selected)}")
        ax.legend(loc="upper left")
        paths.append(out_dir / f"weighted_qq_{fam.lower()}.png")
        fig.savefig(paths[-1])
        plt.close(fig)

    # 4. Tail zooms: lower-tail CDF and upper-tail survival, log scale
    fig, (lo_ax, hi_ax) = plt.subplots(1, 2, figsize=(10, 4.2))
    lo_mask = ecdf.p <= 0.10
    hi_mask = ecdf.p >= 0.90
    lo_grid = np.linspace(ecdf.x[lo_mask].min() * 0.9, ecdf.x[lo_mask].max(), 400)
    hi_grid = np.linspace(ecdf.x[hi_mask].min(), ecdf.x.max() * 1.05, 400)
    lo_ax.scatter(ecdf.x[lo_mask], ecdf.p[lo_mask], s=10, color=INK_SECONDARY, label="Weighted empirical", zorder=3)
    hi_ax.scatter(ecdf.x[hi_mask], 1 - ecdf.p[hi_mask], s=10, color=INK_SECONDARY, label="Weighted empirical", zorder=3)
    for fam in FAMILIES:
        dist, ok = fits[fam]
        style = _style(fam, ok, selected)
        lo_ax.plot(lo_grid, np.clip(dist.cdf(lo_grid), 1e-6, None), label=_label(fam, ok, selected), **style)
        hi_ax.plot(hi_grid, np.clip(dist.sf(hi_grid), 1e-6, None), **style)
    lo_ax.set(yscale="log", ylim=(1e-4, 0.15), xlabel="Body weight (kg)", ylabel="F(x)", title="Lower tail (F <= 0.10)")
    hi_ax.set(yscale="log", ylim=(1e-5, 0.15), xlabel="Body weight (kg)", ylabel="1 - F(x)", title="Upper tail (F >= 0.90)")
    lo_ax.legend(loc="lower right")
    fig.suptitle(f"{title_pop}: tail agreement of LS-fitted candidates", fontsize=11)
    paths.append(out_dir / "weighted_tails.png")
    fig.savefig(paths[-1])
    plt.close(fig)
    return paths
