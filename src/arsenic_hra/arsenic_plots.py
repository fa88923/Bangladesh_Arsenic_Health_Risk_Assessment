"""Diagnostic figures for the district arsenic fits.

Four figures per district, named ``results/figures/arsenic_<subject>_<district>.png``:

    arsenic_ecdf_cdf_<district>.png        reverse-KM ECDF vs the three fitted CDFs
    arsenic_histogram_pdf_<district>.png   detected-value histogram vs fitted PDFs
    arsenic_qq_<district>.png              reverse-KM Q-Q for every candidate
    arsenic_tails_<district>.png           lower-tail CDF and upper-tail survival zoom

Styling mirrors ``bodyweight_plots.py``. The histogram deliberately shows only the
detected concentrations: censored wells are drawn as a shaded band at the censoring
bounds, never as substituted concentrations.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .arsenic_fitting import FAMILIES, PARAM_NAMES, frozen, objective_mask, reverse_km_cdf
from .paths import ensure_output_dir

FAMILY_COLORS = {"Normal": "#2a78d6", "Lognormal": "#eb6834", "Gamma": "#1baf7a"}
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
EMPIRICAL_FILL = "#d6d5cf"
GRID = "#e6e5e0"
CENSORED_BAND = "#b9b7ae"

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
        if not np.all(np.isfinite(list(params.values()))):
            continue
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


def plot_district(district: str, values_mgL, censored, fit_table: pd.DataFrame, selected: str,
                  summary_row: pd.Series | dict | None, out_dir: Path) -> list[Path]:
    """Write all four diagnostic figures for one district and return their paths."""
    ensure_output_dir(out_dir)
    x = np.asarray(values_mgL, dtype=float)
    cens = np.asarray(censored, dtype=bool)
    curve = reverse_km_cdf(x, cens)
    fits = _fits(fit_table)
    detected = x[~cens]
    limits = np.unique(x[cens]) if cens.any() else np.array([])

    pct = float(summary_row["censoring_percent"]) if summary_row is not None else 100.0 * cens.mean()
    subtitle = (f"n={x.size} ({int((~cens).sum())} detected, {int(cens.sum())} left-censored, "
                f"{pct:.2f}%); reverse-KM flat at F~={curve.censoring_fraction:.3f} below "
                f"{curve.smallest_detected:.4g} mg/L")
    hi = float(np.max(x)) * 1.05
    grid = np.linspace(0.0, hi, 1000)

    mask = objective_mask(curve)
    paths: list[Path] = []

    # 1. Reverse-KM ECDF vs fitted CDFs, log x-axis (the scale spans four decades).
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.step(curve.x, curve.F, where="post", color=INK_SECONDARY, linewidth=1.5,
            label="Reverse-KM ECDF (censoring-aware)", zorder=4)
    ax.scatter(curve.x[mask], curve.F[mask], s=16, color=INK_SECONDARY, zorder=5,
               label="objective points (decision 9.1 V2)")
    dropped = curve.x[~mask]
    if dropped.size:
        ax.scatter(dropped, curve.F[~mask], s=18, facecolor="white", edgecolor=INK_SECONDARY,
                   linewidth=0.9, zorder=5, label="plateau points not used")
    for fam in FAMILIES:
        if fam not in fits:
            continue
        dist, ok = fits[fam]
        ax.plot(grid, dist.cdf(grid), label=_label(fam, ok, selected), **_style(fam, ok, selected))
    ax.set(xscale="log", xlabel="Arsenic concentration (mg/L, log scale)", ylabel="Cumulative probability",
           ylim=(0, 1.02), title=f"{district}: reverse-KM ECDF and LS-fitted CDFs\n{subtitle}")
    ax.legend(loc="upper left")
    paths.append(out_dir / f"arsenic_ecdf_cdf_{district}.png")
    fig.savefig(paths[-1])
    plt.close(fig)

    # 2. Detected-value histogram vs fitted PDFs; censored wells as a band, never substituted.
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    bins = np.linspace(0.0, hi, 26)
    ax.hist(detected, bins=bins, density=True, color=EMPIRICAL_FILL, edgecolor="white",
            linewidth=0.4, label=f"detected values (n={detected.size})")
    for limit in limits:
        ax.axvline(limit, color=CENSORED_BAND, linewidth=1.2, linestyle=(0, (2, 3)), zorder=1)
    if limits.size:
        ax.axvspan(0.0, float(limits.max()), color=CENSORED_BAND, alpha=0.22, zorder=0,
                   label=f"{int(cens.sum())} censored wells (bounds {', '.join(f'{v:g}' for v in limits)} mg/L)")
    for fam in FAMILIES:
        if fam not in fits:
            continue
        dist, ok = fits[fam]
        ax.plot(grid, dist.pdf(grid), label=_label(fam, ok, selected), **_style(fam, ok, selected))
    ax.set(xlabel="Arsenic concentration (mg/L)", ylabel="Density",
           title=f"{district}: detected-value histogram and LS-fitted PDFs\n"
                 "censored wells shown as bounds, not substituted concentrations")
    ax.legend(loc="upper right")
    paths.append(out_dir / f"arsenic_histogram_pdf_{district}.png")
    fig.savefig(paths[-1])
    plt.close(fig)

    # 3. Reverse-KM Q-Q, one panel per candidate.
    fig, axes = plt.subplots(1, len(FAMILIES), figsize=(13, 4.6), sharey=True)
    for ax, fam in zip(axes, FAMILIES):
        if fam not in fits:
            ax.set_axis_off()
            continue
        dist, ok = fits[fam]
        theo = np.asarray(dist.ppf(curve.F), dtype=float)
        finite = np.isfinite(theo) & (theo >= 0)
        lim = (0.0, max(float(np.max(theo[finite])) if finite.any() else hi, float(np.max(curve.x))) * 1.05)
        ax.plot(lim, lim, color=INK_SECONDARY, linewidth=1, linestyle=(0, (3, 3)), label="y = x")
        ax.scatter(theo[mask & finite], curve.x[mask & finite], s=20, color=FAMILY_COLORS[fam],
                   edgecolor="white", linewidth=0.6, label="objective points")
        if (~mask & finite).any():
            ax.scatter(theo[~mask & finite], curve.x[~mask & finite], s=20, facecolor="white",
                       edgecolor=FAMILY_COLORS[fam], linewidth=1.0, label="plateau points")
        ax.set(xlim=lim, ylim=lim, xlabel=f"{fam} theoretical quantile (mg/L)",
               title=_label(fam, ok, selected))
    axes[0].set_ylabel("Reverse-KM empirical quantile (mg/L)")
    axes[0].legend(loc="upper left")
    fig.suptitle(f"{district}: reverse-KM Q-Q by candidate family", fontsize=11)
    paths.append(out_dir / f"arsenic_qq_{district}.png")
    fig.savefig(paths[-1])
    plt.close(fig)

    # 4. Tail zoom: lower-tail CDF and upper-tail survival, log x and log y.
    fig, (lo_ax, hi_ax) = plt.subplots(1, 2, figsize=(11, 4.4))
    lo_grid = np.logspace(np.log10(max(curve.smallest_detected, 1e-6) * 0.5), np.log10(hi), 500)
    hi_grid = np.logspace(np.log10(max(curve.smallest_detected, 1e-6)), np.log10(hi), 500)
    lo_ax.step(curve.x[mask], np.clip(curve.F[mask], 1e-6, None), where="post",
               color=INK_SECONDARY, linewidth=1.4, label="Reverse-KM (objective region)")
    hi_ax.step(curve.x[mask], np.clip(1.0 - curve.F[mask], 1e-6, None), where="post",
               color=INK_SECONDARY, linewidth=1.4, label="Reverse-KM (objective region)")
    for fam in FAMILIES:
        if fam not in fits:
            continue
        dist, ok = fits[fam]
        style = _style(fam, ok, selected)
        lo_ax.plot(lo_grid, np.clip(dist.cdf(lo_grid), 1e-6, None), label=_label(fam, ok, selected), **style)
        hi_ax.plot(hi_grid, np.clip(dist.sf(hi_grid), 1e-6, None), **style)
    lo_ax.set(xscale="log", yscale="log", xlabel="Arsenic (mg/L, log)", ylabel="F(x)",
              title="Lower tail: censoring-aware empirical CDF")
    hi_ax.set(xscale="log", yscale="log", xlabel="Arsenic (mg/L, log)", ylabel="1 - F(x)",
              title="Upper tail: exceedance probability")
    lo_ax.legend(loc="upper left", fontsize=8)
    fig.suptitle(f"{district}: tail agreement of the LS-fitted candidates", fontsize=11)
    paths.append(out_dir / f"arsenic_tails_{district}.png")
    fig.savefig(paths[-1])
    plt.close(fig)
    return paths
