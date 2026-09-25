"""Phase 5: survey-weighted body-weight distribution fitting.

Primary fit: survey-weighted CDF least squares against midpoint plotting positions,

    theta_LS = argmin sum_j q_j [F(x_(j); theta) - p_j]^2,

solved with ``scipy.optimize.least_squares`` on residuals ``sqrt(q_j) [F - p_j]``.

Robustness fit: survey-weighted pseudo-maximum likelihood,

    theta_PML = argmax sum_i w~_i log f(x_i; theta),   w~_i = n w_i / sum w,

solved with ``scipy.optimize.minimize``.

Parameterizations (all Lognormal/Gamma fits use loc = 0):

- Normal(mu, sigma)                 -> norm(loc=mu, scale=sigma)
- Lognormal(mu_log, sigma_log)      -> lognorm(s=sigma_log, loc=0, scale=exp(mu_log))
- Gamma(shape k, scale theta)       -> gamma(a=k, loc=0, scale=theta)
- Triangular(a, m, b), 0 < a < m < b -> triang(c=(m-a)/(b-a), loc=a, scale=b-a)
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import optimize, stats
from scipy.special import expit, logit

FAMILIES = ("Normal", "Lognormal", "Gamma", "Triangular")

PARAM_NAMES = {
    "Normal": ("mu", "sigma"),
    "Lognormal": ("mu_log", "sigma_log"),
    "Gamma": ("shape_k", "scale_theta"),
    "Triangular": ("a_min", "m_mode", "b_max"),
}

# Predeclared admissibility thresholds (set before looking at fit results).
# A Normal fit is inadmissible if the expected number of negative draws in the
# principal N = 10,000 run is at least one.
NORMAL_MAX_NEGATIVE_PROB = 1.0 / 10_000
# A Triangular fit is inadmissible if its support excludes more than 0.1% of survey mass.
TRIANGULAR_MAX_EXCLUDED_MASS = 0.001

TAIL_PROBS = (0.01, 0.05, 0.95, 0.99)
REPORT_PROBS = (0.01, 0.05, 0.50, 0.95, 0.99)

LS_OPTIONS = dict(method="trf", ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=20_000)


# ---------------------------------------------------------------------------
# Weighted empirical CDF
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class WeightedECDF:
    x: np.ndarray  # distinct sorted values
    q: np.ndarray  # normalized survey mass at each value, sums to 1
    p: np.ndarray  # midpoint plotting positions
    F: np.ndarray  # ordinary weighted ECDF (upper edge of each jump)


def weighted_ecdf_midpoints(x, w) -> WeightedECDF:
    """Aggregate survey mass at tied values and return weighted midpoint plotting positions."""
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    if x.shape != w.shape or x.size == 0:
        raise ValueError("x and w must be non-empty arrays of equal shape")
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(w)) and np.all(w > 0)):
        raise ValueError("x must be finite and w finite and positive")
    agg = pd.DataFrame({"x": x, "w": w}).groupby("x", sort=True)["w"].sum()
    x_u = agg.index.to_numpy(dtype=float)
    q_u = agg.to_numpy(dtype=float)
    q_u = q_u / q_u.sum()
    F_u = np.cumsum(q_u)
    p_u = F_u - 0.5 * q_u
    return WeightedECDF(x=x_u, q=q_u, p=p_u, F=F_u)


def weighted_moments(x, w) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    q = np.asarray(w, dtype=float) / np.sum(w)
    mean = float(np.sum(q * x))
    return mean, float(np.sqrt(np.sum(q * (x - mean) ** 2)))


# ---------------------------------------------------------------------------
# Family definitions
# ---------------------------------------------------------------------------

def frozen(family: str, params: dict):
    """Return the frozen SciPy distribution for natural parameters ``params``."""
    if family == "Normal":
        return stats.norm(loc=params["mu"], scale=params["sigma"])
    if family == "Lognormal":
        return stats.lognorm(s=params["sigma_log"], loc=0.0, scale=np.exp(params["mu_log"]))
    if family == "Gamma":
        return stats.gamma(a=params["shape_k"], loc=0.0, scale=params["scale_theta"])
    if family == "Triangular":
        a, m, b = params["a_min"], params["m_mode"], params["b_max"]
        return stats.triang(c=(m - a) / (b - a), loc=a, scale=b - a)
    raise ValueError(f"unknown family {family!r}")


def scipy_spec(family: str, params: dict) -> dict:
    """Explicit SciPy constructor arguments, for downstream samplers."""
    if family == "Normal":
        return {"scipy_dist": "norm", "loc": params["mu"], "scale": params["sigma"]}
    if family == "Lognormal":
        return {"scipy_dist": "lognorm", "s": params["sigma_log"], "loc": 0.0, "scale": float(np.exp(params["mu_log"]))}
    if family == "Gamma":
        return {"scipy_dist": "gamma", "a": params["shape_k"], "loc": 0.0, "scale": params["scale_theta"]}
    if family == "Triangular":
        a, m, b = params["a_min"], params["m_mode"], params["b_max"]
        return {"scipy_dist": "triang", "c": (m - a) / (b - a), "loc": a, "scale": b - a}
    raise ValueError(f"unknown family {family!r}")


def initial_params(family: str, x, w) -> dict:
    """Survey-weighted moment starting values."""
    mean, sd = weighted_moments(x, w)
    if family == "Normal":
        return {"mu": mean, "sigma": sd}
    if family == "Lognormal":
        mlog, slog = weighted_moments(np.log(x), w)
        return {"mu_log": mlog, "sigma_log": slog}
    if family == "Gamma":
        return {"shape_k": mean ** 2 / sd ** 2, "scale_theta": sd ** 2 / mean}
    if family == "Triangular":
        x = np.asarray(x, dtype=float)
        lo, hi = float(x.min()), float(x.max())
        mode = float(np.clip(3 * mean - lo - hi, lo + 1e-3 * (hi - lo), hi - 1e-3 * (hi - lo)))
        return {"a_min": 0.95 * lo, "m_mode": mode, "b_max": 1.05 * hi}
    raise ValueError(f"unknown family {family!r}")


# Unconstrained <-> natural transforms. The Triangular LS transform follows the
# plan (a = e^alpha, m = a + e^beta, b = m + e^gamma). The Triangular PML
# transform additionally forces the support to contain every observation, since
# the likelihood is -inf otherwise.

def _to_natural_ls(family: str, t: np.ndarray, x=None) -> dict:
    if family == "Normal":
        return {"mu": t[0], "sigma": np.exp(t[1])}
    if family == "Lognormal":
        return {"mu_log": t[0], "sigma_log": np.exp(t[1])}
    if family == "Gamma":
        return {"shape_k": np.exp(t[0]), "scale_theta": np.exp(t[1])}
    a = np.exp(t[0])
    m = a + np.exp(t[1])
    return {"a_min": a, "m_mode": m, "b_max": m + np.exp(t[2])}


def _from_natural_ls(family: str, p: dict, x=None) -> np.ndarray:
    if family == "Normal":
        return np.array([p["mu"], np.log(p["sigma"])])
    if family == "Lognormal":
        return np.array([p["mu_log"], np.log(p["sigma_log"])])
    if family == "Gamma":
        return np.array([np.log(p["shape_k"]), np.log(p["scale_theta"])])
    a, m, b = p["a_min"], p["m_mode"], p["b_max"]
    return np.array([np.log(a), np.log(m - a), np.log(b - m)])


def _to_natural_pml(family: str, t: np.ndarray, x) -> dict:
    if family != "Triangular":
        return _to_natural_ls(family, t)
    lo, hi = float(np.min(x)), float(np.max(x))
    a = lo * expit(t[0])
    b = hi + np.exp(t[2])
    return {"a_min": a, "m_mode": a + (b - a) * expit(t[1]), "b_max": b}


def _from_natural_pml(family: str, p: dict, x) -> np.ndarray:
    if family != "Triangular":
        return _from_natural_ls(family, p)
    # A start whose support does not cover the data (e.g. an LS fit) is pulled back
    # with a real margin; pinning it at the data edge saturates the sigmoid.
    lo, hi = float(np.min(x)), float(np.max(x))
    margin = 0.02 * (hi - lo)
    a = p["a_min"] if p["a_min"] < lo - 1e-3 * margin else max(lo - margin, 0.5 * lo)
    b = p["b_max"] if p["b_max"] > hi + 1e-3 * margin else hi + margin
    c = np.clip((p["m_mode"] - a) / (b - a), 1e-6, 1 - 1e-6)
    return np.array([logit(a / lo), logit(c), np.log(b - hi)])


def _plain(p: dict) -> dict:
    return {k: float(v) for k, v in p.items()}


# ---------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------

@dataclass
class FitResult:
    family: str
    method: str
    params: dict
    success: bool
    message: str
    objective: float
    n_starts: int = 1
    extra: dict = field(default_factory=dict)

    def dist(self):
        return frozen(self.family, self.params)


def cdf_ls_residuals(t, family: str, ecdf: WeightedECDF) -> np.ndarray:
    """sqrt(q_j) * [F_theta(x_j) - p_j]."""
    F = frozen(family, _to_natural_ls(family, t)).cdf(ecdf.x)
    return np.sqrt(ecdf.q) * (F - ecdf.p)


def weighted_sse(dist, ecdf: WeightedECDF) -> float:
    return float(np.sum(ecdf.q * (dist.cdf(ecdf.x) - ecdf.p) ** 2))


def weighted_max_discrepancy(dist, ecdf: WeightedECDF) -> float:
    return float(np.max(np.abs(dist.cdf(ecdf.x) - ecdf.p)))


def _triangular_ls_starts(x, w) -> list[dict]:
    ecdf = weighted_ecdf_midpoints(x, w)
    lo, hi = float(ecdf.x[0]), float(ecdf.x[-1])
    q01, q50, q99 = np.interp([0.01, 0.5, 0.99], ecdf.p, ecdf.x)
    mode_mass = float(ecdf.x[np.argmax(ecdf.q)])
    starts = [initial_params("Triangular", x, w)]
    for a0 in (0.9 * lo, max(q01 - (q50 - q01), 0.05 * lo)):
        for b0 in (1.05 * hi, q99 + (q99 - q50)):
            for m0 in (q50, mode_mass):
                if 0 < a0 < m0 < b0:
                    starts.append({"a_min": a0, "m_mode": m0, "b_max": b0})
    return starts


def fit_cdf_least_squares(x, w, family: str) -> FitResult:
    """Survey-weighted CDF least-squares fit of one candidate family."""
    ecdf = weighted_ecdf_midpoints(x, w)
    starts = _triangular_ls_starts(x, w) if family == "Triangular" else [initial_params(family, x, w)]
    best = None
    for start in starts:
        res = optimize.least_squares(cdf_ls_residuals, _from_natural_ls(family, start),
                                     args=(family, ecdf), **LS_OPTIONS)
        sse = float(np.sum(res.fun ** 2))
        if best is None or sse < best[1] - 1e-15:
            best = (res, sse)
    res, sse = best
    return FitResult(
        family=family,
        method="weighted_cdf_least_squares",
        params=_plain(_to_natural_ls(family, res.x)),
        success=bool(res.success),
        message=str(res.message),
        objective=sse,
        n_starts=len(starts),
        extra={"nfev": int(res.nfev)},
    )


def negative_weighted_loglik(t, family: str, x: np.ndarray, w_pml: np.ndarray) -> float:
    """Negative survey-weighted pseudo-log-likelihood with weights normalized to sum to n."""
    logf = frozen(family, _to_natural_pml(family, t, x)).logpdf(x)
    if not np.all(np.isfinite(logf)):
        return np.inf
    return float(-np.sum(w_pml * logf))


def fit_weighted_pseudo_mle(x, w, family: str, ls_start: dict | None = None) -> FitResult:
    """Robustness fit by survey-weighted pseudo-maximum likelihood."""
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    w_pml = len(w) * w / w.sum()
    starts = [initial_params(family, x, w)]
    if ls_start is not None:
        starts.append(ls_start)
    best = None
    for start in starts:
        t0 = _from_natural_pml(family, start, x)
        if family == "Triangular":
            # Nelder-Mead on the non-smooth triangular likelihood, restarted once from its own solution.
            nm = {"xatol": 1e-10, "fatol": 1e-10, "maxiter": 20_000, "maxfev": 40_000}
            res = optimize.minimize(negative_weighted_loglik, t0, args=(family, x, w_pml), method="Nelder-Mead", options=nm)
            res = optimize.minimize(negative_weighted_loglik, res.x, args=(family, x, w_pml), method="Nelder-Mead", options=nm)
        else:
            res = optimize.minimize(negative_weighted_loglik, t0, args=(family, x, w_pml), method="L-BFGS-B",
                                    options={"ftol": 1e-15, "gtol": 1e-10, "maxiter": 10_000})
        if np.isfinite(res.fun) and (best is None or res.fun < best.fun):
            best = res
    if best is None:
        return FitResult(family, "weighted_pseudo_mle", {k: np.nan for k in PARAM_NAMES[family]},
                         False, "no finite objective from any start", np.inf, len(starts))
    return FitResult(
        family=family,
        method="weighted_pseudo_mle",
        params=_plain(_to_natural_pml(family, best.x, x)),
        success=bool(best.success),
        message=str(best.message),
        objective=float(best.fun),
        n_starts=len(starts),
    )


# ---------------------------------------------------------------------------
# Evaluation and selection
# ---------------------------------------------------------------------------

def evaluate_distribution_fit(x, w, ls_fit: FitResult, pml_fit: FitResult) -> dict:
    """SSE_w, D_w, tail quantiles, support checks, and LS-vs-PML agreement for one family."""
    family = ls_fit.family
    ecdf = weighted_ecdf_midpoints(x, w)
    d_ls = ls_fit.dist()
    row: dict = {"distribution": family}
    for name in PARAM_NAMES[family]:
        row[f"ls_{name}"] = ls_fit.params[name]
    row["ls_converged"] = ls_fit.success
    row["ls_message"] = ls_fit.message
    row["ls_n_starts"] = ls_fit.n_starts
    row["weighted_cdf_sse"] = weighted_sse(d_ls, ecdf)
    row["weighted_cdf_rmse"] = float(np.sqrt(row["weighted_cdf_sse"]))
    row["weighted_cdf_max_discrepancy"] = weighted_max_discrepancy(d_ls, ecdf)

    for name in PARAM_NAMES[family]:
        row[f"pml_{name}"] = pml_fit.params[name]
    row["pml_converged"] = pml_fit.success
    row["pml_message"] = pml_fit.message
    row["pml_neg_weighted_loglik"] = pml_fit.objective
    if np.isfinite(pml_fit.objective):
        d_pml = pml_fit.dist()
        row["pml_weighted_cdf_sse"] = weighted_sse(d_pml, ecdf)
        row["pml_weighted_cdf_max_discrepancy"] = weighted_max_discrepancy(d_pml, ecdf)
    else:
        row["pml_weighted_cdf_sse"] = row["pml_weighted_cdf_max_discrepancy"] = np.nan
    rel = [abs(ls_fit.params[k] - pml_fit.params[k]) / abs(pml_fit.params[k]) for k in PARAM_NAMES[family]]
    row["ls_pml_max_rel_param_diff"] = float(np.max(rel))

    emp_q = np.interp(REPORT_PROBS, ecdf.p, ecdf.x)
    fit_q = d_ls.ppf(REPORT_PROBS)
    for prob, e, f in zip(REPORT_PROBS, emp_q, fit_q):
        tag = f"P{int(round(prob * 100))}"
        row[f"empirical_{tag}_kg"] = float(e)
        row[f"fitted_{tag}_kg"] = float(f)
    tail_idx = [REPORT_PROBS.index(p) for p in TAIL_PROBS]
    row["tail_max_rel_quantile_error"] = float(np.max(np.abs(fit_q[tail_idx] - emp_q[tail_idx]) / emp_q[tail_idx]))

    # Physical support
    row["prob_bw_nonpositive"] = float(d_ls.cdf(0.0))
    lo, hi = d_ls.support()
    outside = (ecdf.x < lo) | (ecdf.x > hi)
    row["survey_mass_outside_support"] = float(ecdf.q[outside].sum())
    reasons = []
    if not ls_fit.success:
        reasons.append("LS optimizer did not report convergence")
    if family == "Normal" and row["prob_bw_nonpositive"] >= NORMAL_MAX_NEGATIVE_PROB:
        reasons.append(f"P(BW<=0)={row['prob_bw_nonpositive']:.2e} >= {NORMAL_MAX_NEGATIVE_PROB:.0e}")
    if family == "Triangular" and row["survey_mass_outside_support"] > TRIANGULAR_MAX_EXCLUDED_MASS:
        reasons.append(f"support excludes {row['survey_mass_outside_support']:.4%} of survey mass")
    if family in ("Lognormal", "Gamma", "Triangular") and lo < 0:
        reasons.append("support extends below zero")
    row["admissible"] = not reasons
    row["inadmissible_reason"] = "; ".join(reasons)
    return row


def fit_population(x, w, population: str) -> pd.DataFrame:
    """Fit all four candidates by weighted LS and weighted PML and evaluate each."""
    rows = []
    for family in FAMILIES:
        ls = fit_cdf_least_squares(x, w, family)
        pml = fit_weighted_pseudo_mle(x, w, family, ls_start=ls.params)
        row = evaluate_distribution_fit(x, w, ls, pml)
        rows.append({"population": population, "n": int(len(x)), **row})
    table = pd.DataFrame(rows)
    table["sse_rank"] = table["weighted_cdf_sse"].rank(method="min").astype(int)
    table["max_discrepancy_rank"] = table["weighted_cdf_max_discrepancy"].rank(method="min").astype(int)
    table["tail_error_rank"] = table["tail_max_rel_quantile_error"].rank(method="min").astype(int)
    return table


def select_bw_distribution(fit_table: pd.DataFrame) -> tuple[pd.Series, str]:
    """Predeclared rule: lowest weighted CDF SSE among physically admissible candidates.

    Returns the selected row and a rationale noting any criterion that disagrees.
    The result is provisional and requires human review before it is final.
    """
    admissible = fit_table[fit_table["admissible"]]
    if admissible.empty:
        raise RuntimeError(f"no admissible body-weight candidate for {fit_table['population'].iloc[0]}")
    ranked = admissible.sort_values("weighted_cdf_sse")
    chosen = ranked.iloc[0]
    notes = [f"lowest weighted CDF SSE among admissible candidates ({chosen['weighted_cdf_sse']:.3e})"]
    overall_best = fit_table.sort_values("weighted_cdf_sse").iloc[0]
    if overall_best["distribution"] != chosen["distribution"]:
        notes.append(f"{overall_best['distribution']} had lower SSE but was rejected: {overall_best['inadmissible_reason']}")
    for rejected in fit_table[~fit_table["admissible"]].itertuples():
        if rejected.distribution != overall_best["distribution"]:
            notes.append(f"{rejected.distribution} inadmissible: {rejected.inadmissible_reason}")
    best_d = admissible.sort_values("weighted_cdf_max_discrepancy").iloc[0]["distribution"]
    notes.append("also lowest max CDF discrepancy" if best_d == chosen["distribution"]
                 else f"max CDF discrepancy favours {best_d}")
    best_tail = admissible.sort_values("tail_max_rel_quantile_error").iloc[0]["distribution"]
    notes.append("also best tail quantile agreement" if best_tail == chosen["distribution"]
                 else f"tail quantile agreement favours {best_tail}")
    if len(ranked) > 1:
        runner = ranked.iloc[1]
        notes.append(f"runner-up {runner['distribution']} SSE ratio {runner['weighted_cdf_sse'] / chosen['weighted_cdf_sse']:.2f}")
    notes.append(f"LS vs PML max relative parameter difference {chosen['ls_pml_max_rel_param_diff']:.2%}")
    return chosen, "; ".join(notes)


def selected_record(chosen: pd.Series, rationale: str) -> dict:
    family = chosen["distribution"]
    params = {name: float(chosen[f"ls_{name}"]) for name in PARAM_NAMES[family]}
    return {
        "population": chosen["population"],
        "distribution": family,
        "parameterization": "natural parameters from weighted CDF least squares; Lognormal/Gamma loc=0",
        "params": params,
        "scipy": {k: (float(v) if not isinstance(v, str) else v) for k, v in scipy_spec(family, params).items()},
        "weighted_cdf_sse": float(chosen["weighted_cdf_sse"]),
        "weighted_cdf_max_discrepancy": float(chosen["weighted_cdf_max_discrepancy"]),
        "n": int(chosen["n"]),
        "status": "provisional - requires team review before freezing",
        "rationale": rationale,
    }


# ---------------------------------------------------------------------------
# Inverse-CDF sampling for Monte Carlo
# ---------------------------------------------------------------------------

def sample_bw(selected: dict, N: int | None = None, rng: np.random.Generator | None = None,
              u: np.ndarray | None = None) -> np.ndarray:
    """Draw body weights in kg by explicit inverse-CDF sampling: U ~ Uniform(0,1), BW = F^{-1}(U).

    Pass either ``u`` directly or ``N`` with a ``numpy.random.Generator``.
    """
    if u is None:
        if N is None or rng is None:
            raise ValueError("provide u, or both N and rng")
        u = rng.uniform(0.0, 1.0, size=N)
    return frozen(selected["distribution"], selected["params"]).ppf(u)
