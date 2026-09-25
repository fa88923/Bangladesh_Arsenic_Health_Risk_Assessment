"""Phase 7: final sampling distributions for C and BW.

Starts from the provisional Phase 4 (arsenic) and Phase 5 (body-weight) selections
and applies the Member 4 decisions recorded in ``simulation.input_selection`` of
config/run_config.json:

- ``arsenic_family_overrides``: districts whose provisional family is replaced by
  another Phase 4 candidate, using that candidate's saved CDF least-squares fit.
- ``child_bw``: a Normal left-truncated at ``lower_kg``, refitted here by the same
  survey-weighted CDF least squares as Phase 5 (robustness: weighted pseudo-MLE).

Parameterizations (SciPy):

- TruncatedNormal(mu, sigma, lower) -> truncnorm(a=(lower-mu)/sigma, b=inf, loc=mu, scale=sigma)
  (mu and sigma are the parent Normal's parameters, not the truncated moments)
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import optimize, stats

from arsenic_hra import arsenic_fitting as af
from arsenic_hra import bodyweight_fitting as bf
from arsenic_hra.bodyweight_preprocessing import weighted_quantile
from arsenic_hra.paths import RESULTS_ARSENIC, RESULTS_BODYWEIGHT

ARSENIC_SELECTED = RESULTS_ARSENIC / "arsenic_selected_distributions.json"
ARSENIC_FITS = RESULTS_ARSENIC / "arsenic_distribution_fit_results.csv"
ARSENIC_PREPROCESSED = RESULTS_ARSENIC / "arsenic_selected_districts_preprocessed.csv"
BW_SELECTED = RESULTS_BODYWEIGHT / "bodyweight_selected_distributions.json"
CHILD_BW_CLEAN = RESULTS_BODYWEIGHT / "child_bw_clean.csv"

TAIL_PROBS = (0.01, 0.05, 0.50, 0.95, 0.99)


# ---------------------------------------------------------------------------
# Distributions
# ---------------------------------------------------------------------------

class TruncatedNormalMixture:
    """Two-component mixture of lower-truncated Normals sharing one bound (Phase 9 age scenario).

    F(x) = w1 F1(x) + (1 - w1) F2(x). The inverse CDF is evaluated by monotone interpolation
    on a dense grid (``GRID`` points spanning the 1e-12 and 1 - 1e-12 quantiles of the mixture),
    which keeps sampling explicit U -> F^-1(U).
    """
    GRID = 200_001

    def __init__(self, w1: float, mu1: float, sigma1: float, mu2: float, sigma2: float, lower: float):
        self.w1 = float(w1)
        self.c1 = frozen("TruncatedNormal", {"mu": mu1, "sigma": sigma1, "lower": lower})
        self.c2 = frozen("TruncatedNormal", {"mu": mu2, "sigma": sigma2, "lower": lower})
        hi = max(self.c1.ppf(1 - 1e-12), self.c2.ppf(1 - 1e-12))
        self._x = np.linspace(lower, hi, self.GRID)
        self._F = self.cdf(self._x)

    def cdf(self, x):
        return self.w1 * self.c1.cdf(x) + (1 - self.w1) * self.c2.cdf(x)

    def ppf(self, u):
        return np.interp(np.asarray(u, dtype=float), self._F, self._x)

    def mean(self):
        return self.w1 * self.c1.mean() + (1 - self.w1) * self.c2.mean()

    def median(self):
        return float(self.ppf(0.5))

    def support(self):
        return (float(self._x[0]), float("inf"))


def frozen(family: str, params: dict):
    """Frozen SciPy distribution for any family used by the simulation."""
    if family == "TruncatedNormalMixture":
        return TruncatedNormalMixture(**params)
    if family == "TruncatedNormal":
        mu, sigma, lower = params["mu"], params["sigma"], params["lower"]
        return stats.truncnorm(a=(lower - mu) / sigma, b=np.inf, loc=mu, scale=sigma)
    if family == "Triangular":
        a, m, b = params["a_min"], params["m_mode"], params["b_max"]
        return stats.triang(c=(m - a) / (b - a), loc=a, scale=b - a)
    return af.frozen(family, params)  # Normal, Lognormal, Gamma (loc = 0)


def scipy_spec(family: str, params: dict) -> dict:
    if family == "TruncatedNormal":
        mu, sigma, lower = params["mu"], params["sigma"], params["lower"]
        return {"scipy_dist": "truncnorm", "a": (lower - mu) / sigma, "b": float("inf"), "loc": mu, "scale": sigma}
    return af.scipy_spec(family, params)


# ---------------------------------------------------------------------------
# Truncated-Normal child body weight (weighted CDF least squares, as Phase 5)
# ---------------------------------------------------------------------------

def _truncnorm(t, lower: float):
    return frozen("TruncatedNormal", {"mu": t[0], "sigma": np.exp(t[1]), "lower": lower})


def fit_truncated_normal_bw(x, w, lower: float) -> dict:
    """Fit a lower-truncated Normal to survey-weighted body weights.

    Primary: argmin sum_j q_j [F(x_j) - p_j]^2 on weighted midpoint plotting positions.
    Robustness: weighted pseudo-MLE with w~ = n w / sum w. ``lower`` is fixed, not fitted.
    """
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    if np.any(x < lower):
        raise ValueError(f"observations below the truncation bound {lower} kg")
    e = bf.weighted_ecdf_midpoints(x, w)
    mean, sd = bf.weighted_moments(x, w)
    start = np.array([mean, np.log(sd)])
    ls = optimize.least_squares(lambda t: np.sqrt(e.q) * (_truncnorm(t, lower).cdf(e.x) - e.p),
                                start, **bf.LS_OPTIONS)
    w_pml = x.size * w / w.sum()
    pml = optimize.minimize(lambda t: -np.sum(w_pml * _truncnorm(t, lower).logpdf(x)), ls.x,
                            method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 10_000})
    params = {"mu": float(ls.x[0]), "sigma": float(np.exp(ls.x[1])), "lower": float(lower)}
    pml_params = {"mu": float(pml.x[0]), "sigma": float(np.exp(pml.x[1])), "lower": float(lower)}
    return {"params": params, "pml_params": pml_params, "ls_converged": bool(ls.success),
            "ls_message": str(ls.message), "pml_converged": bool(pml.success), "pml_message": str(pml.message)}


def bw_fit_metrics(family: str, params: dict, x, w) -> dict:
    """Weighted SSE, D_w, tail-quantile agreement, and E[1/BW] (dose scales with 1/BW)."""
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    d = frozen(family, params)
    e = bf.weighted_ecdf_midpoints(x, w)
    F = d.cdf(e.x)
    emp_q = weighted_quantile(x, w, TAIL_PROBS)
    fit_q = d.ppf(TAIL_PROBS)
    lo = float(d.support()[0])
    inv_bw = float(d.expect(lambda t: 1.0 / t, lb=max(lo, 1e-12), ub=float(d.ppf(1 - 1e-12)), limit=200))
    row = {"family": family, "weighted_cdf_sse": float(np.sum(e.q * (F - e.p) ** 2)),
           "weighted_cdf_max_discrepancy": float(np.max(np.abs(F - e.p))),
           "tail_max_rel_quantile_error": float(np.max(np.abs(fit_q - emp_q) / emp_q)),
           "P(BW<=0)": float(d.cdf(0.0)), "P(BW<1kg)": float(d.cdf(1.0)),
           "E_inv_BW_per_kg": inv_bw, "empirical_E_inv_BW_per_kg": float(np.sum(w / x) / np.sum(w))}
    for p, eq, fq in zip(TAIL_PROBS, emp_q, fit_q):
        tag = f"P{int(round(p * 100))}"
        row[f"empirical_{tag}_kg"], row[f"fitted_{tag}_kg"] = float(eq), float(fq)
    return row


# ---------------------------------------------------------------------------
# Censoring-aware district mean for the deterministic benchmark
# ---------------------------------------------------------------------------

def reverse_km_mean(values_mgL, censored) -> float:
    """E[C] = integral_0^max (1 - F_RKM(t)) dt, with F = F_hat_0 on the unresolved plateau.

    This is the lower bound of the KM mean (censored mass placed at 0); placing it at the
    smallest detection instead adds F_hat_0 * smallest_detected, < 0.0005 mg/L here.
    """
    c = af.reverse_km_cdf(values_mgL, censored)
    order = np.argsort(c.x)
    x, F = c.x[order], c.F[order]
    keep = x > 0
    grid = np.r_[0.0, x[keep]]
    heights = np.r_[c.censoring_fraction, F[keep]]
    return float(np.sum((1.0 - heights[:-1]) * np.diff(grid)))


# ---------------------------------------------------------------------------
# Assemble the final input records
# ---------------------------------------------------------------------------

def _record(kind: str, key: str, family: str, params: dict, source: str, decision: str, **extra) -> dict:
    return {kind: key, "distribution": family, "params": {k: float(v) for k, v in params.items()},
            "scipy": scipy_spec(family, params), "source": source, "decision": decision, **extra}


def arsenic_inputs(cfg: dict) -> dict[str, dict]:
    """District -> final arsenic sampling record (mg/L)."""
    overrides = cfg["simulation"]["input_selection"]["arsenic_family_overrides"]
    provisional = {r["district"]: r for r in json.loads(ARSENIC_SELECTED.read_text(encoding="utf-8"))}
    fits = pd.read_csv(ARSENIC_FITS)
    pre = pd.read_csv(ARSENIC_PREPROCESSED)
    out = {}
    for district in cfg["arsenic"]["selected_districts"]:
        prov = provisional[district]
        g = pre[pre.DISTRICT == district]
        km_mean = reverse_km_mean(g["As_model_value_mgL"].to_numpy(), g["As_censored"].astype(bool).to_numpy())
        common = {"n": prov["n"], "censoring_percent": prov["censoring_percent"],
                  "reverse_km_mean_mgL": km_mean, "provisional_distribution": prov["distribution"]}
        family = overrides.get(district)
        if family is None or family == prov["distribution"]:
            out[district] = _record("district", district, prov["distribution"], prov["params"],
                                    "Phase 4 provisional selection (CDF least squares)",
                                    "accepted as provisional", **common)
            continue
        row = fits[(fits.DISTRICT == district) & (fits.distribution == family)].iloc[0]
        if not bool(row["admissible"]):
            raise ValueError(f"{district}: override family {family} is inadmissible in Phase 4")
        names = af_param_names(family)
        params = {name: float(row[f"ls_{name}"]) for name in names}
        out[district] = _record("district", district, family, params,
                                "Phase 4 fit table, CDF least-squares candidate",
                                f"override of provisional {prov['distribution']} (see simulation report D7)", **common)
    return out


def af_param_names(family: str) -> tuple[str, ...]:
    return {"Normal": ("mu", "sigma"), "Lognormal": ("mu_log", "sigma_log"),
            "Gamma": ("shape_k", "scale_theta")}[family]


def bodyweight_inputs(cfg: dict) -> tuple[dict[str, dict], pd.DataFrame]:
    """Population -> final BW sampling record (kg), plus the child-BW comparison table."""
    provisional = {r["population"]: r for r in json.loads(BW_SELECTED.read_text(encoding="utf-8"))}
    out = {"adult": _record("population", "adult", provisional["adult"]["distribution"],
                            provisional["adult"]["params"], "Phase 5 provisional selection",
                            "accepted as provisional")}
    spec = cfg["simulation"]["input_selection"]["child_bw"]
    child = pd.read_csv(CHILD_BW_CLEAN)
    x, w = child["BW_kg"].to_numpy(), child["survey_weight"].to_numpy()
    comparison = [dict(bw_fit_metrics(provisional["child"]["distribution"], provisional["child"]["params"], x, w),
                       role="Phase 5 provisional")]
    if spec["family"] == "TruncatedNormal":
        fit = fit_truncated_normal_bw(x, w, spec["lower_kg"])
        comparison.append(dict(bw_fit_metrics("TruncatedNormal", fit["params"], x, w), role="selected",
                               ls_converged=fit["ls_converged"], pml_converged=fit["pml_converged"],
                               pml_mu=fit["pml_params"]["mu"], pml_sigma=fit["pml_params"]["sigma"],
                               ls_pml_max_rel_param_diff=max(
                                   abs(fit["pml_params"][k] - fit["params"][k]) / fit["params"][k]
                                   for k in ("mu", "sigma"))))
        out["child"] = _record("population", "child", "TruncatedNormal", fit["params"],
                               "Phase 7 refit: weighted CDF least squares on child_bw_clean.csv",
                               f"replaces provisional {provisional['child']['distribution']} (see simulation report D8)",
                               ls_converged=fit["ls_converged"])
    else:
        comparison[0]["role"] = "selected"
        out["child"] = _record("population", "child", provisional["child"]["distribution"],
                               provisional["child"]["params"], "Phase 5 provisional selection", "accepted")
    return out, pd.DataFrame(comparison)
