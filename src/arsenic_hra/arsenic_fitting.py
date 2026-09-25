"""Phase 4: censoring-aware district arsenic distribution fitting.

Implements `plan.md` sections 3.A-3.F and
`docs/arsenic_preprocessing_and_distribution_plan.md` sections 8-13.

Step A  reverse Kaplan-Meier empirical CDF          ``reverse_km_cdf``
Step B  primary fit, CDF least squares              ``fit_cdf_least_squares``
Step C  robustness fit, censored maximum likelihood ``fit_censored_mle``
Step D  bootstrap stability                         ``arsenic_bootstrap``
Step E  admissibility and selection                 ``evaluate_district_fit`` / ``select_arsenic_distribution``
Step F  sampling interface                          ``sample_arsenic``

Censored wells are never substituted. A `< L` row enters the reverse-KM curve
through its bound ``L`` only, and the censored likelihood uses ``log F(L)``.

Objective region (frozen decision 9.1, variant V2 "above_initial_mass"): the
reverse-KM curve is flat at height ``F_hat_0`` for every ``x`` below the smallest
detected value and opens with a placeholder point at ``x = 0``. Those points
carry no shape information about the fitted density, and the placeholder pulls
the fit, so ``objective_mask`` drops them.

Parameterizations (Lognormal and Gamma are forced to zero location):

- Normal(mu, sigma)            -> norm(loc=mu, scale=sigma)
- Lognormal(mu_log, sigma_log) -> lognorm(s=sigma_log, loc=0, scale=exp(mu_log))
- Gamma(shape_k, scale_theta)  -> gamma(a=shape_k, loc=0, scale=scale_theta)
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from scipy import optimize, stats

FAMILIES = ("Normal", "Lognormal", "Gamma")

PARAM_NAMES = {
    "Normal": ("mu", "sigma"),
    "Lognormal": ("mu_log", "sigma_log"),
    "Gamma": ("shape_k", "scale_theta"),
}

# Families whose support must start exactly at zero.
ZERO_LOCATION_FAMILIES = ("Lognormal", "Gamma")

# Predeclared constants, mirrored by the `arsenic` block of config/run_config.json.
# A Normal fit is inadmissible if its censored-MLE reproduction puts at least this
# much probability mass at or below zero (this is the body-weight module's gate).
NORMAL_MAX_NONPOSITIVE_PROB = 1e-4
# Fitted parameters must be finite and inside this magnitude band.
PARAM_MIN_MAGNITUDE = 1e-8
PARAM_MAX_MAGNITUDE = 1e6
# Admissibility/selection constants pinned from the frozen config when the runner
# calls ``configure_constants``; the values above are the fallbacks.
CONFIGURED: dict[str, float] = {}

OBJECTIVE_REGION = "above_initial_mass"
OBJECTIVE_REGION_TIE = 1e-12
SELECTION_RULE = "lowest_cdf_rmse_among_admissible"
SELECTION_STATUS = "provisional - requires team review before freezing"

TAIL_PROBS = (0.01, 0.05, 0.95, 0.99)
REPORT_PROBS = (0.01, 0.05, 0.50, 0.95, 0.99)

# Multi-start least squares. Tight tolerances are affordable (the objective is a
# ~20-50 point residual vector) and keep the saved parameters reproducible.
LS_OPTIONS = dict(method="trf", ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=20_000)
# A deliberately looser configuration for the bootstrap inner loop.
BOOTSTRAP_LS_OPTIONS = dict(method="trf", ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=2_000)


def configure_constants(acfg: dict) -> None:
    """Pin the admissibility/selection constants from the run configuration.

    The runner calls this with the `arsenic` block, so the frozen values live in
    config/run_config.json and are auditable; the module-level defaults apply when
    the module is used standalone (tests, notebooks).
    """
    global NORMAL_MAX_NONPOSITIVE_PROB, PARAM_MIN_MAGNITUDE, PARAM_MAX_MAGNITUDE
    global OBJECTIVE_REGION, SELECTION_RULE, SELECTION_STATUS
    NORMAL_MAX_NONPOSITIVE_PROB = float(acfg["normal_max_nonpositive_prob"])
    band = acfg["parameter_magnitude_band"]
    PARAM_MIN_MAGNITUDE, PARAM_MAX_MAGNITUDE = float(band["min"]), float(band["max"])
    OBJECTIVE_REGION = str(acfg["objective_region"])
    SELECTION_RULE = str(acfg["selection_rule"])
    SELECTION_STATUS = str(acfg["selection_status"])
    CONFIGURED.update(
        normal_max_nonpositive_prob=NORMAL_MAX_NONPOSITIVE_PROB,
        parameter_min_magnitude=PARAM_MIN_MAGNITUDE,
        parameter_max_magnitude=PARAM_MAX_MAGNITUDE,
        objective_region=OBJECTIVE_REGION,
        selection_rule=SELECTION_RULE,
    )


# ---------------------------------------------------------------------------
# Step A - reverse Kaplan-Meier empirical CDF
# ---------------------------------------------------------------------------

def _censoring_above_detection(values_mgL, censored) -> bool:
    """True when a censoring bound exceeds the largest detected concentration.

    `lifelines`' reverse-KM estimator is only valid while the censored mass sits at
    or below every detection: it sorts the event table in descending order and takes
    ``exp()`` of a cumulative hazard, so a bound above the largest detection makes
    the running sum decrease, ``exp()`` saturates, and the plateau collapses to 0.
    The BGS data reports ``< 0.006`` and can also report a detection of exactly
    0.006, so this is detectable rather than hypothetical. Flagged rows are surfaced
    in the outputs rather than silently absorbed.
    """
    x = np.asarray(values_mgL, dtype=float)
    cens = np.asarray(censored, dtype=bool)
    if not cens.any() or not (~cens).any():
        return False
    return bool(x[cens].max() > x[~cens].max())


def censoring_plateau_note(curve: RKMCurve) -> str:
    """Describe how the unresolved low region relates to the raw censoring fraction.

    ``F_hat_0`` is generally *below* the raw censored share, because some detections
    sit below the largest censoring bound and those detections carry real
    information that lifts the estimate above a pure censored-only plateau. The
    reverse-KM curve is the correct quantity; the raw share is only a naive upper
    reference. When ``F_hat_0 = 0`` every bound sits at or above every detection, so
    the low tail is fully resolved and no plateau caveat is warranted.
    """
    if curve.censoring_fraction <= 0.0:
        return ("the censoring fraction is fully resolved (F_hat_0 = 0): every censoring bound is at or "
                "above every detection, so detections constrain the low tail and no unresolved region "
                "remains")
    raw = curve.n_censored / curve.n
    if curve.censoring_fraction + 0.005 < raw:
        return (f"the reverse-KM plateau is F_hat_0 = {curve.censoring_fraction:.4f}, below the raw "
                f"censored share {raw:.4f}, because detections below the largest bound still constrain "
                f"the low tail; the region below {curve.smallest_detected:.4g} mg/L remains unresolved")
    return (f"the reverse-KM plateau is F_hat_0 = {curve.censoring_fraction:.4f}, essentially the raw "
            f"censored share {raw:.4f}; the region below {curve.smallest_detected:.4g} mg/L is unresolved")


@dataclass(frozen=True)
class RKMCurve:
    """Reverse-KM step curve: locations ``x`` (mg/L) and heights ``F_hat``."""

    x: np.ndarray
    F: np.ndarray
    censoring_fraction: float  # flat height of the initial plateau, F_hat_0
    n: int
    n_censored: int
    n_detected: int
    smallest_detected: float
    censoring_above_detection: bool = False


def reverse_km_cdf(values_mgL, censored) -> RKMCurve:
    """Reverse (left-censored) Kaplan-Meier empirical CDF of arsenic in mg/L.

    ``censored[i]`` is True when row ``i`` is a `< L` observation, in which case
    ``values_mgL[i]`` is its censoring bound and is used only as an upper bound.

    Raises ValueError when every value is censored, when none is censored, or when
    the arrays disagree; the curve is meaningful only with both kinds of row.
    """
    x = np.asarray(values_mgL, dtype=float)
    cens = np.asarray(censored, dtype=bool)
    if x.shape != cens.shape or x.size == 0:
        raise ValueError("values_mgL and censored must be non-empty arrays of equal shape")
    if not np.all(np.isfinite(x)):
        raise ValueError("values_mgL must be finite")
    if np.any(x < 0):
        raise ValueError("arsenic concentrations must be non-negative")
    n_censored = int(cens.sum())
    n_detected = int(x.size - n_censored)
    if n_censored == x.size:
        raise ValueError("every value is censored: the reverse-KM curve is not identifiable")
    if n_censored == 0:
        raise ValueError("no censored values: use an ordinary empirical CDF instead")
    smallest_detected = float(x[~cens].min())
    censoring_above_detection = _censoring_above_detection(x, cens)

    kmf = KaplanMeierFitter()
    kmf.fit_left_censoring(durations=x, event_observed=~cens)
    xj = kmf.cumulative_density_.index.to_numpy(dtype=float)
    Fj = kmf.cumulative_density_.iloc[:, 0].to_numpy(dtype=float)
    if xj.size == 0:
        raise ValueError("reverse-KM returned an empty curve")
    if not (np.all(np.isfinite(xj)) and np.all(np.isfinite(Fj))):
        raise ValueError("reverse-KM returned non-finite points")
    if np.any(np.diff(Fj) < -1e-12):
        raise ValueError("reverse-KM curve is not monotone non-decreasing")
    return RKMCurve(
        x=xj,
        F=np.clip(Fj, 0.0, 1.0),
        censoring_fraction=float(np.clip(Fj[0], 0.0, 1.0)),
        n=int(x.size),
        n_censored=n_censored,
        n_detected=n_detected,
        smallest_detected=smallest_detected,
        censoring_above_detection=censoring_above_detection,
    )


def objective_mask(curve: RKMCurve, region: str | None = None) -> np.ndarray:
    """Boolean mask of the reverse-KM points used by the CDF least-squares objective.

    ``"above_initial_mass"`` (frozen decision 9.1 V2) keeps only points strictly
    above the initial plateau height. ``"all"`` (V1) keeps everything, including
    the ``x = 0`` placeholder. ``"above_max_limit"`` (V3) keeps only points at or
    above the largest censoring bound, i.e. where every point is a detection.
    """
    region = OBJECTIVE_REGION if region is None else region
    F = curve.F
    if region == "all":
        return np.ones(F.size, dtype=bool)
    if region == "above_initial_mass":
        return F > curve.censoring_fraction + OBJECTIVE_REGION_TIE
    if region == "above_max_limit":
        return curve.x >= curve.smallest_detected - OBJECTIVE_REGION_TIE
    raise ValueError(f"unknown objective region {region!r}")


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
    raise ValueError(f"unknown family {family!r}")


def scipy_spec(family: str, params: dict) -> dict:
    """Explicit SciPy constructor arguments, for downstream samplers."""
    if family == "Normal":
        return {"scipy_dist": "norm", "loc": params["mu"], "scale": params["sigma"]}
    if family == "Lognormal":
        return {"scipy_dist": "lognorm", "s": params["sigma_log"], "loc": 0.0,
                "scale": float(np.exp(params["mu_log"]))}
    if family == "Gamma":
        return {"scipy_dist": "gamma", "a": params["shape_k"], "loc": 0.0, "scale": params["scale_theta"]}
    raise ValueError(f"unknown family {family!r}")


def _to_natural(family: str, t: np.ndarray) -> dict:
    """Unconstrained -> natural parameters, enforcing positivity of the scale/shape."""
    if family == "Normal":
        return {"mu": float(t[0]), "sigma": float(np.exp(t[1]))}
    if family == "Lognormal":
        return {"mu_log": float(t[0]), "sigma_log": float(np.exp(t[1]))}
    if family == "Gamma":
        return {"shape_k": float(np.exp(t[0])), "scale_theta": float(np.exp(t[1]))}
    raise ValueError(f"unknown family {family!r}")


def _from_natural(family: str, params: dict) -> np.ndarray:
    if family == "Normal":
        return np.array([params["mu"], np.log(params["sigma"])])
    if family == "Lognormal":
        return np.array([params["mu_log"], np.log(params["sigma_log"])])
    if family == "Gamma":
        return np.array([np.log(params["shape_k"]), np.log(params["scale_theta"])])
    raise ValueError(f"unknown family {family!r}")


def defective_starts(family: str, values_mgL, censored) -> list[dict]:
    """Extra LS starting values that guard the censored-mass-dominated replicates.

    A resampled heavy-censoring district can push the moment start into a region
    where the CDF least-squares surface is nearly flat, and the censored-MLE start
    may fail outright. These deterministic alternatives (a hard-coded small shape,
    and a fit to the detected values inflated to cover the censored mass) give the
    optimizer a second route. On the ten real districts the multi-start optimum is
    unchanged; they matter only under resampling.
    """
    x = np.asarray(values_mgL, dtype=float)
    cens = np.asarray(censored, dtype=bool)
    det = x[~cens]
    starts: list[dict] = []
    if det.size >= 2 and det.mean() > 0:
        mean, sd = float(det.mean()), float(det.std(ddof=1))
        if family == "Normal":
            for inflation in (1.0, 2.0):
                starts.append({"mu": mean, "sigma": max(sd, 1e-9) * inflation})
        elif family == "Lognormal":
            logs = det[det > 0]
            if logs.size >= 2:
                lo, ls = float(np.log(logs).mean()), max(float(np.log(logs).std(ddof=1)), 1e-6)
                for inflation in (1.0, 2.0):
                    starts.append({"mu_log": lo, "sigma_log": ls * inflation})
        elif family == "Gamma":
            for scale in (sd ** 2 / mean, mean / 4.0, mean):
                if scale > 0:
                    starts.append({"shape_k": max(mean / scale, 1e-6), "scale_theta": scale})
    return starts


def _plain(params: dict) -> dict:
    return {k: float(v) for k, v in params.items()}


# ---------------------------------------------------------------------------
# Starting values
# ---------------------------------------------------------------------------

def moment_start(family: str, detected_mgL) -> dict:
    """Method-of-moments start using *detected values only* (never the bounds)."""
    det = np.asarray(detected_mgL, dtype=float)
    if det.size < 2:
        raise ValueError("at least two detected values are required for a moment start")
    mean = float(det.mean())
    sd = float(det.std(ddof=1))
    if sd <= 0 or mean <= 0:
        raise ValueError("detected values are degenerate: moment start undefined")
    if family == "Normal":
        return {"mu": mean, "sigma": sd}
    if family == "Lognormal":
        logs = np.log(det[det > 0])
        if logs.size < 2:
            raise ValueError("too few positive detected values for a lognormal moment start")
        return {"mu_log": float(logs.mean()), "sigma_log": float(max(logs.std(ddof=1), 1e-6))}
    if family == "Gamma":
        return {"shape_k": mean ** 2 / sd ** 2, "scale_theta": sd ** 2 / mean}
    raise ValueError(f"unknown family {family!r}")


def censored_data(values_mgL, censored) -> stats.CensoredData:
    """SciPy left-censored data object: detected values exact, censored rows bounded."""
    x = np.asarray(values_mgL, dtype=float)
    cens = np.asarray(censored, dtype=bool)
    if x.shape != cens.shape or x.size == 0:
        raise ValueError("values_mgL and censored must be non-empty arrays of equal shape")
    if not np.any(cens):
        raise ValueError("no censored values: censored maximum likelihood is not applicable")
    return stats.CensoredData.left_censored(x, cens)


def _safe_mle_start(family: str, cdata: stats.CensoredData) -> dict | None:
    """Censored-MLE start; None when SciPy's fitter fails on this sample."""
    try:
        if family == "Normal":
            mu, sigma = stats.norm.fit(cdata)
            return {"mu": float(mu), "sigma": float(sigma)}
        if family == "Lognormal":
            s, loc, scale = stats.lognorm.fit(cdata, floc=0)
            if loc != 0:
                return None
            return {"mu_log": float(np.log(scale)), "sigma_log": float(s)}
        if family == "Gamma":
            a, loc, scale = stats.gamma.fit(cdata, floc=0)
            if loc != 0:
                return None
            return {"shape_k": float(a), "scale_theta": float(scale)}
    except Exception:  # SciPy raises a variety of errors on degenerate samples
        return None
    raise ValueError(f"unknown family {family!r}")


# ---------------------------------------------------------------------------
# Step C - censored maximum likelihood
# ---------------------------------------------------------------------------

@dataclass
class MLEFit:
    family: str
    params: dict
    loglikelihood: float
    n_params: int
    n: int
    success: bool
    message: str
    aic: float = float("nan")
    aicc: float = float("nan")


def censored_loglikelihood(family: str, params: dict, values_mgL, censored) -> float:
    """log L = sum_detected log f(x_i) + sum_censored log F(L_i). No substitution."""
    x = np.asarray(values_mgL, dtype=float)
    cens = np.asarray(censored, dtype=bool)
    dist = frozen(family, params)
    with np.errstate(divide="ignore", invalid="ignore"):
        ll = float(np.sum(dist.logpdf(x[~cens])) + np.sum(np.log(dist.cdf(x[cens]))))
    if not np.isfinite(ll):
        raise ValueError("censored log-likelihood is not finite for these parameters")
    return ll


def censored_aic(ll: float, k: int, n: int) -> tuple[float, float]:
    """AIC = 2k - 2l, AICc = AIC + 2k(k+1)/(n-k-1)."""
    aic = 2.0 * k - 2.0 * ll
    denom = n - k - 1
    aicc = aic + (2.0 * k * (k + 1) / denom if denom > 0 else float("nan"))
    return float(aic), float(aicc)


def fit_censored_mle(values_mgL, censored, family: str) -> MLEFit:
    """Robustness fit by censored maximum likelihood, via ``scipy.stats`` fitters."""
    x = np.asarray(values_mgL, dtype=float)
    cens = np.asarray(censored, dtype=bool)
    cdata = censored_data(x, cens)  # raises when there is nothing censored
    try:
        if family == "Normal":
            mu, sigma = stats.norm.fit(cdata)
            params = {"mu": float(mu), "sigma": float(sigma)}
        elif family == "Lognormal":
            s, loc, scale = stats.lognorm.fit(cdata, floc=0)
            params = {"mu_log": float(np.log(scale)), "sigma_log": float(s)}
        elif family == "Gamma":
            a, loc, scale = stats.gamma.fit(cdata, floc=0)
            params = {"shape_k": float(a), "scale_theta": float(scale)}
        else:
            raise ValueError(f"unknown family {family!r}")
    except ValueError:
        raise
    except Exception as exc:  # convergence failure inside SciPy's fitter
        return MLEFit(family, {k: float("nan") for k in PARAM_NAMES[family]}, float("nan"),
                      len(PARAM_NAMES[family]), int(x.size), False, f"scipy fit failed: {exc}")

    try:
        ll = censored_loglikelihood(family, params, x, cens)
        success, message = True, "scipy.stats censored fit converged"
    except ValueError as exc:
        ll, success, message = float("nan"), False, str(exc)
    aic, aicc = censored_aic(ll, len(PARAM_NAMES[family]), int(x.size)) if np.isfinite(ll) else (np.nan, np.nan)
    return MLEFit(family, _plain(params), ll, len(PARAM_NAMES[family]), int(x.size), success, message, aic, aicc)


# ---------------------------------------------------------------------------
# Step B - CDF least squares
# ---------------------------------------------------------------------------

@dataclass
class LSFit:
    family: str
    params: dict
    success: bool
    message: str
    sse: float
    n_points: int
    n_starts: int
    start_labels: list[str] = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    def dist(self):
        return frozen(self.family, self.params)


def ls_residuals(t, family: str, x_points: np.ndarray, F_target: np.ndarray) -> np.ndarray:
    """F(x_j; theta) - F_hat_RKM(x_j): the unweighted CDF least-squares residual."""
    return frozen(family, _to_natural(family, t)).cdf(x_points) - F_target


def cdf_rmse_from_sse(sse: float, n_points: int) -> float:
    return float(np.sqrt(sse / n_points)) if n_points else float("nan")


def fit_cdf_least_squares(values_mgL, censored, family: str, curve: RKMCurve | None = None,
                          region: str | None = None, mle_start: dict | None = None,
                          options: dict | None = None, extra_starts: list[dict] | None = None) -> LSFit:
    """Step B: multi-start CDF least squares against the reverse-KM curve.

    Starts are (a) method of moments on the detected values and (b) the censored-MLE
    solution; ``extra_starts`` appends more (the bootstrap uses them to harden the
    high-censoring replicates). The lowest-SSE result wins.
    """
    x = np.asarray(values_mgL, dtype=float)
    cens = np.asarray(censored, dtype=bool)
    curve = reverse_km_cdf(x, cens) if curve is None else curve
    mask = objective_mask(curve, region)
    if mask.sum() < 2:
        raise ValueError(f"objective region {region or OBJECTIVE_REGION!r} leaves "
                         f"{int(mask.sum())} point(s): at least 2 are required")
    xj, Fj = curve.x[mask], curve.F[mask]

    detected = x[~cens]
    starts: list[tuple[str, dict]] = []
    try:
        starts.append(("moments_detected", moment_start(family, detected)))
    except ValueError:
        pass
    if mle_start is None:
        try:
            cdata = censored_data(x, cens)
            mle_start = _safe_mle_start(family, cdata)
        except ValueError:
            mle_start = None
    if mle_start is not None:
        starts.append(("censored_mle", mle_start))
    for i, extra in enumerate(extra_starts or []):
        starts.append((f"extra_{i}", extra))
    if not starts:
        raise ValueError("no usable starting values for the CDF least-squares fit")

    options = LS_OPTIONS if options is None else options
    best = None
    for label, start in starts:
        try:
            res = optimize.least_squares(ls_residuals, _from_natural(family, start),
                                         args=(family, xj, Fj), **options)
        except (ValueError, FloatingPointError, RuntimeError):
            continue
        if not np.all(np.isfinite(res.x)):
            continue
        params = _to_natural(family, res.x)
        if not np.all(np.isfinite(list(params.values()))):
            continue
        sse = float(np.sum(res.fun ** 2))
        if best is None or sse < best[0]:
            best = (sse, params, res, label)

    if best is None:
        return LSFit(family, {k: float("nan") for k in PARAM_NAMES[family]}, False,
                     "every start failed", float("inf"), int(mask.sum()), len(starts),
                     [s[0] for s in starts])
    sse, params, res, label = best
    return LSFit(
        family=family,
        params=_plain(params),
        success=bool(res.success),
        message=str(res.message),
        sse=sse,
        n_points=int(mask.sum()),
        n_starts=len(starts),
        start_labels=[s[0] for s in starts],
        extra={"best_start": label, "nfev": int(res.nfev),
               "objective_region_points": int(mask.sum()),
               "objective_region_total_points": int(curve.x.size)},
    )


# ---------------------------------------------------------------------------
# Fit metrics and Step E - admissibility and evaluation
# ---------------------------------------------------------------------------

def cdf_sse(dist, x_points, F_target) -> float:
    return float(np.sum((dist.cdf(x_points) - F_target) ** 2))


def cdf_max_discrepancy(dist, x_points, F_target) -> float:
    """Kolmogorov-Smirnov style distance on the reverse-KM points."""
    return float(np.max(np.abs(dist.cdf(x_points) - F_target)))


def empirical_quantiles(curve: RKMCurve, probs) -> np.ndarray:
    """Reverse-KM quantiles; points below the plateau are unresolved, so the curve
    is read as a left-continuous step and clipped at the smallest detected value."""
    F = curve.F
    q = np.interp(np.asarray(probs, dtype=float), F, curve.x)
    return np.maximum(q, curve.smallest_detected * (F[0] > 0))


def scale_parameters(family: str, params: dict) -> dict:
    """The strictly-positive scale-like parameters of a family.

    Only these are band-checked by the magnitude gate; ``Gamma.shape_k`` is exempt
    because ``k < 1`` is legal (an unbounded density at zero) and is exactly what a
    heavily censored district produces.
    """
    if family == "Normal":
        return {"sigma": params["sigma"]}
    if family == "Lognormal":
        return {"sigma_log": params["sigma_log"], "scale=exp(mu_log)": float(np.exp(params["mu_log"]))}
    if family == "Gamma":
        return {"scale_theta": params["scale_theta"]}
    raise ValueError(f"unknown family {family!r}")


def admissibility_reasons(family: str, ls_fit: LSFit, mle_fit: MLEFit | None,
                          params_for_support: dict | None = None) -> list[str]:
    """Reasons a candidate is inadmissible; empty list means admissible."""
    reasons: list[str] = []
    params = ls_fit.params if params_for_support is None else params_for_support
    if not ls_fit.success:
        reasons.append(f"LS optimizer did not report convergence ({ls_fit.message})")
    values = np.array([params[k] for k in PARAM_NAMES[family]], dtype=float)
    if not np.all(np.isfinite(values)):
        reasons.append("fitted parameters are not finite")
    else:
        scales = scale_parameters(family, params)
        if any(v <= 0 for v in scales.values()):
            reasons.append("fitted scale parameters must be strictly positive: "
                           + ", ".join(f"{k}={v:g}" for k, v in scales.items() if v <= 0))
        bad = {k: v for k, v in scales.items()
               if not (PARAM_MIN_MAGNITUDE <= abs(v) <= PARAM_MAX_MAGNITUDE)}
        if bad:
            reasons.append(f"scale parameters outside the sane magnitude band "
                           f"[{PARAM_MIN_MAGNITUDE:g}, {PARAM_MAX_MAGNITUDE:g}]: {bad}")
    if family in ZERO_LOCATION_FAMILIES:
        lo, _ = frozen(family, params).support()
        if lo != 0.0:
            reasons.append(f"{family} support starts at {lo!r}, not 0")
    if family == "Normal":
        p_nonpositive = float(frozen(family, params).cdf(0.0))
        if p_nonpositive >= NORMAL_MAX_NONPOSITIVE_PROB:
            reasons.append(f"P(C<=0)={p_nonpositive:.3e} >= {NORMAL_MAX_NONPOSITIVE_PROB:.0e}")
    return reasons


def evaluate_district_fit(district: str, values_mgL, censored, family: str,
                          curve: RKMCurve | None = None, region: str | None = None,
                          mle_start: dict | None = None, options: dict | None = None,
                          extra_starts: list[dict] | None = None) -> dict:
    """Fit one candidate family both ways for one district and evaluate it (Steps B, C, E).

    Returns a flat row: LS parameters, SSE/RMSE/max discrepancy, tail errors,
    censored-MLE parameters, log-likelihood, AIC/AICc, and the admissibility verdict.
    """
    if family not in FAMILIES:
        raise ValueError(f"unknown family {family!r}")
    x = np.asarray(values_mgL, dtype=float)
    cens = np.asarray(censored, dtype=bool)
    curve = reverse_km_cdf(x, cens) if curve is None else curve
    return _evaluate_row(family, x, cens, curve, region, mle_start, options, extra_starts)


def _evaluate_row(family: str, x: np.ndarray, cens: np.ndarray, curve: RKMCurve,
                  region: str | None, mle_start: dict | None, options: dict | None,
                  extra_starts: list[dict] | None) -> dict:
    """Build one district x family evaluation row."""
    mask = objective_mask(curve, region)
    xj, Fj = curve.x[mask], curve.F[mask]
    ls = fit_cdf_least_squares(x, cens, family, curve=curve, region=region,
                               mle_start=mle_start, options=options, extra_starts=extra_starts)
    mle = fit_censored_mle(x, cens, family)
    values = np.array([ls.params[k] for k in PARAM_NAMES[family]], dtype=float)

    if np.all(np.isfinite(values)):
        dist = ls.dist()
        lo, _ = dist.support()
        zero_loc = bool(lo == 0.0)
        row: dict = {"distribution": family}
        for name in PARAM_NAMES[family]:
            row[f"ls_{name}"] = ls.params[name]
        row["ls_converged"] = ls.success
        row["ls_message"] = ls.message
        row["ls_n_starts"] = ls.n_starts
        row["ls_best_start"] = ls.extra.get("best_start")
        row["objective_region"] = region or OBJECTIVE_REGION
        row["objective_points"] = int(mask.sum())
        row["objective_points_total"] = int(curve.x.size)

        row["cdf_sse"] = ls.sse
        row["cdf_rmse"] = cdf_rmse_from_sse(ls.sse, int(mask.sum()))
        row["cdf_max_discrepancy"] = cdf_max_discrepancy(dist, xj, Fj)
        row["cdf_max_discrepancy_all_points"] = cdf_max_discrepancy(dist, curve.x, curve.F)
        row["support_lower"] = float(lo)
        row["support_starts_at_zero"] = zero_loc
        row["prob_arsenic_nonpositive"] = float(dist.cdf(0.0))
    else:
        row = {"distribution": family}
        for name in PARAM_NAMES[family]:
            row[f"ls_{name}"] = np.nan
        row.update(ls_converged=ls.success, ls_message=ls.message, ls_n_starts=ls.n_starts,
                   ls_best_start=ls.extra.get("best_start"), objective_region=region or OBJECTIVE_REGION,
                   objective_points=int(mask.sum()), objective_points_total=int(curve.x.size),
                   cdf_sse=np.nan, cdf_rmse=np.nan, cdf_max_discrepancy=np.nan,
                   cdf_max_discrepancy_all_points=np.nan, support_lower=np.nan,
                   support_starts_at_zero=False, prob_arsenic_nonpositive=np.nan)
        zero_loc = False

    for name in PARAM_NAMES[family]:
        row[f"mle_{name}"] = mle.params[name]
    row["mle_converged"] = mle.success
    row["mle_message"] = mle.message
    row["mle_loglikelihood"] = mle.loglikelihood
    row["mle_k"] = mle.n_params
    row["aic"] = mle.aic
    row["aicc"] = mle.aicc

    if dist is not None and mle.success:
        try:
            recomputed = censored_loglikelihood(family, mle.params, x, cens)
            row["mle_censored_loglikelihood"] = recomputed
            row["mle_loglik_matches_recompute"] = bool(np.isclose(recomputed, mle.loglikelihood, rtol=1e-9))
        except ValueError:
            row["mle_censored_loglikelihood"] = np.nan
            row["mle_loglik_matches_recompute"] = False
        rel = [abs(ls.params[k] - mle.params[k]) / abs(mle.params[k])
               for k in PARAM_NAMES[family] if mle.params[k] != 0]
        row["ls_mle_max_rel_param_diff"] = float(np.max(rel)) if rel else np.nan
    else:
        row["mle_censored_loglikelihood"] = np.nan
        row["mle_loglik_matches_recompute"] = False
        row["ls_mle_max_rel_param_diff"] = np.nan

    # Tail agreement against the reverse-KM curve.
    emp_q = empirical_quantiles(curve, REPORT_PROBS)
    fit_q = np.asarray(dist.ppf(REPORT_PROBS), dtype=float) if dist is not None else np.full(len(REPORT_PROBS), np.nan)
    for prob, e, f in zip(REPORT_PROBS, emp_q, fit_q):
        tag = f"P{int(round(prob * 100))}"
        row[f"empirical_{tag}_mgL"] = float(e)
        row[f"fitted_{tag}_mgL"] = float(f)
    if dist is not None:
        idx = [REPORT_PROBS.index(p) for p in TAIL_PROBS]
        with np.errstate(divide="ignore", invalid="ignore"):
            rel_err = np.abs(fit_q[idx] - emp_q[idx]) / emp_q[idx]
        row["tail_max_rel_quantile_error"] = float(np.nanmax(rel_err)) if np.any(np.isfinite(rel_err)) else np.nan
    else:
        row["tail_max_rel_quantile_error"] = np.nan

    reasons = admissibility_reasons(family, ls, mle)
    row["admissible"] = not reasons
    row["inadmissible_reason"] = "; ".join(reasons)
    row["zero_location_enforced"] = bool(zero_loc and family in ZERO_LOCATION_FAMILIES)
    return row


def fit_district(district: str, values_mgL, censored, region: str | None = None,
                 families=FAMILIES, options: dict | None = None) -> pd.DataFrame:
    """Fit all candidate families for one district and rank them by CDF RMSE."""
    x = np.asarray(values_mgL, dtype=float)
    cens = np.asarray(censored, dtype=bool)
    curve = reverse_km_cdf(x, cens)
    rows = []
    for family in families:
        row = _evaluate_row(family, x, cens, curve, region, None, options, None)
        rows.append({"DISTRICT": district, "n": int(x.size),
                     "n_detected": curve.n_detected, "n_censored": curve.n_censored,
                     "censoring_percent": round(100.0 * curve.n_censored / x.size, 2),
                     "rkm_points": int(curve.x.size),
                     "rkm_censoring_fraction": curve.censoring_fraction,
                     "censoring_above_detection": curve.censoring_above_detection, **row})
    table = pd.DataFrame(rows)
    table["rmse_rank"] = table["cdf_rmse"].rank(method="min").astype(int)
    table["aicc_rank"] = table["aicc"].rank(method="min").astype(int)
    table["admissible_rmse_rank"] = table.loc[table["admissible"], "cdf_rmse"].rank(method="min")
    return table


# ---------------------------------------------------------------------------
# Step E - selection
# ---------------------------------------------------------------------------

def high_censoring_note(row: pd.Series | dict, warning_percent: float = 50.0,
                        min_detected_mgL: float | None = None,
                        plateau_note: str | None = None) -> str:
    """Caveat text for districts whose censoring burden makes the low tail unresolved.

    ``row`` is a fit-table row (it carries ``censoring_percent`` and
    ``rkm_censoring_fraction``). ``min_detected_mgL`` overrides the smallest detected
    concentration when it is known from the preprocessing summary, and
    ``plateau_note`` supplies the resolved/unresolved detail from ``censoring_plateau_note``.
    """
    get = row.get
    pct = float(get("censoring_percent"))
    if pct < 30.0:
        return ""
    frac = float(get("rkm_censoring_fraction"))
    smallest = min_detected_mgL
    if smallest is None:
        smallest = float(get("minimum_detected_As_mgL")) if get("minimum_detected_As_mgL") is not None else float("nan")
    severity = "HIGH-CENSORING DISTRICT" if pct >= warning_percent else "notable censoring"
    detail = f" ({plateau_note})" if plateau_note else ""
    return (f"{severity}: {pct:.2f}% of wells are left-censored and the reverse-KM curve is flat at "
            f"F~={frac:.3f} for every concentration below the smallest detected value "
            f"({smallest:.4g} mg/L){detail}. Nothing in that region constrains the fitted density and it "
            f"is excluded from the objective (decision 9.1 V2), so this fit is inferred from the "
            f"resolvable region only and is less certain than a low-censoring district's.")


def select_arsenic_distribution(fit_table: pd.DataFrame, warning_percent: float = 50.0,
                                min_detected_mgL: float | None = None,
                                curve: RKMCurve | None = None) -> tuple[pd.Series, str]:
    """Predeclared rule: lowest CDF RMSE among admissible candidates.

    AICc and the tail-quantile error are reported as corroboration or dissent, never
    as a tie-break. The outcome is provisional until a reviewer signs it off.
    """
    district = fit_table["DISTRICT"].iloc[0] if "DISTRICT" in fit_table else "district"
    admissible = fit_table[fit_table["admissible"]]
    if admissible.empty:
        raise RuntimeError(f"no admissible arsenic candidate for {district}")
    ranked = admissible.sort_values("cdf_rmse", kind="stable")
    chosen = ranked.iloc[0]
    runner = ranked.iloc[1] if len(ranked) > 1 else None

    notes = [f"rule [{SELECTION_RULE}]: lowest CDF RMSE among admissible candidates "
             f"({chosen['cdf_rmse']:.4f} on {int(chosen['objective_points'])} of "
             f"{int(chosen['objective_points_total'])} reverse-KM points, region {chosen['objective_region']})"]

    rejected = fit_table[~fit_table["admissible"]]
    for row in rejected.itertuples():
        if row.distribution == "Normal":
            notes.append(f"Normal was rejected by the admissibility gate "
                         f"(P(C<=0)={row.prob_arsenic_nonpositive:.3f}) despite a lower raw SSE; "
                         f"the realistic contest is Lognormal vs Gamma")
        else:
            notes.append(f"{row.distribution} inadmissible: {row.inadmissible_reason}")

    # Robustness corroboration.
    best_aicc = admissible.sort_values("aicc", kind="stable").iloc[0]["distribution"]
    if best_aicc == chosen["distribution"]:
        notes.append("censored-MLE AICc corroborates the LS winner")
    else:
        delta = float(admissible.loc[admissible["distribution"] == best_aicc, "aicc"].iloc[0] - chosen["aicc"])
        notes.append(f"DISSENT: censored-MLE AICc favours {best_aicc} (delta AICc = {delta:+.2f} "
                     f"against the LS winner); recorded, not used as a tie-break")
    best_tail = admissible.sort_values("tail_max_rel_quantile_error", kind="stable").iloc[0]["distribution"]
    notes.append("tail quantile agreement also favours the selected family" if best_tail == chosen["distribution"]
                 else f"tail quantile agreement favours {best_tail}")
    best_ks = admissible.sort_values("cdf_max_discrepancy", kind="stable").iloc[0]["distribution"]
    notes.append("smallest max CDF discrepancy too" if best_ks == chosen["distribution"]
                 else f"max CDF discrepancy favours {best_ks}")

    if runner is not None:
        ratio = runner["cdf_rmse"] / chosen["cdf_rmse"] if chosen["cdf_rmse"] > 0 else np.inf
        notes.append(f"runner-up {runner['distribution']} RMSE ratio {ratio:.2f}")
    if np.isfinite(chosen["ls_mle_max_rel_param_diff"]):
        notes.append(f"LS vs censored-MLE max relative parameter difference "
                     f"{chosen['ls_mle_max_rel_param_diff']:.1%}")

    if chosen["distribution"] == "Lognormal":
        notes.append("parameterization note: mu_log is the mean of log(C) (negative here because "
                     "arsenic is measured in mg/L) and sigma_log is the log-scale standard deviation; "
                     "neither is a support bound")
    note = high_censoring_note(chosen, warning_percent, min_detected_mgL,
                               plateau_note=censoring_plateau_note(curve) if curve is not None else None)
    if note:
        notes.append(note)
    return chosen, "; ".join(notes)


def selected_record(chosen: pd.Series, rationale: str) -> dict:
    """Machine-readable record for downstream (Monte Carlo) consumers."""
    family = chosen["distribution"]
    params = {name: float(chosen[f"ls_{name}"]) for name in PARAM_NAMES[family]}
    return {
        "district": chosen["DISTRICT"],
        "distribution": family,
        "parameterization": "natural parameters from censoring-aware CDF least squares; "
                            "Lognormal/Gamma forced to loc=0",
        "params": params,
        "scipy": {k: (float(v) if not isinstance(v, str) else v) for k, v in scipy_spec(family, params).items()},
        "n": int(chosen["n"]),
        "n_detected": int(chosen["n_detected"]),
        "n_censored": int(chosen["n_censored"]),
        "censoring_percent": float(chosen["censoring_percent"]),
        "rkm_censoring_fraction": float(chosen["rkm_censoring_fraction"]),
        "cdf_sse": float(chosen["cdf_sse"]),
        "cdf_rmse": float(chosen["cdf_rmse"]),
        "cdf_max_discrepancy": float(chosen["cdf_max_discrepancy"]),
        "aicc": float(chosen["aicc"]),
        "objective_region": str(chosen["objective_region"]),
        "objective_points": int(chosen["objective_points"]),
        "status": SELECTION_STATUS,
        "rationale": rationale,
    }


def selection_table(selected: list[dict]) -> pd.DataFrame:
    """The selected-record list flattened into the project's CSV convention."""
    rows = []
    for s in selected:
        rows.append({
            "DISTRICT": s["district"],
            "distribution": s["distribution"],
            "parameterization": s["parameterization"],
            **{f"param_{k}": v for k, v in s["params"].items()},
            **{f"scipy_{k}": v for k, v in s["scipy"].items()},
            "n": s["n"], "n_detected": s["n_detected"], "n_censored": s["n_censored"],
            "censoring_percent": s["censoring_percent"],
            "rkm_censoring_fraction": s["rkm_censoring_fraction"],
            "cdf_sse": s["cdf_sse"], "cdf_rmse": s["cdf_rmse"],
            "cdf_max_discrepancy": s["cdf_max_discrepancy"], "aicc": s["aicc"],
            "objective_region": s["objective_region"], "objective_points": s["objective_points"],
            "status": s["status"], "rationale": s["rationale"],
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Step F - sampling interface
# ---------------------------------------------------------------------------

def sample_arsenic(selected: dict, N: int | None = None, rng: np.random.Generator | None = None,
                   u: np.ndarray | None = None) -> np.ndarray:
    """Draw arsenic concentrations (mg/L) by explicit inverse-CDF sampling.

    ``u = rng.uniform(0, 1, size=N)`` then ``dist.ppf(u)``. Pass either ``u``
    directly or ``N`` together with a ``numpy.random.Generator``. The result is
    finite and >= 0: the Normal branch clips the (inadmissible but supported)
    negative tail to zero rather than returning a negative concentration.
    """
    if u is None:
        if N is None or rng is None:
            raise ValueError("provide u, or both N and rng")
        u = rng.uniform(0.0, 1.0, size=N)
    u = np.asarray(u, dtype=float)
    if not np.all(np.isfinite(u)):
        raise ValueError("u must be finite")
    if np.any((u < 0.0) | (u > 1.0)):
        raise ValueError("u must lie in [0, 1]")
    family = selected["distribution"]
    values = np.asarray(frozen(family, selected["params"]).ppf(u), dtype=float)
    if family == "Normal":
        values = np.maximum(values, 0.0)
    if not np.all(np.isfinite(values)):
        raise ValueError("sampling produced non-finite arsenic concentrations")
    return np.maximum(values, 0.0)
