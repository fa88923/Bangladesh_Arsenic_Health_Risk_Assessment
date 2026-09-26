"""Phase 6: adult and child risk-model parameter sets.

Reads ``risk_model`` from config/run_config.json and separates, for each
population, the fixed constants (ED, AT, Kp, CF, RfDs, CSF) from the sampled
inputs (C, IR, BW, EF, ET, SA). C and BW come from the Phase 4 and Phase 5
selected fits; this module supplies the base-paper inputs IR, EF, ET, and SA.

Base-paper "mean +/- SD" Lognormal values (Yadav & Kalkal 2024, Table 2) are
read as arithmetic mean and SD by default (decision 2026-09-25) and converted to
log-space with

    sigma_log^2 = ln(1 + sd^2 / mean^2),   mu_log = ln(mean) - sigma_log^2 / 2.

Setting ``plus_minus_interpretation`` to ``log_space`` instead uses the reported
numbers directly as (mu_log, sigma_log), the pre-2026-09-25 reading.

Parameterizations (SciPy):

- Lognormal(mu_log, sigma_log)       -> lognorm(s=sigma_log, loc=0, scale=exp(mu_log))
- Triangular(min, mode, max)         -> triang(c=(mode-min)/(max-min), loc=min, scale=max-min)
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats

from arsenic_hra.paths import load_run_config

POPULATIONS = ("adult", "child")

# Sampled inputs, in the order used for iteration storage and sensitivity analysis.
STOCHASTIC_INPUTS = ("C_mg_per_L", "IR_L_per_day", "BW_kg", "EF_days_per_year", "ET_h_per_day", "SA_m2")
# Inputs taken from this module's config section; C and BW come from Phases 4 and 5.
BASE_PAPER_STOCHASTIC = ("IR_L_per_day", "EF_days_per_year", "ET_h_per_day", "SA_m2")
FIXED_INPUTS = ("ED_years", "AT_noncancer_days", "AT_cancer_days", "Kp_cm_per_h", "CF_L_per_cm_m2",
                "RfD_ing_mg_per_kg_day", "RfD_dermal_mg_per_kg_day", "CSF_per_mg_per_kg_day")

UNITS = {
    "C_mg_per_L": "mg/L",
    "IR_L_per_day": "L/day",
    "BW_kg": "kg",
    "EF_days_per_year": "days/year",
    "ET_h_per_day": "h/day",
    "SA_m2": "m^2",
    "ED_years": "years",
    "AT_noncancer_days": "days",
    "AT_cancer_days": "days",
    "Kp_cm_per_h": "cm/h",
    "CF_L_per_cm_m2": "L/(cm*m^2)  [= L*m/(m^3*cm) as printed in the base paper]",
    "RfD_ing_mg_per_kg_day": "mg/(kg*day)",
    "RfD_dermal_mg_per_kg_day": "mg/(kg*day)",
    "CSF_per_mg_per_kg_day": "(mg/(kg*day))^-1",
    "ADD_ing_mg_per_kg_day": "mg/(kg*day)",
    "ADD_dermal_mg_per_kg_day": "mg/(kg*day)",
    "ADD_ing_cancer_mg_per_kg_day": "mg/(kg*day)",
    "HQ_ing": "dimensionless",
    "HQ_dermal": "dimensionless",
    "HI": "dimensionless",
    "ELCR": "dimensionless (excess lifetime probability)",
}


def lognormal_from_arithmetic(mean: float, sd: float) -> tuple[float, float]:
    """Return (mu_log, sigma_log) of the Lognormal with arithmetic ``mean`` and ``sd``."""
    if not (mean > 0 and sd > 0):
        raise ValueError("Lognormal arithmetic mean and SD must be positive")
    sigma2 = np.log1p((sd / mean) ** 2)
    return float(np.log(mean) - sigma2 / 2.0), float(np.sqrt(sigma2))


@dataclass(frozen=True)
class InputDistribution:
    """One sampled base-paper input, with explicit natural and SciPy parameters."""
    name: str
    family: str
    params: dict
    scipy: dict
    reported: dict
    interpretation: str

    def frozen(self):
        spec = {k: v for k, v in self.scipy.items() if k != "scipy_dist"}
        return getattr(stats, self.scipy["scipy_dist"])(**spec)

    def ppf(self, u):
        return self.frozen().ppf(u)


def input_distribution(name: str, spec: dict, interpretation: str) -> InputDistribution:
    family = spec["family"]
    reported = {k: v for k, v in spec.items() if k not in ("family", "source")}
    if family == "Lognormal":
        if interpretation == "arithmetic_mean_sd":
            mu_log, sigma_log = lognormal_from_arithmetic(spec["reported_mean"], spec["reported_sd"])
        elif interpretation == "log_space":
            mu_log, sigma_log = float(spec["reported_mean"]), float(spec["reported_sd"])
        else:
            raise ValueError(f"unknown plus_minus_interpretation {interpretation!r}")
        params = {"mu_log": mu_log, "sigma_log": sigma_log}
        scipy_args = {"scipy_dist": "lognorm", "s": sigma_log, "loc": 0.0, "scale": float(np.exp(mu_log))}
    elif family == "Triangular":
        a, m, b = float(spec["min"]), float(spec["mode"]), float(spec["max"])
        if not (0 < a <= m <= b and a < b):
            raise ValueError(f"{name}: Triangular needs 0 < min <= mode <= max, min < max")
        params = {"min": a, "mode": m, "max": b}
        scipy_args = {"scipy_dist": "triang", "c": (m - a) / (b - a), "loc": a, "scale": b - a}
        interpretation = "min/mode/max"
    else:
        raise ValueError(f"{name}: unsupported family {family!r}")
    return InputDistribution(name, family, params, scipy_args, reported, interpretation)


@dataclass(frozen=True)
class PopulationParameters:
    """Everything the risk equations need for one population, except C and BW."""
    population: str
    fixed: dict
    deterministic: dict
    stochastic: dict = field(default_factory=dict)  # name -> InputDistribution
    cancer_averaging_time: str = "lifetime"
    plus_minus_interpretation: str = "arithmetic_mean_sd"


def population_parameters(population: str, cfg: dict | None = None) -> PopulationParameters:
    if population not in POPULATIONS:
        raise ValueError(f"unknown population {population!r}")
    rcfg = (cfg or load_run_config())["risk_model"]
    pcfg = rcfg["populations"][population]
    const = rcfg["constants"]
    ed = float(pcfg["ED_years"])
    days = float(const["days_per_year"])
    at_mode = rcfg["cancer_averaging_time"]
    if at_mode == "lifetime":
        at_cancer = float(const["lifetime_years"]) * days
    elif at_mode == "exposure_duration":
        at_cancer = ed * days
    else:
        raise ValueError(f"unknown cancer_averaging_time {at_mode!r}")
    fixed = {
        "ED_years": ed,
        "AT_noncancer_days": ed * days,
        "AT_cancer_days": at_cancer,
        "Kp_cm_per_h": float(const["Kp_cm_per_h"]),
        "CF_L_per_cm_m2": float(const["CF_L_per_cm_m2"]),
        "RfD_ing_mg_per_kg_day": float(const["RfD_ing_mg_per_kg_day"]),
        "RfD_dermal_mg_per_kg_day": float(const["RfD_dermal_mg_per_kg_day"]),
        "CSF_per_mg_per_kg_day": float(const["CSF_per_mg_per_kg_day"]),
    }
    interp = rcfg["plus_minus_interpretation"]
    stochastic = {name: input_distribution(name, spec, interp) for name, spec in pcfg["stochastic"].items()}
    if set(stochastic) != set(BASE_PAPER_STOCHASTIC):
        raise ValueError(f"{population}: stochastic inputs must be exactly {BASE_PAPER_STOCHASTIC}")
    return PopulationParameters(population, fixed, dict(pcfg["deterministic"]), stochastic, at_mode, interp)


def parameter_table(cfg: dict | None = None) -> list[dict]:
    """Flat rows describing every model input for both populations (for results/tables)."""
    rows = []
    for population in POPULATIONS:
        p = population_parameters(population, cfg)
        for name, value in p.fixed.items():
            rows.append({"population": population, "parameter": name, "role": "fixed", "unit": UNITS[name],
                         "value": value, "family": "Point", "parameterization": "", "deterministic_value": value})
        for name, dist in p.stochastic.items():
            d = dist.frozen()
            rows.append({
                "population": population, "parameter": name, "role": "sampled", "unit": UNITS[name],
                "value": "", "family": dist.family,
                "parameterization": "; ".join(f"{k}={v:.6g}" for k, v in dist.params.items())
                                    + f" [{dist.interpretation}]",
                "deterministic_value": p.deterministic[name],
                "reported": "; ".join(f"{k}={v}" for k, v in dist.reported.items()),
                "dist_mean": float(d.mean()), "dist_median": float(d.median()),
                "dist_P5": float(d.ppf(0.05)), "dist_P95": float(d.ppf(0.95)),
            })
        for name, family, source in (
                ("C_mg_per_L", "Phase 4 selected fit, per district",
                 "results/arsenic/arsenic_selected_distributions.json"),
                ("BW_kg", "Phase 5 selected fit", "results/bodyweight/bodyweight_selected_distributions.json")):
            rows.append({"population": population, "parameter": name, "role": "sampled", "unit": UNITS[name],
                         "value": "", "family": family, "parameterization": source, "deterministic_value": ""})
    return rows
