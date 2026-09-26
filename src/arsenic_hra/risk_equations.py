"""Phase 6: base-paper risk equations (Yadav & Kalkal 2024, equations 1-6).

Every function takes scalars or NumPy arrays (broadcast element-wise) and keeps
the full published form; nothing is algebraically simplified, so EF * ED and AT
appear explicitly even where they cancel.

    ADD_ing    = C * IR * EF * ED / (BW * AT)                          (1)  mg/(kg*day)
    ADD_dermal = C * Kp * EF * ED * ET * SA * CF / (BW * AT)           (2)  mg/(kg*day)
    HQ_ing     = ADD_ing / RfD_ing                                     (3)  dimensionless
    HQ_dermal  = ADD_dermal / RfD_dermal                               (4)  dimensionless
    HI         = HQ_ing + HQ_dermal                                    (5)  dimensionless
    ELCR       = ADD_ing,cancer * CSF                                  (6)  dimensionless

Units: C mg/L, IR L/day, EF days/year, ED years, BW kg, AT days, Kp cm/h,
ET h/day, SA m^2, CF L/(cm*m^2). CF = 10 because 1 cm * 1 m^2 = 0.01 m^3 = 10 L,
so C*Kp*ET*SA*CF is in mg/day.

AT: the non-cancer doses (1)-(2) use AT = ED * 365. ELCR uses ADD_ing,cancer,
equation (1) evaluated with the cancer averaging time: a 70-year lifetime
(70 * 365 days) under the 2026-09-25 decision, or ED * 365 as in the base paper
(set ``risk_model.cancer_averaging_time``). For adults (ED = 70) the two agree.
"""
from __future__ import annotations

import numpy as np

from arsenic_hra.risk_parameters import PopulationParameters

OUTPUTS = ("ADD_ing_mg_per_kg_day", "ADD_dermal_mg_per_kg_day", "HQ_ing", "HQ_dermal", "HI",
           "ADD_ing_cancer_mg_per_kg_day", "ELCR")

EQUATIONS = [
    {"output": "ADD_ing_mg_per_kg_day", "paper_eq": "(1)", "equation": "C * IR * EF * ED / (BW * AT_noncancer)",
     "inputs": "C [mg/L], IR [L/day], EF [days/year], ED [years], BW [kg], AT_noncancer = ED*365 [days]",
     "unit": "mg/(kg*day)"},
    {"output": "ADD_dermal_mg_per_kg_day", "paper_eq": "(2)",
     "equation": "C * Kp * EF * ED * ET * SA * CF / (BW * AT_noncancer)",
     "inputs": "C [mg/L], Kp [cm/h], EF [days/year], ED [years], ET [h/day], SA [m^2], CF [L/(cm*m^2)], "
               "BW [kg], AT_noncancer [days]",
     "unit": "mg/(kg*day)"},
    {"output": "HQ_ing", "paper_eq": "(3)", "equation": "ADD_ing / RfD_ing",
     "inputs": "ADD_ing [mg/(kg*day)], RfD_ing = 0.0003 [mg/(kg*day)]", "unit": "dimensionless"},
    {"output": "HQ_dermal", "paper_eq": "(4)", "equation": "ADD_dermal / RfD_dermal",
     "inputs": "ADD_dermal [mg/(kg*day)], RfD_dermal = 0.000285 [mg/(kg*day)]", "unit": "dimensionless"},
    {"output": "HI", "paper_eq": "(5)", "equation": "HQ_ing + HQ_dermal", "inputs": "HQ_ing, HQ_dermal",
     "unit": "dimensionless"},
    {"output": "ADD_ing_cancer_mg_per_kg_day", "paper_eq": "(1) with AT_cancer",
     "equation": "C * IR * EF * ED / (BW * AT_cancer)",
     "inputs": "as ADD_ing; AT_cancer = 70*365 [days] (lifetime) or ED*365 (base paper)", "unit": "mg/(kg*day)"},
    {"output": "ELCR", "paper_eq": "(6)", "equation": "ADD_ing_cancer * CSF",
     "inputs": "ADD_ing_cancer [mg/(kg*day)], CSF = 1.5 [(mg/(kg*day))^-1]",
     "unit": "dimensionless (excess lifetime probability)"},
]


def add_ingestion(C, IR, EF, ED, BW, AT):
    """Equation (1): ingestion average daily dose, mg/(kg*day)."""
    return (C * IR * EF * ED) / (BW * AT)


def add_dermal(C, Kp, EF, ED, ET, SA, CF, BW, AT):
    """Equation (2): dermal average daily dose, mg/(kg*day)."""
    return (C * Kp * EF * ED * ET * SA * CF) / (BW * AT)


def hazard_quotient(ADD, RfD):
    """Equations (3)-(4): hazard quotient, dimensionless."""
    return ADD / RfD


def hazard_index(HQ_ing, HQ_dermal):
    """Equation (5): hazard index, dimensionless."""
    return HQ_ing + HQ_dermal


def elcr(ADD_ing_cancer, CSF):
    """Equation (6): excess lifetime cancer risk, dimensionless."""
    return ADD_ing_cancer * CSF


def _validate_inputs(inputs: dict) -> dict:
    required = ("C_mg_per_L", "IR_L_per_day", "BW_kg", "EF_days_per_year", "ET_h_per_day", "SA_m2")
    missing = [k for k in required if k not in inputs]
    if missing:
        raise KeyError(f"missing risk inputs: {missing}")
    arrays = {k: np.asarray(inputs[k], dtype=float) for k in required}
    for name, a in arrays.items():
        if not np.all(np.isfinite(a)):
            raise ValueError(f"{name} contains non-finite values")
        if name == "C_mg_per_L":
            if np.any(a < 0):
                raise ValueError("C_mg_per_L must be >= 0")
        elif np.any(a <= 0):
            raise ValueError(f"{name} must be > 0")
    if np.any(arrays["EF_days_per_year"] > 366):
        raise ValueError("EF_days_per_year cannot exceed 366")
    if np.any(arrays["ET_h_per_day"] > 24):
        raise ValueError("ET_h_per_day cannot exceed 24")
    return arrays


def evaluate_risk(inputs: dict, params: PopulationParameters) -> dict:
    """Apply equations (1)-(6) to one population.

    ``inputs`` maps the six sampled-input names (see ``risk_parameters.STOCHASTIC_INPUTS``)
    to scalars or equal-length arrays; the fixed constants come from ``params``.
    Returns the seven outputs in ``OUTPUTS`` as float arrays.
    """
    x = _validate_inputs(inputs)
    f = params.fixed
    C, IR, BW = x["C_mg_per_L"], x["IR_L_per_day"], x["BW_kg"]
    EF, ET, SA = x["EF_days_per_year"], x["ET_h_per_day"], x["SA_m2"]
    ED = f["ED_years"]

    add_ing = add_ingestion(C, IR, EF, ED, BW, f["AT_noncancer_days"])
    add_derm = add_dermal(C, f["Kp_cm_per_h"], EF, ED, ET, SA, f["CF_L_per_cm_m2"], BW, f["AT_noncancer_days"])
    hq_ing = hazard_quotient(add_ing, f["RfD_ing_mg_per_kg_day"])
    hq_derm = hazard_quotient(add_derm, f["RfD_dermal_mg_per_kg_day"])
    hi = hazard_index(hq_ing, hq_derm)
    add_ing_cancer = add_ingestion(C, IR, EF, ED, BW, f["AT_cancer_days"])
    risk = elcr(add_ing_cancer, f["CSF_per_mg_per_kg_day"])

    out = dict(zip(OUTPUTS, (add_ing, add_derm, hq_ing, hq_derm, hi, add_ing_cancer, risk)))
    out = {k: np.asarray(v, dtype=float) for k, v in out.items()}
    for name, a in out.items():
        if not np.all(np.isfinite(a)) or np.any(a < 0):
            raise FloatingPointError(f"{name} produced non-finite or negative values")
    return out


def deterministic_inputs(params: PopulationParameters, C_mg_per_L: float, BW_kg: float) -> dict:
    """Point inputs for the deterministic benchmark: base-paper Table 2 values plus C and BW."""
    return {"C_mg_per_L": C_mg_per_L, "BW_kg": BW_kg, **params.deterministic}
