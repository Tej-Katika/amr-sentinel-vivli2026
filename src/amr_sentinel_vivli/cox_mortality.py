"""Step 1 — SPIDAAR Cox proportional-hazards mortality model (pre-reg §7).

Tests pre-registered hypothesis **H1**: resistant infections carry a higher
30-day mortality hazard than susceptible infections of the same species, after
propensity-score IPTW adjustment. Reports HR, 95% CI, Schoenfeld PH test, and
the VanderWeele e-value for unmeasured confounding.

Pre-data scaffold — implement after Vivli data access (pre-reg §15).
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from . import config

_NOT_YET = "Pre-data scaffold — implement Step 1 after Vivli data access (pre-reg §15)."


@dataclass
class CoxResult:
    """Result of the pre-registered Step 1 Cox model."""

    hazard_ratio: float
    ci_low: float
    ci_high: float
    e_value: float
    e_value_ci: float
    ph_test_pvalue: float
    n_isolates: int
    h1_supported: bool  # HR > 1 at 95% after IPTW (pre-reg §10 decision rule)


def estimate_propensity_iptw(spidaar: pd.DataFrame) -> pd.Series:
    """Logistic propensity for resistance over ``config.PROPENSITY_COVARIATES``,
    returned as IPTW weights truncated at ``config.IPTW_TRUNCATION``.
    """
    raise NotImplementedError(_NOT_YET)


def fit_cox_mortality(spidaar: pd.DataFrame) -> CoxResult:
    """Fit the IPTW-weighted Cox PH model and compute the e-value.

    On a non-significant HR (≤ 1 at 95% after IPTW), the pre-registered escape
    hatch (pre-reg §10) pivots to a descriptive resistance-burden framing.
    """
    raise NotImplementedError(_NOT_YET)
