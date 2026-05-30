"""Step 2 — Bayesian hierarchical resistance projection (pre-reg §7).

Tests pre-registered hypothesis **H2**. PyMC model on logit(SIR rate) with
country/pathogen/drug fixed effects, country random intercepts + slopes on time,
and a regularized-horseshoe prior on interaction terms. NUTS sampling with
pre-specified convergence criteria; posterior projection to 2025-2030.

Pre-data scaffold — implement after Vivli data access (pre-reg §15).
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from . import config

_NOT_YET = "Pre-data scaffold — implement Step 2 after Vivli data access (pre-reg §15)."


@dataclass
class ProjectionResult:
    """Posterior resistance-rate projection per country × pathogen × drug × year."""

    posterior: object  # arviz.InferenceData
    converged: bool     # R-hat < RHAT_MAX and ESS thresholds met (pre-reg §7)
    max_rhat: float
    min_bulk_ess: float
    min_tail_ess: float
    n_divergences: int


def build_model(atlas_smart: pd.DataFrame):
    """Construct the PyMC hierarchical model (pre-reg §7 Step 2)."""
    raise NotImplementedError(_NOT_YET)


def fit_projection(atlas_smart: pd.DataFrame) -> ProjectionResult:
    """Sample the model (``config.NUTS_CHAINS`` chains, warmup/draws per pre-reg)
    and check convergence. On non-convergence after the pre-registered re-run,
    fall back to pooled country×pathogen models (pre-reg §10).
    """
    raise NotImplementedError(_NOT_YET)
