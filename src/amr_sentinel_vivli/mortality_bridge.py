"""Step 3 — Bridge from projected resistance to attributable mortality (pre-reg §7).

Applies the Step 1 Cox hazard ratio (sampled from its full distribution) to the
Step 2 projected resistance rates through the population-attributable-fraction
formulation, scaled by WHO Global Health Estimates population denominators.
Uncertainty is propagated by a ``config.MONTE_CARLO_DRAWS``-draw Monte Carlo.

SUPERSEDED by the excess-LOS pivot (see docs/strategy_2026.md): the resistance→mortality
bridge is replaced by the excess-bed-day burden (Component 1) and the empiric-adequacy
g-formula (Component 5); retained for provenance. Pre-data scaffold otherwise.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from .bayesian_projection import ProjectionResult
from .cox_mortality import CoxResult

_NOT_YET = "Pre-data scaffold — implement Step 3 after Vivli data access (pre-reg §15)."


@dataclass
class AttributableMortality:
    """Attributable-death estimates per country × pathogen × year (median + 95% PI)."""

    table: pd.DataFrame  # columns: country, pathogen, year, deaths_median, pi_low, pi_high


def population_attributable_fraction(p_resistant: float, hazard_ratio: float) -> float:
    """PAF = p(R)·(HR-1) / (1 + p(R)·(HR-1)) — pre-reg §7 Step 3."""
    num = p_resistant * (hazard_ratio - 1.0)
    return num / (1.0 + num)


def bridge_to_mortality(
    cox: CoxResult,
    projection: ProjectionResult,
    population_denominators: pd.DataFrame,
) -> AttributableMortality:
    """Monte-Carlo propagate (Cox HR distribution × projection posterior × WHO GHE
    intervals) into attributable deaths per country × pathogen × year.
    """
    raise NotImplementedError(_NOT_YET)
