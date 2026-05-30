"""Step 5 — Stewardship g-formula simulator (pre-reg §7).

Tests pre-registered hypothesis **H4**. A Bayesian g-formula (Snowden et al.,
AJE 2011) over counterfactual susceptibility for empiric-therapy switches at
facility level, surfaced as a Streamlit "what-if" tool that catchment-region
users re-run with their own antibiogram, breakpoints, and availability params.

Pre-data scaffold — implement after Vivli data access (pre-reg §15).
"""

from __future__ import annotations

from dataclasses import dataclass

from .mortality_bridge import AttributableMortality

_NOT_YET = "Pre-data scaffold — implement Step 5 after Vivli data access (pre-reg §15)."


@dataclass
class InterventionEffect:
    """Projected facility-level mortality change for one empiric-switch scenario."""

    organism: str
    site: str
    switch: str                 # e.g. "ceftriaxone -> meropenem"
    mortality_reduction: float  # posterior mean
    ci_low: float
    ci_high: float
    credible_excludes_zero: bool  # H4 success criterion (pre-reg §7)


def simulate_switch(
    burden: AttributableMortality,
    organism: str,
    site: str,
    switch: str,
) -> InterventionEffect:
    """Run the Bayesian g-formula for a single empiric-switch scenario."""
    raise NotImplementedError(_NOT_YET)
