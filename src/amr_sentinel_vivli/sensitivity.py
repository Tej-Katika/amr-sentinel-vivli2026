"""Pre-registered sensitivity analyses (pre-reg §8).

A registry of the robustness checks committed before data access, so the final
report can confirm each was run. Implementations land alongside the steps they
probe (pre-data scaffold — pre-reg §15).
"""

from __future__ import annotations

# Pre-registered sensitivity analyses (pre-reg §8). Keep in sync with the report.
SENSITIVITY_ANALYSES: tuple[str, ...] = (
    "Step 1 e-value (VanderWeele & Ding 2017) unmeasured-confounding bound",
    "Step 1 IPTW truncation: untruncated vs 1st/99th vs 5th/95th",
    "Step 1 breakpoint regime: EUCAST v15.0 vs CLSI M100-Ed35",
    "Step 2 prior sensitivity: horseshoe tau 0.1x/10x; country RE Student-t(3)",
    "Step 3 attributable-fraction: PAF vs Cassini et al. 2019 counterfactual framing",
    "Step 4 investment normalization: per DALY vs per death vs per case; log(1+x)",
    "Subgroups: age (<18 vs >=18), sex, infection site, calendar period (2021 vs 2023)",
)


def manifest() -> tuple[str, ...]:
    """Return the locked list of pre-registered sensitivity analyses."""
    return SENSITIVITY_ANALYSES
