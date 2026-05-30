"""Step 4 — Public R&D landscape alignment analysis (pre-reg §7; reframed).

Tests pre-registered hypothesis **H3**. Characterizes how cumulative *public and
philanthropic* R&D investment (Global AMR R&D Hub, 2017-2024, single locked
snapshot) aligns with projected mortality burden, via a Spearman rank correlation
with bootstrap CI.

IMPORTANT (pre-reg §6/§7): the Hub captures public + philanthropic funding ONLY;
private-sector pipelines are excluded. Consistent with the Hub's own published
guidance, this step characterizes *alignment* and does NOT claim definitive R&D
"gaps". Every figure must cite ``config.RD_HUB_SNAPSHOT_DATE``.

Pre-data scaffold — implement after Vivli data access (pre-reg §15).
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from . import config
from .mortality_bridge import AttributableMortality

_NOT_YET = "Pre-data scaffold — implement Step 4 after Vivli data access (pre-reg §15)."


@dataclass
class AlignmentResult:
    """Burden-vs-public-investment alignment for WHO BPPL pathogens."""

    spearman_rho: float
    ci_low: float
    ci_high: float
    snapshot_date: str             # config.RD_HUB_SNAPSHOT_DATE, for figure captions
    scope_note: str                # config.RD_HUB_SCOPE
    per_pathogen: pd.DataFrame     # descriptive table (burden vs investment); no per-row test


def alignment_caption() -> str:
    """Standard figure-caption suffix required on every Step 4 figure."""
    return (
        f"Global AMR R&D Hub snapshot {config.RD_HUB_SNAPSHOT_DATE}; "
        f"{config.RD_HUB_SCOPE}; window {config.RD_HUB_WINDOW[0]}-{config.RD_HUB_WINDOW[1]}."
    )


def analyze_alignment(
    burden: AttributableMortality,
    rd_hub: pd.DataFrame,
) -> AlignmentResult:
    """Rank burden against public+philanthropic investment and report Spearman ρ
    with a bootstrap CI (the single pre-registered inferential test for Step 4).
    """
    if config.RD_HUB_SNAPSHOT_DATE is None:
        raise ValueError("Lock config.RD_HUB_SNAPSHOT_DATE before running Step 4 (pre-reg §6/§7).")
    raise NotImplementedError(_NOT_YET)
