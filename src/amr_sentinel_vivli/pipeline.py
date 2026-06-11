"""End-to-end orchestration of the pivoted analysis (see docs/strategy_2026.md).

The headline moved from a resistance→mortality bridge to resistance-attributable
**excess bed-days** (competing risks), with the empiric-adequacy g-formula as the
systemic-leverage centerpiece; SMART is excluded so external surveillance is ATLAS-only
(a catchment nowcast + SPIDAAR frame-contrast, not a 2030 projection). Components run on
the delivered data inside the Vivli secure environment.
"""

from __future__ import annotations

from . import data_loading
from .bayesian_excess_los import bayesian_excess_los, prior_sensitivity
from .bayesian_projection import frame_contrast, run_nowcast
from .evidence_synthesis import run_evidence_synthesis
from .excess_los import bootstrap_excess_los_ci, cif_decomposition, standardized_excess_los
from .excess_los_sensitivity import (
    ascertainment_weighted_excess_los,
    exposure_assignment_bounds,
    simulate_rmst_precision,
)
from .rd_alignment import catchment_alignment, gram_panel_alignment
from .stewardship_gformula import run_stewardship_gformula


def run() -> dict:
    """Run the pivoted components in order on the delivered Vivli data."""
    spidaar = data_loading.load_spidaar()
    spidaar_isolates = data_loading.load_spidaar_isolates()
    atlas = data_loading.load_atlas()

    # Component 1 (primary) + 1b (co-primary honesty analyses)
    excess = bootstrap_excess_los_ci(spidaar)
    excess_std = standardized_excess_los(spidaar)
    cif = cif_decomposition(spidaar)
    power = simulate_rmst_precision(excess_std["standardized_rmst_susceptible"],
                                    n_resistant=135, n_susceptible=21)
    ascertainment = {
        "weighted": ascertainment_weighted_excess_los(spidaar),
        "bounds": exposure_assignment_bounds(spidaar),
    }
    # Component 1c (co-primary honesty analysis): Bayesian evidence synthesis of the
    # adjusted SSA/LMIC resistance->mortality literature, with our cohort placed alongside.
    evidence = run_evidence_synthesis(spidaar)

    # Component 2 (secondary): Bayesian partial-pooled excess bed-days
    bayesian = {
        "posterior": bayesian_excess_los(spidaar),
        "prior_sensitivity": prior_sensitivity(spidaar),
    }

    # Component 3: ATLAS-only catchment nowcast + SPIDAAR frame-contrast
    nowcast = run_nowcast(atlas)
    contrast = frame_contrast(atlas, spidaar_isolates)

    # Component 5: empiric-adequacy stewardship g-formula (centerpiece)
    stewardship = run_stewardship_gformula(spidaar)

    # Component 4 (R&D mismatch, Cross-Domain): verified GRAM burden vs the locked Hub
    # funding snapshot (config.RD_HUB_SNAPSHOT_DATE); the one unfetched funding split is
    # propagated as Monte-Carlo uncertainty — see rd_alignment.gram_panel_alignment.
    rd_alignment = gram_panel_alignment()
    # Catchment-specific Cross-Domain finding: re-weight the mismatch by the local severe-HAI
    # pathogen mix (SPIDAAR isolates) — the global ranking made regional with our own data.
    rd_alignment_catchment = catchment_alignment(spidaar_isolates)

    return {
        "excess_los": excess,
        "excess_los_standardized": excess_std,
        "cif": cif,
        "power_simulation": power,
        "ascertainment_sensitivity": ascertainment,
        "evidence_synthesis": evidence,
        "bayesian_excess_los": bayesian,
        "nowcast": nowcast,
        "frame_contrast": contrast,
        "stewardship": stewardship,
        "rd_alignment": rd_alignment,
        "rd_alignment_catchment": rd_alignment_catchment,
    }


if __name__ == "__main__":
    run()
