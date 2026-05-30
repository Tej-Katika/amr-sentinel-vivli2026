"""End-to-end orchestration of the five pre-registered steps (pre-reg §7).

Documents the intended data flow. Each call delegates to a step module that is a
pre-data stub until Vivli data access (pre-reg §15), so running this now raises
``NotImplementedError`` at the first step — by design.
"""

from __future__ import annotations

from . import data_loading
from .bayesian_projection import fit_projection
from .cox_mortality import fit_cox_mortality
from .mortality_bridge import bridge_to_mortality
from .rd_alignment import analyze_alignment
from .stewardship_gformula import simulate_switch


def run() -> None:
    """Run Steps 1-5 in order. Intended to execute inside the Vivli environment."""
    spidaar = data_loading.load_spidaar()
    atlas = data_loading.load_atlas()
    smart = data_loading.load_smart()
    rd_hub = data_loading.load_rd_hub_snapshot()

    cox = fit_cox_mortality(spidaar)                       # Step 1 (H1)
    projection = fit_projection(_combine(atlas, smart))    # Step 2 (H2)
    burden = bridge_to_mortality(cox, projection, _who_ghe_denominators())  # Step 3
    _alignment = analyze_alignment(burden, rd_hub)         # Step 4 (H3)
    _intervention = simulate_switch(burden, organism="", site="", switch="")  # Step 5 (H4)


def _combine(atlas, smart):
    raise NotImplementedError("Pre-data scaffold (pre-reg §15).")


def _who_ghe_denominators():
    raise NotImplementedError("Pre-data scaffold (pre-reg §15).")


if __name__ == "__main__":
    run()
