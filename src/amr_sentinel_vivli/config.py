"""Locked, pre-registered configuration for the AMR Sentinel Vivli 2026 analysis.

Every constant here is fixed by the OSF pre-registration
(https://doi.org/10.17605/OSF.IO/BFQDP, deposited 2026-05-29, before data access).
Changing a value after data access is a *deviation* and MUST be logged per
pre-registration §11. Treat this module as read-only.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

# --- Provenance ------------------------------------------------------------
OSF_DOI = "10.17605/OSF.IO/BFQDP"
PREREGISTRATION_URL = f"https://doi.org/{OSF_DOI}"

# --- Reproducibility (pre-reg §12) -----------------------------------------
MASTER_SEED = 20260526  # per-step seeds derived via numpy.random.SeedSequence
PYTHON_VERSION = "3.12"

# --- Geographic scope (SPIDAAR catchment; pre-reg §5) ----------------------
CATCHMENT_COUNTRIES: tuple[str, ...] = ("Ghana", "Kenya", "Malawi", "Uganda")

# --- Breakpoint / agent regimes (pre-reg §6) -------------------------------
PRIMARY_BREAKPOINT_REGIME = "EUCAST v15.0 (2025)"
SENSITIVITY_BREAKPOINT_REGIME = "CLSI M100-Ed35 (2025)"
EMPIRIC_AGENT_LIST = "WHO Model List of Essential Medicines, 23rd ed. (2023)"

# --- Step 1: SPIDAAR Cox mortality model (pre-reg §6-7) ---------------------
OUTCOME = "mortality_30d"  # 30-day all-cause mortality (binary)
PROPENSITY_COVARIATES: tuple[str, ...] = (
    "age", "sex", "infection_site", "country", "year", "organism",
)
IPTW_TRUNCATION = (0.01, 0.99)        # 1st / 99th percentile (primary)
IPTW_TRUNCATION_SENSITIVITY = ((None, None), (0.05, 0.95))  # pre-reg §8
PREDICTED_HR_RANGE = (1.5, 3.5)       # H1 predicted hazard ratio

# --- Step 2: Bayesian hierarchical projection (pre-reg §7) ------------------
PROJECTION_YEARS: tuple[int, ...] = tuple(range(2025, 2031))  # 2025-2030
NUTS_CHAINS = 4
NUTS_WARMUP = 2000
NUTS_DRAWS = 2000                     # per chain -> 8000 posterior draws
RHAT_MAX = 1.01
BULK_ESS_MIN = 400
TAIL_ESS_MIN = 400

# --- Step 3: mortality bridge (pre-reg §7) ---------------------------------
MONTE_CARLO_DRAWS = 10_000

# --- Step 4: Public R&D landscape alignment (pre-reg §6-7; reframed) --------
# Global AMR R&D Hub captures PUBLIC + PHILANTHROPIC funders only; private-sector
# pipelines are NOT captured. The Hub dashboard is subject to retrospective
# revision, so the snapshot is LOCKED to a frozen, dated, peer-reviewed extract and
# cited in every figure caption. We characterize alignment, not "gaps".
#
# SNAPSHOT LOCKED to Czaplewski et al., "An overview of global public and
# philanthropic investments into antibacterial therapeutics (2017-23)", Lancet
# Microbe 2026 (epub 2026-01-09; doi:10.1016/S2666-5247(25)00216-2) — a frozen,
# peer-reviewed extract of the Hub Dynamic Dashboard. Preferred over a hand-pulled
# live-dashboard snapshot because it is reproducible and citable. Window is the
# extract's analysed span (2017-2023). See docs/reference_rd_alignment_2026-06-09.md.
RD_HUB_WINDOW = (2017, 2023)
RD_HUB_SNAPSHOT_DATE: str | None = "2026-01-09"  # Czaplewski et al. epub; LOCKED
RD_HUB_SOURCE = (
    "Czaplewski et al., Lancet Microbe 2026 (epub 2026-01-09), "
    "Global AMR R&D Hub Dynamic Dashboard extract 2017-2023"
)
RD_HUB_SCOPE = "public + philanthropic only (private-sector R&D not captured)"

# --- Multiple-comparisons policy (pre-reg §8.6) ----------------------------
FDR_Q = 0.05  # Benjamini-Hochberg within H1 per-pathogen breakdown

# --- Paths -----------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"        # gitignored; Vivli DUA-restricted
FIGURES_DIR = PROJECT_ROOT / "figures"
RESULTS_DIR = PROJECT_ROOT / "results"

# Number of pipeline steps that draw randomness (for seed spawning).
_N_SEED_STREAMS = 16


def step_seed(step: int) -> int:
    """Deterministic per-step seed derived from MASTER_SEED (pre-reg §12).

    Args:
        step: pipeline step index (1-5) or any stable integer key.

    Returns:
        A 32-bit integer seed unique and reproducible for that step.
    """
    children = np.random.SeedSequence(MASTER_SEED).spawn(_N_SEED_STREAMS)
    return int(children[step % _N_SEED_STREAMS].generate_state(1)[0])
