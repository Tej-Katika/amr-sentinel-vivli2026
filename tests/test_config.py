"""Guardrail tests: the locked, pre-registered constants must not silently drift.

If a value here legitimately changes after data access, that is a *deviation*
(pre-reg §11): update the pre-registration record and this test together, and
document it in the final report. A failing test here is a feature, not a bug.
"""

from amr_sentinel_vivli import config


def test_osf_doi_matches_deposit():
    assert config.OSF_DOI == "10.17605/OSF.IO/BFQDP"


def test_master_seed_locked():
    assert config.MASTER_SEED == 20260526


def test_catchment_is_spidaar_only():
    assert config.CATCHMENT_COUNTRIES == ("Ghana", "Kenya", "Malawi", "Uganda")


def test_rd_hub_window_is_full_collection_span():
    # Reframed Step 4: 2017-2024 (Hub's actual collection start), not 2020-2024.
    assert config.RD_HUB_WINDOW == (2017, 2024)


def test_primary_breakpoint_regime():
    assert config.PRIMARY_BREAKPOINT_REGIME.startswith("EUCAST v15.0")


def test_rd_hub_snapshot_must_be_locked_before_use():
    # Intentionally unset in the scaffold; loaders/Step 4 refuse to run until set.
    assert config.RD_HUB_SNAPSHOT_DATE is None


def test_step_seed_deterministic_and_distinct():
    assert config.step_seed(1) == config.step_seed(1)
    assert config.step_seed(1) != config.step_seed(2)
