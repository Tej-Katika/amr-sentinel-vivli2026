"""Tests for the excess-length-of-stay (excess bed-days) competing-risks core.

Synthetic rows only (no real SPIDAAR record). The frame mirrors the ``load_spidaar``
columns the estimator reads; expected restricted-mean values are hand-computed from
the Kaplan-Meier of time-to-exit (discharge OR death) under enrolment time origin.
"""

import numpy as np
import pandas as pd
import pytest

from amr_sentinel_vivli.excess_los import (
    DEFAULT_TAU,
    EXIT_CENSORED,
    EXIT_DEATH,
    EXIT_DISCHARGE,
    EXIT_FRAME_COLUMNS,
    bootstrap_excess_los_ci,
    build_exit_frame,
    cif_at_tau,
    cif_decomposition,
    excess_los,
    km_rmst,
    standardized_excess_los,
)


def _frame(rows):
    cols = ["pid", "country", "resistant", "los", "dead", "days_to_death",
            "days_observed", "days_to_enrolment", "severity"]
    return pd.DataFrame(rows)[cols]


def _synthetic():
    return _frame([
        # resistant: discharge @ los-enr = 18
        dict(pid=1, country="Kenya", resistant=1, los=20, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=2, severity=3),
        # resistant: death @ dtpta-enr = 8 (competing event)
        dict(pid=2, country="Kenya", resistant=1, los=None, dead=1, days_to_death=10,
             days_observed=30, days_to_enrolment=2, severity=2),
        # susceptible: discharge @ 10
        dict(pid=3, country="Ghana", resistant=0, los=12, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=2, severity=2),
        # susceptible: still admitted -> censored @ nobsd = 15
        dict(pid=4, country="Ghana", resistant=0, los=None, dead=0, days_to_death=None,
             days_observed=15, days_to_enrolment=0, severity=3),
        # resistant: los-enr = 35 > min(nobsd,tau)=28 -> censored @ 28
        dict(pid=5, country="Uganda", resistant=1, los=40, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=5, severity=3),
        # susceptible: los-enr = -2 -> clipped to discharge @ 0
        dict(pid=6, country="Uganda", resistant=0, los=3, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=5, severity=2),
        # unascertained exposure -> excluded from the contrast
        dict(pid=7, country="Kenya", resistant=None, los=9, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=1, severity=2),
    ])


def test_exit_frame_columns_and_event_construction():
    ef = build_exit_frame(_synthetic()).set_index("pid")
    assert list(build_exit_frame(_synthetic()).columns) == list(EXIT_FRAME_COLUMNS)
    assert ef.loc[1, "time"] == 18 and ef.loc[1, "exit_type"] == EXIT_DISCHARGE
    assert ef.loc[2, "time"] == 8 and ef.loc[2, "exit_type"] == EXIT_DEATH
    assert ef.loc[2, "exited"]


def test_exit_frame_admin_censor_and_clip():
    ef = build_exit_frame(_synthetic()).set_index("pid")
    # los-enr exceeds min(nobsd, tau) -> censored at tau, still admitted
    assert ef.loc[5, "time"] == DEFAULT_TAU
    assert ef.loc[5, "exit_type"] == EXIT_CENSORED and not ef.loc[5, "exited"]
    # still-admitted censored at observation window
    assert ef.loc[4, "time"] == 15 and ef.loc[4, "exit_type"] == EXIT_CENSORED
    # negative enrolment-origin time clipped to 0 (discharge at enrolment)
    assert ef.loc[6, "time"] == 0 and ef.loc[6, "exit_type"] == EXIT_DISCHARGE


def test_km_rmst_simple_cases():
    # everyone exits at day 10 -> expected bed-days = 10
    t = np.array([10.0, 10.0, 10.0])
    assert km_rmst(t, np.array([True, True, True]), tau=28) == pytest.approx(10.0)
    # no events observed -> stays admitted whole horizon -> RMST = tau
    assert km_rmst(t, np.array([False, False, False]), tau=28) == pytest.approx(28.0)


def test_excess_los_point_estimate():
    out = excess_los(_synthetic())
    assert out["n_resistant"] == 3 and out["n_susceptible"] == 3   # pid 7 excluded
    assert out["rmst_resistant"] == pytest.approx(18.0)
    assert out["rmst_susceptible"] == pytest.approx(12.6667, abs=1e-3)
    assert out["excess_los"] == pytest.approx(5.3333, abs=1e-3)


def test_bootstrap_ci_brackets_point_and_is_reproducible():
    a = bootstrap_excess_los_ci(_synthetic(), n_boot=200, seed=123)
    b = bootstrap_excess_los_ci(_synthetic(), n_boot=200, seed=123)
    assert a["ci_lower"] == b["ci_lower"] and a["ci_upper"] == b["ci_upper"]  # seeded
    assert a["ci_lower"] <= a["excess_los"] <= a["ci_upper"]


def test_cif_at_tau_competing_risks():
    # 1 death @5, 2 discharges @10 -> death CIF 1/3, discharge CIF 2/3, none left
    time = np.array([5.0, 10.0, 10.0])
    et = np.array([EXIT_DEATH, EXIT_DISCHARGE, EXIT_DISCHARGE])
    cif = cif_at_tau(time, et, tau=28)
    assert cif["death_cif"] == pytest.approx(1 / 3)
    assert cif["discharge_cif"] == pytest.approx(2 / 3)
    assert cif["still_admitted"] == pytest.approx(0.0, abs=1e-9)


def test_standardized_excess_los_uses_both_strata():
    out = standardized_excess_los(_synthetic())
    assert set(out["strata_used"]) == {"low", "high"}    # both arms populated in each
    assert out["strata_dropped"] == []
    lo = out["per_stratum"]["low"]
    assert lo["n_resistant"] >= 1 and lo["n_susceptible"] >= 1
    assert np.isfinite(out["standardized_excess_los"])


def test_standardized_drops_unsupported_stratum():
    # susceptible only in 'low' (severity 2); 'high' has no susceptible -> dropped
    frame = _frame([
        dict(pid=1, country="Kenya", resistant=1, los=20, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=2, severity=3),
        dict(pid=2, country="Kenya", resistant=1, los=12, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=2, severity=2),
        dict(pid=3, country="Ghana", resistant=0, los=15, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=2, severity=2),
    ])
    out = standardized_excess_los(frame)
    assert out["strata_used"] == ["low"] and out["strata_dropped"] == ["high"]


def test_cif_decomposition_arms():
    out = cif_decomposition(_synthetic())
    assert out["resistant"]["n"] == 3 and out["susceptible"]["n"] == 3
    for arm in ("resistant", "susceptible"):
        total = (out[arm]["discharge_cif"] + out[arm]["death_cif"]
                 + out[arm]["still_admitted"])
        assert total == pytest.approx(1.0, abs=1e-9)
