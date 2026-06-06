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
    excess_los,
    km_rmst,
)


def _frame(rows):
    cols = ["pid", "country", "resistant", "los", "dead", "days_to_death",
            "days_observed", "days_to_enrolment"]
    return pd.DataFrame(rows)[cols]


def _synthetic():
    return _frame([
        # resistant: discharge @ los-enr = 18
        dict(pid=1, country="Kenya", resistant=1, los=20, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=2),
        # resistant: death @ dtpta-enr = 8 (competing event)
        dict(pid=2, country="Kenya", resistant=1, los=None, dead=1, days_to_death=10,
             days_observed=30, days_to_enrolment=2),
        # susceptible: discharge @ 10
        dict(pid=3, country="Ghana", resistant=0, los=12, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=2),
        # susceptible: still admitted -> censored @ nobsd = 15
        dict(pid=4, country="Ghana", resistant=0, los=None, dead=0, days_to_death=None,
             days_observed=15, days_to_enrolment=0),
        # resistant: los-enr = 35 > min(nobsd,tau)=28 -> censored @ 28
        dict(pid=5, country="Uganda", resistant=1, los=40, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=5),
        # susceptible: los-enr = -2 -> clipped to discharge @ 0
        dict(pid=6, country="Uganda", resistant=0, los=3, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=5),
        # unascertained exposure -> excluded from the contrast
        dict(pid=7, country="Kenya", resistant=None, los=9, dead=0, days_to_death=None,
             days_observed=30, days_to_enrolment=1),
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
