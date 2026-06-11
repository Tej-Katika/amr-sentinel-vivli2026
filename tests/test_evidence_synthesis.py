"""Tests for the Bayesian evidence synthesis (Component 1c). Synthetic inputs only.

The grid posterior is checked against analytic/symmetry cases; the cohort crude-HR helper
is checked against a hand-computed rate ratio.
"""

import numpy as np
import pandas as pd
import pytest

from amr_sentinel_vivli.evidence_synthesis import (
    EXTERNAL_MORTALITY_STUDIES,
    bayes_random_effects,
    cohort_death_cause_specific_hr,
    fixed_effect_pool,
    loghr_se_from_ci,
    run_evidence_synthesis,
)

_COLS = ["pid", "country", "resistant", "organism", "severity",
         "los", "dead", "days_to_death", "days_observed", "days_to_enrolment"]


def _row(pid, resistant, *, los=None, dead=0, dtpta=None, nobsd=30, enr=0):
    return dict(pid=pid, country="Ghana", resistant=resistant, organism="Escherichia coli",
                severity=2, los=los, dead=dead, days_to_death=dtpta,
                days_observed=nobsd, days_to_enrolment=enr)


def test_loghr_se_from_ci_roundtrip():
    loghr, se = loghr_se_from_ci(0.74, 0.42, 1.30)
    assert loghr == pytest.approx(np.log(0.74))
    assert se == pytest.approx((np.log(1.30) - np.log(0.42)) / (2 * 1.959963985))
    # SE is taken from the CI *width*; a perfectly symmetric CI round-trips exactly.
    loghr2, se2 = loghr_se_from_ci(1.0, 0.5, 2.0)
    z = 1.959963985
    assert np.exp(loghr2 - z * se2) == pytest.approx(0.5, rel=1e-9)
    assert np.exp(loghr2 + z * se2) == pytest.approx(2.0, rel=1e-9)
    # Published (rounded) point estimates need not be the geometric CI centre, so the
    # reconstructed endpoints are only approximately the reported ones.
    assert np.exp(loghr - z * se) == pytest.approx(0.42, rel=0.02)


def test_fixed_effect_pool_equal_weight_is_mean():
    studies = {
        "A": {"loghr": np.log(0.5), "se": 0.2},
        "B": {"loghr": np.log(2.0), "se": 0.2},
    }
    out = fixed_effect_pool(studies)
    assert out["pooled_loghr"] == pytest.approx(0.0, abs=1e-9)   # symmetric on log scale
    assert out["pooled_hr"] == pytest.approx(1.0)
    assert out["se"] == pytest.approx(0.2 / np.sqrt(2))


def test_bayes_random_effects_symmetry_null():
    studies = {"A": {"loghr": 0.0, "se": 0.3}, "B": {"loghr": 0.0, "se": 0.3}}
    out = bayes_random_effects(studies)
    assert out["pooled_hr"] == pytest.approx(1.0, abs=0.02)
    assert out["p_hr_gt_1"] == pytest.approx(0.5, abs=0.02)
    lo, hi = out["hr_ci"]
    assert lo < 1.0 < hi


def test_bayes_random_effects_concordant_positive():
    studies = {"A": {"loghr": 1.0, "se": 0.2}, "B": {"loghr": 1.1, "se": 0.2}}
    out = bayes_random_effects(studies)
    assert out["pooled_loghr_median"] > 0.5
    assert out["p_hr_gt_1"] > 0.95


def test_prediction_interval_wider_than_credible():
    # With genuine between-study spread, the prediction interval must exceed the CrI for mu.
    studies = {"A": {"loghr": -0.7, "se": 0.2}, "B": {"loghr": 0.7, "se": 0.2}}
    out = bayes_random_effects(studies)
    cr_lo, cr_hi = out["hr_ci"]
    pr_lo, pr_hi = out["prediction_hr_ci"]
    assert pr_lo < cr_lo and pr_hi > cr_hi


def test_cohort_crude_hr_matches_rate_ratio():
    rows = [
        _row(1, 1, dead=1, dtpta=4), _row(2, 1, dead=1, dtpta=4), _row(3, 1, los=10),
        _row(4, 0, dead=1, dtpta=5), _row(5, 0, los=8),
    ]
    df = pd.DataFrame(rows)[_COLS]
    out = cohort_death_cause_specific_hr(df)
    assert out["deaths_resistant"] == 2
    assert out["deaths_susceptible"] == 1
    assert out["person_days_resistant"] == pytest.approx(4 + 4 + 10)
    assert out["person_days_susceptible"] == pytest.approx(5 + 8)
    rate_r = 2 / 18
    rate_s = 1 / 13
    assert out["hr"] == pytest.approx(rate_r / rate_s)
    assert out["se"] == pytest.approx(np.sqrt(1 / 2 + 1 / 1))
    assert out["adjusted"] is False


def test_run_evidence_synthesis_keys_and_null_primary():
    rows = [_row(i, 1, dead=1, dtpta=5) for i in range(3)] + \
           [_row(10 + i, 1, los=12) for i in range(5)] + \
           [_row(20, 0, dead=1, dtpta=6)] + [_row(21 + i, 0, los=9) for i in range(4)]
    df = pd.DataFrame(rows)[_COLS]
    out = run_evidence_synthesis(df)
    for k in ("cohort_crude", "primary_pool_adjusted", "sensitivity_pool_with_cohort",
              "fixed_effect_adjusted", "tau_prior_sensitivity"):
        assert k in out
    primary = out["primary_pool_adjusted"]
    assert primary["n_studies"] == len(EXTERNAL_MORTALITY_STUDIES)
    # the two adjusted external studies pool to a null (CrI spans HR=1)
    lo, hi = primary["hr_ci"]
    assert lo < 1.0 < hi
