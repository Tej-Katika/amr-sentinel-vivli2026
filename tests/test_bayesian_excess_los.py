"""Tests for the Bayesian partial-pooled excess-LOS companion (Component 2).

Synthetic rows only. The Laplace/penalized-Poisson pieces are checked against closed-form
or known-coefficient cases; the end-to-end estimator is checked for valid posteriors.
"""

import numpy as np
import pandas as pd
import pytest

from amr_sentinel_vivli.bayesian_excess_los import (
    BIN_EDGES,
    _fit_penalized_poisson,
    _gaussian_overlap,
    _gram_negative,
    _rmst_admitted,
    bayesian_excess_los,
    expand_person_intervals,
    prior_sensitivity,
)

_COLS = ["pid", "country", "resistant", "organism", "severity",
         "los", "dead", "days_to_death", "days_observed", "days_to_enrolment"]


def _row(pid, country, resistant, organism, severity, *, los=None, dead=0, dtpta=None,
         nobsd=30, enr=0):
    return dict(pid=pid, country=country, resistant=resistant, organism=organism,
                severity=severity, los=los, dead=dead, days_to_death=dtpta,
                days_observed=nobsd, days_to_enrolment=enr)


def test_gram_negative_classifier():
    assert _gram_negative("Escherichia coli") == 1
    assert _gram_negative("Klebsiella pneumoniae, Acinetobacter spp") == 1
    assert _gram_negative("Staphylococcus aureus") == 0
    assert _gram_negative("Coagulase negative Staphylococcus") == 0


def test_expand_person_intervals_exposure_and_events():
    df = pd.DataFrame([_row(1, "Kenya", 1, "Escherichia coli", 3, los=10)])[_COLS]
    long = expand_person_intervals(df)
    # bins (0-2,2-7,7-14): exposures 2,5,3 = 10 total; discharge in the 3rd bin (7-14)
    assert list(long["bin"]) == [0, 1, 2]
    assert long["exposure"].sum() == pytest.approx(10.0)
    assert long.iloc[-1]["discharge"] == 1 and long["death"].sum() == 0
    assert long["gram_neg"].iloc[0] == 1 and long["sev_high"].iloc[0] == 1


def test_expand_person_intervals_death_event():
    df = pd.DataFrame([_row(1, "Kenya", 1, "Klebsiella pneumoniae", 2, dead=1, dtpta=5)])[_COLS]
    long = expand_person_intervals(df)
    assert long["death"].sum() == 1 and long["discharge"].sum() == 0
    assert long["exposure"].sum() == pytest.approx(5.0)


def test_fit_penalized_poisson_recovers_coefficients():
    rng = np.random.default_rng(0)
    n = 4000
    x = np.column_stack([np.ones(n), rng.normal(size=n), rng.normal(size=n)])
    true = np.array([0.2, 0.5, -0.3])
    offset = np.zeros(n)
    y = rng.poisson(np.exp(x @ true))
    weak_prec = np.array([1e-6, 1e-6, 1e-6])
    beta, cov = _fit_penalized_poisson(x, y, offset, np.zeros(3), weak_prec)
    assert beta == pytest.approx(true, abs=0.08)
    assert cov.shape == (3, 3) and np.all(np.diag(cov) > 0)


def test_fit_penalized_poisson_prior_shrinks_toward_mean():
    rng = np.random.default_rng(1)
    n = 60
    x = np.column_stack([np.ones(n), rng.normal(size=n)])
    y = rng.poisson(np.exp(x @ np.array([0.0, 1.5])))
    weak = _fit_penalized_poisson(x, y, np.zeros(n), np.zeros(2), np.array([1e-6, 1e-6]))[0]
    strong = _fit_penalized_poisson(x, y, np.zeros(n), np.zeros(2), np.array([1e-6, 100.0]))[0]
    # a tight prior centered at 0 pulls the slope meaningfully toward 0
    assert 0 < strong[1] < weak[1]


def test_rmst_admitted_constant_hazard_matches_closed_form():
    widths = np.diff(BIN_EDGES).astype(float)
    lam = np.full(len(widths), 0.05)
    tau = float(sum(widths))
    expected = (1 - np.exp(-0.05 * tau)) / 0.05
    assert _rmst_admitted(lam, widths) == pytest.approx(expected, abs=1e-6)


def test_rmst_admitted_zero_hazard_is_full_horizon():
    widths = np.diff(BIN_EDGES).astype(float)
    assert _rmst_admitted(np.zeros(len(widths)), widths) == pytest.approx(sum(widths))


def test_gaussian_overlap_bounds():
    assert _gaussian_overlap(0, 1, 0, 1) == pytest.approx(1.0, abs=1e-2)
    assert _gaussian_overlap(0, 0.3, 10, 0.3) < 0.01


def _synthetic_cohort(n_per=12):
    rng = np.random.default_rng(3)
    rows, pid = [], 0
    for country in ("Ghana", "Kenya", "Malawi", "Uganda"):
        for resistant in (0, 1):
            for _ in range(n_per):
                pid += 1
                organism = "Escherichia coli" if rng.random() < 0.8 else "Staphylococcus aureus"
                sev = 3 if rng.random() < 0.5 else 2
                if rng.random() < 0.1:
                    rows.append(_row(pid, country, resistant, organism, sev,
                                     dead=1, dtpta=int(rng.integers(2, 14))))
                else:
                    rows.append(_row(pid, country, resistant, organism, sev,
                                     los=int(rng.integers(4, 24))))
    return pd.DataFrame(rows)[_COLS]


def test_bayesian_excess_los_valid_posterior():
    out = bayesian_excess_los(_synthetic_cohort(), n_draws=800, seed=1)
    assert 0.0 <= out["p_excess_gt_0"] <= 1.0
    assert 0.0 <= out["p_excess_gt_1"] <= out["p_excess_gt_0"] + 1e-9
    lo, hi = out["excess_los_hdi"]
    assert lo <= out["excess_los_mean"] <= hi
    assert np.isfinite(out["discharge_log_hr_mean"])


def test_bayesian_excess_los_reproducible():
    a = bayesian_excess_los(_synthetic_cohort(), n_draws=500, seed=7)
    b = bayesian_excess_los(_synthetic_cohort(), n_draws=500, seed=7)
    assert a["excess_los_mean"] == b["excess_los_mean"]


def test_prior_sensitivity_separates_prior_weight():
    ps = prior_sensitivity(_synthetic_cohort(), n_draws=600, seed=2)
    assert set(ps) == {"null_centered", "flat", "informative_mbira"}
    # the flat prior contributes almost no posterior precision; the tight one contributes more
    assert (ps["flat"]["prior_precision_fraction"]
            < ps["informative_mbira"]["prior_precision_fraction"])
    for res in ps.values():
        assert 0.0 <= res["prior_posterior_overlap"] <= 1.0
