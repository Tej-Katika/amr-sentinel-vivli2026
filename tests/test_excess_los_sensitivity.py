"""Tests for the excess-LOS co-primary sensitivity analyses (Component 1b) and figure data.

Synthetic rows only (no real SPIDAAR record). The frame mirrors the ``load_spidaar``
columns the analyses read; unascertained patients carry ``resistant = None``.
"""

import numpy as np
import pandas as pd
import pytest

from amr_sentinel_vivli.excess_los import excess_los, km_rmst
from amr_sentinel_vivli.excess_los_figures import state_occupation_curves
from amr_sentinel_vivli.excess_los_sensitivity import (
    _rate_for_rmst,
    ascertainment_comparison,
    ascertainment_weighted_excess_los,
    exposure_assignment_bounds,
    relinkage_precision_sweep,
    simulate_rmst_precision,
    weighted_km_rmst,
)

_COLS = ["pid", "country", "resistant", "amrp", "severity",
         "los", "dead", "days_to_death", "days_observed", "days_to_enrolment"]


def _row(pid, country, resistant, severity, *, los=None, dead=0, dtpta=None, nobsd=30, enr=0):
    amrp = {1: 2, 0: 0}.get(resistant, -1)   # unascertained (resistant None) -> amrp -1
    return dict(pid=pid, country=country, resistant=resistant, amrp=amrp, severity=severity,
                los=los, dead=dead, days_to_death=dtpta, days_observed=nobsd, days_to_enrolment=enr)


def _frame(rows):
    return pd.DataFrame(rows)[_COLS]


def _synthetic():
    return _frame([
        _row(1, "Kenya", 1, 3, los=20),
        _row(2, "Kenya", 1, 2, los=18),
        _row(3, "Kenya", 1, 3, dead=1, dtpta=8),
        _row(4, "Ghana", 0, 2, los=12),
        _row(5, "Ghana", 0, 3, los=15),
        _row(6, "Kenya", None, 3, los=6),      # unascertained, short stay
        _row(7, "Ghana", None, 2, los=9),      # unascertained
        _row(8, "Kenya", None, 3, dead=1, dtpta=4),  # unascertained death
    ])


# --- (a) power / precision ---------------------------------------------------

def test_rate_for_rmst_inverts():
    tau = 28
    for target in (5.0, 12.0, 22.0):
        lam = _rate_for_rmst(target, tau)
        rmst = (1.0 - np.exp(-lam * tau)) / lam
        assert rmst == pytest.approx(target, abs=1e-3)


def test_rate_for_rmst_rejects_out_of_range():
    with pytest.raises(ValueError):
        _rate_for_rmst(30.0, 28)


def test_power_sim_calibrated_and_monotone():
    sim = simulate_rmst_precision(18.0, 135, 21,
                                  true_effects=(0.0, 3.0), tau=28, n_sim=600, seed=1)
    null, alt = sim["by_effect"]
    assert null["effect"] == 0.0 and alt["effect"] == 3.0
    assert null["mean_delta"] == pytest.approx(0.0, abs=0.4)   # unbiased under the null
    assert null["power"] < 0.15                                # ~type-I control (5% + noise)
    assert alt["power"] > null["power"]                        # power rises with the effect
    assert null["ci_width"] > 0


def test_power_sim_reproducible():
    a = simulate_rmst_precision(18.0, 135, 21, true_effects=(2.0,), n_sim=300, seed=7)
    b = simulate_rmst_precision(18.0, 135, 21, true_effects=(2.0,), n_sim=300, seed=7)
    assert a["by_effect"][0]["se"] == b["by_effect"][0]["se"]


def test_relinkage_sweep_tightens_with_n():
    sw = relinkage_precision_sweep(18.0, 135, susceptible_grid=(21, 100), n_sim=600, seed=3)
    widths = [r["ci_width"] for r in sw["sweep"]]
    assert widths[1] < widths[0]                              # larger susceptible arm -> tighter CI


# --- weighted RMST -----------------------------------------------------------

def test_weighted_km_rmst_equals_unweighted_for_equal_weights():
    t = np.array([5.0, 10.0, 28.0])
    e = np.array([True, True, False])
    w = np.ones(3)
    assert weighted_km_rmst(t, e, w, tau=28) == pytest.approx(km_rmst(t, e, tau=28))


def test_weighted_km_rmst_upweights_a_patient():
    # duplicating a short-stay patient via weight 2 must lower the restricted mean
    t = np.array([5.0, 20.0])
    e = np.array([True, True])
    base = weighted_km_rmst(t, e, np.array([1.0, 1.0]), tau=28)
    heavier_short = weighted_km_rmst(t, e, np.array([2.0, 1.0]), tau=28)
    assert heavier_short < base


# --- (b) ascertainment-selection sensitivity --------------------------------

def test_ascertainment_comparison_counts_and_fields():
    c = ascertainment_comparison(_synthetic())
    assert c["ascertained"]["n"] == 5 and c["unascertained"]["n"] == 3
    assert "mean_bed_days_to_tau" in c["ascertained"]
    assert set(c["country_ascertainment_rate"]) == {"Ghana", "Kenya"}


def test_ascertainment_weighted_matches_complete_case_value():
    df = _synthetic()
    out = ascertainment_weighted_excess_los(df)
    assert out["complete_case_excess_los"] == pytest.approx(excess_los(df)["excess_los"])
    assert np.isfinite(out["ascertainment_weighted_excess_los"])
    assert out["n_strata_no_ascertained"] == 0


def test_exposure_bounds_bracket_and_order():
    df = _synthetic()
    b = exposure_assignment_bounds(df)
    cc = b["complete_case_excess_los"]
    assert b["delta_min"] <= cc <= b["delta_max"]
    for key in ("all_unascertained_resistant", "all_unascertained_susceptible", "mar_anchor"):
        assert b["delta_min"] - 1e-9 <= b[key] <= b["delta_max"] + 1e-9
    assert b["n_unascertained"] == 3
    assert 0.0 <= b["ascertained_resistant_fraction"] <= 1.0


def test_exposure_bounds_reproducible():
    a = exposure_assignment_bounds(_synthetic())
    b = exposure_assignment_bounds(_synthetic())
    assert a["mar_anchor"] == b["mar_anchor"]            # seeded MAR split


# --- figure data -------------------------------------------------------------

def test_state_occupation_curves_are_valid_probabilities():
    cur = state_occupation_curves(_synthetic())
    assert cur["days"][0] == 0 and cur["days"][-1] == 28
    for arm in ("resistant", "susceptible"):
        p = np.array(cur[arm]["p_admitted"])
        assert p[0] == pytest.approx(1.0)               # everyone admitted at t=0
        assert np.all(np.diff(p) <= 1e-9)               # occupancy is non-increasing
        total = (np.array(cur[arm]["cif_discharge"]) + np.array(cur[arm]["cif_death"]) + p)
        assert np.allclose(total, 1.0)                  # states partition probability


def test_plot_excess_los_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    from amr_sentinel_vivli.excess_los_figures import plot_excess_los
    out = plot_excess_los(_synthetic(), tmp_path / "fig.png")
    assert out.exists() and out.stat().st_size > 0
