"""Tests for the empiric-adequacy stewardship g-formula (Component 5).

Synthetic rows only (no real SPIDAAR record). The frame mirrors the ``load_spidaar``
columns the estimator reads; the treatment node is ``treatment_adequacy`` (txadp).
Expected standardized values are hand-computed from the Kaplan-Meier restricted means
of time-to-exit within each confounder stratum.
"""

import pandas as pd
import pytest

from amr_sentinel_vivli.stewardship_gformula import (
    ADEQUATE,
    BED_DAY_COST_USD_2010,
    INADEQUATE,
    TXADP_ADEQUATE_CODE,
    TXADP_INADEQUATE_CODE,
    adequacy_scenario,
    binarize_adequacy,
    build_calibration_artifact,
    cost_of_bed_days,
    gformula_bed_days,
    gformula_by_resistance,
    monte_carlo_cost,
    positivity_diagnostic,
    run_stewardship_gformula,
    scale_to_population,
    synthetic_cohort,
)

_COLS = ["pid", "country", "resistant", "treatment_adequacy", "severity",
         "los", "dead", "days_to_death", "days_observed", "days_to_enrolment"]


def _row(pid, country, resistant, adeq, severity, *, los=None, dead=0,
         dtpta=None, nobsd=30, enr=0):
    return dict(pid=pid, country=country, resistant=resistant, treatment_adequacy=adeq,
                severity=severity, los=los, dead=dead, days_to_death=dtpta,
                days_observed=nobsd, days_to_enrolment=enr)


def _balanced_frame():
    """One stratum (Kenya/high severity) with both arms; clean integer exit times.

    Adequate arm discharges at days {10, 10}; inadequate arm discharges at {20, 20}
    (enrolment origin, enr=0). KM restricted mean to tau is then exactly the discharge
    day for each arm (everyone exits at one time). avertable = 20 - 10 = 10 bed-days.
    """
    AQ, IQ = TXADP_ADEQUATE_CODE, TXADP_INADEQUATE_CODE
    return pd.DataFrame([
        _row(1, "Kenya", 1, AQ,3, los=10),
        _row(2, "Kenya", 1, AQ,3, los=10),
        _row(3, "Kenya", 1, IQ,3, los=20),
        _row(4, "Kenya", 1, IQ,3, los=20),
    ])[_COLS]


def test_binarize_adequacy_codes_and_unknown():
    s = binarize_adequacy(pd.Series([TXADP_ADEQUATE_CODE, TXADP_INADEQUATE_CODE, 9, -1, None]))
    assert s.iloc[0] == 1.0 and s.iloc[1] == 0.0
    assert s.iloc[2:].isna().all()          # unknown sentinels -> NaN (excluded)


def test_gformula_point_estimate_single_stratum():
    out = gformula_bed_days(_balanced_frame())
    assert out["n_ascertained"] == 4 and out["n_on_support"] == 4
    assert out["strata_dropped"] == []
    assert out["bed_days_set_adequate"] == pytest.approx(10.0)
    assert out["bed_days_set_inadequate"] == pytest.approx(20.0)
    # avertable = inadequate - adequate (bed-days saved by guaranteeing adequacy)
    assert out["avertable_bed_days"] == pytest.approx(10.0)
    # natural course is the observed 50/50 mix -> midpoint
    assert out["bed_days_natural"] == pytest.approx(15.0)
    assert out["avertable_vs_natural"] == pytest.approx(5.0)


def test_gformula_standardizes_over_confounder_distribution():
    """Two strata with different effects standardize to the cohort's stratum weights."""
    AQ, IQ = TXADP_ADEQUATE_CODE, TXADP_INADEQUATE_CODE
    frame = pd.DataFrame([
        # Kenya/high (weight 2/6): adequate 10, inadequate 20 -> within-effect 10
        _row(1, "Kenya", 1, AQ,3, los=10),
        _row(2, "Kenya", 1, IQ,3, los=20),
        # Ghana/low (weight 4/6): adequate 5, inadequate 9 -> within-effect 4
        _row(3, "Ghana", 0, AQ,2, los=5),
        _row(4, "Ghana", 0, AQ,2, los=5),
        _row(5, "Ghana", 0, IQ,2, los=9),
        _row(6, "Ghana", 0, IQ,2, los=9),
    ])[_COLS]
    out = gformula_bed_days(frame)
    # standardized adequate = (2/6)*10 + (4/6)*5 = 6.6667
    assert out["bed_days_set_adequate"] == pytest.approx(20 / 3, abs=1e-6)
    # standardized inadequate = (2/6)*20 + (4/6)*9 = 12.6667
    assert out["bed_days_set_inadequate"] == pytest.approx(38 / 3, abs=1e-6)
    assert out["avertable_bed_days"] == pytest.approx(6.0, abs=1e-6)


def test_gformula_drops_unsupported_stratum():
    """A stratum with only one arm is dropped (positivity) and reported, not imputed."""
    AQ, IQ = TXADP_ADEQUATE_CODE, TXADP_INADEQUATE_CODE
    frame = pd.DataFrame([
        _row(1, "Kenya", 1, AQ,3, los=10),
        _row(2, "Kenya", 1, IQ,3, los=20),
        _row(3, "Ghana", 0, AQ,2, los=5),   # Ghana/low has NO inadequate -> dropped
    ])[_COLS]
    out = gformula_bed_days(frame)
    assert ("low", "Ghana") in out["strata_dropped"]       # key order = (severity_bin, country)
    assert ("high", "Kenya") in out["strata_supported"]
    assert out["n_on_support"] == 2


def test_gformula_raises_when_no_overlap():
    AQ, IQ = TXADP_ADEQUATE_CODE, TXADP_INADEQUATE_CODE
    frame = pd.DataFrame([
        _row(1, "Kenya", 1, AQ,3, los=10),
        _row(2, "Ghana", 0, IQ,2, los=20),  # different strata, no overlap anywhere
    ])[_COLS]
    with pytest.raises(ValueError, match="not identified"):
        gformula_bed_days(frame)


def test_unknown_adequacy_excluded():
    frame = _balanced_frame()
    extra = pd.DataFrame([_row(99, "Kenya", 1, 9, 3, los=99)])[_COLS]  # unknown txadp
    out = gformula_bed_days(pd.concat([frame, extra], ignore_index=True))
    assert out["n_ascertained"] == 4            # the unknown-adequacy row is dropped


def test_death_competing_contrast_present():
    """Inadequate arm has a death; the death CIF contrast is surfaced alongside bed-days."""
    AQ, IQ = TXADP_ADEQUATE_CODE, TXADP_INADEQUATE_CODE
    frame = pd.DataFrame([
        _row(1, "Kenya", 1, AQ,3, los=10),
        _row(2, "Kenya", 1, AQ,3, los=10),
        _row(3, "Kenya", 1, IQ,3, los=20),
        _row(4, "Kenya", 1, IQ,3, dead=1, dtpta=8),   # death competes
    ])[_COLS]
    out = gformula_bed_days(frame)
    assert out["death_cif_set_inadequate"] > out["death_cif_set_adequate"]
    assert out["averted_death_fraction"] == pytest.approx(
        out["death_cif_set_inadequate"] - out["death_cif_set_adequate"]
    )


def test_by_resistance_controlled_direct_effect():
    AQ, IQ = TXADP_ADEQUATE_CODE, TXADP_INADEQUATE_CODE
    frame = pd.DataFrame([
        _row(1, "Kenya", 1, AQ,3, los=10),
        _row(2, "Kenya", 1, IQ,3, los=20),
        _row(3, "Kenya", 0, AQ,3, los=6),
        _row(4, "Kenya", 0, IQ,3, los=8),
    ])[_COLS]
    out = gformula_by_resistance(frame)
    assert out["resistant"]["avertable_bed_days"] == pytest.approx(10.0)
    assert out["susceptible"]["avertable_bed_days"] == pytest.approx(2.0)


def test_positivity_diagnostic_flags_off_support():
    AQ, IQ = TXADP_ADEQUATE_CODE, TXADP_INADEQUATE_CODE
    frame = pd.DataFrame([
        _row(1, "Kenya", 1, AQ,3, los=10),
        _row(2, "Kenya", 1, IQ,3, los=20),  # Kenya/high on support
        _row(3, "Ghana", 0, AQ,2, los=5),
        _row(4, "Ghana", 0, AQ,2, los=5),   # Ghana/low only adequate -> off support
    ])[_COLS]
    diag = positivity_diagnostic(frame)
    assert diag["n_strata"] == 2 and diag["n_strata_off_support"] == 1
    assert diag["frac_off_support"] == pytest.approx(0.5)   # 2 of 4 patients off-support


def test_cost_point_estimate_usd2010():
    cost = cost_of_bed_days({"Ghana": 100.0, "Kenya": 200.0}, currency="usd2010")
    assert cost["by_country"]["Ghana"] == pytest.approx(100.0 * BED_DAY_COST_USD_2010["Ghana"])
    assert cost["total"] == pytest.approx(
        100.0 * BED_DAY_COST_USD_2010["Ghana"] + 200.0 * BED_DAY_COST_USD_2010["Kenya"]
    )


def test_cost_unknown_country_raises():
    with pytest.raises(KeyError):
        cost_of_bed_days({"Nowhere": 10.0})


def test_scale_to_population():
    out = scale_to_population(2.5, {"Ghana": 40, "Kenya": 10})
    assert out == {"Ghana": 100.0, "Kenya": 25.0}


def test_monte_carlo_cost_reproducible_and_brackets_point():
    a = monte_carlo_cost({"Malawi": 1000.0}, draws=4000, seed=7)
    b = monte_carlo_cost({"Malawi": 1000.0}, draws=4000, seed=7)
    assert a["mean"] == b["mean"] and a["ci_lower"] == b["ci_lower"]   # seeded
    assert a["ci_lower"] < a["median"] < a["ci_upper"]
    # lognormal median = 1000 * PPP median (6.53) for Malawi
    assert a["median"] == pytest.approx(1000.0 * 6.53, rel=0.05)


def test_run_entrypoint_bundles_components():
    out = run_stewardship_gformula(_balanced_frame())
    assert set(out) == {"positivity", "pooled", "by_resistance"}
    assert out["pooled"]["avertable_bed_days"] == pytest.approx(10.0)


def test_adequate_inadequate_constants():
    assert (ADEQUATE, INADEQUATE) == (1, 0)


# --- Streamlit calibration artifact + scenario calculator -------------------

def test_synthetic_cohort_has_positivity_and_no_real_columns():
    df = synthetic_cohort(seed=1, per_cell=6)
    a = binarize_adequacy(df["treatment_adequacy"])
    assert (a == 1).sum() > 0 and (a == 0).sum() > 0
    # every (severity_bin, country) cell carries both arms -> g-formula identifies, no drops
    res = gformula_bed_days(df)
    assert res["strata_dropped"] == []
    # synthetic story matches reality: adequacy averts deaths (+) and adds bed-days
    assert res["averted_death_fraction"] > 0
    assert res["bed_days_set_adequate"] > res["bed_days_set_inadequate"]


def test_synthetic_cohort_is_deterministic():
    a = synthetic_cohort(seed=42, per_cell=5)
    b = synthetic_cohort(seed=42, per_cell=5)
    assert a.equals(b)


def test_build_calibration_artifact_is_deidentified_and_json_safe():
    art = build_calibration_artifact(synthetic_cohort(seed=2, per_cell=8), source="synthetic-demo")
    assert art["schema_version"] == 1
    assert art["provenance"]["source"] == "synthetic-demo"
    assert set(art["per_patient"]) == {"pooled", "resistant", "susceptible"}
    pp = art["per_patient"]["pooled"]
    assert pp["delta_bed_days_per_upgrade"] == pytest.approx(
        pp["bed_days_set_adequate"] - pp["bed_days_set_inadequate"]
    )
    # no patient rows / stratum cells / tuple keys leak -> must JSON round-trip
    import json
    assert json.loads(json.dumps(art)) == art


def test_build_calibration_artifact_suppresses_small_cells():
    # a stratum with a 1-patient arm is on-support but below min_cell_n -> counted suppressed
    AQ, IQ = TXADP_ADEQUATE_CODE, TXADP_INADEQUATE_CODE
    frame = pd.DataFrame([
        _row(1, "Kenya", 1, AQ, 3, los=10),
        _row(2, "Kenya", 1, IQ, 3, los=20),   # Kenya/high: 1 vs 1 -> below min_cell_n=5
    ] + [
        _row(10 + i, "Ghana", 0, AQ, 2, los=5) for i in range(6)
    ] + [
        _row(20 + i, "Ghana", 0, IQ, 2, los=9) for i in range(6)
    ])[_COLS]
    art = build_calibration_artifact(frame, min_cell_n=5)
    assert art["provenance"]["n_strata_supported"] == 2
    assert art["provenance"]["n_cells_below_min_n_suppressed"] == 1   # Kenya/high


def _demo_artifact():
    return build_calibration_artifact(
        synthetic_cohort(seed=7, per_cell=10), source="synthetic-demo")


def test_adequacy_scenario_directions_and_scaling():
    art = _demo_artifact()
    sc = adequacy_scenario(art, n_patients=200, current_adequacy=0.5, target_adequacy=0.8,
                           country="Kenya", arm="pooled")
    assert sc["patients_upgraded"] == pytest.approx(60.0)               # 0.3 * 200
    assert sc["averted_deaths"] > 0                                     # adequacy saves lives
    assert sc["added_bed_days"] > 0                                     # ...and adds bed-days
    assert sc["unit_cost_per_bed_day"] == pytest.approx(BED_DAY_COST_USD_2010["Kenya"])
    assert sc["added_cost"] == pytest.approx(sc["added_bed_days"] * sc["unit_cost_per_bed_day"])
    assert "ecological-calibration" in sc["firewall"].lower()


def test_adequacy_scenario_linear_in_n_and_gap():
    art = _demo_artifact()
    base = adequacy_scenario(art, n_patients=100, current_adequacy=0.4, target_adequacy=0.6,
                             country="Ghana")
    doubled = adequacy_scenario(art, n_patients=200, current_adequacy=0.4, target_adequacy=0.6,
                                country="Ghana")
    assert doubled["averted_deaths"] == pytest.approx(2 * base["averted_deaths"])
    assert doubled["added_bed_days"] == pytest.approx(2 * base["added_bed_days"])


def test_adequacy_scenario_target_below_current_is_zero():
    art = _demo_artifact()
    sc = adequacy_scenario(art, n_patients=500, current_adequacy=0.8, target_adequacy=0.5,
                           country="Malawi")
    assert sc["patients_upgraded"] == 0.0
    assert sc["averted_deaths"] == 0.0 and sc["added_bed_days"] == 0.0


def test_adequacy_scenario_validates_inputs():
    art = _demo_artifact()
    with pytest.raises(ValueError):
        adequacy_scenario(art, n_patients=10, current_adequacy=1.5, target_adequacy=0.5,
                          country="Kenya")
    with pytest.raises(KeyError):
        adequacy_scenario(art, n_patients=10, current_adequacy=0.1, target_adequacy=0.5,
                          country="Nowhere")
    with pytest.raises(ValueError):
        adequacy_scenario(art, n_patients=10, current_adequacy=0.1, target_adequacy=0.5,
                          country="Kenya", currency="euro")
