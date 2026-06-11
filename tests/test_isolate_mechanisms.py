"""Tests for the isolate-level mechanism breakdown + relinkage pre-staging (Component 1d).

Synthetic isolates only (no real SPIDAAR record).
"""

import pandas as pd
import pytest

from amr_sentinel_vivli.isolate_mechanisms import (
    MECHANISMS,
    mechanism_cooccurrence,
    mechanism_prevalence,
    relink_isolates_to_patients,
    relinkage_recovery_projection,
    run_isolate_mechanisms,
)

_ISO_COLS = ["iid", "country", "organism", "organism_group", "specimen",
             "infection_site", "clinical_relevance", "sample_month", "sample_year",
             "resistant", "c3r", "mdr", "mrsa", "amrtx_resistant"]


def _iso(iid, country, group, *, c3r=None, mdr=None, mrsa=None):
    return dict(iid=iid, country=country, organism="Escherichia coli", organism_group=group,
                specimen="blood", infection_site="BSI", clinical_relevance="P",
                sample_month=1, sample_year=2022, resistant=1,
                c3r=c3r, mdr=mdr, mrsa=mrsa, amrtx_resistant=1)


def _frame(rows):
    return pd.DataFrame(rows)[_ISO_COLS]


def test_mechanism_prevalence_overall_and_sentinel_exclusion():
    iso = _frame([
        _iso(1, "Ghana", 3, c3r=1), _iso(2, "Ghana", 3, c3r=1), _iso(3, "Kenya", 3, c3r=0),
        _iso(4, "Kenya", 3, c3r=9),    # sentinel -> excluded
        _iso(5, "Kenya", 3, c3r=None),  # NaN -> excluded
    ])
    out = mechanism_prevalence(iso, "c3r")
    assert out["n_ascertained"] == 3        # sentinels dropped from the denominator
    assert out["n_resistant"] == 2
    assert out["prevalence"] == pytest.approx(2 / 3)
    assert 0.0 <= out["ci_lower"] < out["prevalence"] < out["ci_upper"] <= 1.0


def test_mechanism_prevalence_by_stratifier():
    iso = _frame([
        _iso(1, "Ghana", 3, mdr=1), _iso(2, "Ghana", 3, mdr=0),
        _iso(3, "Kenya", 3, mdr=1), _iso(4, "Kenya", 3, mdr=1),
    ])
    out = mechanism_prevalence(iso, "mdr", by="country")
    assert out["strata"]["Ghana"]["prevalence"] == pytest.approx(0.5)
    assert out["strata"]["Kenya"]["prevalence"] == pytest.approx(1.0)
    assert out["strata"]["Kenya"]["n_ascertained"] == 2


def test_mechanism_prevalence_invalid_mechanism():
    with pytest.raises(ValueError, match="mechanism must be"):
        mechanism_prevalence(_frame([_iso(1, "Ghana", 3, c3r=1)]), "esbl")


def test_mechanism_cooccurrence_counts():
    iso = _frame([
        _iso(1, "Ghana", 3, c3r=1, mdr=1), _iso(2, "Ghana", 3, c3r=1, mdr=1),
        _iso(3, "Ghana", 3, c3r=1, mdr=0), _iso(4, "Ghana", 3, c3r=0, mdr=1),
        _iso(5, "Ghana", 3, c3r=0, mdr=0),
        _iso(6, "Ghana", 3, c3r=9, mdr=1),   # not doubly ascertained -> excluded
    ])
    co = mechanism_cooccurrence(iso, "c3r", "mdr")
    assert co["n_both_ascertained"] == 5
    assert co["counts"] == {"a1_b1": 2, "a1_b0": 1, "a0_b1": 1, "a0_b0": 1}
    assert co["p_b_given_a"] == pytest.approx(2 / 3)   # MDR among 3GC-R = 2 of 3
    assert co["p_a_given_b"] == pytest.approx(2 / 3)   # 3GC-R among MDR = 2 of 3


def test_relink_merges_outcome_onto_isolate():
    iso = _frame([_iso(10, "Ghana", 3, c3r=1), _iso(11, "Ghana", 3, c3r=0),
                  _iso(12, "Ghana", 3, c3r=1)])
    patients = pd.DataFrame({"pid": [100, 101], "mortality_30d": [1, 0], "age": [40, 55]})
    crosswalk = pd.DataFrame({"iid": [10, 11], "pid": [100, 101]})  # 12 has no link
    linked = relink_isolates_to_patients(iso, patients, crosswalk)
    assert len(linked) == 2                              # inner join drops the unlinked isolate
    assert set(linked["iid"]) == {10, 11}
    row10 = linked[linked["iid"] == 10].iloc[0]
    assert row10["pid"] == 100 and row10["mortality_30d"] == 1   # outcome attached


def test_relink_missing_key_raises():
    iso = _frame([_iso(1, "Ghana", 3, c3r=1)])
    patients = pd.DataFrame({"pid": [1], "mortality_30d": [0]})
    with pytest.raises(KeyError, match="pid"):
        relink_isolates_to_patients(iso, patients, pd.DataFrame({"iid": [1]}))


def test_relinkage_recovery_projection_arms_and_precision():
    # 4 ascertained MDR isolates: 2 resistant, 2 susceptible
    iso = _frame([_iso(1, "Ghana", 3, mdr=1), _iso(2, "Ghana", 3, mdr=1),
                  _iso(3, "Ghana", 3, mdr=0), _iso(4, "Ghana", 3, mdr=0)])
    proj = relinkage_recovery_projection(iso, patient_arms=(135, 21))
    mdr = proj["by_mechanism"]["mdr"]
    assert mdr["n_resistant"] == 2 and mdr["n_susceptible"] == 2
    assert mdr["susceptible_gain_vs_patient"] == 2 - 21
    base_var = 1 / 135 + 1 / 21
    iso_var = 1 / 2 + 1 / 2
    assert mdr["precision_gain_factor"] == pytest.approx((base_var / iso_var) ** 0.5)


def test_run_isolate_mechanisms_keys():
    iso = _frame([_iso(1, "Ghana", 3, c3r=1, mdr=1, mrsa=0),
                  _iso(2, "Kenya", 5, c3r=0, mdr=1, mrsa=1)])
    out = run_isolate_mechanisms(iso)
    assert set(out) == {"panel", "cooccurrence_c3r_mdr", "by_organism", "relinkage_projection"}
    assert set(out["panel"]) == set(MECHANISMS)
