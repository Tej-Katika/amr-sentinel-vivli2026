"""Tests for the SPIDAAR patient-level loader transforms.

These exercise the pure raw→analysis mapping on **synthetic** rows only — the
real SPIDAAR workbooks are Vivli DUA-restricted and must never enter the repo.
The synthetic frame mirrors the codebook coding (e.g. ``dead`` in {0,1,9},
``amrp`` in {-1,0,1,2}, ``agegr`` 0..14) without any real patient record.
"""

import pandas as pd
import pytest

from amr_sentinel_vivli import config
from amr_sentinel_vivli.data_loading import (
    MORTALITY_HORIZON_DAYS,
    SPIDAAR_COLUMNS,
    _build_spidaar_analysis_frame,
    _restrict_to_catchment,
    _validate_spidaar,
)


def _raw(rows):
    """Build a raw patientdata-'data'-shaped frame from dict rows."""
    cols = ["pid", "ctry", "agegr", "sex", "chaicat", "isol", "amrp", "dead", "dtpta", "nobsd"]
    return pd.DataFrame(rows)[cols]


def _synthetic():
    return _raw([
        # death within 30d -> event, time = dtpta
        dict(pid=1, ctry="Kenya", agegr=14, sex=0, chaicat=1, isol="Escherichia coli",
             amrp=2, dead=1, dtpta=8, nobsd=8),
        # death after 30d -> survivor at horizon (event 0, time 30)
        dict(pid=2, ctry="Ghana", agegr=6, sex=1, chaicat=4, isol="Klebsiella pneumoniae",
             amrp=0, dead=1, dtpta=40, nobsd=40),
        # alive, observed past horizon -> event 0, time 30
        dict(pid=3, ctry="Uganda", agegr=0, sex=0, chaicat=2, isol="Pseudomonas aeruginosa",
             amrp=1, dead=0, dtpta=None, nobsd=45),
        # alive, censored early (obs < horizon) -> event 0, time = nobsd
        dict(pid=4, ctry="Malawi", agegr=9, sex=3, chaicat=10, isol=None,
             amrp=-1, dead=0, dtpta=None, nobsd=12),
        # deceased-censored (dead==9) within 30d counts as all-cause event
        dict(pid=5, ctry="Kenya", agegr=2, sex=1, chaicat=3, isol="Staphylococcus aureus",
             amrp=2, dead=9, dtpta=6, nobsd=30),
    ])


def test_returns_declared_columns_in_order():
    out = _build_spidaar_analysis_frame(_synthetic())
    assert list(out.columns) == list(SPIDAAR_COLUMNS)
    assert len(out) == 5


def test_mortality_event_within_horizon():
    out = _build_spidaar_analysis_frame(_synthetic()).set_index("pid")
    assert out.loc[1, "mortality_30d"] == 1          # died day 8
    assert out.loc[1, "time_at_risk"] == 8
    assert out.loc[5, "mortality_30d"] == 1          # dead==9, died day 6, all-cause


def test_death_after_horizon_is_censored_survivor():
    out = _build_spidaar_analysis_frame(_synthetic()).set_index("pid")
    assert out.loc[2, "mortality_30d"] == 0
    assert out.loc[2, "time_at_risk"] == MORTALITY_HORIZON_DAYS


def test_survivor_time_is_capped_and_early_censor_kept():
    out = _build_spidaar_analysis_frame(_synthetic()).set_index("pid")
    assert out.loc[3, "mortality_30d"] == 0
    assert out.loc[3, "time_at_risk"] == MORTALITY_HORIZON_DAYS   # min(45, 30)
    assert out.loc[4, "time_at_risk"] == 12                       # min(12, 30)


def test_resistance_exposure_mapping():
    out = _build_spidaar_analysis_frame(_synthetic()).set_index("pid")
    assert out.loc[1, "resistant"] == 1              # amrp 2
    assert out.loc[2, "resistant"] == 0              # amrp 0
    assert out.loc[3, "resistant"] == 0              # amrp 1 (mixed S-untested)
    assert pd.isna(out.loc[4, "resistant"])          # amrp -1 unascertainable


def test_covariate_decoding():
    out = _build_spidaar_analysis_frame(_synthetic()).set_index("pid")
    assert out.loc[1, "infection_site"] == "BSI"
    assert out.loc[1, "sex"] == "Male"
    assert out.loc[2, "sex"] == "Female"
    assert pd.isna(out.loc[4, "sex"])                # code 3 = Missing
    assert out.loc[1, "age_group"] == "65+"
    assert out.loc[1, "age"] == 14
    assert pd.isna(out.loc[4, "infection_site"])     # chaicat 10 = No HAI confirmed


def test_restrict_to_catchment_drops_out_of_scope():
    raw = _raw([
        dict(pid=1, ctry="Kenya", agegr=1, sex=0, chaicat=1, isol="x",
             amrp=0, dead=0, dtpta=None, nobsd=10),
        dict(pid=2, ctry="Nigeria", agegr=1, sex=0, chaicat=1, isol="x",
             amrp=0, dead=0, dtpta=None, nobsd=10),
    ])
    kept = _restrict_to_catchment(raw)
    assert list(kept["ctry"]) == ["Kenya"]
    assert set(kept["ctry"]) <= set(config.CATCHMENT_COUNTRIES)


def test_validate_passes_on_clean_frame():
    _validate_spidaar(_build_spidaar_analysis_frame(_synthetic()))


def test_validate_rejects_out_of_catchment():
    bad = _build_spidaar_analysis_frame(_synthetic())
    bad.loc[0, "country"] = "Nigeria"
    with pytest.raises(ValueError, match="catchment"):
        _validate_spidaar(bad)
