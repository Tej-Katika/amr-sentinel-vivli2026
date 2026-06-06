"""Tests for the SPIDAAR isolate-level loader transforms.

Like the patient-level tests, these exercise the pure raw→analysis mapping on
**synthetic** rows only — the real SPIDAAR workbooks are Vivli DUA-restricted and
must never enter the repo. The synthetic frame mirrors the isolate codebook coding
(``rx``/``c3r`` in {0,1,9}, ``mdr``/``mrsa`` in {0,1,99}, ``amrtx`` in {-1,0,1,2},
``sampdat`` as ``MM/YYYY``) without any real isolate record.
"""

import pandas as pd
import pytest

from amr_sentinel_vivli import config
from amr_sentinel_vivli.data_loading import (
    ISOLATE_COLUMNS,
    _build_spidaar_isolate_frame,
    _restrict_to_catchment,
    _validate_spidaar_isolates,
)


def _raw(rows):
    """Build a raw isolatedata-'data'-shaped frame from dict rows."""
    cols = ["iid", "isolate", "group", "chaicat", "clinrel", "stype", "sampdat",
            "c3r", "mdr", "mrsa", "amrtx", "rx", "ctry"]
    return pd.DataFrame(rows)[cols]


def _synthetic():
    return _raw([
        # 3GC-R + MDR Gram-negative, resistant composite, treated-class resistant
        dict(iid=1, isolate="Escherichia coli", group=3, chaicat=1, clinrel="P",
             stype="blood", sampdat="03/2022", c3r=1, mdr=1, mrsa=99, amrtx=2, rx=1, ctry="Kenya"),
        # fully susceptible Gram-negative
        dict(iid=2, isolate="Klebsiella pneumoniae", group=3, chaicat=4, clinrel="Unclear",
             stype="urine", sampdat="10/2022", c3r=0, mdr=0, mrsa=99, amrtx=0, rx=0, ctry="Ghana"),
        # MRSA: mrsa positive, 3GC-R not applicable (sentinel 9), MDR positive
        dict(iid=3, isolate="Staphylococcus aureus", group=7, chaicat=3, clinrel="P",
             stype="pus", sampdat="06/2022", c3r=9, mdr=1, mrsa=1, amrtx=1, rx=1, ctry="Uganda"),
        # not-ascertained across the board (sentinels) + no-HAI site + bad date
        dict(iid=4, isolate="Acinetobacter spp", group=5, chaicat=10, clinrel="Unclear",
             stype="sputum", sampdat=None, c3r=9, mdr=99, mrsa=99, amrtx=-1, rx=9, ctry="Malawi"),
    ])


def test_returns_declared_columns_in_order():
    out = _build_spidaar_isolate_frame(_synthetic())
    assert list(out.columns) == list(ISOLATE_COLUMNS)
    assert len(out) == 4


def test_composite_resistance_exposure_from_rx():
    out = _build_spidaar_isolate_frame(_synthetic()).set_index("iid")
    assert out.loc[1, "resistant"] == 1            # rx 1
    assert out.loc[2, "resistant"] == 0            # rx 0
    assert pd.isna(out.loc[4, "resistant"])        # rx 9 -> not ascertained


def test_per_mechanism_flags_and_sentinels():
    out = _build_spidaar_isolate_frame(_synthetic()).set_index("iid")
    assert out.loc[1, "c3r"] == 1 and out.loc[1, "mdr"] == 1
    assert pd.isna(out.loc[1, "mrsa"])             # 99 -> not applicable
    assert pd.isna(out.loc[3, "c3r"])              # 9 -> not tested on Gram-positive
    assert out.loc[3, "mrsa"] == 1
    assert pd.isna(out.loc[4, "mdr"]) and pd.isna(out.loc[4, "mrsa"])


def test_amrtx_resistance_mapping():
    out = _build_spidaar_isolate_frame(_synthetic()).set_index("iid")
    assert out.loc[1, "amrtx_resistant"] == 1      # amrtx 2
    assert out.loc[2, "amrtx_resistant"] == 0      # amrtx 0
    assert out.loc[3, "amrtx_resistant"] == 0      # amrtx 1 (mixed S-untested)
    assert pd.isna(out.loc[4, "amrtx_resistant"])  # amrtx -1


def test_sample_date_parsing():
    out = _build_spidaar_isolate_frame(_synthetic()).set_index("iid")
    assert out.loc[1, "sample_month"] == 3 and out.loc[1, "sample_year"] == 2022
    assert out.loc[2, "sample_month"] == 10
    assert pd.isna(out.loc[4, "sample_month"])     # unparseable date -> NaT -> NA


def test_covariate_passthrough_and_decoding():
    out = _build_spidaar_isolate_frame(_synthetic()).set_index("iid")
    assert out.loc[1, "organism"] == "Escherichia coli"
    assert out.loc[1, "specimen"] == "blood"
    assert out.loc[1, "infection_site"] == "BSI"
    assert out.loc[1, "clinical_relevance"] == "P"
    assert out.loc[1, "organism_group"] == 3
    assert pd.isna(out.loc[4, "infection_site"])   # chaicat 10 = No HAI confirmed


def test_restrict_to_catchment_drops_out_of_scope():
    raw = _raw([
        dict(iid=1, isolate="x", group=3, chaicat=1, clinrel="P", stype="blood",
             sampdat="01/2022", c3r=1, mdr=1, mrsa=99, amrtx=2, rx=1, ctry="Kenya"),
        dict(iid=2, isolate="x", group=3, chaicat=1, clinrel="P", stype="blood",
             sampdat="01/2022", c3r=1, mdr=1, mrsa=99, amrtx=2, rx=1, ctry="Nigeria"),
    ])
    kept = _restrict_to_catchment(raw)
    assert list(kept["ctry"]) == ["Kenya"]
    assert set(kept["ctry"]) <= set(config.CATCHMENT_COUNTRIES)


def test_validate_passes_on_clean_frame():
    _validate_spidaar_isolates(_build_spidaar_isolate_frame(_synthetic()))


def test_validate_rejects_out_of_catchment():
    bad = _build_spidaar_isolate_frame(_synthetic())
    bad.loc[0, "country"] = "Nigeria"
    with pytest.raises(ValueError, match="catchment"):
        _validate_spidaar_isolates(bad)


def test_validate_rejects_nonbinary_flag():
    bad = _build_spidaar_isolate_frame(_synthetic())
    bad.loc[0, "resistant"] = 2
    with pytest.raises(ValueError, match="binary"):
        _validate_spidaar_isolates(bad)
