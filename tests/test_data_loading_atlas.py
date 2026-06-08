"""Tests for the ATLAS loader's pure raw->panel mapping (no real ATLAS record)."""

import pandas as pd
import pytest

from amr_sentinel_vivli.data_loading import ATLAS_COLUMNS, _build_atlas_frame, _validate_atlas


def _raw(rows):
    cols = ["Isolate Id", "Species", "Country", "Year", "Source", "Speciality",
            "Ceftazidime", "Ceftazidime_I"]
    return pd.DataFrame(rows)[cols]


def _row(iid, species, country, year, interp, mic="8"):
    return dict(**{"Isolate Id": iid, "Species": species, "Country": country, "Year": year,
                   "Source": "Blood", "Speciality": "Medicine", "Ceftazidime": mic,
                   "Ceftazidime_I": interp})


def test_build_atlas_frame_columns_and_resistance_mapping():
    raw = _raw([
        _row(1, "Escherichia coli", "Kenya", 2022, "Resistant", mic=">128"),
        _row(2, "Klebsiella pneumoniae", "Ghana", 2021, "Susceptible", mic="<=1"),
        _row(3, "Escherichia coli", "Uganda", 2023, "Intermediate", mic="4"),
        _row(4, "Escherichia coli", "Malawi", 2022, None),          # not ascertained
    ])
    out = _build_atlas_frame(raw)
    assert list(out.columns) == list(ATLAS_COLUMNS)
    assert out.loc[0, "resistant"] == 1                              # Resistant -> 1
    assert out.loc[1, "resistant"] == 0                              # Susceptible -> 0
    assert out.loc[2, "resistant"] == 0                              # Intermediate -> non-resistant
    assert pd.isna(out.loc[3, "resistant"])                         # missing -> NaN
    assert out.loc[0, "mic"] == ">128" and out.loc[0, "drug"] == "Ceftazidime"


def test_validate_atlas_rejects_non_catchment():
    raw = _raw([_row(1, "Escherichia coli", "Narnia", 2022, "Resistant")])
    frame = _build_atlas_frame(raw)
    with pytest.raises(ValueError, match="outside the catchment"):
        _validate_atlas(frame, catchment_only=True)
    _validate_atlas(frame, catchment_only=False)                    # allowed when not restricting


def test_validate_atlas_rejects_bad_year():
    raw = _raw([_row(1, "Escherichia coli", "Kenya", 1850, "Resistant")])
    with pytest.raises(ValueError, match="year range"):
        _validate_atlas(_build_atlas_frame(raw), catchment_only=True)
