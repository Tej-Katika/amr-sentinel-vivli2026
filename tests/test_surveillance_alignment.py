"""Tests for Component 6: the surveillance blind-spot axis (surveillance_alignment).

Pure functions exercised on synthetic ATLAS-shaped frames — no real-file dependency.
"""

from __future__ import annotations

import pandas as pd
import pytest

from amr_sentinel_vivli.rd_alignment import GRAM_BURDEN_2019
from amr_sentinel_vivli.surveillance_alignment import (
    _gini,
    geographic_concentration,
    panel_surveillance_counts,
    surveillance_burden_mismatch,
    three_axis_alignment,
)

_SPECIES = list(GRAM_BURDEN_2019)


def _atlas(rows: list[tuple[str, str]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["Species", "Country"])


def test_panel_counts_map_species_case_insensitive_and_genus():
    atlas = _atlas([
        ("Escherichia coli", "Kenya"),
        ("escherichia coli", "Kenya"),          # case-insensitive
        ("Klebsiella pneumoniae", "Ghana"),
        ("Klebsiella oxytoca", "Ghana"),        # NOT K. pneumoniae -> excluded
        ("Acinetobacter baumannii", "Uganda"),
        ("Acinetobacter spp.", "Uganda"),       # genus-level aggregation -> A. baumannii
        ("Staphylococcus aureus", "United States"),
        ("Streptococcus pneumoniae", "United States"),
        ("Pseudomonas aeruginosa", "United States"),
        ("Candida albicans", "United States"),  # non-panel -> ignored
    ])
    c = panel_surveillance_counts(atlas)
    assert c["Escherichia coli"] == 2
    assert c["Klebsiella pneumoniae"] == 1       # oxytoca excluded
    assert c["Acinetobacter baumannii"] == 2     # genus-aggregated
    assert c["Staphylococcus aureus"] == 1
    assert set(c) == set(_SPECIES)               # exactly the six panel species


def test_panel_counts_country_filter():
    atlas = _atlas([("Escherichia coli", "Kenya"), ("Escherichia coli", "United States")])
    c = panel_surveillance_counts(atlas, countries=frozenset({"Kenya"}))
    assert c["Escherichia coli"] == 1


def test_surveillance_mismatch_flags_undersurveilled_high_burden():
    # Give the highest-burden species (E. coli, assoc deaths) very few isolates and a
    # low-burden species (P. aeruginosa) many -> E. coli should be under-surveilled (>0),
    # P. aeruginosa over-surveilled (<0), and the ranking sorted by descending log2.
    counts = {p: 1000 + 10 * i for i, p in enumerate(_SPECIES)}  # distinct baselines
    counts["Escherichia coli"] = 10
    counts["Pseudomonas aeruginosa"] = 100000
    out = surveillance_burden_mismatch(counts=counts, draws=2000, seed=4)
    pp = out["per_pathogen"]
    assert pp["Escherichia coli"]["log2_mismatch_median"] > 0
    assert pp["Pseudomonas aeruginosa"]["log2_mismatch_median"] < 0
    meds = [pp[p]["log2_mismatch_median"] for p in out["undersurveilled_ranking"]]
    assert meds == sorted(meds, reverse=True)
    assert out["undersurveilled_ranking"][0] == "Escherichia coli"


def test_surveillance_mismatch_reproducible_and_rejects_zero():
    counts = {p: 100 * (i + 1) for i, p in enumerate(_SPECIES)}  # non-degenerate
    a = surveillance_burden_mismatch(counts=counts, draws=1500, seed=7)
    b = surveillance_burden_mismatch(counts=counts, draws=1500, seed=7)
    assert a["per_pathogen"] == b["per_pathogen"]
    bad = dict(counts, **{"Streptococcus pneumoniae": 0})
    with pytest.raises(ValueError):
        surveillance_burden_mismatch(counts=bad)


def test_three_axis_shares_sum_to_one_and_neglect_flag():
    # Starve S. pneumoniae of surveillance; its funding share is already the smallest in the
    # scoped Hub extract, so it should flag as neglected on all axes (high burden, low both).
    counts = dict.fromkeys(_SPECIES, 5000)
    counts["Streptococcus pneumoniae"] = 50
    out = three_axis_alignment(counts=counts)
    pp = out["per_pathogen"]
    for axis in ("burden_share", "funding_share", "surveillance_share"):
        assert sum(pp[p][axis] for p in _SPECIES) == pytest.approx(1.0)
    assert "Streptococcus pneumoniae" in out["neglected_on_all_axes"]
    assert set(out["neglected_on_all_axes"]) <= set(_SPECIES)


def test_geographic_concentration_ssa_blind_spot():
    rows = ([("Escherichia coli", "United States")] * 70
            + [("Escherichia coli", "Kenya")] * 2
            + [("Escherichia coli", "Nigeria")] * 3
            + [("Escherichia coli", "France")] * 25)
    g = geographic_concentration(_atlas(rows))
    assert g["total_isolates"] == 100
    assert g["ssa"]["isolates"] == 5                 # Kenya + Nigeria
    assert g["ssa"]["share"] == pytest.approx(0.05)
    assert g["catchment"]["isolates"] == 2           # Kenya only (of Gh/Ke/Ma/Ug)
    assert g["top_country"]["country"] == "United States"
    assert g["top_country_to_ssa_ratio"] == pytest.approx(70 / 5)
    assert 0.0 <= g["gini_country_concentration"] <= 1.0


def test_gini_bounds():
    assert _gini([1, 1, 1, 1]) == pytest.approx(0.0, abs=1e-9)   # perfectly even
    assert _gini([0, 0, 0, 100]) > 0.6                            # concentrated
