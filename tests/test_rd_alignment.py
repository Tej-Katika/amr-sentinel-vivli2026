"""Tests for the Cross-Domain R&D mismatch index (Component 4). Synthetic inputs only."""

import numpy as np
import pytest

from amr_sentinel_vivli import config
from amr_sentinel_vivli.rd_alignment import (
    GRAM_PANEL,
    RD_HUB_SNAPSHOT_2026,
    alignment_caption,
    analyze_alignment,
    cross_cutting_headline,
    cross_cutting_share,
    mismatch_index,
    monte_carlo_mismatch_ranking,
    spearman_burden_funding,
)


def test_cross_cutting_share_and_flag():
    out = cross_cutting_share({"A": 30, "B": 10}, cross_cutting_funding=60)
    assert out["pathogen_specific_funding"] == 40
    assert out["total_funding"] == 100
    assert out["cross_cutting_fraction"] == pytest.approx(0.6)
    assert out["flagged"] is True                       # majority is cross-cutting
    assert cross_cutting_share({"A": 80}, 20)["flagged"] is False


def test_mismatch_index_values_and_ranking():
    out = mismatch_index({"A": 60, "B": 40}, {"A": 30, "B": 70})
    by = {r["pathogen"]: r for r in out["ranking"]}
    assert by["A"]["mismatch"] == pytest.approx(0.6 / 0.3)      # 2.0
    assert by["A"]["log2_mismatch"] == pytest.approx(1.0)
    assert by["B"]["log2_mismatch"] == pytest.approx(np.log2(0.4 / 0.7))
    assert out["ranking"][0]["pathogen"] == "A"                 # most under-funded first


def test_mismatch_index_floor_bounds_small_denominator():
    burden = {"A": 50, "B": 50}
    funding = {"A": 1, "B": 99}        # A is barely funded -> A's share 0.01
    no_floor = mismatch_index(burden, funding)
    with_floor = mismatch_index(burden, funding, floor=0.1)
    a_no = next(r for r in no_floor["ranking"] if r["pathogen"] == "A")["mismatch"]
    a_fl = next(r for r in with_floor["ranking"] if r["pathogen"] == "A")["mismatch"]
    assert a_fl < a_no                                          # floor tempers the blow-up


def test_mismatch_index_requires_two_pathogens():
    with pytest.raises(ValueError):
        mismatch_index({"A": 1}, {"A": 1})


def test_spearman_perfect_negative_and_caveat():
    out = spearman_burden_funding({"A": 3, "B": 2, "C": 1}, {"A": 1, "B": 2, "C": 3},
                                  n_boot=200, seed=1)
    assert out["spearman_rho"] == pytest.approx(-1.0)
    assert out["n_pathogens"] == 3 and "n=3" in out["caveat"]


def test_spearman_reproducible():
    a = spearman_burden_funding({"A": 1, "B": 2, "C": 3, "D": 5},
                                {"A": 2, "B": 1, "C": 4, "D": 3}, n_boot=300, seed=7)
    b = spearman_burden_funding({"A": 1, "B": 2, "C": 3, "D": 5},
                                {"A": 2, "B": 1, "C": 4, "D": 3}, n_boot=300, seed=7)
    assert a["ci_lower"] == b["ci_lower"] and a["ci_upper"] == b["ci_upper"]


def test_monte_carlo_ranking_probabilities_sum_to_one():
    ui = {"A": (60, 50, 72), "B": (40, 33, 48)}
    funding = {"A": 30, "B": 70}
    mc = monte_carlo_mismatch_ranking(ui, funding, draws=500, seed=2)
    probs = [v["p_most_underfunded"] for v in mc["per_pathogen"].values()]
    assert sum(probs) == pytest.approx(1.0)
    assert all("log2_mismatch_median" in v for v in mc["per_pathogen"].values())


def test_analyze_alignment_gated_on_snapshot(monkeypatch):
    burden = {"A": 60, "B": 40}
    funding = {"A": 30, "B": 70}
    monkeypatch.setattr(config, "RD_HUB_SNAPSHOT_DATE", None)
    with pytest.raises(ValueError, match="RD_HUB_SNAPSHOT_DATE"):
        analyze_alignment(burden, funding)

    monkeypatch.setattr(config, "RD_HUB_SNAPSHOT_DATE", "2026-06-08")
    out = analyze_alignment(burden, funding, cross_cutting_funding=20, floor=0.05)
    assert set(out) == {"snapshot_caption", "cross_cutting", "ranking_no_floor",
                        "ranking_with_floor", "spearman"}
    assert "2026-06-08" in out["snapshot_caption"]


def test_gram_panel_documents_verified_facts():
    assert len(GRAM_PANEL["pathogens"]) == 6
    assert GRAM_PANEL["combined_attributable_deaths_2019"] == 929_000
    assert "PMC8841637" in GRAM_PANEL["source"]


def test_alignment_caption_carries_scope(monkeypatch):
    monkeypatch.setattr(config, "RD_HUB_SNAPSHOT_DATE", "2026-06-08")
    cap = alignment_caption()
    assert "2026-06-08" in cap and "public + philanthropic" in cap


def test_rd_hub_snapshot_named_funding_sums_to_species_specific():
    # Named top-five + the folded "other" bucket must reconstruct the species-specific
    # total ($1058M); guards the verified figures against silent edits.
    named = sum(RD_HUB_SNAPSHOT_2026["named_species_funding_musd"].values())
    total = named + RD_HUB_SNAPSHOT_2026["other_species_specific_musd"]
    assert total == pytest.approx(RD_HUB_SNAPSHOT_2026["species_specific_total_musd"])
    assert "Czaplewski" in RD_HUB_SNAPSHOT_2026["source"]


def test_cross_cutting_headline_is_flagged_majority():
    # The de-gated headline: pathogen-specific R&D is a MINORITY of the $2.51B total
    # ($1058M / $2510M = 42%), so cross-cutting is the majority and the flag fires.
    head = cross_cutting_headline()
    assert head["pathogen_specific_funding"] == pytest.approx(1058.0)
    assert head["total_funding"] == pytest.approx(2510.0)
    assert head["cross_cutting_fraction"] == pytest.approx(1452.0 / 2510.0, rel=1e-6)
    assert head["flagged"] is True
    assert "Czaplewski" in head["source"]
