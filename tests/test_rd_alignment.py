"""Tests for the Cross-Domain R&D mismatch index (Component 4). Synthetic inputs only."""

import numpy as np
import pandas as pd
import pytest

from amr_sentinel_vivli import config
from amr_sentinel_vivli.rd_alignment import (
    GRAM_BURDEN_2019,
    GRAM_PANEL,
    RD_HUB_SNAPSHOT_2026,
    alignment_caption,
    analyze_alignment,
    catchment_alignment,
    catchment_pathogen_counts,
    cross_cutting_headline,
    cross_cutting_share,
    gram_panel_alignment,
    mismatch_index,
    monte_carlo_mismatch_ranking,
    reconcile_hub_snapshot,
    snapshot_from_hub_export,
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


def _hub_export_frame():
    """Tidy frame shaped like data_loading.load_rd_hub_snapshot output (synthetic)."""
    return pd.DataFrame([
        # species-specific: E.coli 15M (two slices), K.pneu 4M
        dict(pathogen="Escherichia coli", year=2018, funder_type="public",
             investment_usd=10e6, is_cross_cutting=False, funder="A"),
        dict(pathogen="Escherichia coli", year=2020, funder_type="philanthropic",
             investment_usd=5e6, is_cross_cutting=False, funder="B"),
        dict(pathogen="Klebsiella pneumoniae", year=2021, funder_type="public",
             investment_usd=4e6, is_cross_cutting=False, funder="A"),
        # cross-cutting: 20M -> total only
        dict(pathogen="Cross-cutting", year=2021, funder_type="public",
             investment_usd=20e6, is_cross_cutting=True, funder="C"),
    ])


def test_snapshot_from_hub_export_shape_and_sums():
    snap = snapshot_from_hub_export(_hub_export_frame())
    assert snap["total_funding_musd"] == pytest.approx(39.0)        # (15+4+20)
    assert snap["species_specific_total_musd"] == pytest.approx(19.0)  # 15+4
    assert snap["named_species_funding_musd"]["Escherichia coli"] == pytest.approx(15.0)
    assert snap["named_species_funding_musd"]["Klebsiella pneumoniae"] == pytest.approx(4.0)
    # cross-cutting label is never folded onto a pathogen
    assert "Cross-cutting" not in snap["named_species_funding_musd"]
    assert snap["n_funders"] == 3
    # the export carries every species-specific pathogen, so no residual "other"
    assert snap["other_species_specific_musd"] == 0.0
    # GRAM-panel species absent from the export are reported as unfunded
    assert "Staphylococcus aureus" in snap["unfunded_gram_panel_species"]
    assert "Escherichia coli" not in snap["unfunded_gram_panel_species"]


def test_snapshot_from_hub_export_drops_into_headline():
    # The export-derived snapshot is a drop-in for cross_cutting_headline.
    snap = snapshot_from_hub_export(_hub_export_frame())
    head = cross_cutting_headline(snapshot=snap)
    assert head["total_funding"] == pytest.approx(39.0)
    assert head["pathogen_specific_funding"] == pytest.approx(19.0)
    assert head["cross_cutting_fraction"] == pytest.approx(20.0 / 39.0)
    assert head["flagged"] is True


def test_reconcile_hub_snapshot_flags_drift():
    # Identical-to-published reconciles ok; a 50% inflated total breaches the 15% tolerance.
    ok = reconcile_hub_snapshot(dict(RD_HUB_SNAPSHOT_2026))
    assert ok["ok"] is True
    assert ok["fields"]["total_funding_musd"]["rel_diff"] == pytest.approx(0.0)

    drifted = dict(RD_HUB_SNAPSHOT_2026)
    drifted["total_funding_musd"] = RD_HUB_SNAPSHOT_2026["total_funding_musd"] * 1.5
    rec = reconcile_hub_snapshot(drifted)
    assert rec["ok"] is False
    assert rec["fields"]["total_funding_musd"]["within_tol"] is False
    assert rec["fields"]["total_funding_musd"]["rel_diff"] == pytest.approx(0.5)


def test_gram_burden_self_consistency_with_headline_totals():
    # Verified Table S22 medians must reconstruct the GRAM headline totals (a strong
    # check that the per-pathogen rows were transcribed correctly): six leading pathogens
    # sum to 3.57M associated and ~929k attributable deaths (counts in thousands).
    assoc = sum(GRAM_BURDEN_2019[p]["assoc_deaths_k"][0] for p in GRAM_BURDEN_2019)
    attrib = sum(GRAM_BURDEN_2019[p]["attrib_deaths_k"][0] for p in GRAM_BURDEN_2019)
    assert assoc == pytest.approx(3572, abs=1)        # 3.57 million
    assert attrib == pytest.approx(929, abs=1)        # 929,000
    # UIs are ordered lo <= median <= hi for every metric.
    for rec in GRAM_BURDEN_2019.values():
        for med, lo, hi in rec.values():
            assert lo <= med <= hi


def test_gram_burden_associated_rank_matches_published_order():
    order = sorted(GRAM_BURDEN_2019, key=lambda p: GRAM_BURDEN_2019[p]["assoc_deaths_k"][0],
                   reverse=True)
    assert order == ["Escherichia coli", "Staphylococcus aureus", "Klebsiella pneumoniae",
                     "Streptococcus pneumoniae", "Acinetobacter baumannii",
                     "Pseudomonas aeruginosa"]


def test_gram_panel_alignment_robust_ranking_and_determinism():
    a = gram_panel_alignment(draws=3000, seed=4)
    b = gram_panel_alignment(draws=3000, seed=4)
    assert a["per_pathogen"] == b["per_pathogen"]          # reproducible
    assert a["cross_cutting"]["flagged"] is True
    # S. pneumoniae and K. pneumoniae carry essentially all the "most under-funded" mass;
    # S. aureus and P. aeruginosa are over-funded (negative log2 mismatch). This ranking
    # is invariant to the unfetched funding split — the point of propagating it.
    top2 = set(a["underfunded_ranking"][:2])
    assert top2 == {"Streptococcus pneumoniae", "Klebsiella pneumoniae"}
    p_top = a["per_pathogen"]
    assert p_top["Streptococcus pneumoniae"]["p_most_underfunded"] \
        + p_top["Klebsiella pneumoniae"]["p_most_underfunded"] > 0.9
    assert p_top["Pseudomonas aeruginosa"]["log2_mismatch_median"] < 0
    assert p_top["Staphylococcus aureus"]["log2_mismatch_median"] < 0


def test_gram_panel_alignment_gated_on_snapshot(monkeypatch):
    monkeypatch.setattr(config, "RD_HUB_SNAPSHOT_DATE", None)
    with pytest.raises(ValueError, match="RD_HUB_SNAPSHOT_DATE"):
        gram_panel_alignment(draws=100)


# --- Catchment-specific alignment -------------------------------------------------

def _iso(organism, resistant):
    return {"organism": organism, "resistant": resistant}


def test_catchment_pathogen_counts_mapping():
    iso = pd.DataFrame([
        _iso("Escherichia coli", 1), _iso("escherichia coli", 0),
        _iso("Klebsiella pneumoniae", 1), _iso("Klebsiella oxytoca", 1),   # oxytoca NOT K.pneu
        _iso("Acinetobacter spp", 1), _iso("Acinetobacter baumannii", 1),  # both -> A.baumannii
        _iso("Staphylococcus aureus", 0),
        _iso("Streptococcus pneumoniae", 1),
        _iso("Pseudomonas aeruginosa", 1),
    ])
    c = catchment_pathogen_counts(iso)
    assert c["Escherichia coli"] == {"isolates": 2, "resistant": 1}
    assert c["Klebsiella pneumoniae"] == {"isolates": 1, "resistant": 1}   # oxytoca excluded
    assert c["Acinetobacter baumannii"] == {"isolates": 2, "resistant": 2}  # genus-aggregated
    assert c["Staphylococcus aureus"] == {"isolates": 1, "resistant": 0}
    assert c["Streptococcus pneumoniae"]["isolates"] == 1


def test_catchment_alignment_shares_and_ranking():
    # High-burden, zero-named-funding pathogen (K. pneumoniae) should top the ranking;
    # a low-burden, high-funding one (S. aureus) should sit at the bottom.
    counts = {
        "Escherichia coli": {"isolates": 50, "resistant": 40},
        "Staphylococcus aureus": {"isolates": 14, "resistant": 11},
        "Klebsiella pneumoniae": {"isolates": 39, "resistant": 34},
        "Streptococcus pneumoniae": {"isolates": 3, "resistant": 2},
        "Acinetobacter baumannii": {"isolates": 22, "resistant": 16},
        "Pseudomonas aeruginosa": {"isolates": 16, "resistant": 10},
    }
    out = catchment_alignment(counts=counts, weight="resistant", draws=3000, seed=5)
    shares = out["burden_shares"]
    assert sum(shares.values()) == pytest.approx(1.0)
    # P(most under-funded) is a proper distribution over the panel
    assert sum(p["p_most_underfunded"] for p in out["per_pathogen"].values()) == pytest.approx(1.0)
    # ranking sorted by descending median log2 mismatch
    meds = [out["per_pathogen"][p]["log2_mismatch_median"] for p in out["underfunded_ranking"]]
    assert meds == sorted(meds, reverse=True)
    # the highest-funded species is over-funded relative to its small catchment burden
    assert out["per_pathogen"]["Staphylococcus aureus"]["log2_mismatch_median"] < 0
    # E. coli / K. pneumoniae (the dominant catchment Gram-negatives) are under-funded
    assert out["per_pathogen"]["Klebsiella pneumoniae"]["log2_mismatch_median"] > 0
    assert out["per_pathogen"]["Escherichia coli"]["log2_mismatch_median"] > 0


def test_catchment_alignment_reproducible_and_weight_validated():
    counts = {p: {"isolates": 5 + 3 * i, "resistant": 4 + 2 * i}
              for i, p in enumerate(GRAM_BURDEN_2019)}
    a = catchment_alignment(counts=counts, draws=1000, seed=5)
    b = catchment_alignment(counts=counts, draws=1000, seed=5)
    assert a["per_pathogen"] == b["per_pathogen"]
    with pytest.raises(ValueError, match="weight"):
        catchment_alignment(counts=counts, weight="bogus", draws=10)


def test_catchment_alignment_gated_on_snapshot(monkeypatch):
    counts = {p: {"isolates": 5, "resistant": 4} for p in GRAM_BURDEN_2019}
    monkeypatch.setattr(config, "RD_HUB_SNAPSHOT_DATE", None)
    with pytest.raises(ValueError, match="RD_HUB_SNAPSHOT_DATE"):
        catchment_alignment(counts=counts, draws=100)
