"""Tests for the ATLAS catchment nowcast + SPIDAAR frame-contrast (Component 3).

Synthetic ATLAS-tidy and SPIDAAR-isolate frames only (no real record).
"""

import numpy as np
import pandas as pd
import pytest

from amr_sentinel_vivli.bayesian_projection import (
    _eb_beta_prior,
    _wilson_ci,
    cell_counts,
    data_availability_matrix,
    frame_contrast,
    nowcast,
    project_forward,
    regional_trend,
    resistant_fraction_multiplier,
)


def _atlas(rows):
    cols = ["isolate_id", "species", "country", "year", "resistant"]
    return pd.DataFrame(rows)[cols]


def _make_atlas(per_country):
    """per_country: {country: (n_resistant, n_total, year)} -> a tidy ATLAS frame."""
    rows, iid = [], 0
    for country, (r, n, year) in per_country.items():
        for i in range(n):
            iid += 1
            rows.append(dict(isolate_id=iid, species="Escherichia coli", country=country,
                             year=year, resistant=1 if i < r else 0))
    return _atlas(rows)


def test_cell_counts_and_availability():
    df = _make_atlas({"Kenya": (6, 10, 2022), "Ghana": (2, 8, 2021)})
    cc = cell_counts(df).set_index("country")
    assert cc.loc["Kenya", "n_tested"] == 10 and cc.loc["Kenya", "n_resistant"] == 6
    assert cc.loc["Kenya", "prop"] == pytest.approx(0.6)
    mat = data_availability_matrix(df)
    assert mat.loc["Kenya", 2022] == 10 and mat.loc["Ghana", 2021] == 8


def test_eb_prior_mean_matches_pooled_rate():
    r = np.array([6, 2, 8])
    n = np.array([10, 8, 10])
    prior = _eb_beta_prior(r, n)
    assert prior["m"] == pytest.approx(r.sum() / n.sum())
    assert prior["a0"] / (prior["a0"] + prior["b0"]) == pytest.approx(prior["m"])
    assert prior["kappa"] >= 1.0


def test_nowcast_shrinks_between_raw_and_pooled():
    df = _make_atlas({"Kenya": (9, 10, 2022), "Ghana": (1, 10, 2022),
                      "Malawi": (5, 10, 2022), "Uganda": (5, 10, 2022)})
    nc = nowcast(df, n_mc=5000, seed=1)
    m = nc["pooled_rate"]
    for c in nc["cells"]:
        lo, hi = sorted([c["raw_prop"], m])
        assert lo - 1e-9 <= c["posterior_mean"] <= hi + 1e-9     # shrinks toward the pool
        assert c["ci_lower"] <= c["posterior_mean"] <= c["ci_upper"]
        assert c["no_pooling_mean"] == pytest.approx(
            (c["n_resistant"] + 0.5) / (c["n_tested"] + 1.0))


def test_nowcast_reproducible_and_since_year_filters():
    df = _make_atlas({"Kenya": (9, 10, 2014), "Ghana": (1, 10, 2022)})
    a = nowcast(df, since_year=2021, n_mc=3000, seed=5)
    b = nowcast(df, since_year=2021, n_mc=3000, seed=5)
    assert a["cells"][0]["ci_lower"] == b["cells"][0]["ci_lower"]
    assert {c["country"] for c in a["cells"]} == {"Ghana"}        # 2014 Kenya dropped


def test_regional_trend_detects_increase():
    df = pd.concat([
        _make_atlas({"Kenya": (2, 10, 2021)}),
        _make_atlas({"Kenya": (5, 10, 2022)}),
        _make_atlas({"Kenya": (8, 10, 2023)}),
    ], ignore_index=True)
    tr = regional_trend(df)
    assert tr["slope_logit_per_year"] > 0


def test_project_forward_flat_and_trend():
    df = _make_atlas({"Kenya": (6, 10, 2022), "Ghana": (4, 10, 2022)})
    nc = nowcast(df, n_mc=2000, seed=2)
    flat = project_forward(nc, scenario="flat")
    assert all(c["current_level"] == c["projected_level"] for c in flat["cells"])
    up = project_forward(nc, scenario="trend", horizon_years=3,
                         trend={"slope_logit_per_year": 0.5})
    assert all(c["projected_level"] > c["current_level"] for c in up["cells"])
    with pytest.raises(ValueError):
        project_forward(nc, scenario="bogus")


def test_resistant_fraction_multiplier_shape():
    df = _make_atlas({"Kenya": (6, 10, 2022)})
    mult = resistant_fraction_multiplier(nowcast(df, n_mc=2000, seed=3))
    assert "Kenya" in mult and "fraction" in mult["Kenya"] and "ci_lower" in mult["Kenya"]


def test_wilson_ci_known_and_empty():
    lo, hi = _wilson_ci(5, 10)
    assert 0.23 < lo < 0.27 and 0.73 < hi < 0.77
    assert all(np.isnan(x) for x in _wilson_ci(0, 0))


def _spidaar_iso(rows):
    cols = ["iid", "country", "organism", "c3r"]
    return pd.DataFrame(rows)[cols]


def test_frame_contrast_overall_and_filter():
    atlas = _make_atlas({"Kenya": (6, 10, 2022), "Ghana": (3, 10, 2022)})
    iso = _spidaar_iso([
        dict(iid=1, country="Kenya", organism="Escherichia coli", c3r=1),
        dict(iid=2, country="Kenya", organism="Escherichia coli", c3r=1),
        dict(iid=3, country="Ghana", organism="Klebsiella pneumoniae", c3r=0),
        dict(iid=4, country="Ghana", organism="Pseudomonas aeruginosa", c3r=0),  # filtered out
    ])
    fc = frame_contrast(atlas, iso, enterobacterales_only=True)
    assert fc["atlas"]["overall"]["prevalence"] == pytest.approx(9 / 20)
    # Pseudomonas dropped -> 3 Enterobacterales isolates, 2 resistant
    assert fc["spidaar"]["overall"]["n_tested"] == 3
    assert fc["spidaar"]["overall"]["prevalence"] == pytest.approx(2 / 3)
    assert "caveat" in fc and "ceftazidime" in fc["atlas_drug"]


def test_plot_frame_contrast_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    from amr_sentinel_vivli.projection_figures import plot_frame_contrast
    atlas = _make_atlas({"Kenya": (6, 10, 2022)})
    iso = _spidaar_iso([dict(iid=1, country="Kenya", organism="Escherichia coli", c3r=1)])
    out = plot_frame_contrast(frame_contrast(atlas, iso), tmp_path / "fc.png")
    assert out.exists() and out.stat().st_size > 0
