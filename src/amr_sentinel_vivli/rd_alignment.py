"""Component 4 — Cross-Domain R&D mismatch index (descriptive; pre-reg H3, reframed).

A single descriptive **Axis-A** index aligning global AMR burden with global public +
philanthropic R&D funding, per pathogen:

    M_p = (global GRAM burden share of pathogen p) / (global Hub funding share of p),
    reported as log2(M_p)  — >0 means under-funded relative to burden.

plus the pre-registered Spearman rank correlation (burden vs funding) with a bootstrap CI.
Descriptive / ecological only — n = 5-6 pathogens, so NO line is fitted and no causal or
powered claim is made (``docs/analysis_plan_2026.md`` Component 4). This contests the
Global AMR R&D Hub **Cross-Domain Award** (Hub dataset + a Vivli Register dataset).

Two hard honesty rules baked in here:
1. **Cross-cutting funding is reported FIRST.** Much Hub-tracked R&D is not pathogen-
   specific; ``cross_cutting_share`` is the headline magnitude and the index is gated on it
   (if pathogen-specific spend is a minority, every downstream caveat widens). We never
   silently redistribute cross-cutting funding onto pathogens.
2. **Low-denominator floor.** A pathogen with near-zero tracked funding makes M_p explode;
   a pre-specified floor on the funding share bounds it, and the ranking is reported BOTH
   with and without the floor.

Data status (both GATED, not shipped as final):
- **GRAM burden** (numerator): the verified panel and combined totals are documented in
  ``GRAM_PANEL`` (Murray et al., Lancet 2022). Per-pathogen DALYs (associated + attributable,
  both reported, DALYs primary) must be lifted from the GRAM appendix tables in the secure
  workflow and passed in — they are NOT hard-coded here to avoid quoting unverified figures.
- **Hub funding** (denominator): requires a FROZEN, dated Global AMR R&D Hub snapshot with
  the exact filter state recorded; ``config.RD_HUB_SNAPSHOT_DATE`` must be locked first
  (enforced by ``data_loading.load_rd_hub_snapshot``). Title everything "public+philanthropic".

The functions below are pure and metric-agnostic (burden can be DALYs or deaths); they are
unit-tested on synthetic inputs.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config

# Verified facts (Murray et al., "Global burden of bacterial antimicrobial resistance in
# 2019", Lancet 2022; PMC8841637). Used for documentation/sanity bounds — NOT per-pathogen
# values (those are appendix-sourced and supplied at call time).
GRAM_PANEL: dict = {
    "pathogens": [
        "Escherichia coli", "Staphylococcus aureus", "Klebsiella pneumoniae",
        "Streptococcus pneumoniae", "Acinetobacter baumannii", "Pseudomonas aeruginosa",
    ],  # the six leading pathogens, ordered by deaths associated with resistance
    "combined_attributable_deaths_2019": 929_000,   # (660 000-1 270 000)
    "combined_associated_deaths_2019": 3_570_000,    # (2.62-4.78 million)
    "all_pathogens_total_dalys_2019": 47_900_000,    # across all 23 pathogens
    # Verified rank orderings (Murray 2022 main text, p638). Exact per-pathogen
    # death/DALY MAGNITUDES live in the GRAM appendix / figure 4 (not the main text)
    # and are supplied at call time; do NOT substitute the all-cause bacterial-
    # infection deaths from Ikuta et al. 2022 ("one in eight deaths") — a different
    # paper and a common conflation (see docs/reference_rd_alignment_2026-06-09.md).
    "rank_by_associated_deaths": [
        "Escherichia coli", "Staphylococcus aureus", "Klebsiella pneumoniae",
        "Streptococcus pneumoniae", "Acinetobacter baumannii", "Pseudomonas aeruginosa",
    ],
    "rank_by_attributable_deaths": [
        "Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus",
        "Acinetobacter baumannii", "Streptococcus pneumoniae", "Mycobacterium tuberculosis",
    ],
    "source": "Murray et al., Lancet 2022 (GRAM 2019); PMC8841637",
    "note": "Per-pathogen DALYs/deaths are appendix-sourced and passed to the index "
            "functions; not hard-coded here.",
}

# Locked funding snapshot (the index denominator). Verified figures from Czaplewski
# et al., "An overview of global public and philanthropic investments into
# antibacterial therapeutics (2017-23)", Lancet Microbe 2026 (epub 2026-01-09;
# doi:10.1016/S2666-5247(25)00216-2) — a frozen, peer-reviewed extract of the Global
# AMR R&D Hub Dynamic Dashboard. All values are US$ millions, public + philanthropic,
# 2017-2023. Named per-pathogen amounts are quoted verbatim from the paper; the three
# lowest-funded named GRAM species (E. coli, A. baumannii, K. pneumoniae) have exact
# magnitudes only in appendix 1 p18 (below P. aeruginosa in figure 3B), so they are
# folded into ``other_species_specific_musd`` here and supplied per-pathogen at call
# time for the numeric index. See docs/reference_rd_alignment_2026-06-09.md.
RD_HUB_SNAPSHOT_2026: dict = {
    "total_funding_musd": 2510.0,        # "US$2.51 billion ... by 130 funders"
    "n_funders": 130,
    "species_specific_total_musd": 1058.0,   # "$1058 million of species-specific funds"
    "named_species_funding_musd": {
        "Mycobacterium tuberculosis": 474.0,   # "a fifth of the overall funds ... $474M"
        "Staphylococcus aureus": 142.0,        # "13% of $1058M ...; n=87"
        "Clostridioides difficile": 141.0,     # "13%, n=47"
        "Neisseria gonorrhoeae": 101.0,        # "10%, n=27"
        "Pseudomonas aeruginosa": 87.0,        # "8%, n=73"
    },
    # species-specific funding not in the named top-five (incl. the GRAM Gram-negatives
    # E. coli, A. baumannii, K. pneumoniae, which rank BELOW P. aeruginosa in figure 3B):
    "other_species_specific_musd": 113.0,    # 1058 - (474+142+141+101+87)
    # S. pneumoniae does not appear among species-specific targets (figure 3B): ~0.
    "unfunded_gram_panel_species": ["Streptococcus pneumoniae"],
    "source": "Czaplewski et al., Lancet Microbe 2026 (epub 2026-01-09)",
}

# Robustness denominator: a live, directly-extracted Global AMR R&D Hub "by Genus" Data
# Table (Dynamic Dashboard / Investment Gallery, last updated 2026-06-30). US$ millions,
# mapped genus -> GRAM panel species. Used ONLY as a cross-check (``genus_robustness_alignment``),
# never as the primary index denominator (which stays the frozen, published Czaplewski
# extract above). Scope was set on the dashboard BEFORE export to MATCH the locked
# public+philanthropic therapeutics 2017-2023 denominator; recorded verbatim from the
# export's "Applied filters" header:
#   * Infectious Agent (Group) = Bacteria, CategoryType (Sector) = Human, Currency = USD,
#     Status = Active or Closed
#   * CategoryType (Research Area) = Therapeutics (NOT basic research / vaccines / diagnostics)
#   * Year in 2017..2023 (7 years; excludes the future-dated 2024-2033 commitments)
#   * FunderType in {Public-Government, Public-Other, Private-Non Profit} (the public+
#     philanthropic set; EXCLUDES Public-Private partnerships)
#   -> total $2,436M, 128 funders ~= the $2.51bn / 130-funder Czaplewski denominator (the
#      small residual is the Hub's documented retrospective revision: live 2026-06-30 vs
#      Czaplewski's 2026-01-09 freeze).
#   * GENUS level: Escherichia/Klebsiella/Acinetobacter/Pseudomonas ~= the panel species;
#     "Staphylococcus spp." and "Streptococcus spp." are broader than S. aureus /
#     S. pneumoniae, but the Therapeutics filter strips out pneumococcal-vaccine R&D, so
#     (unlike the old all-areas pull) the genus row no longer masks the pneumococcus gap.
RD_HUB_GENUS_SNAPSHOT: dict = {
    "funding_musd": {
        "Escherichia coli": 42.872256,          # Escherichia spp.
        "Staphylococcus aureus": 219.914081,    # Staphylococcus spp. (genus > S. aureus)
        "Klebsiella pneumoniae": 53.920828,     # Klebsiella spp.
        "Streptococcus pneumoniae": 27.236003,  # Streptococcus spp. (therapeutics only)
        "Acinetobacter baumannii": 71.619879,   # Acinetobacter spp.
        "Pseudomonas aeruginosa": 146.874858,   # Pseudomonas spp.
    },
    "extract_date": "2026-06-30",
    "scope": ("Hub 'by Genus' Data Table: Bacteria/Human/USD/Active-or-Closed; Research Area "
              "= Therapeutics; Year 2017-2023; FunderType = public+philanthropic (Public-"
              "Government + Public-Other + Private-Non Profit, excl. Public-Private); genus-"
              "level. Scoped to match the locked Czaplewski denominator (total $2,436M / 128 "
              "funders ~= $2.51bn / 130) — live cross-check, not a replacement."),
    "source": "Global AMR R&D Hub Dynamic Dashboard, by-Genus export (last updated 2026-06-30)",
}

# Verified GRAM-2019 per-pathogen burden numerator, lifted from Murray et al. appendix 1
# Table S22 ("Global deaths and DALYs ... by pathogen-drug combination, 2019", the
# 'Resistance to one or more antibiotics' aggregate row per pathogen). Counts are in
# THOUSANDS with 95% UI (median, lo, hi). Self-consistency check (enforced by a test):
# associated-death medians sum to 3,572k = 3.57M and attributable to 928.6k ~= 929k,
# matching the headline totals; both rank orders match the main text exactly. These are
# the AMR-burden figures (NOT the Ikuta 2022 all-cause infection deaths). Source:
# PMC8841637 appendix 1 Table S22. See docs/reference_rd_alignment_2026-06-09.md.
def _b(ad, dd, adaly, ddaly) -> dict:
    """Pack one pathogen's (median, lo, hi) tuples: assoc/attrib deaths and DALYs (k)."""
    return {"assoc_deaths_k": ad, "attrib_deaths_k": dd,
            "assoc_dalys_k": adaly, "attrib_dalys_k": ddaly}


GRAM_BURDEN_2019: dict = {
    # _b(assoc_deaths, attrib_deaths, assoc_DALYs, attrib_DALYs); each (median, lo, hi), thousands
    "Escherichia coli": _b((829, 601, 1120), (219, 152, 316),
                           (28000, 21000, 36900), (7520, 5270, 10500)),
    "Staphylococcus aureus":    _b((748, 554, 1000), (178, 104, 280),
                                   (24900, 18600, 32700), (5870, 3550, 9220)),
    "Klebsiella pneumoniae":    _b((642, 465, 863), (193, 130, 272),
                                   (27400, 20300, 36100), (8200, 5550, 11400)),
    "Streptococcus pneumoniae": _b((596, 490, 727), (122, 82.4, 166),
                                   (29800, 24400, 36700), (6110, 4050, 8330)),
    "Acinetobacter baumannii":  _b((423, 252, 647), (132, 75.7, 213),
                                   (11800, 7290, 17800), (3670, 2150, 5760)),
    "Pseudomonas aeruginosa":   _b((334, 234, 457), (84.6, 53, 127),
                                   (12000, 8630, 16100), (3050, 1980, 4530)),
}


def cross_cutting_share(funding_by_pathogen: dict, cross_cutting_funding: float) -> dict:
    """Headline magnitude: how much tracked funding is NOT pathogen-specific.

    Reported BEFORE any index. ``flagged`` is True when pathogen-specific funding is a
    minority of the total — the signal to widen every downstream caveat.
    """
    specific = float(sum(funding_by_pathogen.values()))
    total = specific + float(cross_cutting_funding)
    frac = (cross_cutting_funding / total) if total > 0 else float("nan")
    return {
        "pathogen_specific_funding": specific,
        "cross_cutting_funding": float(cross_cutting_funding),
        "total_funding": total,
        "cross_cutting_fraction": frac,
        "flagged": bool(frac > 0.5),
    }


def cross_cutting_headline(snapshot: dict | None = None) -> dict:
    """The de-gated Component-4 headline: how little antibacterial R&D is pathogen-specific.

    Computed straight from the locked Czaplewski/Hub snapshot (``RD_HUB_SNAPSHOT_2026``),
    so it needs no per-pathogen burden and no appendix magnitudes. Reuses
    :func:`cross_cutting_share`. The headline finding: pathogen-specific funding is a
    MINORITY of total antibacterial R&D, so ``flagged`` is True and every per-pathogen
    caveat widens accordingly.
    """
    s = RD_HUB_SNAPSHOT_2026 if snapshot is None else snapshot
    funding_by_pathogen = dict(s["named_species_funding_musd"])
    funding_by_pathogen["Other species-specific"] = s["other_species_specific_musd"]
    cross_cutting = s["total_funding_musd"] - s["species_specific_total_musd"]
    head = cross_cutting_share(funding_by_pathogen, cross_cutting)
    head["source"] = s["source"]
    head["total_funding_musd"] = s["total_funding_musd"]
    return head


def snapshot_from_hub_export(frame, source: str | None = None) -> dict:
    """Collapse a tidy Hub investment export into a ``RD_HUB_SNAPSHOT_2026``-shaped dict.

    ``frame`` is the output of :func:`data_loading.load_rd_hub_snapshot` (columns
    ``pathogen, year, funder_type, investment_usd, is_cross_cutting`` [+ ``funder``]). The
    returned dict is a drop-in denominator for :func:`cross_cutting_headline` and the
    per-pathogen index: species-specific rows are summed per pathogen (US$ millions),
    cross-cutting rows go to the total only, and any GRAM-panel species absent from the
    export is reported under ``unfunded_gram_panel_species``.

    This sources the index denominator from a direct dated Hub export rather than the
    published Czaplewski transcription; keep the latter as the reconciliation cross-check
    (:func:`reconcile_hub_snapshot`).
    """
    total_usd = float(frame["investment_usd"].sum())
    specific = frame[~frame["is_cross_cutting"]]
    per_pathogen = (
        specific.groupby("pathogen")["investment_usd"].sum().div(1e6).round(6)
    )
    named = {str(k): float(v) for k, v in per_pathogen.items()}
    species_specific_usd = float(specific["investment_usd"].sum())

    panel = set(GRAM_BURDEN_2019)
    unfunded = sorted(panel - {p for p, v in named.items() if v > 0})
    n_funders = int(frame["funder"].nunique()) if "funder" in frame.columns else None

    return {
        "total_funding_musd": round(total_usd / 1e6, 6),
        "n_funders": n_funders,
        "species_specific_total_musd": round(species_specific_usd / 1e6, 6),
        "named_species_funding_musd": named,
        # The export carries every species-specific pathogen explicitly, so there is no
        # residual "other" bucket (kept for shape-compatibility with the published snapshot).
        "other_species_specific_musd": 0.0,
        "unfunded_gram_panel_species": unfunded,
        "source": source or (
            f"Global AMR R&D Hub Dynamic Dashboard export {config.RD_HUB_SNAPSHOT_DATE}"
        ),
    }


def reconcile_hub_snapshot(
    export_snapshot: dict,
    published: dict | None = None,
    rel_tol: float = 0.15,
) -> dict:
    """Cross-check a Hub-export snapshot against the published Czaplewski extract.

    Compares the two headline magnitudes (total and species-specific funding) and returns a
    per-field record (export value, published value, relative difference, within-tolerance
    flag). The Hub dashboard is retrospectively revised, so some drift from the frozen
    Czaplewski extract is EXPECTED — this surfaces its size rather than asserting equality.
    ``ok`` is True when every compared field is within ``rel_tol``.
    """
    ref = RD_HUB_SNAPSHOT_2026 if published is None else published
    fields = ("total_funding_musd", "species_specific_total_musd")
    rows = {}
    for f in fields:
        exp = float(export_snapshot[f])
        pub = float(ref[f])
        rel = abs(exp - pub) / pub if pub else float("nan")
        rows[f] = {
            "export": exp,
            "published": pub,
            "rel_diff": rel,
            "within_tol": bool(rel <= rel_tol),
        }
    return {
        "fields": rows,
        "rel_tol": rel_tol,
        "ok": all(r["within_tol"] for r in rows.values()),
        "export_source": export_snapshot.get("source"),
        "published_source": ref.get("source"),
    }


def mismatch_index(
    burden_by_pathogen: dict,
    funding_by_pathogen: dict,
    floor: float | None = None,
) -> dict:
    """Per-pathogen Axis-A mismatch M_p = burden_share / funding_share (+ log2).

    Shares are computed over the pathogens present in BOTH inputs. ``floor`` (if given) is a
    minimum funding share applied before the ratio, bounding the divide-by-small instability
    for sparsely-funded pathogens; call twice (with and without) to report both rankings.
    Returns per-pathogen records sorted by log2(M_p) descending (most under-funded first).
    """
    pathogens = sorted(set(burden_by_pathogen) & set(funding_by_pathogen))
    if len(pathogens) < 2:
        raise ValueError("Need >= 2 pathogens present in both burden and funding inputs.")

    b = np.array([float(burden_by_pathogen[p]) for p in pathogens])
    f = np.array([float(funding_by_pathogen[p]) for p in pathogens])
    if b.sum() <= 0 or f.sum() <= 0:
        raise ValueError("Burden and funding totals must be positive.")

    burden_share = b / b.sum()
    funding_share = f / f.sum()
    fs = funding_share if floor is None else np.maximum(funding_share, floor)

    m = burden_share / fs
    records = [
        {"pathogen": p, "burden_share": float(bs), "funding_share": float(fsh),
         "mismatch": float(mi), "log2_mismatch": float(np.log2(mi))}
        for p, bs, fsh, mi in zip(pathogens, burden_share, fs, m, strict=True)
    ]
    records.sort(key=lambda r: r["log2_mismatch"], reverse=True)
    return {"floor": floor, "n_pathogens": len(pathogens), "ranking": records}


def _spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rank correlation (Pearson on average-tie ranks); no scipy."""
    rx = pd.Series(x).rank(method="average").to_numpy()
    ry = pd.Series(y).rank(method="average").to_numpy()
    return float(np.corrcoef(rx, ry)[0, 1])


def spearman_burden_funding(
    burden_by_pathogen: dict,
    funding_by_pathogen: dict,
    n_boot: int = 2000,
    seed: int | None = None,
) -> dict:
    """Spearman ρ between burden and funding across pathogens, with a bootstrap CI.

    The single pre-registered inferential summary for Component 4 — DESCRIPTIVE only. With
    n = 5-6 pathogens the bootstrap CI is wide and is reported with that caveat; no p-value,
    no fitted line.
    """
    if seed is None:
        seed = config.step_seed(4)
    rng = np.random.default_rng(seed)
    pathogens = sorted(set(burden_by_pathogen) & set(funding_by_pathogen))
    x = np.array([float(burden_by_pathogen[p]) for p in pathogens])
    y = np.array([float(funding_by_pathogen[p]) for p in pathogens])
    n = len(pathogens)

    rho = _spearman_rho(x, y)
    draws = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        draws[b] = _spearman_rho(x[idx], y[idx]) if np.unique(idx).size > 1 else np.nan
    draws = draws[np.isfinite(draws)]
    lo, hi = np.quantile(draws, [0.025, 0.975]) if draws.size else (float("nan"), float("nan"))
    return {
        "spearman_rho": rho,
        "ci_lower": float(lo),
        "ci_upper": float(hi),
        "n_pathogens": n,
        "caveat": (f"Descriptive rank correlation over n={n} pathogens; the bootstrap CI is "
                   "wide and no line is fitted. Not a powered or causal test."),
    }


def monte_carlo_mismatch_ranking(
    burden_ui_by_pathogen: dict,
    funding_by_pathogen: dict,
    floor: float | None = None,
    draws: int | None = None,
    seed: int | None = None,
) -> dict:
    """Propagate GRAM burden uncertainty intervals into the mismatch ranking's stability.

    ``burden_ui_by_pathogen`` maps pathogen -> (median, lo95, hi95). Each draw samples a
    lognormal burden per pathogen (matched to the UI), recomputes M_p, and records each
    pathogen's rank. Returns, per pathogen, the median log2(M_p) and the share of draws in
    which it is the single most under-funded — i.e. how robust the headline ranking is.
    """
    if draws is None:
        draws = config.MONTE_CARLO_DRAWS
    if seed is None:
        seed = config.step_seed(4)
    rng = np.random.default_rng(seed)

    pathogens = sorted(set(burden_ui_by_pathogen) & set(funding_by_pathogen))
    z = 1.959963985
    mu = {}
    sigma = {}
    for p in pathogens:
        med, lo, hi = burden_ui_by_pathogen[p]
        mu[p] = np.log(med)
        sigma[p] = (np.log(hi) - np.log(lo)) / (2 * z)

    log2m = {p: np.empty(int(draws)) for p in pathogens}
    top_count = {p: 0 for p in pathogens}
    for d in range(int(draws)):
        burden = {p: float(rng.lognormal(mu[p], sigma[p])) for p in pathogens}
        rank = mismatch_index(burden, funding_by_pathogen, floor=floor)["ranking"]
        for r in rank:
            log2m[r["pathogen"]][d] = r["log2_mismatch"]
        top_count[rank[0]["pathogen"]] += 1

    return {
        "floor": floor,
        "draws": int(draws),
        "per_pathogen": {
            p: {"log2_mismatch_median": float(np.median(log2m[p])),
                "log2_mismatch_ci": [float(np.quantile(log2m[p], 0.025)),
                                     float(np.quantile(log2m[p], 0.975))],
                "p_most_underfunded": top_count[p] / int(draws)}
            for p in pathogens
        },
    }


def _sample_funding_split(rng) -> dict:
    """Sample the unfetched $113M E.coli/A.baumannii/K.pneumoniae split honestly.

    The only verified constraints (Czaplewski figure 3B / appendix 1 p18): the three sum
    to $113M, are ordered E.coli > A.baumannii > K.pneumoniae, and each is below
    P.aeruginosa ($87M). S.pneumoniae is absent from the species-specific figure, so its
    funding is below the smallest shown bar (K.pneumoniae). We sample uniformly over this
    order-constrained region rather than guess point values; the index is then reported as
    a distribution. Replace with the appendix exact values when available.
    """
    while True:
        parts = np.sort(rng.random(3))[::-1]
        e, a, k = (parts / parts.sum()) * 113.0
        if e < 87.0:  # below P. aeruginosa
            break
    return {
        "Staphylococcus aureus": 142.0,
        "Pseudomonas aeruginosa": 87.0,
        "Escherichia coli": float(e),
        "Acinetobacter baumannii": float(a),
        "Klebsiella pneumoniae": float(k),
        "Streptococcus pneumoniae": float(rng.uniform(0.0, k)),  # below smallest shown bar
    }


def gram_panel_alignment(
    burden_metric: str = "assoc_deaths_k",
    floor: float = 0.02,
    draws: int | None = None,
    seed: int | None = None,
) -> dict:
    """Closed-out Component-4 index over the six GRAM leading pathogens.

    Numerator = verified ``GRAM_BURDEN_2019`` (Table S22, with 95% UIs); denominator =
    the locked ``RD_HUB_SNAPSHOT_2026`` funding, with the unfetched $113M Gram-negative
    split propagated as uncertainty via :func:`_sample_funding_split`. Each Monte-Carlo
    draw samples a lognormal burden per pathogen (matched to its UI) AND a feasible funding
    split, recomputes the floored log2 mismatch, and records the ranking. Returns the
    cross-cutting headline, the Spearman summary (rank-determined, robust to the split),
    and per-pathogen log2 medians + CIs + P(most under-funded). Descriptive only (n=6).
    """
    if config.RD_HUB_SNAPSHOT_DATE is None:
        raise ValueError("Lock config.RD_HUB_SNAPSHOT_DATE before running Component 4.")
    if draws is None:
        draws = config.MONTE_CARLO_DRAWS
    if seed is None:
        seed = config.step_seed(4)
    rng = np.random.default_rng(seed)
    z = 1.959963985

    pathogens = list(GRAM_BURDEN_2019)
    mu, sigma = {}, {}
    for p in pathogens:
        med, lo, hi = GRAM_BURDEN_2019[p][burden_metric]
        mu[p] = np.log(med)
        sigma[p] = (np.log(hi) - np.log(lo)) / (2 * z)

    log2m = {p: np.empty(int(draws)) for p in pathogens}
    top_count = {p: 0 for p in pathogens}
    for d in range(int(draws)):
        burden = {p: float(rng.lognormal(mu[p], sigma[p])) for p in pathogens}
        funding = _sample_funding_split(rng)
        rank = mismatch_index(burden, funding, floor=floor)["ranking"]
        for r in rank:
            log2m[r["pathogen"]][d] = r["log2_mismatch"]
        top_count[rank[0]["pathogen"]] += 1

    # Spearman on the median (point) configuration — funding ranks are fixed by figure 3B,
    # so the point estimate is robust to the split; the bootstrap CI is wide at n=6.
    median_burden = {p: GRAM_BURDEN_2019[p][burden_metric][0] for p in pathogens}
    median_funding = _sample_funding_split(np.random.default_rng(seed))
    spearman = spearman_burden_funding(median_burden, median_funding, seed=seed)

    per_pathogen = {
        p: {"log2_mismatch_median": float(np.median(log2m[p])),
            "log2_mismatch_ci": [float(np.quantile(log2m[p], 0.025)),
                                 float(np.quantile(log2m[p], 0.975))],
            "p_most_underfunded": top_count[p] / int(draws)}
        for p in pathogens
    }
    ranked = sorted(per_pathogen, key=lambda p: per_pathogen[p]["log2_mismatch_median"],
                    reverse=True)
    return {
        "burden_metric": burden_metric,
        "floor": floor,
        "draws": int(draws),
        "cross_cutting": cross_cutting_headline(),
        "spearman": spearman,
        "per_pathogen": per_pathogen,
        "underfunded_ranking": ranked,
        "caption": alignment_caption(),
        "note": ("Numerator: GRAM-2019 Table S22 (verified, 95% UI). Denominator: Czaplewski "
                 "et al. 2026; the $113M E.coli/A.baumannii/K.pneumoniae split is unfetched "
                 "(appendix 1 p18) and propagated as uncertainty, not fabricated. Descriptive, "
                 "n=6, no fitted line."),
    }


def genus_robustness_alignment(
    burden_metric: str = "assoc_deaths_k",
    funding_musd: dict | None = None,
    floor: float = 0.02,
    draws: int | None = None,
    seed: int | None = None,
) -> dict:
    """Robustness cross-check of the Component-4 ranking against a live Hub genus extract.

    Re-runs the per-pathogen mismatch index with the same verified ``GRAM_BURDEN_2019``
    numerator (burden UIs propagated by :func:`monte_carlo_mismatch_ranking`) but swaps the
    frozen, published Czaplewski denominator for a live, directly-extracted Hub "by Genus"
    Data Table (``RD_HUB_GENUS_SNAPSHOT``) scoped on the dashboard to MATCH that denominator
    (public+philanthropic therapeutics, 2017-2023, genus-level; see the constant's note) — so
    it is an independent-source corroboration, NOT a replacement for the primary index.

    Returns the per-pathogen log2 medians/CIs + P(most under-funded), the under-funded
    ranking, the Spearman summary, and a direct comparison to the primary Czaplewski ranking
    (``rank_match`` per pathogen + the qualitative direction agreement). The finding is
    corroborated: the community Gram-negatives stay under-funded and S. aureus / P. aeruginosa
    stay over-funded; and because the Therapeutics filter strips pneumococcal-vaccine R&D, the
    live extract now RECOVERS the S. pneumoniae under-funding signal that an all-areas genus
    pull had masked.
    """
    funding = RD_HUB_GENUS_SNAPSHOT["funding_musd"] if funding_musd is None else funding_musd
    if draws is None:
        draws = config.MONTE_CARLO_DRAWS
    if seed is None:
        seed = config.step_seed(4)

    burden_ui = {p: GRAM_BURDEN_2019[p][burden_metric] for p in GRAM_BURDEN_2019}
    mc = monte_carlo_mismatch_ranking(burden_ui, funding, floor=floor, draws=draws, seed=seed)
    per_pathogen = mc["per_pathogen"]
    ranked = sorted(per_pathogen, key=lambda p: per_pathogen[p]["log2_mismatch_median"],
                    reverse=True)

    median_burden = {p: GRAM_BURDEN_2019[p][burden_metric][0] for p in GRAM_BURDEN_2019}
    spearman = spearman_burden_funding(median_burden, funding, seed=seed)

    # Compare against the primary (Czaplewski) index run on the same metric/floor/seed.
    primary = gram_panel_alignment(burden_metric=burden_metric, floor=floor,
                                   draws=draws, seed=seed)
    primary_rank = primary["underfunded_ranking"]
    comparison = {
        "primary_ranking": primary_rank,
        "genus_ranking": ranked,
        "rank_match": {p: (primary_rank.index(p) == ranked.index(p)) for p in ranked},
        "direction_agrees": {
            p: ((per_pathogen[p]["log2_mismatch_median"] > 0)
                == (primary["per_pathogen"][p]["log2_mismatch_median"] > 0))
            for p in ranked
        },
    }
    return {
        "burden_metric": burden_metric,
        "floor": floor,
        "draws": int(draws),
        "funding_source": RD_HUB_GENUS_SNAPSHOT["source"],
        "scope_caveat": RD_HUB_GENUS_SNAPSHOT["scope"],
        "spearman": spearman,
        "per_pathogen": per_pathogen,
        "underfunded_ranking": ranked,
        "comparison": comparison,
        "note": ("Independent-source corroboration: same GRAM-2019 burden numerator, "
                 "denominator swapped to a live Hub by-Genus extract scoped to match the "
                 "Czaplewski denominator (public+phil therapeutics 2017-2023, genus-level). "
                 "The Gram-negative under-funding and S. aureus / P. aeruginosa over-funding "
                 "are corroborated; the Therapeutics filter strips pneumococcal-vaccine R&D so "
                 "the live extract recovers the S. pneumoniae under-funding signal."),
    }


def analyze_alignment(
    burden_by_pathogen: dict,
    funding_by_pathogen: dict,
    cross_cutting_funding: float = 0.0,
    floor: float | None = None,
) -> dict:
    """Component-4 entrypoint: cross-cutting headline FIRST, then index + Spearman.

    Requires ``config.RD_HUB_SNAPSHOT_DATE`` to be locked (the funding inputs come from a
    frozen, dated Hub snapshot). Returns the cross-cutting magnitude, the mismatch ranking
    with and without the low-denominator floor, and the Spearman summary — all descriptive.
    """
    if config.RD_HUB_SNAPSHOT_DATE is None:
        raise ValueError(
            "Lock config.RD_HUB_SNAPSHOT_DATE (a frozen, dated Hub snapshot) before running "
            "Component 4 (pre-reg §6/§7; Hub data is retrospectively revised)."
        )
    return {
        "snapshot_caption": alignment_caption(),
        "cross_cutting": cross_cutting_share(funding_by_pathogen, cross_cutting_funding),
        "ranking_no_floor": mismatch_index(burden_by_pathogen, funding_by_pathogen),
        "ranking_with_floor": mismatch_index(burden_by_pathogen, funding_by_pathogen, floor=floor),
        "spearman": spearman_burden_funding(burden_by_pathogen, funding_by_pathogen),
    }


def alignment_caption() -> str:
    """Standard figure-caption suffix required on every Component-4 figure."""
    return (
        f"Global AMR R&D Hub snapshot {config.RD_HUB_SNAPSHOT_DATE} "
        f"({config.RD_HUB_SOURCE}); {config.RD_HUB_SCOPE}; "
        f"window {config.RD_HUB_WINDOW[0]}-{config.RD_HUB_WINDOW[1]}. "
        f"Burden: {GRAM_PANEL['source']}."
    )


# --- Catchment-specific alignment (the Cross-Domain finding made local) ----------
#
# The global GRAM ranking is dominated by E. coli / S. aureus / S. pneumoniae and answers a
# global question. The Cross-Domain Award asks what the *catchment* (SSA hospital) data say.
# We re-weight the burden by the pathogen mix actually observed in the SPIDAAR severe-HAI
# isolates — the only delivered Register source covering all six GRAM panel species (the ATLAS
# catchment cell is E. coli / K. pneumoniae only) — and re-run the same mismatch index against
# the same global funding shares. This is an *isolation-frequency* burden proxy, NOT a mortality
# burden: it says which pathogens dominate the local resistant-HAI caseload, then asks whether
# global R&D money follows them. Labelled as such throughout.

# Genus/species substrings that map an organism string onto a GRAM panel pathogen. Acinetobacter
# is aggregated at genus level onto A. baumannii (the dominant clinical Acinetobacter and the
# GRAM/ATLAS panel representative); Klebsiella is matched species-specifically to K. pneumoniae.
_PANEL_ALIASES: dict = {
    "Escherichia coli": ("escherichia coli",),
    "Staphylococcus aureus": ("staphylococcus aureus",),
    "Klebsiella pneumoniae": ("klebsiella pneumoniae",),
    "Streptococcus pneumoniae": ("streptococcus pneumoniae",),
    "Acinetobacter baumannii": ("acinetobacter",),  # genus-level aggregation (documented)
    "Pseudomonas aeruginosa": ("pseudomonas aeruginosa",),
}


def catchment_pathogen_counts(isolates: pd.DataFrame) -> dict:
    """Count catchment isolates and resistant isolates per GRAM panel pathogen.

    Maps the (possibly polymicrobial) ``organism`` string onto the six panel species via
    :data:`_PANEL_ALIASES`. Returns ``{species: {"isolates": n, "resistant": n}}`` over the
    six panel species only (non-panel organisms are out of scope for this index).
    """
    org = isolates["organism"].astype(str).str.lower()
    res = pd.to_numeric(isolates["resistant"], errors="coerce").fillna(0).to_numpy()
    counts = {}
    for species, aliases in _PANEL_ALIASES.items():
        mask = np.zeros(len(isolates), dtype=bool)
        for a in aliases:
            mask |= org.str.contains(a, regex=False).to_numpy()
        counts[species] = {"isolates": int(mask.sum()), "resistant": int(res[mask].sum())}
    return counts


def catchment_alignment(
    isolates: pd.DataFrame | None = None,
    counts: dict | None = None,
    weight: str = "resistant",
    floor: float = 0.02,
    draws: int | None = None,
    seed: int | None = None,
) -> dict:
    """Catchment-specific Axis-A mismatch: local pathogen mix vs global funding share.

    Burden share per pathogen is the catchment isolate frequency (``weight="frequency"``) or
    the resistant-isolate frequency (``weight="resistant"``, the AMR-burden proxy and default),
    taken from the SPIDAAR isolates. Each Monte-Carlo draw samples the burden shares from a
    Dirichlet (Jeffreys prior on the observed counts — propagating small-cell uncertainty such
    as the ~3 S. pneumoniae isolates) and the funding split from :func:`_sample_funding_split`,
    recomputes the floored log2 mismatch, and records the ranking. Returns per-pathogen log2
    medians + CIs + P(most under-funded), the Spearman summary, and the point burden shares —
    mirroring :func:`gram_panel_alignment` so the catchment and global results are comparable.
    """
    if config.RD_HUB_SNAPSHOT_DATE is None:
        raise ValueError("Lock config.RD_HUB_SNAPSHOT_DATE before running Component 4.")
    if weight not in ("resistant", "frequency"):
        raise ValueError("weight must be 'resistant' or 'frequency'.")
    if counts is None:
        if isolates is None:
            raise ValueError("Provide either an isolates frame or a precomputed counts dict.")
        counts = catchment_pathogen_counts(isolates)
    if draws is None:
        draws = config.MONTE_CARLO_DRAWS
    if seed is None:
        seed = config.step_seed(5)
    rng = np.random.default_rng(seed)

    key = "resistant" if weight == "resistant" else "isolates"
    pathogens = list(counts)
    n = np.array([counts[p][key] for p in pathogens], dtype=float)
    if n.sum() <= 0:
        raise ValueError("No catchment isolates in the panel under the chosen weight.")
    alpha = n + 0.5  # Jeffreys-prior Dirichlet over the observed composition

    log2m = {p: np.empty(int(draws)) for p in pathogens}
    top_count = {p: 0 for p in pathogens}
    for d in range(int(draws)):
        shares = rng.dirichlet(alpha)
        burden = {p: float(s) for p, s in zip(pathogens, shares, strict=True)}
        funding = _sample_funding_split(rng)
        rank = mismatch_index(burden, funding, floor=floor)["ranking"]
        for r in rank:
            log2m[r["pathogen"]][d] = r["log2_mismatch"]
        top_count[rank[0]["pathogen"]] += 1

    point_shares = {p: float(n[i] / n.sum()) for i, p in enumerate(pathogens)}
    median_funding = _sample_funding_split(np.random.default_rng(seed))
    spearman = spearman_burden_funding(point_shares, median_funding, seed=seed)

    per_pathogen = {
        p: {"burden_share": point_shares[p],
            "log2_mismatch_median": float(np.median(log2m[p])),
            "log2_mismatch_ci": [float(np.quantile(log2m[p], 0.025)),
                                 float(np.quantile(log2m[p], 0.975))],
            "p_most_underfunded": top_count[p] / int(draws)}
        for p in pathogens
    }
    ranked = sorted(per_pathogen, key=lambda p: per_pathogen[p]["log2_mismatch_median"],
                    reverse=True)
    return {
        "weight": weight,
        "floor": floor,
        "draws": int(draws),
        "catchment_counts": counts,
        "burden_shares": point_shares,
        "spearman": spearman,
        "per_pathogen": per_pathogen,
        "underfunded_ranking": ranked,
        "source": ("SPIDAAR severe-HAI isolates (catchment); funding "
                   + RD_HUB_SNAPSHOT_2026["source"]),
        "note": ("Burden is an isolation-frequency proxy from the catchment HAI isolates (NOT a "
                 "mortality burden); funding shares are the global Hub snapshot with the $113M "
                 "Gram-negative split propagated as uncertainty. Descriptive, n=6, no fitted "
                 "line."),
    }
