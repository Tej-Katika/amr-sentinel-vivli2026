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
