"""Component 6: the surveillance blind-spot axis (the third mismatch).

Component 4 (``rd_alignment``) shows AMR R&D *funding* is misaligned with mortality burden.
This module adds a third axis from the ATLAS surveillance register itself: where the
*surveillance* is — i.e. how many isolates are actually observed — relative to the same
burden, both per pathogen and per region. The finding the challenge data make unavoidable:
burden, funding, AND surveillance are all misaligned, and all three neglect the same
sub-Saharan-African community Gram-negatives. The richest dataset in the catalogue
(~1M isolates, 83 countries, 2004-2024) is used at full scale, not just the catchment cell.

Design mirrors ``rd_alignment`` so the three axes are directly comparable: the same verified
``GRAM_BURDEN_2019`` numerator and the same Monte-Carlo / Spearman mismatch machinery, with
the funding denominator swapped for ATLAS isolate counts. Pure/deterministic; the heavy ATLAS
read lives in ``data_loading.load_atlas_backbone``. Descriptive only (n=6 pathogens).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config
from .rd_alignment import (
    _PANEL_ALIASES,
    GRAM_BURDEN_2019,
    RD_HUB_GENUS_SNAPSHOT,
    monte_carlo_mismatch_ranking,
    spearman_burden_funding,
)

# Sub-Saharan Africa country labels as ATLAS spells them (10 are actually present;
# the rest are carried so a future extract is classified correctly). Used only to size the
# geographic surveillance blind spot. "Ivory Coast" is ATLAS's spelling of Cote d'Ivoire.
SSA_COUNTRIES: frozenset = frozenset({
    "South Africa", "Nigeria", "Ivory Coast", "Cote d'Ivoire", "Kenya", "Cameroon",
    "Uganda", "Malawi", "Ghana", "Namibia", "Mauritius", "Tanzania", "Ethiopia",
    "Senegal", "Zimbabwe", "Zambia", "Mozambique", "Rwanda", "Botswana", "Mali",
    "Burkina Faso", "Angola", "Gabon", "Madagascar", "Democratic Republic of the Congo",
    "Congo", "Niger", "Benin", "Togo", "Eswatini", "Gambia",
})

# Extended GRAM-2019 burden panel for the n=16 surveillance-mismatch (gives the descriptive
# n=6 index real inferential weight). Verified verbatim from Murray et al. (Lancet 2022,
# PMC8841637) appendix Table S22, the per-pathogen "Resistance to one or more antibiotics"
# aggregate row, ASSOCIATED deaths (median, lo95, hi95) in thousands. Scope = pathogens with
# both a Table-S22 associated-death estimate AND meaningful ATLAS surveillance (TB/gonorrhoea/
# Salmonella/Shigella are in S22 but out of ATLAS scope, so excluded). The six panel species
# reuse the exact values already in ``rd_alignment.GRAM_BURDEN_2019`` (a test enforces the match).
GRAM_ASSOC_DEATHS_EXTENDED: dict = {
    "Escherichia coli": (829.0, 601.0, 1120.0),
    "Staphylococcus aureus": (748.0, 554.0, 1000.0),
    "Klebsiella pneumoniae": (642.0, 465.0, 863.0),
    "Streptococcus pneumoniae": (596.0, 490.0, 727.0),
    "Acinetobacter baumannii": (423.0, 252.0, 647.0),
    "Pseudomonas aeruginosa": (334.0, 234.0, 457.0),
    "Enterococcus faecium": (200.0, 123.0, 303.0),
    "Enterobacter spp.": (185.0, 122.0, 264.0),
    "Group B Streptococcus": (173.0, 125.0, 232.0),
    "Enterococcus faecalis": (112.0, 69.2, 167.0),
    "Proteus spp.": (80.0, 52.7, 115.0),
    "Serratia spp.": (42.7, 27.0, 65.3),
    "Group A Streptococcus": (39.0, 18.3, 77.1),
    "Citrobacter spp.": (35.5, 21.7, 52.9),
    "Haemophilus influenzae": (31.5, 25.5, 39.0),
    "Morganella spp.": (3.0, 1.93, 4.59),
}

# ATLAS ``Species`` -> extended-panel pathogen. Genus-level entries ("... spp.") aggregate all
# species of the genus (matching Table S22's genus rows); Group A/B Streptococcus map to
# S. pyogenes / S. agalactiae; the six panel species stay species-specific (Acinetobacter and
# the "spp." rows are genus-aggregated, as documented for the core panel).
EXTENDED_PANEL_ALIASES: dict = {
    "Escherichia coli": ("escherichia coli",),
    "Staphylococcus aureus": ("staphylococcus aureus",),
    "Klebsiella pneumoniae": ("klebsiella pneumoniae",),
    "Streptococcus pneumoniae": ("streptococcus pneumoniae",),
    "Acinetobacter baumannii": ("acinetobacter",),
    "Pseudomonas aeruginosa": ("pseudomonas aeruginosa",),
    "Enterococcus faecium": ("enterococcus faecium",),
    "Enterobacter spp.": ("enterobacter",),
    "Group B Streptococcus": ("streptococcus agalactiae",),
    "Enterococcus faecalis": ("enterococcus faecalis",),
    "Proteus spp.": ("proteus",),
    "Serratia spp.": ("serratia",),
    "Group A Streptococcus": ("streptococcus pyogenes",),
    "Citrobacter spp.": ("citrobacter",),
    "Haemophilus influenzae": ("haemophilus influenzae",),
    "Morganella spp.": ("morganella",),
}


def panel_surveillance_counts(
    atlas: pd.DataFrame,
    countries: frozenset | None = None,
    aliases: dict | None = None,
) -> dict:
    """ATLAS isolate counts per panel pathogen (global by default).

    Maps the ``Species`` string onto the panel via ``aliases`` (default the six-species
    :data:`rd_alignment._PANEL_ALIASES`; pass :data:`EXTENDED_PANEL_ALIASES` for the n=16
    panel) so ATLAS, SPIDAAR-catchment and the burden/funding axes share one taxonomy.
    ``countries`` optionally restricts to a country set.
    """
    aliases = aliases or _PANEL_ALIASES
    df = atlas if countries is None else atlas[atlas["Country"].astype(str).isin(countries)]
    sp = df["Species"].astype(str).str.lower()
    counts: dict = {}
    for pathogen, keys in aliases.items():
        mask = np.zeros(len(df), dtype=bool)
        for a in keys:
            mask |= sp.str.contains(a, regex=False).to_numpy()
        counts[pathogen] = int(mask.sum())
    return counts


def expanded_surveillance_mismatch(
    atlas: pd.DataFrame | None = None,
    counts: dict | None = None,
    floor: float = 0.02,
    draws: int | None = None,
    seed: int | None = None,
) -> dict:
    """Burden<->surveillance mismatch over the FULL n=16 GRAM/ATLAS pathogen panel.

    Same machinery as :func:`surveillance_burden_mismatch` but over 16 pathogens
    (:data:`GRAM_ASSOC_DEATHS_EXTENDED` associated-death burden vs ATLAS isolate counts) — so
    the burden<->surveillance Spearman correlation carries a real bootstrap CI instead of the
    n=6 point estimate. Returns the correlation, per-pathogen log2 medians/CIs, and the
    under-surveilled ranking. Still descriptive, but powered enough to state a direction.
    """
    if counts is None:
        if atlas is None:
            raise ValueError("Pass either an ATLAS frame or precomputed extended-panel counts.")
        counts = panel_surveillance_counts(atlas, aliases=EXTENDED_PANEL_ALIASES)
    if draws is None:
        draws = config.MONTE_CARLO_DRAWS
    if seed is None:
        seed = config.step_seed(4)

    surveillance = {p: float(counts[p]) for p in GRAM_ASSOC_DEATHS_EXTENDED}
    if min(surveillance.values()) <= 0:
        raise ValueError(f"Every panel pathogen needs >0 surveillance isolates: {surveillance}")

    mc = monte_carlo_mismatch_ranking(GRAM_ASSOC_DEATHS_EXTENDED, surveillance,
                                      floor=floor, draws=draws, seed=seed)
    per_pathogen = mc["per_pathogen"]
    ranked = sorted(per_pathogen, key=lambda p: per_pathogen[p]["log2_mismatch_median"],
                    reverse=True)

    median_burden = {p: GRAM_ASSOC_DEATHS_EXTENDED[p][0] for p in GRAM_ASSOC_DEATHS_EXTENDED}
    spearman = spearman_burden_funding(median_burden, surveillance, seed=seed)

    return {
        "n_pathogens": len(GRAM_ASSOC_DEATHS_EXTENDED),
        "floor": floor,
        "draws": int(draws),
        "surveillance_counts": {p: int(counts[p]) for p in GRAM_ASSOC_DEATHS_EXTENDED},
        "spearman": spearman,
        "per_pathogen": per_pathogen,
        "undersurveilled_ranking": ranked,
        "note": ("Burden<->surveillance mismatch over the full n=16 GRAM/ATLAS panel "
                 "(associated-death burden vs ATLAS isolate counts). Powers the n=6 index up: "
                 "the Spearman correlation now carries a bootstrap CI."),
    }


def surveillance_burden_mismatch(
    atlas: pd.DataFrame | None = None,
    counts: dict | None = None,
    burden_metric: str = "assoc_deaths_k",
    floor: float = 0.02,
    draws: int | None = None,
    seed: int | None = None,
) -> dict:
    """Burden<->surveillance mismatch over the six panel species (mirrors the funding index).

    Numerator = verified ``GRAM_BURDEN_2019`` (with 95% UIs, Monte-Carlo'd); denominator =
    ATLAS isolate counts. Positive log2 = UNDER-surveilled relative to burden. Returns the
    per-pathogen log2 medians/CIs + P(most under-surveilled), the ranking, and the Spearman
    summary. Descriptive only (n=6).
    """
    if counts is None:
        if atlas is None:
            raise ValueError("Pass either an ATLAS frame or precomputed panel counts.")
        counts = panel_surveillance_counts(atlas)
    if draws is None:
        draws = config.MONTE_CARLO_DRAWS
    if seed is None:
        seed = config.step_seed(4)

    surveillance = {p: float(counts[p]) for p in GRAM_BURDEN_2019}
    if min(surveillance.values()) <= 0:
        raise ValueError(f"Every panel species needs >0 surveillance isolates: {surveillance}")

    burden_ui = {p: GRAM_BURDEN_2019[p][burden_metric] for p in GRAM_BURDEN_2019}
    mc = monte_carlo_mismatch_ranking(burden_ui, surveillance, floor=floor, draws=draws, seed=seed)
    per_pathogen = mc["per_pathogen"]
    ranked = sorted(per_pathogen, key=lambda p: per_pathogen[p]["log2_mismatch_median"],
                    reverse=True)

    median_burden = {p: GRAM_BURDEN_2019[p][burden_metric][0] for p in GRAM_BURDEN_2019}
    spearman = spearman_burden_funding(median_burden, surveillance, seed=seed)

    return {
        "burden_metric": burden_metric,
        "floor": floor,
        "draws": int(draws),
        "surveillance_counts": {p: int(counts[p]) for p in GRAM_BURDEN_2019},
        "per_pathogen": per_pathogen,
        "undersurveilled_ranking": ranked,
        "spearman": spearman,
        "note": ("Third mismatch axis: same GRAM-2019 burden numerator as Component 4, "
                 "denominator = ATLAS isolate counts. Positive log2 = under-surveilled "
                 "relative to mortality burden. Descriptive (n=6)."),
    }


def _shares(values: dict) -> dict:
    total = float(sum(values.values()))
    return {k: (float(v) / total if total > 0 else float("nan")) for k, v in values.items()}


def three_axis_alignment(
    atlas: pd.DataFrame | None = None,
    counts: dict | None = None,
    burden_metric: str = "assoc_deaths_k",
) -> dict:
    """Burden / funding / surveillance share per pathogen, and who is neglected on all three.

    Burden share = median ``GRAM_BURDEN_2019``; funding share = the in-scope live Hub by-genus
    extract (``RD_HUB_GENUS_SNAPSHOT``, public+phil therapeutics 2017-2023); surveillance
    share = ATLAS isolate counts. Each normalised over the six panel species. A pathogen is
    flagged ``neglected_on_all_axes`` when its burden share is above the panel mean while BOTH
    its funding and surveillance shares are below it — the burden the system under-funds AND
    under-watches.
    """
    if counts is None:
        if atlas is None:
            raise ValueError("Pass either an ATLAS frame or precomputed panel counts.")
        counts = panel_surveillance_counts(atlas)

    species = list(GRAM_BURDEN_2019)
    burden_share = _shares({p: GRAM_BURDEN_2019[p][burden_metric][0] for p in species})
    funding_share = _shares({p: RD_HUB_GENUS_SNAPSHOT["funding_musd"][p] for p in species})
    surveillance_share = _shares({p: float(counts[p]) for p in species})

    mean_share = 1.0 / len(species)
    per_pathogen = {}
    neglected = []
    for p in species:
        flag = (burden_share[p] > mean_share
                and funding_share[p] < mean_share
                and surveillance_share[p] < mean_share)
        per_pathogen[p] = {
            "burden_share": burden_share[p],
            "funding_share": funding_share[p],
            "surveillance_share": surveillance_share[p],
            "neglected_on_all_axes": flag,
        }
        if flag:
            neglected.append(p)

    return {
        "burden_metric": burden_metric,
        "funding_source": RD_HUB_GENUS_SNAPSHOT["source"],
        "per_pathogen": per_pathogen,
        "neglected_on_all_axes": neglected,
        "note": ("Shares normalised over the six panel species on each axis. 'Neglected on "
                 "all axes' = above-mean burden share but below-mean funding AND surveillance "
                 "share."),
    }


def _gini(values: np.ndarray) -> float:
    """Gini concentration of a non-negative vector (0 = even, ->1 = concentrated)."""
    x = np.sort(np.asarray(values, dtype=float))
    n = x.size
    if n == 0 or x.sum() == 0:
        return float("nan")
    cum = np.cumsum(x)
    return float((n + 1 - 2 * np.sum(cum) / cum[-1]) / n)


def geographic_concentration(
    atlas: pd.DataFrame,
    ssa: frozenset = SSA_COUNTRIES,
    catchment: tuple = config.CATCHMENT_COUNTRIES,
) -> dict:
    """Where the surveillance physically is: the SSA blind spot, quantified from ATLAS alone.

    Sub-Saharan Africa carries the highest AMR-attributable death rate of any world region
    (Murray et al., Lancet 2022) yet, as this returns, contributes a tiny share of the global
    surveillance isolates. Reports the SSA and catchment isolate shares, the single largest
    contributor, the top-contributor-to-SSA ratio, and the Gini concentration across
    countries. Self-contained (no external denominator needed).
    """
    country = atlas["Country"].astype(str)
    total = int(len(country))
    vc = country.value_counts()
    ssa_n = int(country.isin(ssa).sum())
    catch_n = int(country.isin(set(catchment)).sum())
    top_country = str(vc.index[0])
    top_n = int(vc.iloc[0])

    return {
        "total_isolates": total,
        "n_countries": int(country.nunique()),
        "ssa": {"isolates": ssa_n, "share": ssa_n / total if total else float("nan")},
        "catchment": {"isolates": catch_n,
                      "share": catch_n / total if total else float("nan")},
        "top_country": {"country": top_country, "isolates": top_n,
                        "share": top_n / total if total else float("nan")},
        "top_country_to_ssa_ratio": (top_n / ssa_n) if ssa_n else float("nan"),
        "gini_country_concentration": _gini(vc.to_numpy()),
        "note": ("SSA bears the highest AMR-attributable death rate of any region (Murray "
                 "2022) yet contributes the listed minority share of global surveillance "
                 "isolates — the surveillance blind spot."),
    }


def run_surveillance_alignment(
    atlas: pd.DataFrame | None = None,
    burden_metric: str = "assoc_deaths_k",
    floor: float = 0.02,
    draws: int | None = None,
    seed: int | None = None,
) -> dict:
    """Bundle the surveillance blind-spot axis for the pipeline (Component 6)."""
    if atlas is None:
        from . import data_loading
        atlas = data_loading.load_atlas_backbone()
    counts = panel_surveillance_counts(atlas)
    counts_ext = panel_surveillance_counts(atlas, aliases=EXTENDED_PANEL_ALIASES)
    return {
        "panel_counts": {p: int(counts[p]) for p in GRAM_BURDEN_2019},
        "burden_surveillance_mismatch": surveillance_burden_mismatch(
            counts=counts, burden_metric=burden_metric, floor=floor, draws=draws, seed=seed),
        "expanded_mismatch": expanded_surveillance_mismatch(
            counts=counts_ext, floor=floor, draws=draws, seed=seed),
        "three_axis": three_axis_alignment(counts=counts, burden_metric=burden_metric),
        "geographic": geographic_concentration(atlas),
    }
