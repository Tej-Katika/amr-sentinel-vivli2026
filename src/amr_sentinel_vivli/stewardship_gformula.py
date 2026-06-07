"""Component 5 (centerpiece, EXPLORATORY): empiric-adequacy stewardship g-formula.

The pivoted thesis (``docs/strategy_2026.md``, ``docs/analysis_plan_2026.md``) is the
**AMR burden paradox**: resistance is *not* a clean per-patient killer (triangulated
null — SPIDAAR + MBIRA + Fiji), so the actionable leverage is **systemic** —
specifically whether the patient received *adequate empiric therapy*, not the
resistance phenotype itself. This module makes that leverage quantitative.

Estimand
--------
A point-treatment **g-formula** (g-computation / standardization; Robins 1986) for the
counterfactual effect of empiric-therapy adequacy ``A`` (adequate=1 / inadequate=0) on
restricted-mean time-in-hospital (**bed-days**) to ``tau = 28`` days, with in-hospital
death a competing event for discharge (so the bed-day contrast is read alongside the
death contrast — adequacy that averts death can *raise* occupancy by keeping survivors
to discharge; we never report bed-days without the competing death CIF):

    E[Y | set A=a] = Σ_ℓ P(L=ℓ) · E[Y | A=a, L=ℓ],     (standardization over L)

    avertable_bed_days       = E[Y | set A=inadequate] − E[Y | set A=adequate]
    avertable_vs_natural     = E[Y | natural course]   − E[Y | set A=adequate]

estimated **nonparametrically by stratification** so the positivity assumption is
*visible*: a confounder stratum that lacks either arm is dropped (never silently
imputed) and reported in ``strata_dropped``. The within-stratum bed-day mean reuses the
Component-1 Kaplan-Meier restricted-mean (``excess_los.km_rmst``) and the death CIF
reuses ``excess_los.cif_at_tau`` — the same tested machinery, so this estimator inherits
their competing-risks correctness. A parametric/Bayesian pooled-hazard g-formula (partial
pooling across countries) is the secure-environment upgrade where PyMC is available;
stratification is the honest, dependency-free, hard-to-misread floor.

Collider / mediator discipline (the single subtlest flaw)
---------------------------------------------------------
Resistance largely *determines* adequacy (a resistant isolate is more likely to defeat
empiric therapy), so resistance is a **determinant of the treatment**, not a generic
confounder of the adequacy→bed-days effect; conditioning on it as a confounder would bias
the estimate. Resistance is therefore **excluded from the confounder set L**. The
principled stratum-specific contrast is the controlled direct effect *within* resistance
strata (``gformula_by_resistance``); the pooled estimate is the population effect.

Status: EXPLORATORY, gated on Gate A (``docs/analysis_plan_2026.md`` §3). The treatment is
``treatment_adequacy`` (raw ``txadp``: ~106 adequate / 52 inadequate / 178 unknown). The
**integer codes** of ``txadp`` are assumed here (``TXADP_ADEQUATE_CODE`` /
``TXADP_INADEQUATE_CODE``) and MUST be verified against the secure-environment codebook
before any result is quoted; a one-line constant change fixes them. With 52 inadequate
patients the per-stratum cells are thin and the contrast is wide / prior-sensitive — this
is a what-if calibration tool, not a confirmatory individual-effect claim.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config
from .excess_los import (
    DEFAULT_TAU,
    bin_severity,
    build_exit_frame,
    cif_at_tau,
    km_rmst,
)

# --- Treatment node: empiric-therapy adequacy from ``txadp`` -----------------
# Coding aligned to the documented interpretation (data-utilization-roadmap /
# analysis_plan §3 Gate A): the realized raw ``txadp`` counts are {0: 106, 1: 52,
# 9: 178}, documented as "106 adequate / 52 inadequate / 178 unknown" -> therefore
# adequate = 0, inadequate = 1. Code 9 (and any other value / NaN) is the unknown
# sentinel and is excluded from the contrast, mirroring the ``amrp`` convention in
# ``data_loading._RESISTANT_FROM_AMRP``.
#
# >>> GATE A (the single most important codebook check) <<< : this code->label
# DIRECTION is not yet confirmed against the official SPIDAAR codebook, and it
# FLIPS THE SIGN of every result here. The clinical coherence check supports it
# (adequate=0 makes adequate therapy AVERT deaths, death CIF 0.071 vs 0.097, the
# expected direction), but verify in the secure env before quoting. Flipping these
# two constants is the one-line switch if the codebook disagrees.
TXADP_ADEQUATE_CODE = 0
TXADP_INADEQUATE_CODE = 1
_ADEQUACY_FROM_TXADP = {TXADP_ADEQUATE_CODE: 1.0, TXADP_INADEQUATE_CODE: 0.0}

ADEQUATE = 1
INADEQUATE = 0

# Parsimonious DAG-justified confounder set L (NOT resistance — see module docstring).
# Default is deliberately small: with 52 inadequate patients, every extra covariate
# empties cells and drops strata. ``severity`` is the dominant confounder of both
# adequacy and recovery; ``country`` carries case-mix / facility differences.
DEFAULT_GFORMULA_COVARIATES = ("severity_bin", "country")

# --- WHO-CHOICE inpatient bed-day unit costs (docs/reference_verified_2026-06-06.md) ---
# Primary-hospital tier (where most excess bed-days accrue). The "hotel" component
# includes personnel/capital/lab/food but EXCLUDES drugs & diagnostics — the correct
# denominator for bed-days (drug costs modelled separately, avoiding double-counting).
# Base case: 2010 USD (market FX). Uncertainty: PPP int-$ median + 95% UI (for the
# Monte-Carlo sensitivity). Source: WHO-CHOICE 2010 I$; Stenberg et al. 2018.
BED_DAY_COST_USD_2010 = {"Ghana": 6.30, "Kenya": 5.45, "Malawi": 3.25, "Uganda": 3.81}
BED_DAY_COST_PPP_INTL = {  # (median, lo95, hi95) PPP international dollars, 2010
    "Ghana": (14.28, 5.94, 32.50),
    "Kenya": (14.03, 5.54, 30.87),
    "Malawi": (6.53, 2.58, 14.13),
    "Uganda": (10.25, 4.08, 23.83),
}

_Z_975 = 1.959963985  # standard-normal 97.5th percentile, for the lognormal 95% UI fit


def binarize_adequacy(txadp) -> pd.Series:
    """Map raw ``txadp`` to the binary treatment {adequate=1, inadequate=0, NaN}.

    Returns a float Series (NaN where unascertained), index-aligned with the input.
    The integer codes are assumed (``TXADP_*_CODE``) and must be codebook-verified
    (Gate A) before any result is quoted — see the module docstring.
    """
    s = pd.to_numeric(pd.Series(txadp), errors="coerce")
    return s.map(_ADEQUACY_FROM_TXADP)


def _covariate_values(sub: pd.DataFrame, name: str) -> np.ndarray:
    """Resolve one confounder column to its standardization key values.

    ``severity_bin`` collapses ``disev`` to low/high (reusing ``excess_los.bin_severity``,
    which keeps both arms populated where a 3-level split would empty a cell); any other
    name is read straight from the frame as a categorical key.
    """
    if name == "severity_bin":
        return bin_severity(sub["severity"])
    col = sub[name]
    return col.astype("string").where(col.notna(), other=None).to_numpy()


def _build_gformula_frame(
    df: pd.DataFrame,
    covariates: tuple[str, ...],
    tau: int,
) -> pd.DataFrame:
    """Adequacy-ascertained patients as a stratification-ready frame.

    Restricts to patients with ascertained adequacy (``binarize_adequacy`` in {0,1}),
    builds the enrolment-origin time-to-exit frame (``excess_los.build_exit_frame``), and
    attaches the binary ``adequacy`` treatment, the ``resistant`` flag, and each
    standardization covariate. Rows missing any covariate key are dropped (cannot be
    standardized) — the count is recoverable as ``n_ascertained - len(returned)``.
    """
    adeq_all = binarize_adequacy(df["treatment_adequacy"])
    asc = df[adeq_all.isin([0.0, 1.0])].reset_index(drop=True)

    ef = build_exit_frame(asc, tau)  # pid, country, resistant, time, exit_type, exited
    ef = ef.assign(
        adequacy=binarize_adequacy(asc["treatment_adequacy"]).astype(int).to_numpy(),
        resistant=pd.to_numeric(asc["resistant"], errors="coerce").to_numpy(),
    )
    for name in covariates:
        ef[name] = _covariate_values(asc, name)

    return ef.dropna(subset=list(covariates))


def gformula_bed_days(
    df: pd.DataFrame,
    covariates: tuple[str, ...] = DEFAULT_GFORMULA_COVARIATES,
    tau: int = DEFAULT_TAU,
) -> dict:
    """Standardization g-formula: counterfactual bed-days under set-adequacy regimes.

    Computes, by stratifying on ``covariates`` and standardizing to the cohort's
    confounder distribution, the counterfactual restricted-mean bed-days under
    "set adequacy = adequate for all" vs "set inadequate for all" vs the natural course,
    plus the competing death CIF at ``tau`` under each set-regime. Strata lacking either
    arm are dropped (positivity made visible) and listed in ``strata_dropped``.

    Returns a dict with the regime means, ``avertable_bed_days`` (= inadequate − adequate,
    bed-days saved by guaranteeing adequacy), ``avertable_vs_natural``, the death-CIF
    contrast, and the per-stratum table. The bed-day contrast MUST be read with the death
    contrast (competing risks): averting death can raise occupancy.
    """
    w = _build_gformula_frame(df, covariates, tau)
    n_total = int(len(w))

    per_stratum: dict = {}
    supported: list = []
    dropped: list = []
    for key, g in w.groupby(list(covariates)):
        key = key if isinstance(key, tuple) else (key,)
        adeq = g[g["adequacy"] == ADEQUATE]
        inad = g[g["adequacy"] == INADEQUATE]
        rec = {
            "n": int(len(g)),
            "n_adequate": int(len(adeq)),
            "n_inadequate": int(len(inad)),
        }
        if len(adeq) >= 1 and len(inad) >= 1:
            rec["rmst_adequate"] = km_rmst(adeq["time"], adeq["exited"], tau)
            rec["rmst_inadequate"] = km_rmst(inad["time"], inad["exited"], tau)
            rec["death_cif_adequate"] = cif_at_tau(
                adeq["time"], adeq["exit_type"], tau)["death_cif"]
            rec["death_cif_inadequate"] = cif_at_tau(
                inad["time"], inad["exit_type"], tau)["death_cif"]
            supported.append(key)
        else:
            dropped.append(key)
        per_stratum[key] = rec

    support_n = sum(per_stratum[k]["n"] for k in supported)
    if support_n == 0:
        raise ValueError(
            "No confounder stratum has both an adequate and an inadequate patient; the "
            "adequacy contrast is not identified at this covariate resolution. Reduce "
            "`covariates` (positivity) or obtain more data (see docs/analysis_plan_2026.md)."
        )

    def _std(field: str) -> float:
        return sum(per_stratum[k]["n"] / support_n * per_stratum[k][field] for k in supported)

    def _natural(field_a: str, field_i: str) -> float:
        total = 0.0
        for k in supported:
            r = per_stratum[k]
            within = (r["n_adequate"] * r[field_a] + r["n_inadequate"] * r[field_i]) / (
                r["n_adequate"] + r["n_inadequate"]
            )
            total += r["n"] / support_n * within
        return total

    e_y_adequate = _std("rmst_adequate")
    e_y_inadequate = _std("rmst_inadequate")
    e_y_natural = _natural("rmst_adequate", "rmst_inadequate")
    e_death_adequate = _std("death_cif_adequate")
    e_death_inadequate = _std("death_cif_inadequate")

    return {
        "tau": tau,
        "covariates": tuple(covariates),
        "n_ascertained": n_total,
        "n_on_support": int(support_n),
        "n_strata": len(per_stratum),
        "strata_supported": supported,
        "strata_dropped": dropped,
        "bed_days_set_adequate": e_y_adequate,
        "bed_days_set_inadequate": e_y_inadequate,
        "bed_days_natural": e_y_natural,
        "avertable_bed_days": e_y_inadequate - e_y_adequate,
        "avertable_vs_natural": e_y_natural - e_y_adequate,
        "death_cif_set_adequate": e_death_adequate,
        "death_cif_set_inadequate": e_death_inadequate,
        "averted_death_fraction": e_death_inadequate - e_death_adequate,
        "per_stratum": per_stratum,
    }


def gformula_by_resistance(
    df: pd.DataFrame,
    covariates: tuple[str, ...] = DEFAULT_GFORMULA_COVARIATES,
    tau: int = DEFAULT_TAU,
) -> dict:
    """Controlled-direct-effect view: the adequacy g-formula within each resistance arm.

    Because resistance determines adequacy (module docstring), the principled
    stratum-specific contrast is computed *within* resistant and within susceptible
    patients separately. Returns ``{"resistant": <result|None>, "susceptible": <result|None>}``;
    an arm whose contrast is not identified at this resolution returns ``None`` with the
    reason in the companion ``*_error`` key.
    """
    out: dict = {}
    for arm, name in ((1, "resistant"), (0, "susceptible")):
        sub = df[pd.to_numeric(df["resistant"], errors="coerce") == arm]
        try:
            out[name] = gformula_bed_days(sub, covariates, tau)
        except ValueError as exc:
            out[name] = None
            out[f"{name}_error"] = str(exc)
    return out


def positivity_diagnostic(
    df: pd.DataFrame,
    covariates: tuple[str, ...] = DEFAULT_GFORMULA_COVARIATES,
    tau: int = DEFAULT_TAU,
) -> dict:
    """Quantify treatment positivity / overlap before trusting any standardized contrast.

    Reports, per confounder stratum, the adequate/inadequate counts and the empirical
    adequacy propensity P(A=adequate | L); and overall: how many strata lack a positivity
    region (either arm empty), the fraction of ascertained patients sitting in those
    off-support strata, and the propensity range across supported strata. A large
    off-support fraction means the standardized estimate rests on a non-overlapping
    minority — say so rather than reporting a point estimate.
    """
    w = _build_gformula_frame(df, covariates, tau)
    n_total = int(len(w))

    rows: list = []
    n_off_support = 0
    for key, g in w.groupby(list(covariates)):
        key = key if isinstance(key, tuple) else (key,)
        n_a = int((g["adequacy"] == ADEQUATE).sum())
        n_i = int((g["adequacy"] == INADEQUATE).sum())
        on_support = n_a >= 1 and n_i >= 1
        if not on_support:
            n_off_support += n_a + n_i
        rows.append({
            "stratum": key,
            "n": n_a + n_i,
            "n_adequate": n_a,
            "n_inadequate": n_i,
            "propensity": (n_a / (n_a + n_i)) if (n_a + n_i) else float("nan"),
            "on_support": on_support,
        })

    props = [r["propensity"] for r in rows if r["on_support"]]
    return {
        "covariates": tuple(covariates),
        "n_ascertained": n_total,
        "n_strata": len(rows),
        "n_strata_off_support": sum(1 for r in rows if not r["on_support"]),
        "frac_off_support": (n_off_support / n_total) if n_total else float("nan"),
        "propensity_min": float(min(props)) if props else float("nan"),
        "propensity_max": float(max(props)) if props else float("nan"),
        "per_stratum": rows,
    }


# --- Bed-day -> cost translation (WHO-CHOICE) --------------------------------

def scale_to_population(per_patient_bed_days: float, n_by_country: dict) -> dict:
    """Scale a per-patient avertable bed-day effect to per-country cohort totals.

    Applies the pooled per-patient effect uniformly to each country's patient count —
    appropriate for the cohort itself; extrapolation to the wider catchment additionally
    requires catchment admission counts (not in the SPIDAAR data) and is the user's to
    supply via ``n_by_country``. Negative effects (adequacy raises occupancy via averted
    death) propagate as negative totals.
    """
    return {c: per_patient_bed_days * float(n) for c, n in n_by_country.items()}


def cost_of_bed_days(bed_days_by_country: dict, currency: str = "usd2010") -> dict:
    """Point cost of a per-country bed-day count at WHO-CHOICE unit costs.

    ``currency``: ``"usd2010"`` (2010 market-FX USD, base case) or ``"ppp"`` (PPP
    international-$ median). Returns the total plus the per-country breakdown. No
    uncertainty — use ``monte_carlo_cost`` for the sourced range.
    """
    if currency == "usd2010":
        table = BED_DAY_COST_USD_2010
    elif currency == "ppp":
        table = {c: v[0] for c, v in BED_DAY_COST_PPP_INTL.items()}
    else:
        raise ValueError("currency must be 'usd2010' or 'ppp'")

    by_country = {}
    for c, bd in bed_days_by_country.items():
        if c not in table:
            raise KeyError(f"No WHO-CHOICE unit cost for country {c!r}")
        by_country[c] = float(bd) * table[c]
    return {
        "currency": currency,
        "total": float(sum(by_country.values())),
        "by_country": by_country,
    }


def monte_carlo_cost(
    bed_days_by_country: dict,
    draws: int | None = None,
    seed: int | None = None,
) -> dict:
    """Monte-Carlo cost of avertable bed-days, propagating the WHO-CHOICE 95% UI.

    Each country's unit cost is drawn from a lognormal calibrated to the PPP int-$ median
    and 95% uncertainty interval (``BED_DAY_COST_PPP_INTL``); per-country costs are summed
    per draw. Returns the mean/median and the 95% interval of the total. Seeded from
    ``config.step_seed(5)`` for reproducibility (pre-reg §12). Currency = PPP int-$ (2010);
    the 2010-USD base-case point is ``cost_of_bed_days(..., 'usd2010')``.
    """
    if draws is None:
        draws = config.MONTE_CARLO_DRAWS
    if seed is None:
        seed = config.step_seed(5)
    rng = np.random.default_rng(seed)

    total = np.zeros(int(draws), dtype=float)
    for c, bd in bed_days_by_country.items():
        if c not in BED_DAY_COST_PPP_INTL:
            raise KeyError(f"No WHO-CHOICE PPP unit cost for country {c!r}")
        median, lo, hi = BED_DAY_COST_PPP_INTL[c]
        mu = np.log(median)
        sigma = (np.log(hi) - np.log(lo)) / (2.0 * _Z_975)
        unit = rng.lognormal(mu, sigma, size=int(draws))
        total += float(bd) * unit

    lo95, hi95 = np.quantile(total, [0.025, 0.975])
    return {
        "currency": "ppp_intl_2010",
        "draws": int(draws),
        "mean": float(total.mean()),
        "median": float(np.median(total)),
        "ci_lower": float(lo95),
        "ci_upper": float(hi95),
    }


def run_stewardship_gformula(
    df: pd.DataFrame,
    covariates: tuple[str, ...] = DEFAULT_GFORMULA_COVARIATES,
    tau: int = DEFAULT_TAU,
) -> dict:
    """Step-5 entrypoint: positivity diagnostic + pooled and within-resistance g-formula.

    Convenience wrapper that runs the positivity check first (so the contrast is never
    read without its overlap context) and then the pooled and controlled-direct-effect
    (by-resistance) adequacy g-formulas on the SPIDAAR patient cohort. EXPLORATORY and
    Gate-A-gated (see module docstring).
    """
    return {
        "positivity": positivity_diagnostic(df, covariates, tau),
        "pooled": gformula_bed_days(df, covariates, tau),
        "by_resistance": gformula_by_resistance(df, covariates, tau),
    }
