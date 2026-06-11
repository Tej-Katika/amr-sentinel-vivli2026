"""Component 1d: isolate-level resistance-mechanism breakdown + relinkage pre-staging.

The patient frame carries only a single collapsed resistance summary (``amrp``); the 244
delivered isolates carry the *mechanism-resolved* phenotype the summary throws away —
third-generation-cephalosporin resistance (``c3r``), multidrug resistance (``mdr``) and
methicillin resistance (``mrsa``) — plus specimen, HAI site and clinical relevance. This
module wrings that detail dry (per-mechanism Wilson-interval prevalence, co-resistance, and
mechanism-by-organism), and pre-stages the **patient↔isolate relinkage** we have requested
from Vivli (``docs/vivli_relinkage_request.md``): the merge machinery that, once a pid↔iid
crosswalk arrives, restores the registered isolate-level same-species contrast, and a
quantified projection of what that link buys.

Why mechanisms matter to the thesis: the patient-level summary ascertains resistance for only
156/336 patients (135 resistant / 21 susceptible). At the isolate level the *same* phenotype is
ascertained far more completely and with real susceptible representation — 3GC-R in 204 isolates
(23 susceptible) and **MDR in 231 (51 susceptible)** — so the relinkage is the single highest-
value route back to a powered contrast. This module makes that case with numbers rather than an
assertion. Pure NumPy/pandas; unit-tested on synthetic isolates.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .bayesian_projection import _wilson_ci

# The three ascertained resistance mechanisms in the isolate file, with readable labels.
MECHANISMS: dict = {
    "c3r": "Third-generation cephalosporin resistance (3GC-R)",
    "mdr": "Multidrug resistance (MDR)",
    "mrsa": "Methicillin resistance (MRSA)",
}

# Default categorical stratifiers for the mechanism panel (all present in ISOLATE_COLUMNS).
DEFAULT_STRATIFIERS: tuple[str, ...] = (
    "country", "infection_site", "specimen", "clinical_relevance", "organism_group",
)

# Patient-level ascertained arm sizes (the deviation cohort) the relinkage would improve on.
# Source: load_spidaar() realized arms (resistant 135 / susceptible 21); see project status.
PATIENT_ARMS_DEFAULT: tuple[int, int] = (135, 21)


def _ascertained(isolates: pd.DataFrame, mechanism: str) -> np.ndarray:
    """The 0/1 ascertained values of one mechanism flag (NaN sentinels dropped)."""
    if mechanism not in MECHANISMS:
        raise ValueError(f"mechanism must be one of {tuple(MECHANISMS)}; got {mechanism!r}")
    s = pd.to_numeric(isolates[mechanism], errors="coerce")
    return s[s.isin([0, 1])].to_numpy()


def mechanism_prevalence(isolates: pd.DataFrame, mechanism: str, by: str | None = None) -> dict:
    """Wilson-interval prevalence of one mechanism, overall or by a stratifier.

    Restricts to isolates *ascertained* for that mechanism (flag in {0,1}); the sentinel
    "not tested / not applicable" codes (e.g. 3GC-R on a Gram-positive) are excluded from the
    denominator, never counted as susceptible. Returns prevalence + 95% Wilson CI + counts,
    overall (``by=None``) or one record per stratum.
    """
    if by is None:
        v = _ascertained(isolates, mechanism)
        r, n = int(v.sum()), int(v.size)
        lo, hi = _wilson_ci(r, n)
        return {"mechanism": mechanism, "n_ascertained": n, "n_resistant": r,
                "prevalence": (r / n) if n else float("nan"), "ci_lower": lo, "ci_upper": hi}

    flag = pd.to_numeric(isolates[mechanism], errors="coerce")
    asc = isolates.assign(_f=flag)[flag.isin([0, 1])]
    out: dict = {}
    for key, g in asc.groupby(by, dropna=False):
        r, n = int(g["_f"].sum()), int(len(g))
        lo, hi = _wilson_ci(r, n)
        out[str(key)] = {"n_ascertained": n, "n_resistant": r,
                         "prevalence": r / n, "ci_lower": lo, "ci_upper": hi}
    return {"mechanism": mechanism, "by": by, "strata": out}


def mechanism_panel(
    isolates: pd.DataFrame,
    stratifiers: tuple[str, ...] = DEFAULT_STRATIFIERS,
) -> dict:
    """All three mechanisms: overall Wilson-interval prevalence + each stratification."""
    panel: dict = {}
    for mech in MECHANISMS:
        rec = {"label": MECHANISMS[mech], "overall": mechanism_prevalence(isolates, mech)}
        rec["by"] = {
            s: mechanism_prevalence(isolates, mech, by=s)["strata"]
            for s in stratifiers if s in isolates.columns
        }
        panel[mech] = rec
    return panel


def mechanism_cooccurrence(isolates: pd.DataFrame, a: str = "c3r", b: str = "mdr") -> dict:
    """Joint distribution of two mechanisms among isolates ascertained for BOTH.

    Returns the 2x2 counts and the conditional co-resistance fractions (P(b|a) and P(a|b)) —
    e.g. how much 3GC-R travels with MDR. Restricted to the doubly-ascertained subset so the
    denominators are honest.
    """
    fa = pd.to_numeric(isolates[a], errors="coerce")
    fb = pd.to_numeric(isolates[b], errors="coerce")
    both = fa.isin([0, 1]) & fb.isin([0, 1])
    xa, xb = fa[both].to_numpy().astype(int), fb[both].to_numpy().astype(int)
    n11 = int(((xa == 1) & (xb == 1)).sum())
    n10 = int(((xa == 1) & (xb == 0)).sum())
    n01 = int(((xa == 0) & (xb == 1)).sum())
    n00 = int(((xa == 0) & (xb == 0)).sum())
    n_a, n_b = n11 + n10, n11 + n01
    return {
        "a": a, "b": b, "n_both_ascertained": int(both.sum()),
        "counts": {"a1_b1": n11, "a1_b0": n10, "a0_b1": n01, "a0_b0": n00},
        "p_b_given_a": (n11 / n_a) if n_a else float("nan"),
        "p_a_given_b": (n11 / n_b) if n_b else float("nan"),
    }


def mechanism_by_organism(
    isolates: pd.DataFrame,
    group_col: str = "organism_group",
) -> dict:
    """Per-organism-group prevalence of each mechanism (which bugs carry which resistance)."""
    return {
        mech: mechanism_prevalence(isolates, mech, by=group_col)["strata"]
        for mech in MECHANISMS
    }


# --- Relinkage pre-staging ------------------------------------------------------

def relink_isolates_to_patients(
    isolates: pd.DataFrame,
    patients: pd.DataFrame,
    crosswalk: pd.DataFrame,
    isolate_key: str = "iid",
    patient_key: str = "pid",
) -> pd.DataFrame:
    """Join isolate phenotypes to patient outcomes via a pid↔iid crosswalk.

    The machinery the registered isolate-level same-species contrast needs, ready for the
    moment Vivli supplies the link. ``crosswalk`` must carry both keys (one row per isolate);
    the result is one row per *linked* isolate carrying its mechanism flags AND its patient's
    outcome/covariates — the confirmatory analysis frame. Pure; raises if a key is missing.
    Until the crosswalk exists this is exercised only on synthetic data (no positional merge is
    ever attempted — the counts do not reconcile; see ``docs/vivli_relinkage_request.md``).
    """
    for name, frame, key in (("crosswalk", crosswalk, isolate_key),
                             ("crosswalk", crosswalk, patient_key),
                             ("isolates", isolates, isolate_key),
                             ("patients", patients, patient_key)):
        if key not in frame.columns:
            raise KeyError(f"{name} is missing the join key {key!r}")
    cw = crosswalk[[isolate_key, patient_key]].drop_duplicates(subset=isolate_key)
    linked = isolates.merge(cw, on=isolate_key, how="inner")
    return linked.merge(patients, on=patient_key, how="left", suffixes=("_iso", "_pat"))


def relinkage_recovery_projection(
    isolates: pd.DataFrame,
    patient_arms: tuple[int, int] = PATIENT_ARMS_DEFAULT,
) -> dict:
    """Quantify what the patient↔isolate link buys, per mechanism — arm sizes + precision.

    For each mechanism, reports the isolate-level ascertained resistant/susceptible counts
    (the arms the registered same-species contrast would use) against the current patient-level
    arms. ``precision_gain_factor`` is the ratio of contrast standard errors,
    ``sqrt[(1/n_res^pat + 1/n_sus^pat) / (1/n_res^iso + 1/n_sus^iso)]`` — how much narrower the
    interval becomes from the better-ascertained arms alone (>1 means tighter). The headline:
    MDR ascertainment recovers the susceptible arm the patient summary cannot.
    """
    n_res_pat, n_sus_pat = patient_arms
    base_var = 1.0 / n_res_pat + 1.0 / n_sus_pat
    out: dict = {}
    for mech in MECHANISMS:
        v = _ascertained(isolates, mech)
        n_res, n_sus = int(v.sum()), int((v == 0).sum())
        iso_var = (1.0 / n_res + 1.0 / n_sus) if (n_res and n_sus) else float("inf")
        out[mech] = {
            "label": MECHANISMS[mech],
            "n_ascertained": int(v.size),
            "n_resistant": n_res,
            "n_susceptible": n_sus,
            "susceptible_gain_vs_patient": n_sus - n_sus_pat,
            "precision_gain_factor": float(np.sqrt(base_var / iso_var)),
        }
    return {
        "patient_arms": {"n_resistant": n_res_pat, "n_susceptible": n_sus_pat},
        "by_mechanism": out,
        "note": ("Isolate-level ascertainment of the SAME phenotype is far more complete than "
                 "the patient summary; MDR in particular recovers a susceptible arm "
                 f"({out['mdr']['n_susceptible']}) the patient frame cannot ({n_sus_pat}). "
                 "Requires the pid↔iid crosswalk (docs/vivli_relinkage_request.md)."),
    }


def run_isolate_mechanisms(isolates: pd.DataFrame) -> dict:
    """Component-1d entrypoint: mechanism panel + co-resistance + relinkage projection."""
    return {
        "panel": mechanism_panel(isolates),
        "cooccurrence_c3r_mdr": mechanism_cooccurrence(isolates, "c3r", "mdr"),
        "by_organism": mechanism_by_organism(isolates),
        "relinkage_projection": relinkage_recovery_projection(isolates),
    }
