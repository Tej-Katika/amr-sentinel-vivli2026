"""Component 1b (co-primary honesty analyses) for the excess-LOS null.

At ~21 confirmed-susceptible patients the headline Δ-RMST (excess bed-days) is
precision-starved, and the dominant validity threat is *exposure-ascertainment
selection* (only 156 of 336 patients have an ascertained resistant/susceptible status;
178 are unascertained, ``amrp == -1``). The plan (``docs/analysis_plan_2026.md`` §1b)
makes two analyses CO-PRIMARY — not deferrable bullets — because at this n they decide
whether the headline exists:

(a) **Power / precision simulation** — under the realized design (n_resistant≈135,
    n_susceptible≈21) what Δ-RMST 95% CI width and power can the cohort actually deliver,
    and how much would the pending re-linkage (more susceptibles) tighten it?

(b) **Ascertainment-selection sensitivity** — does the complete-case contrast survive
    (i) reweighting the ascertained to the full population (MAR on observed covariates),
    and (ii) extreme / tipping-point assignment of the unascertained patients' exposure?
    The unascertained have OBSERVED outcomes (bed-days) — only their *exposure* is
    missing — so this is exposure-MNAR, not outcome-missingness. The key defense: the
    feared mechanism (longer-stayers are cultured → ascertained-resistant → a *spurious
    positive* excess-LOS) biases TOWARD positive; the complete-case result is NEGATIVE,
    so that bias cannot manufacture it.

Hand-rolled in numpy (the local env lacks lifelines/scipy/pymc); reuses the tested
``excess_los`` competing-risks machinery (``build_exit_frame`` / ``km_rmst``).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config
from .excess_los import DEFAULT_TAU, bin_severity, build_exit_frame, excess_los, km_rmst

_Z_975 = 1.959963985


# --- (a) Power / precision simulation ----------------------------------------

def _rate_for_rmst(target: float, tau: int) -> float:
    """Exponential rate λ whose restricted mean to ``tau`` equals ``target``.

    RMST(λ) = (1 − e^{−λτ}) / λ is strictly decreasing from τ (λ→0) to 0 (λ→∞), so a
    unique λ is found by bisection for any ``target`` in (0, τ).
    """
    if not 0.0 < target < tau:
        raise ValueError(f"target RMST must be in (0, {tau}); got {target}")
    lo, hi = 1e-9, 100.0

    def rmst(lam: float) -> float:
        return (1.0 - np.exp(-lam * tau)) / lam

    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if rmst(mid) > target:   # too much time in state -> need more hazard
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _simulate_arm_rmst(rate: float, n: int, tau: int, rng) -> float:
    """One simulated arm: n exponential exit times, admin-censored at τ, KM restricted mean."""
    t = rng.exponential(1.0 / rate, size=n)
    exited = t <= tau
    time = np.minimum(t, float(tau))
    return km_rmst(time, exited, tau)


def simulate_rmst_precision(
    rmst_baseline: float,
    n_resistant: int,
    n_susceptible: int,
    true_effects=(0.0, 1.0, 1.65, 2.6, 3.0),
    tau: int = DEFAULT_TAU,
    n_sim: int = 1000,
    seed: int | None = None,
) -> dict:
    """Δ-RMST precision/power under the realized design (single-cause exponential proxy).

    The susceptible arm's time-to-exit is modelled as exponential calibrated to
    ``rmst_baseline``; the resistant arm is calibrated to ``rmst_baseline + δ`` for each
    ``δ`` in ``true_effects``. For each δ we simulate ``n_sim`` two-arm datasets, compute
    Δ̂-RMST per sim, and report the sampling SE, the achievable 95% CI width
    (2·1.96·SE), the power to exclude 0, and the mean Δ̂ (a bias check). A competing-risks
    DGP would change details but not the binding lesson — that the **susceptible arm size
    drives precision**; this proxy is transparent and dependency-free.

    Returns ``{"design": ..., "by_effect": [ {effect, se, ci_width, power, mean_delta}, ... ]}``.
    """
    if seed is None:
        seed = config.step_seed(1)
    rng = np.random.default_rng(seed)
    rate_s = _rate_for_rmst(rmst_baseline, tau)

    by_effect = []
    for delta in true_effects:
        rate_r = _rate_for_rmst(rmst_baseline + delta, tau)
        deltas = np.empty(n_sim, dtype=float)
        for i in range(n_sim):
            deltas[i] = (_simulate_arm_rmst(rate_r, n_resistant, tau, rng)
                         - _simulate_arm_rmst(rate_s, n_susceptible, tau, rng))
        se = float(deltas.std(ddof=1))
        power = float(np.mean(np.abs(deltas) > _Z_975 * se)) if se > 0 else float("nan")
        by_effect.append({
            "effect": float(delta),
            "se": se,
            "ci_width": 2.0 * _Z_975 * se,
            "power": power,
            "mean_delta": float(deltas.mean()),
        })
    return {
        "design": {
            "rmst_baseline": float(rmst_baseline),
            "n_resistant": int(n_resistant),
            "n_susceptible": int(n_susceptible),
            "tau": tau,
            "n_sim": int(n_sim),
        },
        "by_effect": by_effect,
    }


def relinkage_precision_sweep(
    rmst_baseline: float,
    n_resistant: int,
    susceptible_grid=(21, 50, 100, 156, 204),
    true_effect: float = 2.6,
    tau: int = DEFAULT_TAU,
    n_sim: int = 1000,
    seed: int | None = None,
) -> dict:
    """How much would re-linkage (a larger susceptible arm) tighten the CI?

    Fixes a plausible true effect (default 2.6 d, the Fiji excess-LOS) and reports the
    achievable Δ-RMST CI width and power at each candidate susceptible-arm size — the
    addendum's "explicit quantification of how much the pending re-linkage would tighten
    it." ``n_resistant`` is held fixed (re-linkage mainly recovers the scarce arm).
    """
    if seed is None:
        seed = config.step_seed(1)
    rows = []
    for i, n_s in enumerate(susceptible_grid):
        sim = simulate_rmst_precision(
            rmst_baseline, n_resistant, int(n_s), true_effects=(true_effect,),
            tau=tau, n_sim=n_sim, seed=seed + i,
        )
        e = sim["by_effect"][0]
        rows.append({"n_susceptible": int(n_s), "ci_width": e["ci_width"], "power": e["power"]})
    return {"true_effect": float(true_effect), "n_resistant": int(n_resistant), "sweep": rows}


# --- (b) Ascertainment-selection sensitivity ---------------------------------

def _ascertainment_masks(df: pd.DataFrame):
    """(ascertained, unascertained) boolean masks: ascertained = resistant in {0,1}."""
    res = pd.to_numeric(df["resistant"], errors="coerce")
    ascertained = res.isin([0, 1])
    return ascertained, ~ascertained


def ascertainment_comparison(df: pd.DataFrame) -> dict:
    """Descriptive contrast of ascertained vs unascertained patients on observed covariates.

    Surfaces whether ascertainment (having a usable resistant/susceptible status) is
    associated with length-of-stay, the observation window, severity, mortality, or
    country — i.e. whether the complete-case cohort is a selected subgroup. This is the
    evidence base for the selection concern, reported before any reweight.
    """
    ascertained, unasc = _ascertainment_masks(df)
    nobsd = pd.to_numeric(df["days_observed"], errors="coerce")
    sev = pd.to_numeric(df["severity"], errors="coerce")
    dead = pd.to_numeric(df["dead"], errors="coerce")
    # Bed-days to tau is defined for everyone (deaths/censored included), unlike raw `los`
    # which is NaN for deaths — so it is the cleaner occupancy summary for the contrast.
    bed_days = build_exit_frame(df, tau=DEFAULT_TAU)["time"].to_numpy()

    def _summ(mask) -> dict:
        m = mask.to_numpy()
        return {
            "n": int(m.sum()),
            "mean_bed_days_to_tau": float(bed_days[m].mean()),
            "mean_days_observed": float(nobsd[m].mean()),
            "frac_high_severity": float((sev[m] == 3).mean()),
            "mortality_rate": float(dead[m].isin([1, 9]).mean()),
        }

    return {
        "ascertained": _summ(ascertained),
        "unascertained": _summ(unasc),
        "country_ascertainment_rate": {
            str(c): float(ascertained[df["country"] == c].mean())
            for c in sorted(df["country"].dropna().unique())
        },
    }


def weighted_km_rmst(times, exited, weights, tau: int = DEFAULT_TAU) -> float:
    """Weighted Kaplan-Meier restricted mean (area under the weighted survival curve).

    Weighted analogue of ``excess_los.km_rmst``: at each event time the survival drops by
    (weighted events) / (weighted at-risk). Used for inverse-ascertainment-probability
    weighting. Pure, deterministic.
    """
    times = np.asarray(times, dtype=float)
    exited = np.asarray(exited, dtype=bool)
    w = np.asarray(weights, dtype=float)
    if times.size == 0:
        return float("nan")
    capped = times > tau
    t = np.minimum(times, float(tau))
    e = exited & ~capped

    area, surv, prev = 0.0, 1.0, 0.0
    for ti in np.unique(t[e]):
        area += surv * (ti - prev)
        at_risk = w[t >= ti].sum()
        d = w[(t == ti) & e].sum()
        surv *= 1.0 - d / at_risk
        prev = ti
    area += surv * (float(tau) - prev)
    return float(area)


def ascertainment_weighted_excess_los(
    df: pd.DataFrame,
    covariates: tuple[str, ...] = ("severity_bin", "country"),
    tau: int = DEFAULT_TAU,
) -> dict:
    """Inverse-ascertainment-probability-weighted Δ-RMST (MAR-on-covariates sensitivity).

    Estimates P(ascertained | stratum) empirically for each ``covariates`` stratum
    (denominator = all patients incl. unascertained), weights each ascertained patient by
    1/P, and recomputes the weighted Δ-RMST — reweighting the complete-case cohort up to
    the full population under MAR given the strata. Compared with the unweighted
    complete-case Δ; agreement means observed-covariate ascertainment selection is not
    driving the headline. Strata with no ascertained patients are skipped (reported).
    """
    asc_mask, _ = _ascertainment_masks(df)
    work = df.copy()
    work["_sev_bin"] = bin_severity(work["severity"])
    key_cols = ["_sev_bin" if c == "severity_bin" else c for c in covariates]
    work = work.dropna(subset=key_cols)
    asc_mask = asc_mask.loc[work.index]

    # P(ascertained | stratum)
    prop = work.assign(_asc=asc_mask.astype(float)).groupby(key_cols)["_asc"].mean()

    asc = work[asc_mask].copy()
    p = prop.reindex(
        pd.MultiIndex.from_frame(asc[key_cols]) if len(key_cols) > 1 else asc[key_cols[0]]
    ).to_numpy()
    asc = asc.assign(_w=1.0 / p)
    asc = asc[np.isfinite(asc["_w"])]

    ef = build_exit_frame(asc, tau)
    ef = ef.assign(
        resistant=pd.to_numeric(asc["resistant"], errors="coerce").to_numpy(),
        w=asc["_w"].to_numpy(),
    )
    r = ef[ef["resistant"] == 1]
    s = ef[ef["resistant"] == 0]
    rmst_r = weighted_km_rmst(r["time"], r["exited"], r["w"], tau)
    rmst_s = weighted_km_rmst(s["time"], s["exited"], s["w"], tau)

    cc = excess_los(df, tau)
    return {
        "tau": tau,
        "complete_case_excess_los": cc["excess_los"],
        "ascertainment_weighted_excess_los": rmst_r - rmst_s,
        "weighted_rmst_resistant": rmst_r,
        "weighted_rmst_susceptible": rmst_s,
        "n_strata": int(prop.notna().sum()),
        "n_strata_no_ascertained": int((prop == 0).sum()),
    }


def exposure_assignment_bounds(df: pd.DataFrame, tau: int = DEFAULT_TAU) -> dict:
    """Envelope + extreme/MAR scenarios for Δ-RMST over the unascertained exposure.

    The 178 unascertained patients have OBSERVED exit times but unknown exposure. We
    recompute Δ-RMST under exposure assignments of that pool:

    * ``delta_min`` / ``delta_max`` — the exact achievable envelope over *all* assignments.
      Δ-RMST is maximised by sending the longest-staying unascertained to the resistant arm
      and the shortest to susceptible (and vice-versa for the minimum); the optimum is a
      threshold on exit time (exchange argument), so a single scan over cut points is exact.
    * ``all_unascertained_resistant`` / ``all_unascertained_susceptible`` — the two
      one-arm extremes (interpretable reference points inside the envelope).
    * ``mar_anchor`` — a representative (seeded random) split at the ascertained resistant
      fraction: the value expected if exposure is missing-at-random w.r.t. outcome. This is
      the plausible scenario; the envelope endpoints are mathematical bounds, not forecasts.

    The defence the plan rests on: the feared mechanism (longer-stayers cultured →
    ascertained-resistant → spurious POSITIVE excess) biases toward the envelope's upper
    end, yet the complete-case and MAR-anchor values are NEGATIVE.
    """
    asc_mask, unasc_mask = _ascertainment_masks(df)
    asc = df[asc_mask].reset_index(drop=True)
    unasc = df[unasc_mask].reset_index(drop=True)

    asc_ef = build_exit_frame(asc, tau)
    asc_ef = asc_ef.assign(resistant=pd.to_numeric(asc["resistant"], errors="coerce").to_numpy())
    un_ef = build_exit_frame(unasc, tau)  # resistant col is NaN; assignment is what we vary

    base_r = asc_ef[asc_ef["resistant"] == 1]
    base_s = asc_ef[asc_ef["resistant"] == 0]

    def _delta(r_add: pd.DataFrame, s_add: pd.DataFrame) -> float:
        r = pd.concat([base_r, r_add])
        s = pd.concat([base_s, s_add])
        return km_rmst(r["time"], r["exited"], tau) - km_rmst(s["time"], s["exited"], tau)

    empty = un_ef.iloc[0:0]
    n = len(un_ef)
    asc_sorted = un_ef.sort_values("time").reset_index(drop=True)  # ascending exit time

    all_r = _delta(un_ef, empty)
    all_s = _delta(empty, un_ef)

    # Threshold scan: top-c longest -> resistant (maximises Δ); top-c longest -> susceptible
    # (minimises Δ). The global envelope is the extremum over both families.
    deltas = []
    for c in range(n + 1):
        top_to_r = _delta(asc_sorted.iloc[n - c:], asc_sorted.iloc[:n - c])
        top_to_s = _delta(asc_sorted.iloc[:n - c], asc_sorted.iloc[n - c:])
        deltas.extend([top_to_r, top_to_s])

    frac_r = len(base_r) / (len(base_r) + len(base_s))
    rng = np.random.default_rng(config.step_seed(1))
    perm = rng.permutation(n)
    k = int(round(frac_r * n))
    mar = _delta(un_ef.iloc[perm[:k]], un_ef.iloc[perm[k:]])

    return {
        "tau": tau,
        "complete_case_excess_los": excess_los(df, tau)["excess_los"],
        "n_unascertained": int(n),
        "delta_min": float(min(deltas)),
        "delta_max": float(max(deltas)),
        "all_unascertained_resistant": all_r,
        "all_unascertained_susceptible": all_s,
        "mar_anchor": mar,
        "ascertained_resistant_fraction": frac_r,
    }
