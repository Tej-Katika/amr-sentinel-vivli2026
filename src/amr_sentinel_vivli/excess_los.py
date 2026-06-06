"""Component 1 (primary): resistance-attributable excess length-of-stay (excess bed-days).

Primary estimand (``docs/analysis_plan_2026.md``): the difference in restricted mean
time spent in the **admitted** state to ``tau = 28`` days between resistant and
susceptible HAI patients — i.e. expected **excess bed-days**. A patient stops
accruing bed-days at the FIRST of discharge-alive or in-hospital death (both exit
the admitted state), so expected bed-days to ``tau`` is the restricted mean of the
time-to-exit survival curve, and

    excess_bed_days = RMST_exit(resistant) - RMST_exit(susceptible).

Death is a *competing* event for discharge: the discharge-vs-death split is reported
separately (cumulative incidence, a later build) so we never claim a spurious
"LOS reduction" that is really excess mortality. But the bed-day **occupancy**
number is the time-to-exit restricted mean, and that is the economic headline.

Time origin = ENROLMENT (verified codebook semantics): ``los`` and ``days_to_death``
are admission-origin, ``days_to_enrolment`` (``enrtpt``) is admission->enrolment, and
``days_observed`` (``nobsd``) is the enrolment-origin observation window. We therefore
place everyone at t=0 at enrolment (this sidesteps immortal-time bias: a patient must
survive admission->enrolment to enter the cohort), with
    time-to-discharge = los - enrtpt,   time-to-death = days_to_death - enrtpt,
and administrative censoring at ``min(nobsd, tau)``. A small fraction discharge at or
before enrolment (``los < enrtpt``); their enrolment-origin time is clipped to 0.

CAVEAT (read before quoting any number): the crude contrast here is CONFOUNDED
(sicker-resistant vs less-sick-susceptible) and exposure-ascertainment selected
(longer-stayers are more likely cultured, hence ascertained-resistant — which can
manufacture a positive excess). Severity standardization, the power/precision
simulation, and the ascertainment-selection sensitivity are separate modules
(next build, co-primary per the plan). This module is the crude occupancy estimator
plus its stratified bootstrap interval — the honest starting point, not the final
causal claim.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config

# Exit-state codes for the admitted -> {discharged-alive | died} competing risks.
EXIT_CENSORED = 0   # still admitted at the censoring horizon
EXIT_DISCHARGE = 1  # discharged alive (event of interest for bed-day accrual end)
EXIT_DEATH = 2      # in-hospital death (competing event)

DEFAULT_TAU = 28    # restricted-mean horizon, inside the 28-31d observation window

EXIT_FRAME_COLUMNS = ("pid", "country", "resistant", "time", "exit_type", "exited")


def build_exit_frame(df: pd.DataFrame, tau: int = DEFAULT_TAU) -> pd.DataFrame:
    """Map a loaded SPIDAAR patient frame to the enrolment-origin time-to-exit frame.

    Pure function (no I/O), unit-tested on synthetic rows. Expects the
    ``load_spidaar`` columns ``pid, country, resistant, los, dead, days_to_death,
    days_observed, days_to_enrolment`` and returns ``EXIT_FRAME_COLUMNS``:

    * ``time``      enrolment-origin follow-up, capped at ``min(nobsd, tau)``
    * ``exit_type`` one of ``EXIT_CENSORED / EXIT_DISCHARGE / EXIT_DEATH``
    * ``exited``    True iff the admitted state was left (discharge or death) by ``time``

    Death (``dead`` in {1, 9} with a death day) takes precedence over discharge.
    Patients with neither a death day nor a discharge ``los`` are treated as still
    admitted and censored at the observation horizon. ``los < enrtpt`` is clipped to 0.
    """
    enr = pd.to_numeric(df["days_to_enrolment"], errors="coerce").fillna(0).astype(float)
    dead = pd.to_numeric(df["dead"], errors="coerce")
    dtpta = pd.to_numeric(df["days_to_death"], errors="coerce")
    los = pd.to_numeric(df["los"], errors="coerce")
    nobsd = pd.to_numeric(df["days_observed"], errors="coerce")

    is_death = dead.isin([1, 9]) & dtpta.notna()
    is_discharge = (~is_death) & los.notna()

    # Raw enrolment-origin event time (inf where neither event is recorded).
    t_event = pd.Series(np.inf, index=df.index, dtype=float)
    t_event[is_death] = dtpta[is_death] - enr[is_death]
    t_event[is_discharge] = los[is_discharge] - enr[is_discharge]
    t_event = t_event.clip(lower=0.0)

    raw_type = pd.Series(EXIT_CENSORED, index=df.index, dtype=int)
    raw_type[is_death] = EXIT_DEATH
    raw_type[is_discharge] = EXIT_DISCHARGE

    # Censoring horizon: enrolment-origin observation window, capped at tau.
    horizon = np.minimum(nobsd.fillna(tau).astype(float), float(tau))

    observed = t_event <= horizon
    time = np.where(observed, t_event, horizon)
    exit_type = np.where(observed, raw_type, EXIT_CENSORED)

    out = pd.DataFrame({
        "pid": df["pid"].to_numpy(),
        "country": df["country"].to_numpy(),
        "resistant": pd.to_numeric(df["resistant"], errors="coerce").to_numpy(),
        "time": np.asarray(time, dtype=float),
        "exit_type": np.asarray(exit_type, dtype=int),
    })
    out["exited"] = out["exit_type"] != EXIT_CENSORED
    return out[list(EXIT_FRAME_COLUMNS)]


def km_rmst(times, exited, tau: int = DEFAULT_TAU) -> float:
    """Restricted mean survival time (area under the Kaplan-Meier curve) to ``tau``.

    Here the "event" is leaving the admitted state (discharge OR death), so the
    RMST is the expected number of bed-days accrued to ``tau``. Pure, deterministic.
    """
    times = np.asarray(times, dtype=float)
    exited = np.asarray(exited, dtype=bool)
    if times.size == 0:
        return float("nan")

    # Anything beyond tau is censored at tau (still admitted).
    capped = times > tau
    t = np.minimum(times, float(tau))
    e = exited & ~capped

    area = 0.0
    surv = 1.0
    prev = 0.0
    for ti in np.unique(t[e]):                 # KM only steps at event times
        area += surv * (ti - prev)
        at_risk = np.sum(t >= ti)
        deaths = np.sum((t == ti) & e)
        surv *= 1.0 - deaths / at_risk
        prev = ti
    area += surv * (float(tau) - prev)         # tail at the last survival level
    return float(area)


def excess_los(df: pd.DataFrame, tau: int = DEFAULT_TAU) -> dict:
    """Crude resistance-attributable excess bed-days (Delta restricted-mean to ``tau``).

    Restricts to the ascertained resistant-vs-susceptible cohort (``resistant`` in
    {0, 1}) and returns the per-arm restricted-mean bed-days and their difference.
    The result is the crude (unadjusted) contrast — see the module caveat.
    """
    ascertained = df[pd.to_numeric(df["resistant"], errors="coerce").isin([0, 1])]
    ef = build_exit_frame(ascertained, tau)
    res = ef[ef["resistant"] == 1]
    sus = ef[ef["resistant"] == 0]

    rmst_r = km_rmst(res["time"], res["exited"], tau)
    rmst_s = km_rmst(sus["time"], sus["exited"], tau)
    return {
        "tau": tau,
        "n_resistant": int(len(res)),
        "n_susceptible": int(len(sus)),
        "rmst_resistant": rmst_r,
        "rmst_susceptible": rmst_s,
        "excess_los": rmst_r - rmst_s,
    }


def bootstrap_excess_los_ci(
    df: pd.DataFrame,
    tau: int = DEFAULT_TAU,
    n_boot: int = 2000,
    alpha: float = 0.05,
    seed: int | None = None,
) -> dict:
    """Stratified bootstrap CI for crude excess bed-days.

    Resamples patients with replacement WITHIN each (arm x country) stratum, so both
    the resistant/susceptible arm sizes and the country mix are preserved (important
    with the small ~21-patient susceptible arm). Percentile interval; a BCa
    refinement is a planned enhancement. Seeded for reproducibility (pre-reg §12).
    """
    if seed is None:
        seed = config.step_seed(1)
    rng = np.random.default_rng(seed)

    ascertained = df[pd.to_numeric(df["resistant"], errors="coerce").isin([0, 1])].copy()
    ascertained["resistant"] = pd.to_numeric(ascertained["resistant"], errors="coerce").astype(int)
    strata = [g.index.to_numpy() for _, g in ascertained.groupby(["resistant", "country"])]

    point = excess_los(ascertained, tau)["excess_los"]
    draws = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        picks = np.concatenate([rng.choice(idx, size=idx.size, replace=True) for idx in strata])
        draws[b] = excess_los(ascertained.loc[picks], tau)["excess_los"]

    lo, hi = np.quantile(draws, [alpha / 2, 1 - alpha / 2])
    return {
        "tau": tau,
        "excess_los": point,
        "ci_lower": float(lo),
        "ci_upper": float(hi),
        "n_boot": n_boot,
        "alpha": alpha,
    }


# --- Severity adjustment + competing-risks decomposition (robustness checks) ----

def bin_severity(severity) -> np.ndarray:
    """Collapse disease severity (``disev`` in {1,2,3}) to low {1,2} / high {3}.

    The susceptible arm has no severity-1 patients, so a 3-level direct
    standardization has an empty cell. The binary low/high split keeps every patient
    and leaves both arms populated in each stratum. Values outside {1,2,3} -> NaN.
    """
    s = pd.to_numeric(pd.Series(severity), errors="coerce")
    out = np.where(s.isin([1, 2]), "low", np.where(s == 3, "high", None))
    return out


def standardized_excess_los(df: pd.DataFrame, tau: int = DEFAULT_TAU, strata=None) -> dict:
    """Severity-standardized excess bed-days by direct (g-formula) standardization.

    Computes arm-specific restricted-mean bed-days within each severity stratum, then
    standardizes both arms to the *pooled* stratum distribution of the ascertained
    cohort — removing confounding by the measured severity mix. Strata where either
    arm is empty are dropped (with the weights renormalized) and reported, so the
    standardized contrast is never built on an unsupported cell.
    """
    asc = df[pd.to_numeric(df["resistant"], errors="coerce").isin([0, 1])].reset_index(drop=True)
    ef = build_exit_frame(asc, tau)
    ef = ef.assign(
        resistant=pd.to_numeric(asc["resistant"], errors="coerce").astype(int).to_numpy(),
        stratum=bin_severity(asc["severity"]) if strata is None else np.asarray(strata),
    )
    ef = ef[ef["stratum"].notna()]

    per_stratum: dict = {}
    used: list = []
    for s, g in ef.groupby("stratum"):
        res, sus = g[g["resistant"] == 1], g[g["resistant"] == 0]
        if len(res) >= 1 and len(sus) >= 1:
            per_stratum[s] = {
                "n": int(len(g)), "n_resistant": int(len(res)), "n_susceptible": int(len(sus)),
                "rmst_resistant": km_rmst(res["time"], res["exited"], tau),
                "rmst_susceptible": km_rmst(sus["time"], sus["exited"], tau),
            }
            used.append(s)
    dropped = [s for s in ef["stratum"].unique() if s not in used]

    total = sum(per_stratum[s]["n"] for s in used)
    std_r = sum(per_stratum[s]["n"] / total * per_stratum[s]["rmst_resistant"] for s in used)
    std_s = sum(per_stratum[s]["n"] / total * per_stratum[s]["rmst_susceptible"] for s in used)
    return {
        "tau": tau,
        "standardized_excess_los": std_r - std_s,
        "standardized_rmst_resistant": std_r,
        "standardized_rmst_susceptible": std_s,
        "strata_used": used,
        "strata_dropped": dropped,
        "per_stratum": per_stratum,
    }


def cif_at_tau(time, exit_type, tau: int = DEFAULT_TAU) -> dict:
    """Aalen-Johansen cumulative incidence of discharge and death by ``tau``.

    Proper competing-risks decomposition: each cause's CIF increments by
    ``S(t-) * d_cause(t) / n_at_risk(t)`` where ``S`` is the all-cause
    (leave-admitted) survival. Returns the discharge CIF, death CIF, and the
    residual still-admitted probability at ``tau``.
    """
    time = np.asarray(time, dtype=float)
    et = np.asarray(exit_type, dtype=int)
    surv = 1.0
    cif_disc = cif_death = 0.0
    for ti in np.unique(time[et != EXIT_CENSORED]):
        if ti > tau:
            break
        at_risk = np.sum(time >= ti)
        d_disc = np.sum((time == ti) & (et == EXIT_DISCHARGE))
        d_death = np.sum((time == ti) & (et == EXIT_DEATH))
        cif_disc += surv * d_disc / at_risk
        cif_death += surv * d_death / at_risk
        surv *= 1.0 - (d_disc + d_death) / at_risk
    return {
        "discharge_cif": float(cif_disc),
        "death_cif": float(cif_death),
        "still_admitted": float(1.0 - cif_disc - cif_death),
    }


def cif_decomposition(df: pd.DataFrame, tau: int = DEFAULT_TAU) -> dict:
    """Per-arm competing-risks CIF (discharge vs death) — the mechanism behind the
    bed-day contrast (does resistance shift exits toward death rather than discharge?)."""
    asc = df[pd.to_numeric(df["resistant"], errors="coerce").isin([0, 1])]
    ef = build_exit_frame(asc, tau)
    out = {}
    for arm, name in ((1, "resistant"), (0, "susceptible")):
        a = ef[ef["resistant"] == arm]
        out[name] = {"n": int(len(a)), **cif_at_tau(a["time"], a["exit_type"], tau)}
    return out
