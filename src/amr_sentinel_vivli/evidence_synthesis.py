"""Component 1c (co-primary honesty analysis): Bayesian evidence synthesis.

The Component-1 cohort is underpowered on mortality (21 susceptible patients, one
susceptible death), so a per-patient mortality contrast is not estimable here — that is
*why* the headline moved to excess bed-days. Rather than wave at the literature, this
module performs a **random-effects Bayesian meta-analysis** of the directly comparable,
adjusted SSA/LMIC evidence on the resistance -> in-hospital-death hazard ratio, and places
our own cohort's *crude* estimate alongside it (not inside the pooled adjusted estimate, to
avoid mixing a confounded crude HR with adjusted ones).

What is pooled
--------------
Two external cohorts report an **adjusted** in-hospital-death hazard ratio for 3GC-resistant
vs 3GC-susceptible Enterobacterales infection, on a comparable scale:

* **MBIRA** (Aiken et al., Lancet Infect Dis 2023; 8 SSA hospitals, 878 BSI): ratio of
  cause-specific hazards 0.74 (95% CI 0.42-1.30).
* **Fiji** (Loftus et al., JGAR 2022; Suva): adjusted in-hospital mortality aHR 1.13
  (0.51-2.53).

Our cohort contributes a third, *crude* cause-specific death HR computed from the competing-
risks exit frame. Its 95% interval spans roughly an order of magnitude (one susceptible
death), so it is shown for completeness and as a sensitivity pool — it barely moves the
adjusted pooled estimate — but the **primary synthesis is the two adjusted studies**.

Inference
---------
A standard normal-normal hierarchical model on the log-HR scale,

    y_i ~ Normal(theta_i, s_i^2),   theta_i ~ Normal(mu, tau^2),

with the per-study theta_i integrated out analytically (y_i ~ Normal(mu, s_i^2 + tau^2)).
The posterior over (mu, tau) is evaluated on a deterministic grid (no sampling, fully
reproducible), with weakly-informative priors mu ~ Normal(0, 1) and tau ~ Half-Normal(0.5).
We report the pooled HR with credible interval, P(HR>1), the between-study SD tau, and the
**prediction interval** for a new setting (the honest target for "what does resistance do to
mortality in the next SSA hospital"). A fixed-effect inverse-variance pool is reported as a
classical comparator, and the tau prior scale is varied as a sensitivity.

Pure NumPy/pandas; unit-tested on synthetic inputs. The secure-environment upgrade is a full
MCMC fit (e.g. a Half-Cauchy tau hyperprior in PyMC); with two-to-three studies the posterior
is prior-regularized either way, which is why the prior is stated and varied.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config
from .excess_los import DEFAULT_TAU, build_exit_frame

_EXIT_DEATH = 2

# Verified, adjusted external evidence (in-hospital death, 3GC-R vs 3GC-S, hazard-ratio
# scale). Figures re-checked against the sources in docs/reference_verified_2026-06-06.md.
EXTERNAL_MORTALITY_STUDIES: dict = {
    "MBIRA (Aiken 2023)": {
        "hr": 0.74, "ci_low": 0.42, "ci_high": 1.30, "n": 878, "adjusted": True,
        "estimand": "in-hospital death, ratio of cause-specific hazards (matched/adjusted)",
        "source": "Aiken et al., Lancet Infect Dis 2023; PMC7617135",
    },
    "Fiji (Loftus 2022)": {
        "hr": 1.13, "ci_low": 0.51, "ci_high": 2.53, "n": 162, "adjusted": True,
        "estimand": "in-hospital mortality, adjusted hazard ratio",
        "source": "Loftus et al., JGAR 2022; PMC9452645",
    },
}


def loghr_se_from_ci(
    hr: float, ci_low: float, ci_high: float, z: float = 1.959963985
) -> tuple[float, float]:
    """Back out (log-HR, SE) from a reported HR and its symmetric-on-log-scale 95% CI."""
    loghr = float(np.log(hr))
    se = float((np.log(ci_high) - np.log(ci_low)) / (2 * z))
    return loghr, se


def cohort_death_cause_specific_hr(df: pd.DataFrame, tau: int = DEFAULT_TAU) -> dict:
    """Crude cause-specific death hazard ratio (resistant vs susceptible) from our cohort.

    Uses the competing-risks exit frame: the cause-specific death hazard in each arm is
    deaths / person-time-at-risk to ``tau``, and the log rate-ratio has the standard
    SE = sqrt(1/d_resistant + 1/d_susceptible). This is **crude** (unadjusted for severity),
    so it is reported with that caveat and given essentially zero weight in the adjusted pool
    via its wide interval; it is not a substitute for the adjusted external estimates.
    """
    asc = df[pd.to_numeric(df["resistant"], errors="coerce").isin([0, 1])].copy()
    ef = build_exit_frame(asc, tau)
    ef = ef.assign(resistant=pd.to_numeric(asc["resistant"], errors="coerce").to_numpy())

    out = {}
    rate = {}
    for arm, key in ((1, "resistant"), (0, "susceptible")):
        s = ef[ef["resistant"] == arm]
        deaths = int((s["exit_type"] == _EXIT_DEATH).sum())
        person_time = float(s["time"].sum())
        out[f"deaths_{key}"] = deaths
        out[f"person_days_{key}"] = person_time
        rate[arm] = deaths / person_time if person_time > 0 else np.nan

    d_r = out["deaths_resistant"]
    d_s = out["deaths_susceptible"]
    loghr = float(np.log(rate[1] / rate[0]))
    se = float(np.sqrt(1.0 / d_r + 1.0 / d_s)) if d_r > 0 and d_s > 0 else float("inf")
    z = 1.959963985
    out.update({
        "hr": float(np.exp(loghr)),
        "loghr": loghr,
        "se": se,
        "ci_low": float(np.exp(loghr - z * se)),
        "ci_high": float(np.exp(loghr + z * se)),
        "adjusted": False,
        "estimand": "in-hospital death, CRUDE cause-specific hazard ratio (unadjusted)",
        "caveat": ("Crude and severity-confounded (resistant patients are sicker); the wide "
                   "interval reflects one susceptible death. Shown for completeness, not pooled "
                   "into the primary adjusted estimate."),
    })
    return out


def _as_y_s(studies: dict) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Pack a studies dict into (log-HR, SE, labels), deriving SE from CI where needed."""
    labels, y, s = [], [], []
    for name, rec in studies.items():
        if "loghr" in rec and "se" in rec:
            loghr, se = float(rec["loghr"]), float(rec["se"])
        else:
            loghr, se = loghr_se_from_ci(rec["hr"], rec["ci_low"], rec["ci_high"])
        labels.append(name)
        y.append(loghr)
        s.append(se)
    return np.array(y), np.array(s), labels


def fixed_effect_pool(studies: dict) -> dict:
    """Inverse-variance fixed-effect pooled log-HR (classical comparator)."""
    y, s, _ = _as_y_s(studies)
    w = 1.0 / s**2
    mu = float(np.sum(w * y) / np.sum(w))
    se = float(np.sqrt(1.0 / np.sum(w)))
    z = 1.959963985
    return {
        "pooled_loghr": mu,
        "pooled_hr": float(np.exp(mu)),
        "hr_ci": [float(np.exp(mu - z * se)), float(np.exp(mu + z * se))],
        "se": se,
        "n_studies": len(y),
    }


def bayes_random_effects(
    studies: dict,
    mu_prior_sd: float = 1.0,
    tau_prior_scale: float = 0.5,
    n_mu: int = 601,
    n_tau: int = 400,
    mu_lim: float = 3.0,
    tau_lim: float = 2.5,
) -> dict:
    """Random-effects Bayesian meta-analysis on the log-HR scale (deterministic grid).

    Priors: ``mu ~ Normal(0, mu_prior_sd)``, ``tau ~ Half-Normal(tau_prior_scale)``. The
    per-study random effect is integrated out analytically, so the joint posterior over
    (mu, tau) is evaluated exactly on a grid. Returns the pooled HR + 95% credible interval,
    P(HR>1), the posterior median between-study SD tau, and the **prediction interval** for a
    new study (the mixture predictive over the (mu, tau) posterior).
    """
    y, s, labels = _as_y_s(studies)
    mu_grid = np.linspace(-mu_lim, mu_lim, n_mu)
    tau_grid = np.linspace(1e-3, tau_lim, n_tau)  # start >0 to keep the predictive proper

    # log joint over the grid: prior(mu) + prior(tau) + sum_i Normal(y_i; mu, s_i^2 + tau^2)
    MU, TAU = np.meshgrid(mu_grid, tau_grid, indexing="ij")  # (n_mu, n_tau)
    log_prior_mu = -0.5 * (MU / mu_prior_sd) ** 2
    log_prior_tau = -0.5 * (TAU / tau_prior_scale) ** 2  # half-normal (tau>=0 by grid)
    log_lik = np.zeros_like(MU)
    for yi, si in zip(y, s, strict=True):
        var = si**2 + TAU**2
        log_lik += -0.5 * np.log(2 * np.pi * var) - 0.5 * (yi - MU) ** 2 / var
    log_post = log_prior_mu + log_prior_tau + log_lik
    post = np.exp(log_post - log_post.max())
    post /= post.sum()

    # Marginal posterior of mu
    p_mu = post.sum(axis=1)
    p_mu /= p_mu.sum()
    cdf_mu = np.cumsum(p_mu)

    def q_mu(p):
        return float(np.interp(p, cdf_mu, mu_grid))

    mu_mean = float(np.sum(mu_grid * p_mu))
    p_hr_gt1 = float(p_mu[mu_grid > 0].sum())

    # Posterior median between-study SD tau
    p_tau = post.sum(axis=0)
    p_tau /= p_tau.sum()
    tau_median = float(np.interp(0.5, np.cumsum(p_tau), tau_grid))

    # Mixture predictive for a new study's true log-HR: integrate Normal(theta; mu, tau^2)
    # over the (mu, tau) posterior, evaluated on a fine theta grid.
    theta = np.linspace(-mu_lim - tau_lim, mu_lim + tau_lim, 2 * n_mu)
    var_t = (TAU**2).ravel()
    pred = (post.ravel()[:, None]
            * np.exp(-0.5 * (theta[None, :] - MU.ravel()[:, None]) ** 2 / var_t[:, None])
            / np.sqrt(2 * np.pi * var_t[:, None])).sum(axis=0)
    pred /= pred.sum()
    cdf_pred = np.cumsum(pred)

    def q_pred(p):
        return float(np.interp(p, cdf_pred, theta))

    z = 1.959963985  # noqa: F841  (kept for parity with CI helpers; grid quantiles used here)
    return {
        "labels": labels,
        "n_studies": len(y),
        "pooled_loghr_mean": mu_mean,
        "pooled_loghr_median": q_mu(0.5),
        "pooled_hr": float(np.exp(q_mu(0.5))),
        "hr_ci": [float(np.exp(q_mu(0.025))), float(np.exp(q_mu(0.975)))],
        "p_hr_gt_1": p_hr_gt1,
        "tau_median": tau_median,
        "prediction_hr_ci": [float(np.exp(q_pred(0.025))), float(np.exp(q_pred(0.975)))],
        "priors": {"mu_prior_sd": mu_prior_sd, "tau_prior_scale": tau_prior_scale},
    }


def run_evidence_synthesis(df: pd.DataFrame, tau: int = DEFAULT_TAU) -> dict:
    """Assemble Component-1c: the adjusted SSA/LMIC mortality synthesis + our cohort.

    Primary = Bayesian random-effects pool of the two adjusted external cohorts. Also reports
    our crude cohort HR, a sensitivity pool that adds it, the fixed-effect comparator, and a
    tau-prior sensitivity. Returns a single dict for the pipeline and report.
    """
    cohort = cohort_death_cause_specific_hr(df, tau)
    cohort_study = {
        "loghr": cohort["loghr"], "se": cohort["se"], "hr": cohort["hr"],
        "ci_low": cohort["ci_low"], "ci_high": cohort["ci_high"],
    }

    primary = bayes_random_effects(EXTERNAL_MORTALITY_STUDIES)
    with_cohort_studies = {**EXTERNAL_MORTALITY_STUDIES, "This cohort (crude)": cohort_study}
    with_cohort = bayes_random_effects(with_cohort_studies)

    tau_sensitivity = {
        f"tau_scale={ts}": bayes_random_effects(EXTERNAL_MORTALITY_STUDIES, tau_prior_scale=ts)
        for ts in (0.25, 0.5, 1.0)
    }

    return {
        "estimand": "Resistance -> in-hospital death hazard ratio (3GC-R vs 3GC-S), SSA/LMIC",
        "cohort_crude": cohort,
        "primary_pool_adjusted": primary,
        "sensitivity_pool_with_cohort": with_cohort,
        "fixed_effect_adjusted": fixed_effect_pool(EXTERNAL_MORTALITY_STUDIES),
        "tau_prior_sensitivity": tau_sensitivity,
        "interpretation": (
            "The pooled adjusted in-hospital-death HR for resistance is null (credible interval "
            "spans 1). Our cohort's crude HR is consistent but uninformative (severity-confounded, "
            "order-of-magnitude interval) and barely moves the pool. This synthesis is why the "
            "headline contribution is the competing-risks excess-bed-days estimand, not mortality."
        ),
        "seed_note": "Deterministic grid posterior; no Monte-Carlo seed required.",
        "master_seed": config.MASTER_SEED,
    }
