"""Component 2 (secondary robustness companion): Bayesian partial-pooled excess bed-days.

A secondary, prior-regularized companion to the Component-1 crude excess-LOS — NOT the
headline (``docs/analysis_plan_2026.md`` Component 2). It re-estimates the same
excess-bed-day Δ-RMST through a **competing-risks piecewise-exponential** model with a
partially-pooled resistance effect, and reports it as posterior probabilities —
P(excess > 0) and P(excess > 1 day) — so the small-sample uncertainty is explicit.

Model
-----
Person-time is expanded into intervals on day bins (1-2 / 3-7 / 8-14 / 15-28, respecting
the ≥48h HAI window and the small counts). In each interval the discharge and death
cause-specific event counts are Poisson with an exposure-time offset:

    discharge_count ~ Poisson(exposure · exp(η_dis)),   death_count ~ Poisson(... η_dth),
    η = bin-intercept + β_resistant·R + β_gram·1[Gram-neg] + β_sev·1[high severity] + u_country.

* **Partial pooling on country only** (u_c ~ Normal(0, σ_country)); **pathogen enters as a
  FIXED Gram-negative-vs-other stratifier**, never a random effect (so the thin Gram-positive
  cell never borrows a Gram-negative-anchored prior).
* **Primary prior is NULL-CENTERED**: β_resistant on the discharge log-hazard ~ Normal(0, 0.5).
  The cited discharge-hazard evidence is null-to-reversed (MBIRA discharge CSH ratio 1.16;
  Fiji 0.99), so an optimistic stay-prolonging prior is unjustified. Informative LMIC priors
  are a pre-registered SENSITIVITY only (``prior_sensitivity``).

Inference here is a **Laplace approximation**: the penalized-Poisson MAP (Newton; the Gaussian
priors act as the penalty) plus the posterior covariance from the inverse Hessian. Coefficient
draws are pushed through the Aalen-Johansen restricted-mean map (g-computation over the cohort's
covariate distribution) to give the Δ-RMST posterior. This is fully reproducible and dependency-
free; the secure-environment upgrade is full NUTS (PyMC) with σ_country as a sampled hyperprior
(here σ_country is fixed and varied as a reported sensitivity, since with four countries it is
prior-dominated regardless). The flat-prior posterior, the prior→posterior overlap, and the
prior's share of the posterior precision are all reported so prior-driven and data-driven
inference stay separable.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config
from .excess_los import DEFAULT_TAU, bin_severity, build_exit_frame

# Day-bin edges -> widths (1-2, 3-7, 8-14, 15-28). Sum of widths = tau = 28.
BIN_EDGES = (0, 2, 7, 14, 28)

# Gram-negative genera (any present -> Gram-negative stratum); else "other" (Gram-positive).
_GRAM_NEGATIVE = (
    "escherichia", "klebsiella", "pseudomonas", "acinetobacter", "enterobacter",
    "proteus", "providencia", "morganella", "serratia", "citrobacter", "pantoea",
    "raoultella", "salmonella", "shigella", "stenotrophomonas", "haemophilus",
    "neisseria", "burkholderia",
)

# Exit-state codes mirror excess_los.build_exit_frame.
_EXIT_DISCHARGE = 1
_EXIT_DEATH = 2


def _gram_negative(organism) -> int:
    """1 if any Gram-negative genus appears in the (possibly polymicrobial) organism string."""
    s = str(organism).lower()
    return int(any(g in s for g in _GRAM_NEGATIVE))


def expand_person_intervals(
    df: pd.DataFrame,
    tau: int = DEFAULT_TAU,
    bin_edges: tuple[int, ...] = BIN_EDGES,
) -> pd.DataFrame:
    """Expand the ascertained cohort to a piecewise-exponential person-interval frame.

    Each ascertained patient (``resistant`` in {0,1}) contributes one row per day-bin they
    are at risk in, with the time-at-risk (``exposure``) in that bin and indicators for a
    discharge or death event occurring in it. Carries the covariates the hazard models use
    (``resistant``, ``gram_neg``, ``sev_high``, ``country``, ``bin``). Pure/deterministic.
    """
    asc = df[pd.to_numeric(df["resistant"], errors="coerce").isin([0, 1])].reset_index(drop=True)
    ef = build_exit_frame(asc, tau)
    resistant = pd.to_numeric(asc["resistant"], errors="coerce").to_numpy()
    sev_high = (bin_severity(asc["severity"]) == "high").astype(int)
    gram = asc["organism"].map(_gram_negative).to_numpy()
    country = asc["country"].to_numpy()

    edges = list(bin_edges)
    rows = []
    for i in range(len(ef)):
        t = float(ef["time"].iloc[i])
        et = int(ef["exit_type"].iloc[i])
        for k in range(len(edges) - 1):
            lo, hi = edges[k], edges[k + 1]
            if t <= lo:
                break                      # already left the admitted state before this bin
            exposure = min(t, hi) - lo
            in_bin_exit = lo < t <= hi
            rows.append({
                "bin": k,
                "exposure": exposure,
                "discharge": int(in_bin_exit and et == _EXIT_DISCHARGE),
                "death": int(in_bin_exit and et == _EXIT_DEATH),
                "resistant": int(resistant[i]),
                "gram_neg": int(gram[i]),
                "sev_high": int(sev_high[i]),
                "country": country[i],
                "patient": i,
            })
    return pd.DataFrame(rows)


def _design(frame: pd.DataFrame, countries: list, n_bins: int):
    """Build the design matrix and column index map for the hazard GLMs."""
    n = len(frame)
    bin_cols = [f"bin{k}" for k in range(n_bins)]
    cov_cols = ["resistant", "gram_neg", "sev_high"]
    country_cols = [f"country::{c}" for c in countries]
    cols = bin_cols + cov_cols + country_cols
    x = np.zeros((n, len(cols)), dtype=float)
    bins = frame["bin"].to_numpy()
    for k in range(n_bins):
        x[:, k] = (bins == k).astype(float)
    base = n_bins
    x[:, base] = frame["resistant"].to_numpy()
    x[:, base + 1] = frame["gram_neg"].to_numpy()
    x[:, base + 2] = frame["sev_high"].to_numpy()
    cbase = base + 3
    fc = frame["country"].to_numpy()
    for j, c in enumerate(countries):
        x[:, cbase + j] = (fc == c).astype(float)
    return x, cols


def _fit_penalized_poisson(x, y, offset, prior_mean, prior_prec, max_iter=200, tol=1e-9):
    """Penalized-Poisson MAP (Newton) with Gaussian prior -> (beta, covariance).

    Maximizes Σ[y·η − e^η] − ½(β−μ0)'Λ(β−μ0) with η = Xβ + offset and Λ = diag(prior_prec).
    Returns the MAP coefficients and the Laplace covariance (X'diag(μ)X + Λ)^{-1}.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    offset = np.asarray(offset, float)
    beta = np.array(prior_mean, dtype=float).copy()
    for _ in range(max_iter):
        eta = np.clip(x @ beta + offset, -30, 30)
        mu = np.exp(eta)
        grad = x.T @ (y - mu) - prior_prec * (beta - prior_mean)
        info = (x.T * mu) @ x + np.diag(prior_prec)
        step = np.linalg.solve(info, grad)
        beta += step
        if np.max(np.abs(step)) < tol:
            break
    eta = np.clip(x @ beta + offset, -30, 30)
    info = (x.T * np.exp(eta)) @ x + np.diag(prior_prec)
    return beta, np.linalg.inv(info)


def _rmst_admitted(lam, widths):
    """Restricted mean time in the admitted state for piecewise-constant exit hazards.

    ``lam`` is the per-bin all-cause exit hazard (discharge + death), ``widths`` the bin
    widths; returns ∫S(t)dt over the horizon with S piecewise-exponential. Vectorized over
    a leading axis (patients/draws): ``lam`` shape (..., n_bins).
    """
    lam = np.asarray(lam, float)
    w = np.asarray(widths, float)
    surv = np.ones(lam.shape[:-1])
    rmst = np.zeros(lam.shape[:-1])
    for k in range(lam.shape[-1]):
        lk = lam[..., k]
        safe = lk > 1e-12
        contrib = np.where(safe, surv * (1 - np.exp(-lk * w[k])) / np.where(safe, lk, 1.0),
                           surv * w[k])
        rmst = rmst + contrib
        surv = surv * np.exp(-lk * w[k])
    return rmst


def _prior_vectors(cols, prior_mu, prior_sd, sigma_country, cov_sd=1.0, bin_sd=10.0,
                   resistant_for_discharge=True):
    """Prior mean + precision per coefficient (null-centered resistant on the discharge model)."""
    mean = np.zeros(len(cols))
    sd = np.empty(len(cols))
    for j, c in enumerate(cols):
        if c.startswith("bin"):
            sd[j] = bin_sd
        elif c == "resistant":
            sd[j] = prior_sd if resistant_for_discharge else 1.0
            mean[j] = prior_mu if resistant_for_discharge else 0.0
        elif c.startswith("country::"):
            sd[j] = sigma_country
        else:
            sd[j] = cov_sd
    return mean, 1.0 / sd**2


def _gaussian_overlap(m0, s0, m1, s1) -> float:
    """Overlapping coefficient of two normals (1 = identical, 0 = disjoint)."""
    # closed form via the crossing point(s); robust numeric version by quadrature
    lo = min(m0 - 6 * s0, m1 - 6 * s1)
    hi = max(m0 + 6 * s0, m1 + 6 * s1)
    xs = np.linspace(lo, hi, 4000)
    p0 = np.exp(-0.5 * ((xs - m0) / s0) ** 2) / (s0 * np.sqrt(2 * np.pi))
    p1 = np.exp(-0.5 * ((xs - m1) / s1) ** 2) / (s1 * np.sqrt(2 * np.pi))
    return float(np.trapezoid(np.minimum(p0, p1), xs))


def bayesian_excess_los(
    df: pd.DataFrame,
    prior_mu: float = 0.0,
    prior_sd: float = 0.5,
    sigma_country: float = 0.25,
    n_draws: int = 4000,
    tau: int = DEFAULT_TAU,
    seed: int | None = None,
) -> dict:
    """Partial-pooled Bayesian excess bed-days via Laplace + g-computation.

    Fits cause-specific discharge and death hazard GLMs (penalized Poisson MAP + Laplace
    covariance), draws coefficients, and for each draw g-computes the population Δ-RMST
    (set resistant=1 vs 0 over the cohort's covariates) through the competing-risks RMST
    map. Returns the posterior mean/HDI of the excess bed-days, P(excess>0) and P(excess>1d),
    and the discharge log-HR posterior. ``prior_mu``/``prior_sd`` set the (null-centered)
    prior on the discharge resistance log-HR.
    """
    if seed is None:
        seed = config.step_seed(2)
    rng = np.random.default_rng(seed)

    long = expand_person_intervals(df, tau=tau)
    countries = sorted(pd.unique(long["country"]))
    n_bins = len(BIN_EDGES) - 1
    widths = np.diff(BIN_EDGES).astype(float)
    x, cols = _design(long, countries, n_bins)
    offset = np.log(long["exposure"].to_numpy())

    mean_dis, prec_dis = _prior_vectors(cols, prior_mu, prior_sd, sigma_country,
                                        resistant_for_discharge=True)
    mean_dth, prec_dth = _prior_vectors(cols, 0.0, 1.0, sigma_country,
                                        resistant_for_discharge=False)
    beta_dis, cov_dis = _fit_penalized_poisson(x, long["discharge"].to_numpy(), offset,
                                               mean_dis, prec_dis)
    beta_dth, cov_dth = _fit_penalized_poisson(x, long["death"].to_numpy(), offset,
                                               mean_dth, prec_dth)

    # Per-patient covariate rows (bins = 0, resistant = 0) for g-computation.
    pat = long.drop_duplicates("patient").sort_values("patient")
    zc, _ = _design(pat.assign(bin=-1), countries, n_bins)  # bin=-1 -> all bin dummies 0
    res_idx = cols.index("resistant")
    bin_idx = [cols.index(f"bin{k}") for k in range(n_bins)]
    zc[:, res_idx] = 0.0  # resistant handled explicitly below

    draws = np.empty(n_draws)
    log_hr = np.empty(n_draws)
    bdis = rng.multivariate_normal(beta_dis, cov_dis, size=n_draws)
    bdth = rng.multivariate_normal(beta_dth, cov_dth, size=n_draws)
    for d in range(n_draws):
        log_hr[d] = bdis[d, res_idx]
        base_dis = zc @ bdis[d]
        base_dth = zc @ bdth[d]
        bin_dis = bdis[d, bin_idx]
        bin_dth = bdth[d, bin_idx]
        delta_patient = np.empty(len(pat))
        for r in (0, 1):
            eta_dis = base_dis[:, None] + bin_dis[None, :] + bdis[d, res_idx] * r
            eta_dth = base_dth[:, None] + bin_dth[None, :] + bdth[d, res_idx] * r
            lam = np.exp(np.clip(eta_dis, -30, 30)) + np.exp(np.clip(eta_dth, -30, 30))
            rmst = _rmst_admitted(lam, widths)
            if r == 0:
                rmst0 = rmst
            else:
                rmst1 = rmst
        delta_patient = rmst1 - rmst0
        draws[d] = delta_patient.mean()

    lo, hi = np.quantile(draws, [0.025, 0.975])
    return {
        "tau": tau,
        "prior_mu": prior_mu,
        "prior_sd": prior_sd,
        "sigma_country": sigma_country,
        "n_draws": n_draws,
        "excess_los_mean": float(draws.mean()),
        "excess_los_hdi": [float(lo), float(hi)],
        "p_excess_gt_0": float(np.mean(draws > 0)),
        "p_excess_gt_1": float(np.mean(draws > 1)),
        "discharge_log_hr_mean": float(log_hr.mean()),
        "discharge_log_hr_sd": float(log_hr.std(ddof=1)),
        "_draws": draws,
    }


def prior_sensitivity(
    df: pd.DataFrame,
    n_draws: int = 4000,
    tau: int = DEFAULT_TAU,
    seed: int | None = None,
) -> dict:
    """Run the null-centered (primary), flat, and informative-prior posteriors side by side.

    Reports each posterior's excess-bed-day summary plus, for the discharge resistance
    log-HR, the prior→posterior overlap and the prior's share of the posterior precision
    (the "effective prior weight") — so prior-driven and data-driven inference are
    separable, as the plan requires. The informative prior (Normal(0.15, 0.30) on the
    discharge log-HR) encodes the MBIRA "resistance speeds discharge" evidence.
    """
    specs = {
        "null_centered": {"prior_mu": 0.0, "prior_sd": 0.5},     # PRIMARY
        "flat": {"prior_mu": 0.0, "prior_sd": 5.0},
        "informative_mbira": {"prior_mu": 0.15, "prior_sd": 0.30},
    }
    out = {}
    for name, spec in specs.items():
        res = bayesian_excess_los(df, n_draws=n_draws, tau=tau, seed=seed, **spec)
        post_m, post_s = res["discharge_log_hr_mean"], res["discharge_log_hr_sd"]
        prior_prec = 1.0 / spec["prior_sd"] ** 2
        post_prec = 1.0 / post_s**2
        res["prior_posterior_overlap"] = _gaussian_overlap(
            spec["prior_mu"], spec["prior_sd"], post_m, post_s)
        res["prior_precision_fraction"] = float(min(1.0, prior_prec / post_prec))
        res.pop("_draws", None)
        out[name] = res
    return out
