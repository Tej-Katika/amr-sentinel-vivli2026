"""Step 2 (re-based) — ATLAS catchment resistance NOWCAST + frame-contrast (Component 3).

The pre-registered 2025-2030 per-country projection is **not data-identified**: the
delivered ATLAS catchment is ~1,519 isolates with no observations after 2023 (Ghana/
Malawi/Uganda 2021-2023 only; Kenya 2013/2014 + 2021-2023), so per-country time slopes
to 2030 would be pure prior. We therefore re-base Step 2 to a **partial-pooled current
resistance-prevalence LEVEL** (a nowcast) for the panel that exists — Enterobacterales
(E. coli + K. pneumoniae) × **ceftazidime** (ceftriaxone interpretation is blank in the
catchment) — with the SPIDAAR severe-HAI cohort used as a deliberate **frame-contrast**
rather than out-of-sample validation. Deviation logged 2026-06-07; see
``docs/analysis_plan_2026.md`` Component 3.

Partial pooling is an **empirical-Bayes beta-binomial**: a Beta(a0, b0) prior whose mean
and concentration are matched to the catchment's grand resistance rate and the
between-country variance (method of moments — no digamma/scipy needed), giving each
country a Beta posterior that shrinks small cells toward the pooled level. Credible
intervals are Monte-Carlo from the posterior (numpy Beta draws; the env lacks scipy/pymc).
With only four countries the prior is influential — per-country results are reported as
"prior-regularized," and a no-pooling (Jeffreys) companion is returned for honesty.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config, data_loading

_JEFFREYS = 0.5  # Beta(0.5, 0.5) reference prior for the no-pooling companion


def cell_counts(df: pd.DataFrame, group=("country",)) -> pd.DataFrame:
    """Collapse the isolate frame to (n_tested, n_resistant, prop) per ``group`` cell.

    Restricts to ascertained isolates (``resistant`` in {0, 1}). Pure/deterministic.
    """
    asc = df[pd.to_numeric(df["resistant"], errors="coerce").isin([0, 1])].copy()
    asc["resistant"] = pd.to_numeric(asc["resistant"], errors="coerce").astype(int)
    g = asc.groupby(list(group), dropna=True)["resistant"].agg(n_tested="count", n_resistant="sum")
    g = g.reset_index()
    g["prop"] = g["n_resistant"] / g["n_tested"]
    return g


def data_availability_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Country × year table of ascertained isolate counts — lead with this.

    Pre-empts the sparsity critique by showing, up front, exactly how thin each cell is
    (many catchment cells are 12-93 isolates; nothing after 2023).
    """
    asc = df[pd.to_numeric(df["resistant"], errors="coerce").isin([0, 1])]
    return (asc.pivot_table(index="country", columns="year", values="resistant",
                            aggfunc="count", fill_value=0)
            .astype(int))


def _eb_beta_prior(r: np.ndarray, n: np.ndarray) -> dict:
    """Empirical-Bayes Beta(a0, b0) prior by moment-matching across cells.

    Mean = pooled grand rate m; concentration κ = m(1−m)/τ² − 1, where τ² is the
    between-cell variance with the within-cell binomial component removed. τ² is floored
    at a small positive value (strong-pooling guard) and κ at 1 pseudo-observation.
    Returns ``a0, b0, m, kappa, tau2``.
    """
    r = np.asarray(r, dtype=float)
    n = np.asarray(n, dtype=float)
    p = r / n
    m = r.sum() / n.sum()
    within = np.mean(p * (1.0 - p) / n)          # mean binomial sampling variance
    total = np.mean((p - m) ** 2)                # observed between-cell variance
    tau2 = max(total - within, 1e-6)             # true between-cell variance (floored)
    kappa = max(m * (1.0 - m) / tau2 - 1.0, 1.0)
    return {"a0": m * kappa, "b0": (1.0 - m) * kappa, "m": float(m), "kappa": float(kappa),
            "tau2": float(tau2)}


def _beta_ci(a: float, b: float, rng, n_mc: int, alpha: float = 0.05) -> tuple[float, float]:
    draws = rng.beta(a, b, size=n_mc)
    lo, hi = np.quantile(draws, [alpha / 2, 1 - alpha / 2])
    return float(lo), float(hi)


def nowcast(
    df: pd.DataFrame,
    since_year: int | None = None,
    group: str = "country",
    n_mc: int = 20000,
    seed: int | None = None,
) -> dict:
    """Partial-pooled resistance-prevalence nowcast per ``group`` cell.

    ``since_year`` restricts to the recent window that defines "current" level (e.g. 2021
    to drop Kenya's 2013/2014 breakpoint-era isolates; see ``regional_trend`` for the
    with/without-era comparison). Returns per-cell posterior mean + 95% credible interval
    (empirical-Bayes Beta posterior, MC interval) alongside the raw proportion, plus the
    fitted prior and a no-pooling (Jeffreys) companion mean for each cell.
    """
    if seed is None:
        seed = config.step_seed(2)
    rng = np.random.default_rng(seed)

    use = df if since_year is None else df[pd.to_numeric(df["year"], errors="coerce") >= since_year]
    counts = cell_counts(use, group=(group,))
    prior = _eb_beta_prior(counts["n_resistant"].to_numpy(), counts["n_tested"].to_numpy())

    cells = []
    for _, row in counts.iterrows():
        r, n = int(row["n_resistant"]), int(row["n_tested"])
        a, b = prior["a0"] + r, prior["b0"] + (n - r)
        lo, hi = _beta_ci(a, b, rng, n_mc)
        cells.append({
            group: row[group],
            "n_tested": n,
            "n_resistant": r,
            "raw_prop": r / n,
            "posterior_mean": a / (a + b),
            "ci_lower": lo,
            "ci_upper": hi,
            "no_pooling_mean": (r + _JEFFREYS) / (n + 2 * _JEFFREYS),
        })
    return {
        "since_year": since_year,
        "group": group,
        "prior": prior,
        "pooled_rate": prior["m"],
        "cells": cells,
    }


def _wls_slope(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    """Weighted simple linear regression -> (slope, intercept)."""
    w = w / w.sum()
    xbar = np.sum(w * x)
    ybar = np.sum(w * y)
    var_x = np.sum(w * (x - xbar) ** 2)
    cov = np.sum(w * (x - xbar) * (y - ybar))
    slope = cov / var_x if var_x > 0 else 0.0
    return float(slope), float(ybar - slope * xbar)


def regional_trend(df: pd.DataFrame, since_year: int | None = None) -> dict:
    """Single pooled (regional) resistance time-slope on the logit scale.

    Pools all catchment countries by year (per-country slopes are not identified), takes
    the empirical-logit of the yearly resistant fraction, and fits a weighted (by n_tested)
    linear trend. Returns the slope per year (logit scale) with and without the option to
    drop pre-2020 (Kenya breakpoint-era) isolates. This is a BORROWED regional trend, not a
    per-country forecast — used only for the optional projection scenario.
    """
    use = df if since_year is None else df[pd.to_numeric(df["year"], errors="coerce") >= since_year]
    by_year = cell_counts(use, group=("year",)).sort_values("year")
    yr = by_year["year"].to_numpy(dtype=float)
    r = by_year["n_resistant"].to_numpy(dtype=float)
    n = by_year["n_tested"].to_numpy(dtype=float)
    logit = np.log((r + 0.5) / (n - r + 0.5))    # Haldane-corrected empirical logit
    slope, intercept = _wls_slope(yr, logit, n)
    return {
        "since_year": since_year,
        "slope_logit_per_year": slope,
        "intercept_logit": intercept,
        "years": yr.astype(int).tolist(),
        "n_by_year": n.astype(int).tolist(),
    }


def project_forward(
    nowcast_result: dict,
    horizon_years: int = 2,
    scenario: str = "flat",
    trend: dict | None = None,
) -> dict:
    """Project each cell's nowcast level forward by ``horizon_years``.

    ``scenario="flat"`` (default, defensible) carries the current level unchanged.
    ``scenario="trend"`` applies the BORROWED regional logit slope (``trend`` from
    ``regional_trend``) equally to every cell — explicitly labelled as borrowed, never a
    per-country data forecast. Never projects to a specific 2030 point.
    """
    if scenario not in ("flat", "trend"):
        raise ValueError("scenario must be 'flat' or 'trend'")
    if scenario == "trend" and trend is None:
        raise ValueError("scenario='trend' requires a `trend` from regional_trend()")

    out = []
    for c in nowcast_result["cells"]:
        p = c["posterior_mean"]
        if scenario == "flat":
            projected = p
        else:
            shifted = np.log(p / (1 - p)) + trend["slope_logit_per_year"] * horizon_years
            projected = float(1.0 / (1.0 + np.exp(-shifted)))
        out.append({nowcast_result["group"]: c[nowcast_result["group"]],
                    "current_level": p, "projected_level": projected})
    return {"scenario": scenario, "horizon_years": horizon_years, "cells": out}


def resistant_fraction_multiplier(nowcast_result: dict) -> dict:
    """Per-country resistant-fraction multiplier (mean + CrI) for the Component-1 bridge.

    Exposes the nowcast posterior as the resistant-fraction input that scales the
    excess-bed-day burden to admissions; the credible interval is propagated into the
    bed-day Monte-Carlo.
    """
    group = nowcast_result["group"]
    return {
        c[group]: {"fraction": c["posterior_mean"], "ci_lower": c["ci_lower"],
                   "ci_upper": c["ci_upper"], "n_tested": c["n_tested"]}
        for c in nowcast_result["cells"]
    }


_SPIDAAR_ENTERO = ("Escherichia coli", "Klebsiella pneumoniae")


def _wilson_ci(r: int, n: int, z: float = 1.959963985) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion (no scipy)."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = r / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def _prev_ci(r: int, n: int) -> dict:
    lo, hi = _wilson_ci(r, n)
    return {"prevalence": (r / n) if n else float("nan"), "ci_lower": lo, "ci_upper": hi,
            "n_tested": int(n), "n_resistant": int(r)}


def frame_contrast(
    atlas_df: pd.DataFrame,
    spidaar_isolates: pd.DataFrame,
    enterobacterales_only: bool = True,
) -> dict:
    """SPIDAAR severe-HAI 3GC-R vs ATLAS mixed-surveillance ceftazidime-R (the frame-shift).

    A deliberate **frame-contrast**, NOT temporal validation: SPIDAAR is a severity-enriched
    severe-HAI *inpatient* cohort while ATLAS is mixed surveillance, so a higher SPIDAAR
    resistance level quantifies the HAI/severity frame-shift in the catchment — a genuine
    triangulation finding. Both are 3rd-generation-cephalosporin indicators but NOT the same
    drug (ATLAS = ceftazidime interpretation; SPIDAAR ``c3r`` = a 3GC-R composite), and the
    sampling frames differ — stated as the ``caveat``. Returns overall + per-country
    prevalence with Wilson intervals for each source, restricted to the comparable
    Enterobacterales panel by default.
    """
    a = atlas_df[pd.to_numeric(atlas_df["resistant"], errors="coerce").isin([0, 1])].copy()
    s = spidaar_isolates.copy()
    s["c3r"] = pd.to_numeric(s["c3r"], errors="coerce")
    s = s[s["c3r"].isin([0, 1])]
    if enterobacterales_only:
        a = a[a["species"].isin(_SPIDAAR_ENTERO)]
        s = s[s["organism"].astype(str).isin(_SPIDAAR_ENTERO)]
    a["resistant"] = pd.to_numeric(a["resistant"], errors="coerce").astype(int)
    s["c3r"] = s["c3r"].astype(int)

    def _by_country(frame, col):
        return {str(c): _prev_ci(int(g[col].sum()), int(g[col].count()))
                for c, g in frame.groupby("country")}

    return {
        "panel": ("Enterobacterales (E. coli + K. pneumoniae) × 3GC"
                  if enterobacterales_only else "all organisms × 3GC"),
        "atlas_drug": "ceftazidime (interpretation; mixed surveillance)",
        "spidaar_marker": "c3r 3GC-R composite (severe-HAI inpatient)",
        "atlas": {"overall": _prev_ci(int(a["resistant"].sum()), int(a["resistant"].count())),
                  "by_country": _by_country(a, "resistant")},
        "spidaar": {"overall": _prev_ci(int(s["c3r"].sum()), int(s["c3r"].count())),
                    "by_country": _by_country(s, "c3r")},
        "caveat": ("ATLAS = ceftazidime interpretation from mixed surveillance; SPIDAAR = "
                   "c3r 3GC-R composite from a severe-HAI inpatient cohort. Different drug "
                   "basis and sampling frame — the contrast measures the frame-shift, not "
                   "an exact like-for-like resistance difference."),
    }


def run_nowcast(df: pd.DataFrame | None = None, since_year: int | None = 2021) -> dict:
    """Step-2 entrypoint: data-availability matrix + partial-pooled catchment nowcast.

    Loads the catchment Enterobacterales × ceftazidime panel if ``df`` is not supplied,
    defaults the nowcast window to ``since_year=2021`` (drops Kenya's breakpoint-era
    isolates), and returns the availability matrix plus the nowcast. EXPLORATORY — see the
    module docstring on the four-country prior influence.
    """
    if df is None:
        df = data_loading.load_atlas()
    return {
        "data_availability": data_availability_matrix(df),
        "nowcast": nowcast(df, since_year=since_year),
        "nowcast_all_years": nowcast(df, since_year=None),
    }
