# OSF Pre-Registration Addendum — AMR Sentinel (Vivli 2026)

**Parent pre-registration:** https://doi.org/10.17605/OSF.IO/BFQDP (deposited 2026-05-29, before data access).
**Addendum date:** 2026-06-06 (committed to the repo 2026-06-07). **Authored after data access.** This addendum is filed under
the parent pre-registration's deviation provision (§11) and is mirrored in
`docs/deviation_log.md`. We state plainly that the changes below were made after
inspecting the shape (not the resistance→outcome associations) of the SPIDAAR data,
and we therefore label the affected analyses as deviations, not original pre-registration.

## 1. Data scope change: SMART excluded → ATLAS-only external surveillance
Merck SMART is excluded from the Vivli 2026 Challenge and cannot be used. All
external isolate-surveillance analyses (Step 2 projection/nowcast) use Pfizer ATLAS
only. Catchment coverage in ATLAS is sparse (~1,519 isolates across the four
catchment countries; Ghana/Malawi/Uganda 2021–2023 only, Kenya 2013/2014/2021–2023,
none after 2023); the wider credible intervals this forces are an irreducible,
honestly-reported consequence of excluding SMART, not a modeling choice.

## 2. Primary estimand change: 30-day mortality → resistance-attributable excess bed-days
The parent pre-registration's primary estimand is 30-day all-cause mortality
(Cox proportional-hazards, predicted HR 1.5–3.5). We change the PRIMARY estimand to
the resistance-attributable difference in restricted mean length-of-stay (excess
bed-days) to a 28-day horizon, estimated from a 3-state competing-risks multistate
model (Admitted → Discharged-alive | In-hospital death; death competes with
discharge for ending bed-day accrual):
  Δ-LOS = ∫₀²⁸ [P_admitted^Resistant(t) − P_admitted^Susceptible(t)] dt.
Rationale: the mortality contrast is underpowered at the delivered sample size and is
null in the larger, fully-ascertained precedent literature; excess bed-days is the
estimand the competing-risks methodology and the Impact case actually support.
30-day mortality is retained as a SECONDARY analysis, reported with wide intervals and
never headlined.

### 2.1 Exposure definition (pinned)
Resistant = `amrp == 2`; Susceptible = `amrp == 0` ONLY (confirmed all-susceptible).
`amrp == 1` (mixed/untested, no resistant isolate observed) and `amrp == -1`
(unascertained) are excluded from the primary contrast. The realized (resistant,
susceptible) counts under this rule are (135, 21), reported before any modeling.

### 2.2 Primary estimator and time origin
Headline: nonparametric Aalen-Johansen state-occupation curves by arm, time origin =
enrolment, integrated over [0, 28] days; 95% CIs by stratified BCa bootstrap (B=2000),
with an explicit small-sample coverage caveat. A single severity-standardized
companion estimate (cause-specific Cox on severity only, g-formula standardization) is
reported alongside. Overlap-weighted / IPTW estimates are bounded sensitivities only;
where shown, the estimand is the overlap population and is labelled as such.
Fine-Gray is reported for literature comparability only (via R `cmprsk`).

## 3. New analyses (declared as co-primary or sensitivity, all post-data)
- **Co-primary (a):** a pre-data-style power/precision simulation under the realized
  (resistant ≈135 / susceptible ≈21) design across a grid of true effects and
  ascertainment mechanisms, reporting the credible-interval width the cohort can deliver.
- **Co-primary (b):** an exposure-ascertainment selection sensitivity for the 178
  unascertained patients — an ascertainment-propensity reweight P(ascertained | severity,
  ward, length-of-stay/observation-window, country) plus MNAR/tipping-point bounds on Δ-LOS.
- **Secondary:** a Bayesian partial-pooling competing-risks layer (country-level pooling;
  pathogen as a fixed Gram-negative-vs-other stratifier). The PRIMARY prior on the
  discharge log-hazard ratio is NULL-CENTERED, Normal(0, 0.5); informative LMIC/SSA
  priors (Fiji, MBIRA, UK) are reported as sensitivity only. We report the flat-prior
  posterior and the effective prior sample size so prior-driven and data-driven inference
  are separable. This prior was NOT pre-registered before data access.
- **Secondary:** ATLAS catchment resistance NOWCAST (last-observed prevalence level,
  partial-pooled) and a SPIDAAR frame-contrast (severe-HAI inpatient vs ATLAS mixed
  surveillance), in place of a 2025–2030 per-country projection.
- **Descriptive:** a Burden-Weighted R&D Mismatch index (Axis A only: global GRAM
  burden share vs global Hub public+philanthropic funding share, DALYs primary,
  log2 ratio + Spearman ρ, n=5–6 pathogens, no inferential claim). This renders the
  parent pre-registration's H3 alignment estimand per-pathogen. A low-denominator
  floor of [VALUE — to be locked by the team] is pre-specified here.
- **Exploratory:** a g-formula empiric-adequacy stewardship what-if tool (Streamlit),
  gated on verifying a patient-side empiric-adequacy field exists; otherwise it ships as
  an ecological/transport calculator, not an estimated individual effect.

## 4. Decision rule on the pending re-linkage
If a `pid`↔`iid` re-linked SPIDAAR extract is obtained before analysis lock, (i) the
within-species resistant-vs-susceptible contrast becomes the primary exposure
definition, (ii) the isolate-ascertained cohort (3GC-R ascertained 204/244) becomes the
primary cohort for the Bayesian layer with all priors reverting to sensitivity, and
(iii) the stewardship adequacy node is refined to an individual-level `amrtx`-based node.

## 5. Honesty statement
At the delivered sample size (≈21 confirmed-susceptible patients) the headline interval
may not exclude zero. We pre-commit to reporting the interval honestly and to framing the
deliverable as the most precise resistance-attributable bed-day estimate a primary SSA-HAI
cohort can currently support, together with an explicit quantification of how much the
pending re-linkage would tighten it. We will not present a prior-driven or selection-biased
point estimate as a detected effect.
