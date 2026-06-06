# Vivli 2026 — Integrated Analysis Plan: Resistance-Attributable Excess Bed-Days (SPIDAAR-Anchored)

**Status:** build-ready plan, post-critique. **Date:** 2026-06-06. **Lead methodologist sign-off pending the gated verifications in §3 and the team decisions in §6.**
**Governance anchor:** OSF pre-registration `10.17605/OSF.IO/BFQDP` (deposited 2026-05-29, before data access). Primary estimand there is 30-day mortality (Cox, predicted HR 1.5–3.5). **This plan is a post-data-access deviation and is logged as such (§5).**

---

> **Local gate resolution (2026-06-06).** Gates the design agents flagged were checked
> against the raw 60-column `patientdata.xls` (they saw only the loader's 10 columns):
> - **Gate L (outcome + confounders): PASS** — `los` (255/336), `disev`, `ward`, `enrtpt`
>   present, plus a rich confounder set (`morst, surgyn, devyn, devices, orgsup,
>   diab_/hiv_/tuber_/prisk*/*_medh`). Component 1 has its outcome.
> - **Gate A (adequacy field): PASS** — `txadp` exists {adequate 106 / inadequate 52 /
>   unknown 178}; Component 5 can be exploratory-causal, not ecological-only.
>   (0/1/9 coding to confirm vs secure-env codebook; matches the data roadmap.)
> - **Gate E (`enrtpt`): PASS = day-count** (range 2–118) → left-truncation sensitivity enabled.
> - **Realized exposure ns (pinned rule):** Resistant(amrp==2)=135, Susceptible(amrp==0)=**21**,
>   not-confirmed(amrp==1)=2, unascertained(-1)=178. Susceptible arm is **21**, not 23.
> - **Still needs secure-env codebook:** Gate T (`amrp` = index isolate / exposure fixed at
>   baseline?) and Gate X (ATLAS panel fields).

## 1. Executive summary

We pivot the headline from resistance→mortality (underpowered: ~14–55 exposure-ascertained deaths; literature null) to **resistance-attributable excess length-of-stay (excess bed-days)** estimated by a **competing-risks multistate model** on the SPIDAAR patient cohort, where in-hospital **death competes with discharge-alive** for ending bed-day accrual.

The honest design, after integrating four adversarial critiques, is **precision-honest and single-headline**, not a triangulated showcase:

- **Headline estimand:** the difference in **restricted mean time in the "admitted" state to τ = 28 days** between resistant and susceptible HAI patients — i.e. excess bed-days = ∫₀²⁸ [P_admitted^R(t) − P_admitted^S(t)] dt — from a 3-state (Admitted → Discharged-alive | Died) **Aalen-Johansen** model, **time origin = enrolment**, **estimated crude plus a single severity-standardized companion**. Overlap/IPTW adjustment is demoted to a bounded sensitivity (non-identified at n=23 susceptible).
- **Two co-primary honesty analyses** that the critiques insisted cannot be deferred: (a) a **pre-data power/precision simulation** under the real (R≈135 / S≈23) design; (b) an **exposure-ascertainment (selection) sensitivity** for the 178 unascertained patients (tipping-point / MNAR bounds + an ascertainment-propensity reweight). At this n, these determine whether the headline exists.
- **The binding limit is reframed as the contribution.** With 23 confirmed-susceptible patients, Fiji-level precision (2.6d, CI 2.5–2.8) will not transfer. We say so, quantify it, and frame the deliverable as *"the most precise resistance-attributable bed-day estimate a primary SSA-HAI cohort can currently support, and an explicit quantification of how much the pending re-linkage would tighten it."* That transparency is the Impact/Methodology story and survives a live grilling; overselling a "detectable" effect does not.

**Supporting components, all demoted from their original ambition:**
- **Bayesian partial pooling** (country only; pathogen as fixed Gram-negative-vs-other stratifier) is a **null-centered, prior-regularized robustness companion**, never the headline; informative LMIC priors appear only as a pre-registered sensitivity. Strategy doc already flags partial pooling as "crowded."
- **ATLAS** is re-based from a 2025–2030 per-country projection (not data-identified: 3 of 4 countries have only 2021–2023; none past 2023) to a **catchment resistance *nowcast* + a SPIDAAR *frame-contrast*** (severe-HAI inpatient vs ATLAS mixed surveillance), with a regional-only pooled time slope offered as an explicitly-borrowed scenario.
- **Cross-domain R&D mismatch** ships as a **single descriptive, pre-registered (H3) like-for-like index** (global public+philanthropic funding share vs global GRAM burden share, DALYs primary), with bed-day weighting as a sensitivity overlay on the 1–2 pathogens that can support it.
- **Stewardship g-formula** is recast from "confirmatory capstone" to a **pre-registered *exploratory* what-if Streamlit tool**, hard-gated on verifying that a patient-side empiric-adequacy field actually exists (it is **not** in the verified loader).

**Innovation thesis (relocated, because the methods themselves are known):** identification of excess bed-days *under a broken patient↔isolate link* in an LMIC HAI cohort, the SSA-HAI frame-contrast against global surveillance, and a re-runnable local-calibration tool — not the multistate estimator, which is established.

---

## 2. Per-component final specs (post-critique), ordered

### Component 1 — Excess-LOS competing-risks primary (`src/amr_sentinel_vivli/excess_los.py`)

**Estimand.** Δ-LOS = E[min(LOS,28) | resistant] − E[min(LOS,28) | susceptible] = ∫₀²⁸ [P_admitted^R(t) − P_admitted^S(t)] dt, where LOS = time-to-discharge-alive and in-hospital death is a competing absorbing event. Restricted mean to τ=**28** days (inside the 28–31d `nobsd` window, off the thin boundary). This is a difference in restricted expected occupancy — causally interpretable, unlike a Fine-Gray subdistribution HR.

**Exposure definition (pinned, the single largest hidden flaw fixed first).** Re-define before any modeling, and change the loader:
- **Resistant = `amrp == 2`** (≥1 resistant isolate).
- **Susceptible = `amrp == 0` ONLY** (confirmed all-susceptible).
- **`amrp == 1` (mixed/untested, no resistant isolate seen) is NOT susceptible** — it moves to a third "not-confirmed" category alongside `amrp == -1`, excluded from the primary contrast and carried only in a labelled sensitivity.
- This corrects the current loader, which maps `amrp==1 → susceptible(0)` (`_RESISTANT_FROM_AMRP = {2:1, 0:0, 1:0}`). **Report the resulting (R, S) n under the new rule up front.** Expect S to *shrink* below 23 — own this number.

**Method (single headline + one companion).**
1. **HEADLINE:** nonparametric **Aalen-Johansen** state-occupation curves, estimated separately by arm, time origin = **enrolment** (t0=0), integrated over [0,28] → crude Δ-RMST. CIs by **stratified BCa bootstrap** (resample within country, B=2000), with an explicit coverage caveat at this n.
2. **COMPANION:** a single **severity-standardized** estimate via cause-specific Cox on **severity only** (g-formula standardization over the empirical severity distribution). Not a 9-dim adjustment; not IPTW.
3. **Cause-specific hazards** (discharge and death) reported as **exploratory mechanistic decomposition only** — with ~56 deaths split by arm, the death-hazard model is uninformative; flag it.
4. **Fine-Gray** reported for literature comparability **only via rpy2→R `cmprsk`** (lifelines has no Fine-Gray fitter — do not claim it in Python), or omitted. Never the causal headline.
5. **Overlap-weighted / IPTW** Δ-RMST: bounded **sensitivity** with effective-sample-size, max-weight, and weight-distribution diagnostics; state explicitly that ATO redefines the estimand to the overlap population.

**Time handling.** Enrolment origin is **primary** (sidesteps both the non-Markov-under-left-truncation tension and immortal-time-on-exposure). Admission origin + left-truncation via `enrtpt` is a **sensitivity**, contingent on `enrtpt` being a verified day-count (§3 gate). `dead==9` (n=1) → death state (matches loader convention); one-patient recode sensitivity. The ~25 patients with neither `los` nor `dtpta` are **retained as censored-still-admitted** at min(`nobsd`,28).

**Burden / cost (estimand-matched scaling — fixes the transportability error).** Scale Δ-LOS to total catchment bed-days using a count **matched to the estimand population**: the crude/severity-standardized (ATE-type) Δ-LOS scales to the full resistant-HAI count; an ATO Δ-LOS, if shown, scales **only** to the overlap-population count, stated in the caption. Cost = bed-days × country-specific bed-day unit cost (WHO-CHOICE / national HTA), Monte-Carlo propagated (`config.MONTE_CARLO_DRAWS=10000`), reported as a sourced **range** with unit-cost sensitivity — never a single point.

**Data inputs (verified-available):** `amrp`, `dead`, `dtpta`, `nobsd`. **Gated-pending-verification:** `los`, `disev`(severity), `ward`, `enrtpt` (referenced in the deviation log as raw fields but **not surfaced in the loader**; §3 Gate L). External: bed-day costs; published excess-LOS anchors.

**Outputs:** headline Δ-RMST (crude + severity-standardized) with BCa CI; stacked state-occupation plot (`figures/excess_los_stateoccupation.png`); CIFs with day-28 contrast; supporting cause-specific HR table (flagged exploratory); total bed-days + cost table with MC CI; E-value on the adjusted discharge HR + tipping-point for Δ-LOS; sensitivity panel (time-origin, dead==9, weights); **honest power statement** (effective n, simulated CI width); `results/excess_los_primary.json`; tests on synthetic frames.

**Key assumptions:** `amrp` reflects the **index** HAI isolate (exposure fixed at baseline → no time-dependent bias) — **verify** (§3 Gate T); `los` = discharge-alive, mutually exclusive with `dtpta`; non-informative administrative censoring within [0,28]; conditional exchangeability (bounded by E-value, not assumed); positivity does **not** fully hold at n=23 (overlap estimand stated explicitly if used); MAR-on-exposure for the 178 unascertained (the dominant threat, probed in Component-2-style sensitivity).

**Libraries:** lifelines (AalenJohansenFitter, CoxPHFitter with `entry=`), numpy/scipy (bootstrap, BCa, integration, MC), pandas, statsmodels (propensity), matplotlib; optional rpy2→R `cmprsk`/`etm` as a cross-validation oracle and the only real Fine-Gray path. **None installed** in the env — add pinned to `pyproject.toml`, log a Tooling deviation (mirror xlrd), ship an import smoke-test.

**Residual risks:** crude Δ-RMST is a confounded (sicker-resistant vs less-sick-susceptible) comparison — severity-standardization is the only feasible mitigation; bootstrap under-coverage at n=23; ascertainment selection bias (Component 1b) can manufacture a positive effect; within-species contrast not identifiable without the link (logged deviation).

---

### Component 1b — Power simulation + ascertainment sensitivity (co-primary)

**Why co-primary:** at n≈23-or-fewer susceptible, these two analyses determine whether the headline is real; the critiques unanimously refused to let them be "output bullets" or "open items."

**(a) Pre-data power/precision simulation.** Simulate Δ-RMST CI width under the actual (R≈135 / S≈23) design across a grid of true effects (0 to ~3 days) and ascertainment mechanisms. Deliverable: the CI width the cohort can deliver. If the CI cannot exclude 0 for plausible effects, the component is reframed honestly as "the precision a primary SSA-HAI cohort currently supports + the value of re-linkage." Output: `results/excess_los_power_sim.json`, a CI-width figure.

**(b) Exposure-ascertainment selection sensitivity.** 178/336 unascertained (`amrp==-1`). Pre-specify an **ascertainment-propensity model** P(ascertained | severity, ward, `los`/`nobsd`, country) → ascertainment-weighted Δ-LOS, **plus** tipping-point / MNAR bounds. Because culturing is plausibly LOS-driven (longer-staying patients more cultured *and* more ascertained-resistant → self-fulfilling positive Δ-LOS), this is the binding validity threat and is reported alongside the headline. Also characterize how the 158 ascertained differ from the 178 unascertained on observed covariates, and qualify the 4-country population claim accordingly.

---

### Component 2 — Bayesian partial pooling (secondary robustness companion) (`bayesian_projection.py` is Step-2 ATLAS; this LOS-pooling layer lives in `excess_los.py` or a sibling module)

**Estimand.** Same excess-bed-day Δ-RMST, with a partially-pooled resistance effect on the **discharge** cause-specific log-hazard. Reported as P(excess > 0) and P(excess > 1 day) on the bed-day scale.

**Method (corrected).** Competing-risks **piecewise-exponential / Poisson** likelihood (patient-day expansion; daily intervals binned 1–2 / 3–7 / 8–14 / 15–28 to respect ≥48h HAI and small counts); discharge and death cause-specific hazards modeled jointly. Hierarchy:
- **Pooling on COUNTRY only** (u_c ~ Normal(0, σ_country)); **pathogen as a FIXED Gram-negative-vs-other stratifier**, not a random effect — so the unmatched *S. aureus* cell never borrows a Gram-negative-anchored prior. State explicitly that *S. aureus* excess-LOS is not estimable here.
- **Primary prior CENTERED AT THE NULL:** mu ~ Normal(0, 0.5) on the discharge log-HR. Justification: the cited discharge-hazard evidence is null-to-reversed (MBIRA discharge CSH ratio 1.16 = resistance *speeds* discharge; Fiji discharge aHR 0.99). The original spec's mu ~ Normal(−0.10, 0.30) is **more optimistic than its own literature** and post-data — removed as primary.
- Informative LMIC/SSA priors (Fiji/MBIRA/UK), translated to the log-HR scale via prior-predictive simulation through the AJ map, appear **only as pre-registered sensitivity**.
- σ_country ~ HalfNormal(0, 0.25); confounder coefficients ~ Normal(0,1); **non-centered** parameterization. Drop the regularized-horseshoe (overkill, fights the funnel).

**Mandatory honesty:** report **effective prior sample size**, the **flat-prior posterior** alongside the informative one, prior-vs-posterior overlap **for σ as well as μ** (σ is prior-dominated with 4 groups — say so; label per-country results "prior-regularized," not a data result). **Drop the "Stage-1/Stage-2 agreement = robustness" framing** — both stages use the same 158 patients and likelihood; concordance is a statement about the prior, not validation.

**Re-linkage decision rule (pre-registered):** if the `pid`↔`iid` crosswalk arrives before lock, switch to the data-dominant isolate-ascertained cohort (3GC-R ascertained 204/244) as primary and revert all priors to sensitivity.

**Libraries:** PyMC≥5, ArviZ, preliz (prior elicitation), lifelines/sksurv (frequentist cross-check), pandas/numpy, matplotlib.

**Residual risks:** flat-prior posterior likely spans roughly −3 to +5 days (cannot exclude 0); P(excess>1d) likely ~0.5–0.7, not headline-grade — **accept and lead with it**; prior-data conflict in the *wrong direction* if SPIDAAR echoes MBIRA; target-population claim limited by ascertainment selection.

**Reconcile before pitch:** the strategy doc says MBIRA found 3GC-R *was* associated with mortality; this spec cites MBIRA death CSH 0.74 (null). **Re-pull MBIRA (PMC7617135) and reconcile** — a judge who knows MBIRA catches this instantly and it undermines all four anchors.

---

### Component 3 — ATLAS-only catchment nowcast + frame-contrast (`bayesian_projection.py`)

**Estimand (re-based).** Partial-pooled **current (last-observed, ~2022–2023) resistance-prevalence LEVEL** per catchment country × cell — which the data support — *not* per-country 2025–2030 trajectories. Verified sparsity: catchment = ~1,519 isolates (0.15% of ATLAS); Ghana/Malawi/Uganda only 2021–2023, Kenya only 2013/2014/2021–2023, **nothing past 2023**. Per-country linear time slopes to 2030 are pure prior/pooling — not data-identified.

**Panel (re-specified to what exists).** **Enterobacterales × CEFTAZIDIME** (E. coli + K. pneumoniae as carrying species; **ceftriaxone interpretation is blank in catchment** — the original spec's anchor cell does not exist), R re-derived from MIC under EUCAST v15.0; secondary MDR composite; *S. aureus*×methicillin down-ranked (thin, year-uneven). Any agent lacking MICs uses ATLAS-supplied interpretation and is excluded from any trend claim.

**Model (lightened).** Single partial-pooling tier (SSA → region → country) with bug×drug fixed effects and **ONE pooled (regional) time slope** (no country-specific slopes — unidentified; Ghana has zero pre-2021 data and is a singleton in "Western Africa"). Non-centered; **beta-binomial** for lab/over-dispersion; **drop the horseshoe**. 8000 NUTS draws is fine but matched to evidence now.

**Projection.** Offered **only as a flat-carried scenario (default)** plus a clearly-labelled **regional-trend (borrowed)** scenario; pre-register a short **1–2yr nowcast** as the defensible horizon. Never a country-specific point forecast to 2030.

**Kenya 2013/2014 handling (pre-registered):** the only long series straddles a decade of breakpoint-vintage change (ceftazidime-R 0.31→0.67 is partly artefact). Use an era indicator or exclude pre-2020; show trend with/without.

**SPIDAAR role (reframed).** SPIDAAR (2021–2022, in-window, severe-HAI severity-enriched frame) is **NOT out-of-sample temporal validation** — it is a **deliberate frame-contrast** (severe-HAI inpatient ~88.7% 3GC-R vs ATLAS mixed-surveillance ceftazidime-R ~0.35–0.92). A genuine, novel triangulation finding for the catchment; the offset-term version quantifies the HAI frame-shift, not robustness.

**Bridge to Impact (made explicit):** the prevalence level (with CrI) feeds Component 1's bed-day burden as the **resistant-fraction multiplier** on admissions, CrI propagated into the bed-day Monte-Carlo. If this bridge cannot be made concrete, demote the whole module to a supporting descriptive figure.

**Outputs:** `results/projection_posterior.nc`; `results/projection_levels.csv` (per cell: median, 50/95% CrI, support_flag, in-cell n_tested); **lead with a data-availability matrix** (country × year × cell n_tested) to pre-empt the sparsity attack; `results/spidaar_framecontrast.csv` + figure; convergence report.

**Libraries:** PyMC≥5, ArviZ, numpy/pandas, pytensor, scipy/statsmodels (Wilson intervals), matplotlib.

**Residual risks:** even levels are noisy (cell n often 12–65); most catchment cells flag "pooling-dominated"; novelty is thin (Waterlow minus its incidence layer) — the **frame-contrast and the explicit sparsity-owning are the defensible contributions**, not the projection.

---

### Component 4 — Cross-domain R&D mismatch (`rd_alignment.py`)

**Estimand.** Single descriptive **Axis A** index: M_p = (global GRAM burden share) / (global Hub public+philanthropic funding share), reported as **log2 M_p**. This *is* the pre-registered H3 alignment estimand rendered per-pathogen, plus the registered Spearman ρ. Descriptive/ecological — not causal.

**Method.** Pathogen panel = intersection of GRAM-22 ∩ Hub "infectious agent" ∩ SPIDAAR organisms (~K. pneumoniae, E. coli, S. aureus/MRSA, S. pneumoniae, Acinetobacter, + "other"). **Burden numerator = DALYs (associated + attributable, both reported)** as PRIMARY. Funding denominator = Hub global public+philanthropic (therapeutics+diagnostics+vaccines primary; basic-research sensitivity). **Extract and report the cross-cutting/non-pathogen-specific slice FIRST as a headline magnitude** — gate the index on it (if pathogen-specific spend is a minority, widen all caveats); never silently redistribute. Monte-Carlo propagate GRAM UIs through b_p; Dirichlet sensitivity on cross-cutting; **pre-specify the low-denominator floor in the OSF addendum**; report ranking with/without floor. **Spearman ρ with CI, descriptive only, explicit n=5-6 caveat; never fit a line.**

**Demotions per critique:** bed-day weighting → **sensitivity overlay on the 1–2 pathogens** (likely K. pneumoniae, E. coli) where SPIDAAR has resistant AND susceptible ascertained cases, with per-cell n printed; the **SSA-burden-vs-global-funding "Axis B"** → supplement only, labelled an equity/geographic contrast (dividing SSA numerator by global denominator is a scope artefact, not a finding); the **GRAM-calibration scatter is CUT**, replaced by a single honest sentence of qualitative concordance with the event count.

**Governance:** any SPIDAAR per-pathogen excess-LOS used here is computed **only after** the OSF addendum + deviation-log entry (§5). Title everything "public+philanthropic"; fixed caveat that private R&D may offset apparent therapeutic under-funding.

**Libraries:** pandas/numpy, scipy.stats (Dirichlet, Spearman, bootstrap), matplotlib/plotly, requests/playwright (Hub has no public API — **frozen dated snapshot mandatory**, exact filter state recorded).

**Residual risks:** n=5-6 pathogens → no powered inference (descriptive only); divide-by-small instability (floor flagged); taxonomy crosswalk losses; live-dashboard drift (snapshot).

---

### Component 5 — Stewardship g-formula what-if (exploratory) (`stewardship_gformula.py`, `app/streamlit_app.py`)

**Status: DOWNGRADED to pre-registered *exploratory* causal what-if, HARD-GATED.** The original "confirmatory Tier-A" rests on a patient-side empiric-adequacy field `txadp`. **`txadp` is NOT in the verified loader** (patient columns are `{pid, ctry, agegr, sex, chaicat, isol, amrp, dead, dtpta, nobsd}`; `SPIDAAR_COLUMNS` has no adequacy field). The only verified adequacy signal is `amrtx` (resistance to *administered* class) on the **isolate** file — **unlinkable** to patients. Until Gate A (§3) confirms a patient-side empiric-window adequacy field exists, **there is no patient-level treatment node and Tier-A cannot be called confirmatory.**

**Estimand (if Gate A passes).** Population counterfactual Δ in restricted-mean time-in-hospital (excess bed-days) to τ=28 under "set empiric adequacy A=1 for all" vs natural course, **defined on the positivity-supported population**. Mortality secondary, wide CI, never headlined.

**Method.** Parametric/Bayesian **g-formula** over competing-risks cause-specific hazards (discrete-day microsimulation, seeded `config.step_seed(5)`), partial-pooled across countries. **Parsimonious DAG-justified L = severity + ward + organism-group + country (random effect)** — NOT a 9-dim set (52 inadequate patients cannot populate it; partial pooling does not manufacture positivity). **Resolve the resistance/adequacy collider explicitly in a DAG:** resistance largely *determines* adequacy, so it is a mediator/determinant, not a generic confounder — either drop `amrp` from L or reframe as a controlled-direct-effect within resistance strata. **Success criterion = sign-consistency across g-formula / IPTW / Fine-Gray cross-check**, not "CrI excludes 0." Report a formal positivity diagnostic (adequacy-score overlap, % off-support). Treat 53%-missing adequacy as plausibly **NMAR** with a tipping-point bound, not a MICE footnote.

**Streamlit tool.** Loads de-identified posterior artifacts + aggregated cell tables only (Vivli k-anonymity / egress review; gate cell display behind a minimum-n). User knobs: local breakpoints/antibiogram → expected adequacy; ward-mix; bed-day cost. **Tier-B `amrtx` cells (ecological) drive the scenario mapping — label every tool output "ecological-calibration-driven scenario, not an estimated individual effect"** (the firewall must be real, not rhetorical). Ships synthetic demo data; Apache-2.0.

**Innovation framing:** the contribution is **identification-under-a-broken-link + a re-runnable LMIC local-calibration tool**, not the g-formula.

**Libraries:** lifelines, numpy/pandas, PyMC+ArviZ, scikit-learn (propensity/missingness), streamlit+plotly, statsmodels (MICE), pytest+ruff.

**Residual risks:** effective n≤158 with 52 in the smaller arm → wide/prior-driven CrIs; NMAR adequacy unfixable by MICE; egress may block real-data artifacts; **whole component gated on Gate A and the re-linkage decision (§6).**

---

## 3. Dependency-ordered build sequence

**Hard gates (verify in the secure env BEFORE the dependent work; not "open items"):**
- **Gate L (loader fields):** confirm `los`, `disev`(severity), `ward`, `enrtpt` exist in the raw `patientdata.xls`. The deviation log line 13 lists `los`/`enrtpt` as day-count fields, but **none are surfaced in the loader and severity/ward are nowhere in code.** If absent, Components 1/2/5 lose their outcome (`los`) and/or confounders — Component 1 then cannot run at all.
- **Gate T (exposure timing):** confirm `amrp` reflects the **index** HAI isolate (exposure fixed at baseline). If later-acquired → time-varying modeling needed.
- **Gate E (`enrtpt` semantics):** day-count (enables left-truncation sensitivity) vs categorical phase (covariate only).
- **Gate A (adequacy field):** confirm whether a patient-side empiric-window adequacy field exists. Drives whether Component 5 is exploratory-causal or ecological-only.
- **Gate X (ATLAS fields):** confirm species / antibiotic / S-I-R vs MIC presence per panel cell.

**Build order:**

| Order | Task | Blocked on re-linkage? | Blocked on a gate? |
|---|---|---|---|
| 0 | **Tooling:** add lifelines, scikit-survival, statsmodels, scipy, PyMC, ArviZ, preliz (+ optional rpy2) to `pyproject.toml` pinned; reconcile `config.PYTHON_VERSION` ("3.12") with the real 3.13 interpreter; import smoke-test; log Tooling deviation. | No | No |
| 1 | **Pin exposure definition** (S = `amrp==0` only); update `_RESISTANT_FROM_AMRP` and `SPIDAAR_COLUMNS`; report (R,S) ns; log deviation. | No | No |
| 2 | **Extend loader** to surface `los`, `disev`, `ward`, `enrtpt`; synthetic-frame tests. | No | **Gate L** |
| 3 | **Power/precision simulation** (Component 1b-a) — can run on synthetic + the (R,S) ns alone. | No | No |
| 4 | **`excess_los.py` core** — crude AJ Δ-RMST (enrolment origin) + severity-standardized companion + bootstrap; tests. | No | **L, T** |
| 5 | **Ascertainment sensitivity** (Component 1b-b). | No | L |
| 6 | **Cost/burden scaling** (estimand-matched) + E-value/tipping-point. | No (costs external) | depends on 4 |
| 7 | **Bayesian LOS pooling** (Component 2), null-centered prior. | **Decision rule fires if link arrives** | depends on 4 |
| 8 | **ATLAS nowcast + frame-contrast** (Component 3); lead with data-availability matrix. | No | **Gate X** |
| 9 | **R&D mismatch index** (Component 4, Axis A); cross-cutting slice first. | bed-day overlay improves with link | OSF addendum first |
| 10 | **Stewardship g-formula + Streamlit** (Component 5). | **Promoted to exploratory-causal only if link arrives** | **Gate A** |

**Not blocked on re-linkage:** Components 1, 1b, 3, 4 (Axis A), and tooling. **Blocked / materially upgraded by re-linkage:** Component 2 (flips prior-dependent→data-dominant), Component 5 (ecological→individual-level), the within-species contrast, and the bed-day overlay in Component 4.

---

## 4. Risk register

| Risk | Severity | Mitigation |
|---|---|---|
| **23 (or fewer) confirmed-susceptible** → wide CIs, likely cannot exclude 0; Fiji precision will not transfer | Binding | Lead with crude AJ + single severity-standardized estimate; ship the power simulation as co-primary; reframe deliverable as "best precision a primary SSA-HAI cohort supports + value of re-linkage"; report effective n / weight diagnostics / simulated CI width; Bayesian HDI as the honest uncertainty. |
| **Exposure mis-definition** (`amrp==1` untested miscoded as susceptible in current loader) | Critical (validity) | Pin S = `amrp==0` only; `amrp==1` is a third category; fix loader; report (R,S) ns up front; pre-register. |
| **No patient↔isolate link** | Structural | Run the whole confirmatory triangle within the patient file (`amrp`/`los`/`dead`); use isolate `amrtx` ecologically only; log within-species deviation; pre-register re-linkage decision rules (Components 2 & 5). |
| **Informative exposure ascertainment** (178/336 unascertained, plausibly LOS-driven) | Critical (selection) | Co-primary ascertainment-propensity reweight + MNAR/tipping-point bounds; characterize ascertained vs unascertained; qualify population claim. |
| **ATLAS thin catchment** (3/4 countries only 2021–2023; none past 2023; Ghana singleton; Kenya breakpoint-era jump) | High | Re-base to nowcast not 2025–30 projection; regional-only pooled slope; flat-carried default scenario; ceftazidime not ceftriaxone; pre-register Kenya era handling; lead with data-availability matrix. |
| **Hub public/philanthropic only; geography = funder location** | High | Title everything "public+philanthropic"; Axis A (global vs global) is the only headline; SSA-vs-global to supplement as equity contrast; caveat private R&D may offset. |
| **n=5-6 pathogens → no powered inference** | High | Index is descriptive only; Spearman ρ with CI + explicit caveat; no line-fitting; rely on rank presentation not p-values. |
| **`txadp` adequacy field may not exist** | High | Gate A; downgrade Component 5 to exploratory/ecological until verified; do not claim confirmatory Tier-A in writing. |
| **Cost/estimand transportability** (scaling ATO effect to full population) | High | Match scaling count to estimand population; ATO scales only to overlap count; cost as sourced range + MC CI + unit-cost sensitivity; no final $ recommendation without local validation. |
| **g-null paradox / collider on resistance-as-confounder** | Medium-High | Parsimonious DAG-justified L; resolve resistance/adequacy structure explicitly (drop or controlled-direct-effect); cross-check IPTW vs g-formula. |
| **Fine-Gray not in lifelines** | Medium | rpy2→R `cmprsk` or omit; never list as Python-native. |
| **Prior dominance / prior-data conflict wrong direction** | Medium | Null-centered primary prior; flat-prior posterior reported; effective-prior-sample-size; informative priors sensitivity-only; reconcile MBIRA citation. |
| **Tooling / Python pin mismatch** (config 3.12 vs 3.13 interpreter; deps uninstalled) | Medium | Reconcile pin; pinned deps + import smoke-test; Tooling deviation. |
| **Pre-registration inversion** (estimand & priors chosen post-data) | Medium (reputational) | Log every change as post-data deviation; OSF addendum; state plainly priors were not pre-registered before data access; let the null-centered, data-honest version carry the pitch. |
| **Vivli egress / k-anonymity** | Medium | Export only de-identified aggregated posteriors; min-n gate on cell display; synthetic demo data public. |
| **Garden of forks / multiplicity** | Medium | Lock one headline; FDR (`config.FDR_Q`) on per-pathogen breakdowns; lean on partial-pooling shrinkage over many stratified tests. |

---

## 5. Verbatim repo text

### (a) `docs/deviation_log.md` — append these rows

```markdown
| 2026-06-06 | Methodological | **Primary estimand changed from 30-day mortality (Cox, predicted HR 1.5–3.5) to resistance-attributable excess length-of-stay (excess bed-days)** via a competing-risks multistate (Aalen-Johansen) model, with in-hospital death as a competing event for discharge-alive. Mortality retained only as an honestly-underpowered secondary. | The mortality contrast is underpowered (~14–55 exposure-ascertained deaths; the resistant-vs-susceptible split is smaller still) and the verified-precedent literature (Fiji aHR 1.13, CI 0.51–2.53) is null at larger, fully-ascertained n. Excess bed-days is the more defensible, more impactful, competing-risks-appropriate estimand. Pre-registered (OSF 10.17605/OSF.IO/BFQDP) before data access; this is a post-data-access deviation, declared as such. See OSF addendum (2026-06-06). |
| 2026-06-06 | Methodological | **Exposure definition pinned: susceptible = `amrp == 0` only (confirmed all-susceptible). `amrp == 1` (mixed/untested, no resistant isolate seen) is NO LONGER mapped to the susceptible arm; it is a third "not-confirmed-susceptible" category excluded from the primary contrast (carried only in a labelled sensitivity).** `_RESISTANT_FROM_AMRP` and `SPIDAAR_COLUMNS` updated accordingly. | The prior mapping `{2:1, 0:0, 1:0}` contaminated the susceptible arm with untested patients, making the causal contrast indefensible ("what is a susceptible patient?"). The corrected definition yields a smaller, cleaner susceptible n, reported up front before any modeling. |
| 2026-06-06 | Methodological | **Primary time origin = enrolment (t0=0), restricted to [0, τ=28 days]. Admission origin + left-truncation via `enrtpt` is a sensitivity only.** τ set to 28 (not 30) to sit inside the 28–31d `nobsd` window. | Enrolment origin sidesteps both the non-Markov-under-left-truncation tension in the Aalen-Johansen estimator and immortal-time bias on the exposure; it is robust whether `enrtpt` is a day-count or a categorical phase. |
| 2026-06-06 | Methodological | **Two co-primary honesty analyses added:** (a) a pre-data power/precision simulation under the real (resistant ≈135 / susceptible ≈23) design; (b) an exposure-ascertainment selection sensitivity for the 178 unascertained patients (ascertainment-propensity reweight + MNAR/tipping-point bounds). | At this sample size these determine whether the headline effect is identifiable; they are analyses, not deferrable open items. |
| 2026-06-06 | Methodological | **Bayesian partial pooling added as a SECONDARY prior-regularized robustness companion** (country-level pooling only; pathogen as a fixed Gram-negative-vs-other stratifier). Primary prior on the discharge log-HR is **null-centered** (Normal(0, 0.5)); informative LMIC/SSA priors are sensitivity only. Decision rule: if a `pid`↔`iid` re-linkage extract arrives before lock, the isolate-ascertained cohort becomes primary and all priors revert to sensitivity. | Cited discharge-hazard evidence (MBIRA discharge CSH ratio 1.16; Fiji 0.99) is null-to-reversed, so an optimistic discharge-slowing prior is unjustified; partial pooling is a precision-stabilizer, not the headline (strategy doc flags it as crowded). The informative prior was NOT pre-registered before data access — declared here. |
| 2026-06-06 | Methodological | **Step 2 (ATLAS) re-based from a 2025–2030 per-country projection to a catchment resistance NOWCAST (last-observed level) + a SPIDAAR frame-contrast.** Per-country time slopes to 2030 dropped (not data-identified: Ghana/Malawi/Uganda only 2021–2023, Kenya 2013/2014/2021–2023, none past 2023). Panel locked to Enterobacterales × ceftazidime (ceftriaxone interpretation blank in catchment) re-derived from MIC under EUCAST v15.0; single regional pooled time slope only; Kenya pre-2020 breakpoint-era handled by indicator/exclusion. | The granted ATLAS extract does not contain enough catchment years to identify per-country trends; honest nowcast + frame-contrast is defensible where a 2030 projection is not. |
| 2026-06-06 | Tooling | Added survival/Bayesian dependencies to `pyproject.toml` (lifelines, scikit-survival, statsmodels, scipy, pymc, arviz, preliz; optional rpy2) with version pins; reconciled `config.PYTHON_VERSION` with the secure-environment interpreter; added an import smoke-test. | None of these were installed; the excess-LOS, Bayesian, and projection modules require them. Mirrors the xlrd precedent (2026-05-31). No change to any analysis method beyond enabling it. |
| 2026-06-06 | Methodological | **Step 5 (stewardship g-formula) downgraded from confirmatory capstone to a pre-registered EXPLORATORY causal what-if tool, gated on verifying that a patient-side empiric-adequacy field exists.** If it does not, the treatment node is the isolate-side `amrtx` (ecological only) and the Streamlit tool ships as a local-calibration calculator, not an estimated individual effect. | The verified patient loader carries no adequacy field; the only adequacy signal (`amrtx`) is isolate-side and unlinkable. Confirmatory status cannot be claimed until the field is verified in the secure environment. |
```

### (b) OSF pre-registration ADDENDUM — write verbatim to `docs/osf_addendum_2026-06-06.md`

```markdown
# OSF Pre-Registration Addendum — AMR Sentinel (Vivli 2026)

**Parent pre-registration:** https://doi.org/10.17605/OSF.IO/BFQDP (deposited 2026-05-29, before data access).
**Addendum date:** 2026-06-06. **Authored after data access.** This addendum is filed under
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
susceptible) counts under this rule are reported before any modeling.

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
  (resistant ≈135 / susceptible ≈23) design across a grid of true effects and
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
At the delivered sample size (≈23 confirmed-susceptible patients) the headline interval
may not exclude zero. We pre-commit to reporting the interval honestly and to framing the
deliverable as the most precise resistance-attributable bed-day estimate a primary SSA-HAI
cohort can currently support, together with an explicit quantification of how much the
pending re-linkage would tighten it. We will not present a prior-driven or selection-biased
point estimate as a detected effect.
```

---

## 6. What still needs the team's decision

1. **Run the gate verifications in the secure env (blocking everything):** Gate L (`los`/`disev`/`ward`/`enrtpt` actually in `patientdata.xls`), Gate T (`amrp` = index isolate, exposure fixed at baseline), Gate E (`enrtpt` day-count vs phase), Gate A (does a patient-side empiric-adequacy field exist?), Gate X (ATLAS MIC/interpretation presence per panel cell). **If Gate L fails, Component 1 has no outcome variable — escalate immediately.**
2. **Confirm the realized (R, S) ns** under S=`amrp==0`-only, and decide whether the susceptible arm is large enough to proceed or whether the entry pivots entirely to "value of re-linkage."
3. **Sign off the overlap-population narrowing** if any IPTW/ATO is shown, and lock the burden/cost scaling rule per estimand.
4. **Lock τ = 28** (vs 30/31) and confirm no clinically-motivated alternative.
5. **Decide the ascertainment-sensitivity strategy** (reweight vs MNAR bounds vs both) and lock it in the addendum.
6. **Lock the Bayesian prior set** (null-centered primary confirmed) and the σ priors.
7. **Reconcile the MBIRA citation** (strategy doc "associated with mortality" vs spec death CSH 0.74) — re-pull PMC7617135 before the pitch.
8. **Source and lock country-specific bed-day unit costs** (WHO-CHOICE vs national HTA) with a citation and uncertainty range.
9. **Lock the R&D index low-denominator floor value** and decide whether Axis B appears at all (supplement vs cut).
10. **Decide Component 5's status** explicitly: exploratory-causal (Gate A passes) vs ecological-only calculator vs gated on re-linkage; and confirm Vivli egress/k-anonymity thresholds for shipping posterior artifacts.
11. **Reconcile `config.PYTHON_VERSION` ("3.12")** with the actual interpreter and lock pinned dependency versions.
12. **Decide the illness-death 4-state ward-escalation model's status** (exploratory exhibit vs cut) once escalation-transition density is seen.

**Relevant files:** `C:\projects\amr-sentinel-vivli2026\src\amr_sentinel_vivli\data_loading.py` (loader/exposure mapping to change), `...\config.py` (Python pin, propensity covariates, panel to add), `...\bayesian_projection.py` / `stewardship_gformula.py` / `rd_alignment.py` (stubs to implement), `docs\deviation_log.md` and the new `docs\osf_addendum_2026-06-06.md` (text in §5), plus new `src\amr_sentinel_vivli\excess_los.py` and `app\streamlit_app.py`.