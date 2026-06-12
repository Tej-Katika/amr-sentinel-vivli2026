# The AMR burden paradox in sub-Saharan African hospitals

### A competing-risks reframing of what antibiotic resistance does to patients — and where the leverage to act actually lies

**Vivli 2026 AMR Surveillance Open Data Re-Use Data Challenge — Final Submission (DRAFT)**
Tejashwar Reddy Katika (Independent Researcher, University of North Texas; Lead) · Akhilesh Reddy Katika (MS Data Science, Flinders University)
Datasets: **SPIDAAR** (primary) · **Pfizer ATLAS** · **Global AMR R&D Hub** · Pre-registration: [OSF 10.17605/OSF.IO/BFQDP](https://doi.org/10.17605/OSF.IO/BFQDP) · Code: GitHub/Zenodo (Apache-2.0, 135 tests) · **Cross-Domain Award eligible**

> **Draft status.** Figures here are from the development run on the delivered files; the confirmatory run is executed in the Vivli secure environment. The R&D-Hub funding snapshot, the GRAM per-pathogen burden, and the SPIDAAR codebook gates (exposure and adequacy coding) are now resolved and the Cross-Domain index is computed; the only un-fetched Component-4 input is the $113M split among three Gram-negatives (appendix-locked), which is propagated as uncertainty rather than guessed and does not affect the ranking. Every deviation from the pre-registration is logged in `docs/deviation_log.md` and a supplementary OSF addendum.
>
> **Figures** are regenerated from the committed code and are gitignored under the Data Use Agreement (the public repository carries synthetic-mirror outputs only); the submitted PDF embeds the egress-reviewed versions produced in the secure environment. Figure paths below resolve when the report is built where those files exist.

---

## 1. The question, and why it is not the obvious one

Global models attribute an enormous mortality burden to antimicrobial resistance (AMR) — 1.27 million deaths *attributable* and 4.95 million *associated* in 2019 (Murray et al., *Lancet* 2022) — and sub-Saharan Africa (SSA) carries the highest AMR-attributable death rate of any world region. SPIDAAR is the **only register in the Vivli AMR catalogue that carries patient outcomes** (length-of-stay and mortality) for this region, so it is uniquely positioned to ask the question the burden models cannot answer from surveillance alone: *what does the resistance phenotype actually do to an individual patient, and what follows for stewardship and R&D priorities?*

Our pre-registered plan anchored on a resistance→mortality bridge. After data access, three facts forced a disciplined pivot — and the pivot **is the contribution**:

1. **The mortality contrast is not estimable here.** Resistance is ascertained for only 156 of 336 patients (135 resistant, **21 susceptible**), ~17 deaths split across arms — and the far larger, directly comparable SSA literature (MBIRA, 8 hospitals, 878 BSI; Fiji) is null after adjustment, which we formalise in a pooled synthesis (§3.1, Figure 2).
2. **SMART is excluded** from the Challenge, so external surveillance is **ATLAS-only**.
3. **SPIDAAR uniquely holds length-of-stay**, detectable where mortality is not (Fiji: null mortality, yet a precise 2.6-day excess-LOS).

So we move the headline onto **resistance-attributable excess bed-days** estimated with **competing risks**, and let the data tell an honest, counter-intuitive story: in this cohort, as in the authoritative SSA literature, **resistance is not a clean per-patient killer**. The actionable burden is **systemic** — whether a patient receives *adequate empiric therapy* — not the phenotype itself. This reframing is the contribution: it turns an underpowered, contradicted mortality question into a tractable, decision-relevant one, and points two deliverables at the real lever — an empiric-adequacy stewardship simulator (the centerpiece) and a competing-risks correction to how AMR burden is costed.

## 2. Data and design

**SPIDAAR** (Ghana, Kenya, Malawi, Uganda): 336 hospitalised patients with healthcare-associated infections; length-of-stay, in-hospital death, severity, ward, empiric-therapy adequacy, and a patient-level resistance summary; plus 244 isolates with per-mechanism resistance. **ATLAS**: 1.0 M isolates, of which 1,519 fall in the catchment; the analysable cell is **Enterobacterales (E. coli + K. pneumoniae) × ceftazidime** (665 isolates; ceftriaxone interpretation is blank in the catchment, and no catchment data exist after 2023). **Global AMR R&D Hub**: public + philanthropic funding, frozen at a dated snapshot.

The exposure is pinned to a defensible contrast — **resistant = `amrp==2`; susceptible = `amrp==0` only** (untested/unascertained excluded). All estimators are competing-risks-aware and reproducible (master seed `20260526`); because the secure environment lacks specialised survival/Bayesian libraries, they are implemented directly in NumPy/pandas and unit-tested on synthetic data (135 tests). Full specification: `docs/analysis_plan_2026.md`.

## 3. Results

### 3.1 Resistance does not prolong stay — a triangulated null (Components 1, 1b, 2)

The headline estimand is the difference in restricted-mean time in the *admitted* state to τ=28 days (excess bed-days), with in-hospital death competing with discharge, time origin at enrolment.

| Estimate | Excess bed-days | Interval |
|---|---|---|
| Crude (Aalen-Johansen Δ-RMST) | **−1.65 d** | 95% bootstrap CI [−4.89, +1.90] |
| Severity-standardised | **−1.87 d** | negative within both severity strata |
| Bayesian, null-centred prior | **−0.94 d** | 95% HDI [−4.3, +2.8]; P(excess>0)=0.30, P(>1d)=0.14 |

The point estimate is **negative** — resistant patients leave the admitted state *faster*, not slower. The competing-risks decomposition shows why (cumulative incidence at day 28): resistant patients both discharge more (63.3% vs 57.1%) **and** die more (11.9% vs 4.8%), while susceptible patients linger (still-admitted 38.1% vs 24.8%). Both exits truncate bed-day accrual, so naïve length-of-stay comparisons that ignore death are biased.

![Bed occupancy by resistance arm and competing outcomes at day 28](../figures/excess_los_stateoccupation.png)

***Figure 1.** Left: probability of remaining admitted over time by arm — the area between the curves to τ=28 is the excess bed-days (resistant occupancy sits below susceptible). Right: competing discharge / death / still-admitted cumulative incidence at day 28. Resistant patients exit faster by both routes; susceptible patients linger.*

**Honesty analyses are co-primary (Component 1b).** With 21 susceptible patients the design is underpowered — ~21% power for even a 2.6-day effect, a ~9-day-wide interval. The fix is the **patient↔isolate re-linkage** requested from Vivli, and the isolate file shows why (Component 1d): the mechanism-resolved phenotype the patient summary collapses is ascertained far more completely — **3GC-R in 204 isolates, MDR in 231, MRSA in 44** (the two Gram-negative mechanisms travel together, P(MDR\|3GC-R)=0.91) — and **MDR ascertainment alone recovers a 51-isolate susceptible arm against the patient frame's 21**, a ~1.5× tightening that would roughly halve the interval (~4.9 d) and lift power to ~55%. We ship the merge machinery and projection now (no positional join — the counts do not reconcile). Crucially, the **exposure-ascertainment selection** that could most plausibly fake a result runs the *wrong way*: ascertained patients are longer-staying and lower-mortality, yet the excess is negative — surviving inverse-ascertainment reweighting (−1.85 d) and the missing-at-random anchor (−0.61 d).

**A formal evidence synthesis, not narrative triangulation (Component 1c).** We pool the *adjusted* SSA/LMIC evidence on the resistance→in-hospital-death hazard ratio in a random-effects Bayesian meta-analysis (Figure 2). The two adjusted cohorts pool to a **null HR of 0.88 (95% credible interval 0.46–1.79; P(HR>1)=0.33)** with a **0.30–2.77 prediction interval** for the next SSA setting; our own cohort adds only a crude, severity-confounded HR (2.76, 0.37–20.8) that barely moves the pool. Across the adjusted evidence, **resistance is not an independent per-patient mortality driver once severity and access are controlled** — which is why the headline is the excess-bed-days estimand SPIDAAR uniquely supports, not another underpowered mortality contrast.

![Random-effects synthesis of resistance→in-hospital-death across SSA/LMIC cohorts](../figures/evidence_forest.png)

***Figure 2.** In-hospital-death hazard ratio (3GC-R vs 3GC-S), log scale. Adjusted external cohorts (MBIRA, Fiji), our crude cohort estimate, the random-effects pooled HR (diamond), and the prediction interval for a new setting. The pooled estimate and every contributing adjusted study straddle the null (HR=1).*

### 3.2 The frame-shift: severe-HAI cohorts run far hotter than surveillance (Component 3)

Because per-country 2025–2030 projections are not data-identified (no catchment ATLAS data after 2023), Step 2 is re-based to a partial-pooled **resistance nowcast** (empirical-Bayes beta-binomial; pooled ceftazidime-R **0.68**; Ghana 0.47, Kenya 0.65, Malawi 0.78, Uganda 0.76) led by a data-availability matrix that states the sparsity up front. Its novel output is a deliberate **frame-contrast**: SPIDAAR's severe-HAI inpatient 3GC-R prevalence is **0.85** against ATLAS mixed-surveillance ceftazidime-R of **0.62**, with the same ordering in every country. The severity/HAI frame shifts measured resistance upward by ~20 percentage points — a transparency caveat for anyone using mixed-surveillance prevalence to reason about inpatient burden.

![SPIDAAR severe-HAI vs ATLAS surveillance 3GC resistance by country](../figures/spidaar_framecontrast.png)

***Figure 3.** 3GC resistance prevalence (Wilson intervals) in the comparable Enterobacterales panel: ATLAS mixed surveillance (ceftazidime-R) vs SPIDAAR severe-HAI inpatients (3GC-R), by country and overall. The severe-HAI frame runs higher in every country.*

### 3.3 The leverage is empiric adequacy, not the phenotype (Component 5 — centerpiece)

If resistance itself is not the per-patient driver, the lever is getting the *right empiric drug* to the patient. A competing-risks **g-formula** on empiric-therapy adequacy (`txadp`; 106 adequate / 52 inadequate / 178 unknown) estimates the counterfactual of raising adequacy, with death competing with discharge. Resistance is treated as a *determinant* of adequacy — a mediator, so it is excluded from the confounder set. The result is the burden paradox made quantitative and clinically coherent: **adequate empiric therapy averts deaths (+2.6 percentage points) while *adding* bed-days**, because patients who would have died now survive to discharge. The death-aversion benefit **concentrates in the resistant arm** (+4.7 pp averted; zero in the susceptible arm) — exactly where first-line therapy is most likely to fail. Positivity is healthy (4% off-support).

**Hardened, not hyped.** A stratified bootstrap (within each confounder×arm cell, preserving positivity) puts honest — and wide — intervals on these contrasts: pooled death-aversion +2.6 pp (95% CI −5.5 to +10.0), resistant arm +4.7 pp (−3.6 to +13.2), bed-days +2.7 (−0.9 to +6.0). The direction holds but the null is not excluded — hence *exploratory*. An **E-value** (VanderWeele & Ding 2017) shows an unmeasured confounder would need a risk-ratio of **2.06 (pooled) / 3.03 (resistant arm)** with both adequacy and death to explain away the point estimate — non-trivial, though weak confounding reaches the near-null CI limit. Conclusive proof needs the larger re-linked cohort, not this one.

This ships as a re-runnable **Streamlit what-if tool**: a catchment-region stewardship programme enters its patient volume, current vs target empiric adequacy, and country bed-day cost (WHO-CHOICE: Ghana $6.30, Kenya $5.45, Malawi $3.25, Uganda $3.81 per bed-day, 2010 USD) and reads off projected deaths averted, bed-days added, and cost. The tool carries only a de-identified calibration artifact (no patient records) and labels every output an *ecological-calibration scenario, not an individual effect*.

> **Gate A — resolved (codebook-confirmed).** The adequacy coding (`txadp` Adequate=0, Inadequate=1, Unknown=9) is confirmed against the official SPIDAAR codebook, so the result sign is settled (adequate therapy averts deaths). The exposure coding (`amrp` resistant=2 / susceptible=0) is likewise codebook-confirmed. The analysis remains exploratory only because of the small inadequate arm (n=52), not the coding.

### 3.4 Burden vs R&D investment (Component 4 — Cross-Domain)

**This component meets the Global AMR R&D Hub Cross-Domain Award brief directly: it joins the Hub funding dataset to two Vivli AMR Register datasets (SPIDAAR and ATLAS).** The Cross-Domain index aligns global GRAM burden share with Hub public+philanthropic funding share per pathogen (log2 mismatch + Spearman ρ; n=5–6 pathogens, descriptive only, no fitted line), reporting the **non-pathogen-specific (cross-cutting) funding magnitude first** and bounding sparsely-funded pathogens with a pre-specified floor. The funding denominator is locked to a frozen, peer-reviewed extract of the Hub Dynamic Dashboard (Czaplewski et al., *Lancet Microbe* 2026; 2017–2023), preferred over a live snapshot for reproducibility.

**Cross-cutting first.** Of US$2.51 billion in public+philanthropic antibacterial-therapeutics R&D (130 funders, 2017–2023), only **$1.06 billion — 42% — is pathogen-specific**; the majority (58%) is cross-cutting, so every per-pathogen claim is appropriately widened. Within the pathogen-specific minority, spend is dominated by *M. tuberculosis* ($474M, ~45% of it) and by Gram-positive/niche targets (*S. aureus* $142M, *C. difficile* $141M, *N. gonorrhoeae* $101M, *P. aeruginosa* $87M).

**Burden–funding alignment (computed).** Pairing the verified GRAM-2019 per-pathogen burden (appendix Table S22, 95% UIs) against this funding, antibacterial R&D shows **no positive alignment with AMR mortality burden** across the six leading pathogens — Spearman ρ from **+0.14** (associated deaths) to **−0.49** (attributable DALYs), CI spanning ±1 at n=6. The floored per-pathogen log2 mismatch (Monte-Carlo over the burden UIs and the unfetched $113M Gram-negative split) is robust: **_S. pneumoniae_ (≈ +2.9) and _K. pneumoniae_ (≈ +2.0) are most under-funded** (together ≈99% of the "most under-funded" probability), while ***P. aeruginosa* (≈ −1.4) and _S. aureus_ (≈ −1.0) are over-funded** — the inverse alignment the systemic-leverage thesis predicts: money flows to a coordinated TB programme and to novel-compound targets, not to the highest-burden community Gram-negatives where access and empiric adequacy are the lever.

**Made local — the Cross-Domain finding in the catchment (computed).** The global ranking answers a global question; the Award asks what *our* data say. Re-weighting the mismatch by the pathogen mix in the SPIDAAR severe-HAI isolates (the only Register source spanning all six panel species; resistant-isolate frequency as a local burden proxy, Dirichlet uncertainty over the counts) **re-orders the priority** (Figure 4): *K. pneumoniae* (≈ +2.5) and *E. coli* (≈ +1.2) become most under-funded, while *S. pneumoniae* — top of the *global* ranking — **flips to over-funded because it is nearly absent from SSA hospital infection (3 of 146 panel isolates)**. The pathogens that actually drive resistant HAI in the catchment are precisely the community Gram-negatives global R&D under-funds, not the pneumococcus. *(An isolation-frequency proxy, not a mortality burden; ρ ≈ 0 under either weighting.)*

![Burden–funding mismatch per pathogen: global GRAM vs SSA catchment](../figures/rd_mismatch_global_vs_catchment.png)

***Figure 4.** Per-pathogen log2 burden/funding mismatch (>0 = under-funded relative to burden) under the global GRAM burden vs the SSA severe-HAI catchment burden. The catchment lens lifts the Gram-negatives (*K. pneumoniae*, *E. coli*) and flips *S. pneumoniae* from most under-funded globally to over-funded locally.*

## 4. Innovation and impact

**Innovation.** The novelty is not the estimator — competing-risks length-of-stay is established — but the *reframing*: identifying excess bed-days under a broken patient↔isolate link in an LMIC HAI cohort, turning the resulting null into a systemic-leverage thesis, and shipping it as a re-runnable local-calibration tool. We treat the small-sample limitation as a quantified result, not a caveat to hide.

**Impact.**

- **A transparency correction to AMR costing.** "Excess bed-days × cost" attributions that ignore competing mortality over-state the bed-day burden of resistance; making the death channel explicit (resistant patients exit *faster*, partly by dying) yields the honest figure and shifts the economic case onto empiric adequacy.
- **An actionable, locally-parameterised tool.** The stewardship simulator converts a defensible counterfactual — adequacy averts deaths, concentrated in the resistant arm — into deaths/bed-days/cost a programme can act on, re-runnable with local antibiograms and WHO-CHOICE costs.
- **R&D realignment (Cross-Domain).** If access and empiric adequacy are the lever, the mismatch index argues for reading R&D alignment against access and diagnostics, not novel compounds alone.
- **Honest inputs to the burden debate.** The frame-contrast (severe-HAI 0.85 vs surveillance 0.62) and the explicit power/precision accounting give modellers and funders calibrated, caveated SSA-HAI inputs where those data are thinnest.

## 5. Limitations (stated plainly)

The susceptible arm (n=21) is small and the headline interval does not exclude zero; we report this as the finding, not around it, and quantify what re-linkage would buy. The patient↔isolate link is broken, so the contrast is patient-level, not within-species. ATLAS catchment coverage is sparse and ends in 2023 (nowcast, not forecast). The Hub captures public+philanthropic funding only. The adequacy analysis is exploratory and Gate-A-dependent. None of these is hidden; each is logged and bounded.

## 6. Reproducibility and deliverables

Open-source pipeline (Apache-2.0, GitHub + Zenodo): SPIDAAR/ATLAS loaders, competing-risks excess-LOS + sensitivity + figure, Bayesian companion, the random-effects Bayesian evidence synthesis, the isolate-mechanism breakdown + relinkage machinery, ATLAS nowcast + frame-contrast, the global and catchment-specific R&D mismatch indices, the g-formula (with stratified-bootstrap CIs + an E-value), and the Streamlit tool; 135 unit tests on synthetic data; pre-registration + deviation log + OSF addendum. No stewardship recommendation is claimed final without local validation; the framework is built as infrastructure for catchment-region researchers, who have been engaged for independent review.

---

*References on a separate page (`docs/references_2026.md`): Murray et al., Lancet 2022 (GRAM); Aiken et al., Lancet Infect Dis 2023 (MBIRA); Loftus et al., JGAR 2022 (Fiji); VanderWeele & Ding, Ann Intern Med 2017 (E-value); Stenberg et al., Cost Eff Resour Alloc 2018 (WHO-CHOICE bed-day costs); Czaplewski et al., Lancet Microbe 2026 (Global AMR R&D Hub funding 2017–23). Full verified figures and verbatim source extracts in `docs/reference_verified_2026-06-06.md` and `docs/reference_rd_alignment_2026-06-09.md`.*
