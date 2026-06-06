# Winning strategy — Vivli 2026 AMR Data Challenge (research-backed pivot)

**Status:** proposed pivot, pending team sign-off. Drafted 2026-06-06 from a
deep-research pass (109 agents, 25/25 claims adversarially verified 2/3+) plus
local validation of the SPIDAAR patient file. Supersedes the EOI's *headline*
while staying inside the SPIDAAR-anchored, multi-dataset spine Vivli granted
access for. Any adopted change is logged in `docs/deviation_log.md` (pre-reg §11).

## TL;DR

Move the headline **off the resistance→mortality bridge** (fragile: ~14
exposure-ascertained deaths, and the literature shows even well-powered cohorts
often can't confirm it) and **onto resistance-attributable excess length-of-stay
(excess bed-days)** estimated with competing-risks methods — an outcome SPIDAAR
uniquely holds, that is far more statistically tractable, and that translates
directly into an economic/stewardship/R&D-prioritization story. Frame the whole
entry as a **primary-cohort calibration/transparency probe** on modelled global
(GRAM) burden estimates in the SSA settings where those models are weakest, and
explicitly contest the **Global AMR R&D Hub Cross-Domain Award**.

## Two new hard constraints discovered (force an EOI re-cast regardless)

1. **Merck SMART is excluded from the Challenge.** All AMR Register datasets are
   available *except* Merck datasets; SMART can only be requested *outside* the
   Challenge. This breaks EOI Method 2 ("ATLAS + SMART projection") as written —
   the projection must be **ATLAS-only** (~917K isolates, 2004–2024).
   *(Vivli 2026 Data Challenge Overview.)*
2. **The registered isolate-level Cox model is impossible** (no patient↔isolate
   key) and the patient-level fallback is power-starved. Independent literature
   confirms the mortality angle is the *least* defensible: the Fiji
   Enterobacterales BSI cohort (n=162, 36 deaths) collapsed to a null adjusted
   HR 1.13 (95% CI 0.51–2.53) — an order of magnitude more deaths than SPIDAAR's
   exposure-ascertained ~14, and still null. A small cohort cannot headline a
   mortality HR. *(Loftus et al., Lancet Reg Health W Pac 2022.)*

## The pivot — primary headline

**"The bed-day burden of antibiotic resistance in sub-Saharan African hospitals:
a primary-cohort estimate, and what it says about global burden models and R&D
priorities."**

- **Primary outcome:** resistance-attributable **excess length-of-stay /
  excess bed-days**, via a **competing-risks multistate model** — Fine-Gray
  subdistribution hazard for cumulative discharge incidence **plus** cause-specific
  Cox — with **death as the competing event** (death precludes discharge, so naive
  KM/Cox overestimate discharge incidence; Fine-Gray is the validated fix, shown
  in a Malawi paediatric cohort, n=15,463). *(PMC10911049.)*
- **Why it wins on defensibility:** in the very precedent where mortality was
  null, the *same data* yielded a precisely-estimated excess-LOS (Fiji: 2.6 days,
  95% CI 2.5–2.8). LOS is detectable where mortality is not. SPIDAAR uniquely has
  LOS among Vivli registers, and is exclusively HAI ≥48h — so excess-bed-days is
  both novel and economically translatable.
- **Why it wins on impact:** excess bed-days × local bed-day cost = a concrete
  cost figure local stewardship programs and funders can act on. AMR is a
  top-tier killer in the catchment (Western SSA has the world's highest
  AMR-attributable death rate; in the WHO African region AMR-associated deaths
  exceed HIV and malaria). *(Murray et al., Lancet 2022; Lancet Glob Health 2023.)*

## Local data validation (SPIDAAR patient file, 336 patients)

| Field | Completeness | Role |
|---|---|---|
| `los` (discharge LOS) | 255/336 | event-of-interest time (discharge) |
| `dtpta` (days to death) | 56 | competing-event time (in-hospital death) |
| `nobsd` (observation days) | 336 | censoring / follow-up window (~28–31d) |
| `dead` (0/1/9) | 336 | 56 deaths total (55 + 1 censored-deceased) |
| `disev` (severity 1/2/3) | 336 | key confounder |
| `ward` (ICU/surg/med) | 336 | confounder |
| `amrp` (exposure) | 336 | **but ascertained only for 158** (−1 for 178) |

**Feasible:** discharge (255) + death (56) competing events give real outcome
power. **Honest limit:** exposure is ascertained for 158 patients with only
**23 susceptible** → wide CIs; Fiji-level precision will NOT transfer. Mitigations:
Bayesian partial pooling across the 4 countries/pathogens; informative priors from
published excess-LOS estimates; and the re-linkage request (would recover exposure
for far more patients via the better-ascertained isolate file — 3GC-R ascertained
in 204/244).

## Re-cast of the five EOI methods

1. **SPIDAAR cohort →** competing-risks **excess-LOS** as the *primary* estimand
   (Fine-Gray + cause-specific Cox; propensity/severity adjustment; e-value for
   unmeasured confounding). Mortality is retained as a **secondary, honestly
   underpowered** result with wide CIs — reported, not headlined.
2. **Projection →** **ATLAS-only** Bayesian hierarchical resistance trajectories
   (SMART excluded). Scope claims to ATLAS coverage.
3. **Bridge →** convert projected resistance prevalence to **excess bed-days /
   cost** per country-pathogen (LOS-based, defensible) as primary; keep a
   mortality bridge **only** as an explicit **calibration comparison against GRAM**
   — framed as transparency/honesty, *not* as a claim that ~14 events can
   calibrate a global model.
4. **R&D Mismatch Index →** keep; **target the Cross-Domain Award** (combinable
   with other awards, open to non-students; requires the Hub DB + a Register
   dataset — we satisfy this). **Scope to public + philanthropic funding only**
   (Hub does not capture private/industry R&D) to avoid overclaiming.
5. **Stewardship g-formula simulator →** keep as the causal/counterfactual
   component + shipped Streamlit tool; recast the "what-if" outcome as
   **avertable bed-days** (mortality secondary). This is the element that matches
   recent winners (causal inference + reproducible open-source code).

## Award-targeting

- **Primary:** Grand Prize / Visionary — judged on **Innovation** ("how creative")
  and **Impact** ("effect on the field if implemented"), plus Methodology/Design;
  finalists give a **live Zoom pitch**. Recent top prizes went to causal-inference
  + reproducible-code entries (2024 Allel spatiotemporal early-warning; 2025
  Vellore causal climate/pollution; 2025 Vania–Gupta counterfactual ML policy).
  Our differentiator: a **real primary outcome cohort** where winners used global
  surveillance data — provided we match them on causal framing + shipped tooling.
- **Cross-Domain (Global AMR R&D Hub):** explicitly contest via Method 4.
- **Rigor layer, not headline:** Bayesian partial pooling and LMIC/equity framing
  are judge-recognized but **crowded** (multiple past winners used them). Use to
  stabilize small per-stratum counts and borrow strength — do not lead with them.

## Positionality (de-risk)

Non-catchment authors → lead with the **open-source-infrastructure-for-local-
researchers** framing the EOI already committed to; keep the senior catchment-
region reviewer in the loop; claim **no stewardship recommendation as final
without local validation**. This converts the positionality risk into the
"reproducible tooling + equity" strength judges reward.

## Key risks / caveats (from adversarial verification)

- **Geographic transfer:** the Fiji null is a *power/method* illustration only —
  the directly comparable SSA study (MBIRA, 8 SSA hospitals) found 3GC-R **was**
  associated with mortality. SPIDAAR's own data must drive any directional claim.
- **Calibration-probe inference** (14 events vs GRAM) is *our* framing; present as
  motivation/transparency, not a validated capability.
- **Exposure imbalance / ascertainment** (23 susceptible; 178 unascertained) is
  the binding limitation — partial pooling + re-linkage are the mitigations.
- Hub DB = public + philanthropic only.

## The path (sequence)

1. **Submit the re-linkage request now** (`docs/vivli_relinkage_request.md`) — it
   gates how strong every angle can be; if a `pid↔iid` extract arrives, exposure
   ascertainment is largely fixed.
2. Team sign-off on this pivot → write the **deviation-log entries** (SMART
   excluded → ATLAS-only; primary estimand mortality → excess-LOS) and a
   **supplementary pre-registration** addendum on OSF before touching outcomes.
3. Implement the **competing-risks LOS module** (Step 1 re-cast) on the validated
   patient frame; partial-pooling layer.
4. ATLAS-only projection; bed-day bridge; R&D Mismatch (Cross-Domain); g-formula
   avertable-bed-days simulator + Streamlit.
5. 5-page report + GitHub/Zenodo; rehearse the live pitch around the bed-day
   headline.

## Sources (verified)

- Vivli 2026 Data Challenge Overview (criteria, SMART exclusion, Cross-Domain):
  https://amr.vivli.org/data-challenge/data-challenge-overview/
- 2024 / 2025 finalist & award pages: https://amr.vivli.org/data-challenge/2024-finalist-and-award-winning-solutions/ ,
  https://amr.vivli.org/data-challenge/2025-finalist-and-award-winning-solutions/
- Murray et al., GRAM, Lancet 2022: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02724-0/fulltext
- WHO African region GRAM, Lancet Glob Health 2023: https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(23)00539-9/fulltext
- Loftus et al., Fiji Enterobacterales cohort (null mortality, precise excess-LOS):
  https://www.sciencedirect.com/science/article/pii/S2213716522001515
- Malawi competing-risks discharge cohort: https://pmc.ncbi.nlm.nih.gov/articles/PMC10911049/
- Global AMR R&D Hub dynamic dashboard / methodology: https://globalamrhub.org/dynamic-dashboard/
