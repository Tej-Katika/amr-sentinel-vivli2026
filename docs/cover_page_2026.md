# Cover Page — Vivli 2026 AMR Data Challenge (MS Form content)

> **How to use this file.** The official Cover Page is a Microsoft Form. This
> document holds paste-ready content for the likely fields plus framing engineered
> to map onto the four *unweighted* judging dimensions (Methodology · Design ·
> Innovation · Impact). Map each block to the corresponding form field; where the
> form's wording differs, keep the content and adjust the label. Nothing here
> restates the 5-page report — it *sells* it. Keep every number identical to
> `docs/final_report_2026.md` (numbers are reconciled against the 91/91 confirmatory
> manifest).

---

## 1. Title

**The AMR burden paradox in sub-Saharan African hospitals** — a competing-risks
reframing of what antibiotic resistance does to patients, and where the leverage to
act actually lies.

## 2. Team

| Role | Name | Affiliation | Status |
|---|---|---|---|
| **Lead Applicant** | Tejashwar Reddy Katika | Independent Researcher; University of North Texas | Non-student |
| Co-applicant | Akhilesh Reddy Katika | MS Data Science, Flinders University | Non-student |

*Team size 2 (within the 2–5 rule). Contact e-mail: tejashwar1029@gmail.com.*

> **Award categories to check on the form (DECIDED).** Select **three**: the
> **Global AMR R&D Hub Cross-Domain Award** (primary path; evaluates student +
> non-student together), the **Impact Award (non-student)**, and the **Innovation
> Award**. Do **not** select any Student award — the Lead is a non-student, which
> bars the Student pools.

## 3. Datasets used (establishes Cross-Domain eligibility)

- **SPIDAAR** — AMR Register, *primary* (Ghana, Kenya, Malawi, Uganda; the only
  Register dataset carrying patient outcomes: length-of-stay + mortality).
- **Pfizer ATLAS** — AMR Register (external surveillance; full 1,011,168-isolate
  register used for the surveillance-mismatch axis).
- **Global AMR R&D Hub** — funding dataset (Cross-Domain).

> **Cross-Domain Award eligibility, stated explicitly:** this entry uses the Global
> AMR R&D Hub dataset **plus two AMR Register datasets** (SPIDAAR + ATLAS).
> Components 4 and 6 join Hub funding to the Registers to build a burden ↔ funding ↔
> surveillance mismatch — the Cross-Domain brief, met directly.

## 4. Plain-language summary (lay abstract)

Global models blame antibiotic resistance for a huge death toll, and sub-Saharan
Africa carries the highest resistance-attributable death *rate* on earth. Using the
only Vivli register that follows patients (SPIDAAR: 336 hospital patients across
four African countries), we asked what resistance actually *does* to an individual
patient — and found the headline is not what the global numbers imply. Once you
account for the fact that death and discharge compete, **resistance did not prolong
hospital stay and was not an independent killer** — a result that matches the
largest, best-adjusted African evidence. The real, actionable burden is
**systemic: whether a patient gets the right antibiotic early** (adequate empiric
therapy), which we quantify and ship as a re-runnable stewardship tool. And the
system's resources point *away* from the burden on **three axes at once** —
research **funding**, disease **surveillance**, and the burden itself all neglect
the same community Gram-negative infections in sub-Saharan Africa (the region holds
just **2.3%** of the world's surveillance isolates; the USA alone holds ~7× more).
Our call is a **reallocation**: point the levers we control — empiric-therapy
access, R&D funding, and surveillance — at the burden they now bypass, and we ship
the instruments to act on it.

## 5. Contribution mapped to the four judging dimensions

> Reviewers score Methodology, Design, Innovation, and Impact separately and
> unweighted. Each block below answers one.

**Methodology (appropriateness + novel application).**
Competing-risks throughout, implemented from scratch in NumPy/pandas because the
secure environment lacks survival/Bayesian libraries, and unit-tested (155 tests):
Aalen-Johansen Δ-RMST for excess bed-days with death competing with discharge; a
random-effects **Bayesian** meta-analysis pooling the adjusted SSA/LMIC
resistance→death evidence; a competing-risks **g-formula** for the empiric-adequacy
counterfactual with stratified-bootstrap CIs and an **E-value** for unmeasured
confounding; an empirical-Bayes beta-binomial resistance nowcast. Every headline
number is re-derived by a self-checking confirmatory harness against a pinned
manifest (**91/91 within tolerance**). This sits squarely in the causal /
Bayesian / counterfactual family that has won prior Challenges.

**Design (alignment with the research questions).**
Pre-registered (OSF) on a resistance→mortality bridge; after data access, three
facts (mortality not estimable at n=21 susceptible; SMART excluded; SPIDAAR
uniquely holds length-of-stay) forced a **disciplined, fully-logged pivot** to an
excess-bed-days estimand the data can actually support — and the pivot *is* the
contribution. Every deviation is in `docs/deviation_log.md` + an OSF addendum. The
exposure and adequacy codings are codebook-confirmed.

**Innovation (creativity).**
The novelty is the **reframing**: turning an underpowered, literature-contradicted
mortality null into a systemic-leverage thesis, and the **triple-mismatch** finding
— burden, funding, *and* surveillance all bending away from the same SSA community
Gram-negatives, with *S. pneumoniae* uniquely neglected on all three axes. We treat
the small-sample limitation as a **quantified result**, not a caveat to hide
(power, HDIs, E-values, inverse-ascertainment reweighting). Shipped as a
re-runnable, locally-parameterised stewardship simulator.

**Impact (scalability, implementation, policy).**
(1) A **transparency correction to AMR costing** — "excess bed-days × cost"
attributions that ignore competing mortality overstate resistance's bed-day burden.
(2) An **actionable tool** — a stewardship simulator taking local patient volume,
empiric adequacy, and WHO-CHOICE bed-day costs to project deaths averted / bed-days
added / cost. (3) A **resource-realignment call (Cross-Domain)** — across 16
pathogens surveillance tracks burden (ρ = +0.84, 95% CI 0.55–0.95, excludes zero)
but funding does not, and surveillance collapses geographically; steer therapeutics
funding and surveillance to the SSA community Gram-negatives, and raise empiric
adequacy where it pays most (+4.7 pp deaths averted in the resistant arm).

## 6. Links (open-access; permitted and encouraged)

- **Code (Apache-2.0, 155 tests, confirmatory harness):** github.com/Tej-Katika/amr-sentinel-vivli2026
- **Pre-registration:** OSF 10.17605/OSF.IO/BFQDP
- **Archived release (Zenodo concept DOI):** 10.5281/zenodo.21186887 — https://doi.org/10.5281/zenodo.21186887

## 7. One-line pitch (for a title/summary field or the live pitch open)

*In the highest-AMR-burden region on earth, the resistance phenotype isn't the
lever the headlines imply — and funding, surveillance, and burden all point away
from where the lever actually is. We quantify the paradox and ship the instruments
to reallocate.*

---

### Open items before submission
- [x] Mint Zenodo DOI → concept DOI 10.5281/zenodo.21186887 wired into §6 + report header + README.
- [ ] Verify the live MS Form's exact field list and map these blocks to it.
- [ ] Confirm the contact e-mail / any ORCID fields the form requests.
