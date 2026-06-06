# SPIDAAR re-linkage data request (draft)

**To:** Vivli AMR Data Sharing team / SPIDAAR data contributor
**From:** AMR Sentinel team (Vivli 2026 AMR Open Data Challenge)
**Re:** Patient↔isolate crosswalk for the SPIDAAR cohort (dataset as delivered in the secure environment)
**Pre-registration:** OSF [10.17605/OSF.IO/BFQDP](https://doi.org/10.17605/OSF.IO/BFQDP) (deposited 2026-05-29, before data access)

## Summary of the request

We request a re-linked SPIDAAR extract that allows each **isolate** record to be
joined to its **patient** record. The minimal acceptable form is **either**:

1. a crosswalk table with one row per isolate carrying both the isolate id
   (`iid`) and the patient id (`pid`); **or**
2. the existing `spidaar_isolatedata` workbook with a `pid` column added; **or**
3. the existing `spidaar_patientdata` workbook with the constituent `iid`(s)
   added per patient.

No additional variables beyond the linking key are required.

## Why this is needed

Our pre-registered primary analysis (Step 1) is an **isolate-level** survival
model: 30-day mortality contrasting resistant vs. susceptible isolates *of the
same species*, with patient covariates (age, sex, infection site, country) and
calendar time. Assembling that model requires the resistance phenotype (isolate
file) and the mortality outcome + demographics (patient file) on the **same
unit**.

As delivered, the two SPIDAAR workbooks cannot be joined:

- `spidaar_isolatedata.xls` — 244 isolates: organism, specimen, HAI site,
  clinical relevance, sample month/year, and resistance (`c3r`, `mdr`, `mrsa`,
  composite `rx`). Keyed by `iid`. **No `pid`.**
- `spidaar_patientdata.xls` — 336 patients: 30-day mortality, age class, sex,
  patient-level resistance summary (`amrp`), isolate count (`nisol`). Keyed by
  `pid`. **No `iid`.**
- There is **no crosswalk sheet**, and the counts do not reconcile for a
  positional/row-order join: 244 isolate rows vs. `sum(nisol) = 235` vs. 207
  patients with ≥1 recorded isolate. Any positional merge would mis-assign
  outcomes to phenotypes and is therefore not defensible.

The mortality outcome and demographics exist **only** patient-side; the
per-isolate resistance phenotype exists **only** isolate-side. Without the link
we cannot run the registered isolate-level model.

## Impact of running without the link (our current deviation)

We have deviated to a **patient-level** Step 1 (logged in
`docs/deviation_log.md`, 2026-06-04), using the pre-derived patient resistance
summary `amrp`. This collapse produces a severely imbalanced and underpowered
comparison — on the order of **~135 resistant vs. ~23 susceptible patients** and
a single-digit event count in the susceptible arm — which the registered
isolate-level design (244 isolates, real phenotype variation, e.g. 3GC-R
ascertained in ~200 isolates) would substantially improve. Restoring the link is
the only path that returns us to the **pre-registered confirmatory analysis**
without a methodological deviation.

## What we will do with it

- Re-run Step 1 at the registered isolate level (resistant vs. susceptible,
  same-species contrast) as the confirmatory analysis; the patient-level result
  becomes a documented robustness check rather than the headline.
- No re-identification is attempted; the link is used only to associate an
  isolate's phenotype with its patient's outcome inside the secure environment.
  All outputs remain aggregate and DUA-compliant.

## Governance

This request seeks **no new variables** and changes none of the locked
pre-registered parameters; it only restores the join the registered design
assumes. If the link cannot be provided, we will retain the patient-level
deviation and report it transparently with the power limitation stated.

---
*Draft prepared for the team to review and submit through the Vivli data-request
channel. Fill in the contact/ticket reference before sending.*
