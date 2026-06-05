# Deviation Log

Records any departure from the pre-registered analysis plan
(OSF: https://doi.org/10.17605/OSF.IO/BFQDP, deposited 2026-05-29, before data access).
Per pre-reg §11, every deviation — methodological or operational — is logged here
with the date, what changed, and why. This file is the audit trail referenced in
the final report's deviation note.

| Date | Type | Change | Rationale |
|------|------|--------|-----------|
| 2026-05-31 | Tooling | Added `xlrd>=2.0.1` to `pyproject.toml` dependencies. | SPIDAAR was delivered as legacy `.xls` files (patient, isolate, definitions); pandas requires the `xlrd` engine to read that format. No change to any analysis method, parameter, or model — purely enables reading the format the data shipped in. |
| 2026-06-04 | Methodological | Step 1 unit of analysis changed from isolate to **patient**; the exposure is the patient-level resistance summary `amrp` (≥1 resistant isolate vs. all susceptible). | SPIDAAR was delivered as two unlinked workbooks — `isolatedata` (244 isolates: organism, site, country, sample month/year, resistance) and `patientdata` (336 patients: 30-day mortality, age, sex, `amrp`). They share no join key (no `pid` in the isolate file, no `iid` in the patient file, no crosswalk sheet) and the counts do not reconcile for a positional join (244 isolate rows vs. `sum(nisol)=235` vs. 207 patients with ≥1 isolate). The mortality outcome and age/sex exist only patient-side, so the pre-registered isolate-level Cox model cannot be assembled. Confirmed with the investigator before implementing. Robustness check (final report): re-run Step 1 at isolate level if a re-linked SPIDAAR extract (`pid`↔`iid`) is later obtained. |
| 2026-06-04 | Methodological | Dropped `year` from the Step 1 propensity covariate set (`config.PROPENSITY_COVARIATES`); remaining covariates: age, sex, infection_site, country, organism. | Calendar year is not recorded in `patientdata` (only day-count fields: `dtpta`, `nobsd`, `enrtpt`, `los`); sample month/year exists only in the unlinkable `isolatedata`. A consequence of the patient-level unit above. The catchment window is a single enrolment period, so residual temporal confounding is expected to be small; noted as a limitation. |
