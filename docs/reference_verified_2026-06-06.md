# Verified evidence base (2026-06-06)

Figures below were fetched from primary sources and adversarially re-verified
(every value re-checked against the source PDF/PMC page). Use these exact numbers
and citations in the report and live pitch.

## MBIRA — the authoritative directly-comparable SSA cohort

**Aiken AM, et al. "Mortality associated with third-generation cephalosporin
resistance in Enterobacterales bloodstream infections at eight sub-Saharan African
hospitals (MBIRA): a prospective cohort study." *Lancet Infect Dis*
2023;23(11):1280–1290.** doi:10.1016/S1473-3099(23)00233-5. PMID 37454672; PMC7617135.
Prospective matched parallel-cohort, 8 SSA hospitals, 878 BSI patients (221 3GC-S,
657 3GC-R) + 1634 matched uninfected controls; enrolled Nov 2020 – Jan 2022.

**Headline (the differential effect of resistance, after matching/adjustment for age,
HIV, indwelling devices, Charlson index, site-stratified):**
- In-hospital death, **ratio of cause-specific HRs (3GC-R vs 3GC-S) = 0.74 (95% CI 0.42–1.30) — NULL.**
- 30-day death, ratio of relative risks = 0.82 (0.53–1.27) — null.
- Hospital discharge, ratio of cause-specific HRs = 1.16 (0.93–1.45) — null (CI crosses 1).
- Median post-enrolment LOS **shorter** in 3GC-R (7 [IQR 3–14] vs 9 [4–18] days).

**Context (infection vs no infection):** any Enterobacterales BSI raised in-hospital
death ~5–7× vs matched uninfected controls (3GC-S CSH 6.79 [4.06–11.37]; 3GC-R CSH
5.01 [3.96–6.32]). Crude in-hospital death was higher in the resistant cohort
(37.1% vs 28.1%) — **but this dissolves after matching/adjustment.**

**Interpretation for our pitch:** resistance *per se* is **not** an independent
per-patient mortality/LOS driver in SSA hospitals; the crude difference is confounded
by infection severity and access to effective therapy. A "resistance kills" headline
would be directly contradicted by MBIRA. ⚠️ This **corrects** an earlier draft claim
in `strategy_2026.md` that "MBIRA found 3GC-R *was* associated with mortality."

## Fiji / Loftus — competing-risks excess-LOS precedent (citation corrected)

**Loftus MJ, et al. "Attributable mortality and excess length of stay associated with
third-generation cephalosporin-resistant Enterobacterales bloodstream infections: a
prospective cohort study in Suva, Fiji." *Journal of Global Antimicrobial Resistance*
2022;30:286–293.** doi:10.1016/j.jgar.2022.06.016. PMID 35738385; PMC9452645.
(Earlier drafts mislabelled the journal as *Lancet Reg Health W Pac* — the ScienceDirect
pii S2213716522001515 is correct; the journal name was wrong.)
- Adjusted in-hospital mortality aHR **1.13 (95% CI 0.51–2.53)** — null.
- Excess LOS attributable to 3GC-R (multistate model) **2.6 days (95% CI 2.5–2.8)**.
- Discharged-alive aHR **0.99 (95% CI 0.65–1.50)**; composite Cox 0.90 (0.63–1.30).

## WHO-CHOICE inpatient bed-day unit costs (hotel component)

Source: WHO-CHOICE "Cost per inpatient bed day by hospital level" dataset (2010 I$,
GBD-2010 aligned); methodology Stenberg et al. 2018, *Cost Eff Resour Alloc* 16:11.
The bed-day "hotel" cost **includes** personnel/capital/equipment/lab/food but
**excludes** drugs and diagnostic tests — the correct denominator for excess bed-days
(drug/diagnostic costs modelled separately, avoiding double-counting).

Base case = **primary-hospital** tier (most excess bed-days occur at district level).
2010 USD (market FX); PPP int-$ with 95% UI in brackets:

| Country | 2010 USD/bed-day | PPP I$ (model) | 95% UI (I$) |
|---|---|---|---|
| Ghana   | 6.30 | 14.28 | 5.94–32.50 |
| Kenya   | 5.45 | 14.03 | 5.54–30.87 |
| Malawi  | 3.25 | 6.53  | 2.58–14.13 |
| Uganda  | 3.81 | 10.25 | 4.08–23.83 |

Secondary/tertiary tiers run higher (e.g. Ghana tertiary US$8.49; Kenya tertiary
US$7.35) — use if the cohort is referral-level. Use the WHO 95% UI as the
Monte-Carlo sensitivity range. Inflate 2010 I$ → 2023 I$ by ~**1.34–1.36×** (US GDP
implicit price deflator; state the base year), or inflate 2010 NCU by each country's
GDP deflator and convert at the 2023 exchange rate for local-currency reporting.

Independent corroboration: Ghana micro-costing (Aboagye et al., *Ghana Med J* 2010;
PMC2996840) US$6.05 mission / 9.95 district / 18.80 referral (2003) — brackets the WHO
range. A Kenya tertiary all-inclusive figure (~US$93/day, Wauye et al., *Global Heart*
2025; PMC12047636) **includes drugs/labs** and is a high-dependency upper bound only.

**Minor caveats from verification:** WHO-CHOICE PDF page references in the source
notes were off-by-one (values correct); the exact BEA deflator index base year should
be stated when used.
