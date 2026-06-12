# References

*Companion reference page for `final_report_2026.md` (Vivli 2026 AMR Surveillance Open Data Re-Use Data Challenge). Every figure quoted in the report is sourced below; full verbatim extracts and verification notes are in `docs/reference_verified_2026-06-06.md` and `docs/reference_rd_alignment_2026-06-09.md`.*

1. **Murray CJL, Ikuta KS, Sharara F, et al.** Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. *The Lancet* 2022;399(10325):629–655. doi:10.1016/S0140-6736(21)02724-0. PMC8841637. *(GRAM 2019 burden numerator: 1.27 M attributable / 4.95 M associated AMR deaths; per-pathogen Table S22.)*

2. **Aiken AM, Barffour R, Akoto AO, et al.** Mortality associated with third-generation cephalosporin resistance in Enterobacterales bloodstream infections at eight sub-Saharan African hospitals (MBIRA): a prospective cohort study. *The Lancet Infectious Diseases* 2023;23(11):1280–1290. doi:10.1016/S1473-3099(23)00233-5. PMID 37454672; PMC7617135. *(Authoritative directly-comparable SSA cohort; resistance→in-hospital-death ratio of cause-specific HRs 0.74 [0.42–1.30], null after adjustment.)*

3. **Loftus MJ, Naidu R, Cheng AC, et al.** Attributable mortality and excess length of stay associated with third-generation cephalosporin-resistant Enterobacterales bloodstream infections: a prospective cohort study in Suva, Fiji. *Journal of Global Antimicrobial Resistance* 2022;30:286–293. doi:10.1016/j.jgar.2022.06.016. PMID 35738385; PMC9452645. *(Multistate excess-LOS precedent: 2.6 days [2.5–2.8]; adjusted mortality aHR 1.13 [0.51–2.53], null.)*

4. **VanderWeele TJ, Ding P.** Sensitivity analysis in observational research: introducing the E-value. *Annals of Internal Medicine* 2017;167(4):268–274. doi:10.7326/M16-2607. PMID 28693043. *(E-value for the g-formula adequacy contrast, §3.3.)*

5. **Stenberg K, Lauer JA, Gkountouras G, Fitzpatrick C, Stanciole A.** Econometric estimation of WHO-CHOICE country-specific costs for inpatient and outpatient health service delivery. *Cost Effectiveness and Resource Allocation* 2018;16:11. doi:10.1186/s12962-018-0095-x. PMC5907387. *(WHO-CHOICE inpatient bed-day unit costs used in the Streamlit stewardship tool: Ghana \$6.30, Kenya \$5.45, Malawi \$3.25, Uganda \$3.81 per bed-day, 2010 USD.)*

6. **Czaplewski L, Hill T, Lienhardt C, et al.** An overview of global public and philanthropic investments into antibacterial therapeutics (2017–23). *The Lancet Microbe* 2026; published online 9 Jan 2026. doi:10.1016/S2666-5247(25)00216-2. PubMed 41525775. *(Frozen, peer-reviewed extract of the Global AMR R&D Hub Dynamic Dashboard; funding denominator for the Cross-Domain index — \$2.51 bn / 130 funders / \$1.06 bn species-specific.)*

---

*Data sources.* **SPIDAAR** (Surveillance Partnership to Improve Data for Action on Antimicrobial Resistance), accessed via the Vivli AMR Surveillance Open Data platform. **Pfizer ATLAS** (Antimicrobial Testing Leadership and Surveillance), `atlas_vivli_2004_2024`. **Global AMR R&D Hub** Dynamic Dashboard (as extracted by Czaplewski et al. 2026, ref. 6). Pre-registration: OSF doi:10.17605/OSF.IO/BFQDP.

*Methods note.* Competing-risks estimands use the Aalen-Johansen cumulative-incidence estimator and restricted-mean (RMST) state-occupation contrasts; all estimators are implemented directly in NumPy/pandas (the secure environment lacks specialised survival/Bayesian libraries) and unit-tested against synthetic data.
