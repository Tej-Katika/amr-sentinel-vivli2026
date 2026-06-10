# Verified evidence base — Component 4 (Cross-Domain R&D alignment)

Pulled and verified 2026-06-09 (web research; verify-as-you-go discipline, cf.
`docs/reference_verified_2026-06-06.md` for MBIRA/Fiji/WHO-CHOICE). Every figure
below is quoted from a primary source; nothing here is recalled or interpolated.
Where exact per-pathogen magnitudes live only in a figure or appendix that could
not be cleanly extracted, that is stated explicitly and the value is left for the
secure-env / manual appendix pull — not guessed.

---

## 1. Funding denominator — LOCKED snapshot

**Czaplewski L, et al. "An overview of global public and philanthropic investments
into antibacterial therapeutics (2017–23)." The Lancet Microbe, 2026; epub 9 Jan
2026. doi:10.1016/S2666-5247(25)00216-2. PubMed 41525775.** A frozen, peer-reviewed
extract of the **Global AMR R&D Hub Dynamic Dashboard** (a GARDP × Hub analysis),
public + philanthropic funders, 2017–2023. Open-access PDF mirrored at
globalamrhub.org (`/wp-content/uploads/2026/01/Czaplewski-et-al-2026.pdf`).

This is `config.RD_HUB_SNAPSHOT_DATE = "2026-01-09"` / `RD_HUB_SOURCE`. It is
preferred over a hand-pulled live-dashboard snapshot because it is reproducible,
dated, and citable — and it already applies the public+philanthropic filter.

### Verified figures (US$ millions, public + philanthropic, 2017–2023)

| Quantity | Value | Verbatim source text |
|---|---|---|
| Total investment | **$2.51 billion** by **130 funders** | "A total of US$2·51 billion was invested in antibacterial R&D by 130 funders" |
| Peak / latest year | $445M (2020) → $363M (2023), −18% | "Funding peaked at $445 million in 2020 but declined by 18% to $363 million in 2023" |
| Species-specific total | **$1058 million** | "$142 million [13%] of $1058 million of species-specific funds" |
| *M. tuberculosis* | **$474M** (≈20% of overall) | "a fifth of the overall funds since 2017 ($474 million …)" |
| *S. aureus* | **$142M** (13% of species-specific), n=87 | "Staphylococcus aureus ($142 million [13%] of $1058 million of species-specific funds; n=87)" |
| *C. difficile* | **$141M** (13%), n=47 | "Clostridioides difficile ($141 million [13%], n=47)" |
| *N. gonorrhoeae* | **$101M** (10%), n=27 | "Neisseria gonorrhoea ($101 million, 10%, n=27)" |
| *P. aeruginosa* | **$87M** (8%), n=73 | "Pseudomonas aeruginosa ($87 million, 8%, n=73; figure 3B…)" |

### Derived (arithmetic on verified figures)
- **Pathogen-specific share = $1058M / $2510M = 42.2%** → **cross-cutting / non-species-
  specific = 57.8%** → trips `cross_cutting_share(...)["flagged"] = True`.
- Named top-five sum to $945M; the remaining **$113M** of species-specific funding is
  the `other_species_specific_musd` bucket (`RD_HUB_SNAPSHOT_2026`).

### Verified RANK order of species-specific funding (figure 3B, descending)
M. tuberculosis › S. aureus › C. difficile › N. gonorrhoeae › P. aeruginosa ›
**Escherichia coli › Acinetobacter baumannii › Klebsiella pneumoniae**.
→ The three highest-burden GRAM Gram-negatives sit at the **bottom** of the funding
ranking. **Streptococcus pneumoniae does not appear among species-specific targets**
(figure 3B) → effectively zero dedicated funding.

### NOT yet extracted (appendix-locked — do NOT guess)
Exact $ for E. coli, A. baumannii, K. pneumoniae (the small bars below P. aeruginosa,
axis maxes at $100M) are in **appendix 1, p18**. Pull these to run the full numeric
log2 mismatch index on the GRAM-6 panel.

---

## 2. Burden numerator — Murray 2022 (GRAM 2019)

**Murray CJL, et al. "Global burden of bacterial antimicrobial resistance in 2019: a
systematic analysis." The Lancet 2022;399:629–655. PMC8841637.** Verified from the
open-access full text (Europe PMC render PDF).

### Verified figures
| Quantity | Value | Verbatim source text |
|---|---|---|
| Six leading pathogens, combined | **929,000 (660,000–1,270,000) deaths attributable**; **3.57M (2.62–4.78M) associated** | "these six pathogens were responsible for 929 000 (95% UI 660 000–1 270 000) … attributable … and 3·57 million (2·62–4·78) … associated …" |
| Rank by **associated** deaths | E. coli › S. aureus › K. pneumoniae › S. pneumoniae › A. baumannii › P. aeruginosa | "E coli, Staphylococcus aureus, K pneumoniae, S pneumoniae, Acinetobacter baumannii, and Pseudomonas aeruginosa, by order of number of deaths" |
| Rank by **attributable** deaths | E. coli › K. pneumoniae › S. aureus › A. baumannii › S. pneumoniae › M. tuberculosis | "E coli … followed by K pneumoniae, S aureus, A baumannii, S pneumoniae, and M tuberculosis" |
| MRSA (pathogen–drug) | >100,000 attributable deaths, 3.5M DALYs | "meticillin-resistant S aureus was the one pathogen–drug combination … with more than 100 000 deaths and 3·5 million DALYs attributable to resistance" |
| E. coli & K. pneumoniae | each ≈200,000 attributable deaths | Oxford GRAM news: "two bacteria each caused close to 200,000 deaths in one year from AMR" |
| All-pathogen DALYs (23) | 47.9M | (documented in `GRAM_PANEL`) |

### NOT yet extracted (appendix-locked — do NOT guess)
Exact per-pathogen global **AMR-associated / attributable death counts and DALYs**
(+95% UI) live in **GRAM figure 4 / appendix**, not the main text. Pull these for the
numeric index Monte-Carlo (`monte_carlo_mismatch_ranking` wants (median, lo95, hi95)).

### ⚠ Conflation guard (caught during this pull)
The figures **S. aureus 1.1M / E. coli 950k / S. pneumoniae 829k / K. pneumoniae 790k
/ P. aeruginosa 559k** are from a **different paper** — Ikuta KS, et al. "Global
mortality associated with 33 bacterial pathogens in 2019" (Lancet 2022, the "1 in 8
deaths" study) — and are **all-cause deaths associated with the pathogen, NOT AMR-
attributable/associated**. Do not use them as the AMR burden numerator. (Same class
of citation error as the MBIRA mix-up; logged so it isn't repeated.)

---

## 3. What is wired vs. still gated

- **Wired now (verified, de-gated):** the cross-cutting headline
  (`cross_cutting_headline()` → 57.8% non-pathogen-specific, flagged) and the
  qualitative inverse-mismatch narrative (highest-burden Gram-negatives + S.
  pneumoniae are the least funded; pathogen-specific spend is dominated by TB and
  Gram-positive/niche targets). Snapshot date + source locked in `config`.
- **Still gated (machinery built + tested, awaiting exact magnitudes):** the numeric
  per-pathogen `mismatch_index` / `spearman_burden_funding` / `monte_carlo_mismatch_
  ranking` on the GRAM-6 panel. Needs: (a) Czaplewski appendix 1 p18 per-species $ for
  E. coli / A. baumannii / K. pneumoniae; (b) GRAM appendix per-pathogen deaths/DALYs
  + 95% UI. Both are public; neither was cleanly extractable from a figure via web
  fetch this session.

## Sources
- Czaplewski et al., Lancet Microbe 2026: https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(25)00216-2/fulltext · PDF: https://globalamrhub.org/wp-content/uploads/2026/01/Czaplewski-et-al-2026.pdf
- Murray et al., Lancet 2022 (GRAM): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8841637/
- GRAM project (Oxford): https://www.tropicalmedicine.ox.ac.uk/gram/news/global-burden-of-bacterial-antimicrobial-resistance
- Ikuta et al., Lancet 2022 (conflation guard): https://www.healthdata.org/news-events/newsroom/news-releases/lancet-one-eight-deaths-2019-linked-bacterial-infections-second
