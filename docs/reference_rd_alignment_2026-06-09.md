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

### Per-pathogen burden — EXTRACTED & VERIFIED (GRAM appendix 1, Table S22)
Pulled 2026-06-09 from the GRAM supplementary `mmc1.pdf` (Europe PMC supplementary bundle
for PMC8841637), Table S22 "Global deaths and DALYs … by pathogen-drug combination, 2019",
the **'Resistance to one or more antibiotics'** aggregate row per pathogen. Counts in
**thousands (median, 95% UI)**. Now shipped verbatim as `rd_alignment.GRAM_BURDEN_2019`.

| Pathogen | Assoc. deaths | Attrib. deaths | Assoc. DALYs | Attrib. DALYs |
|---|---|---|---|---|
| E. coli | 829 (601–1120) | 219 (152–316) | 28,000 (21,000–36,900) | 7,520 (5,270–10,500) |
| S. aureus | 748 (554–1000) | 178 (104–280) | 24,900 (18,600–32,700) | 5,870 (3,550–9,220) |
| K. pneumoniae | 642 (465–863) | 193 (130–272) | 27,400 (20,300–36,100) | 8,200 (5,550–11,400) |
| S. pneumoniae | 596 (490–727) | 122 (82.4–166) | 29,800 (24,400–36,700) | 6,110 (4,050–8,330) |
| A. baumannii | 423 (252–647) | 132 (75.7–213) | 11,800 (7,290–17,800) | 3,670 (2,150–5,760) |
| P. aeruginosa | 334 (234–457) | 84.6 (53–127) | 12,000 (8,630–16,100) | 3,050 (1,980–4,530) |

**Self-consistency (transcription check):** associated-death medians sum to **3,572k = 3.57M**
and attributable to **928.6k ≈ 929k** — both match the headline totals exactly; the
associated and attributable rank orders both match the main text. (Guarded by a unit test.)

### The ONE remaining un-fetched value (appendix-locked, propagated not guessed)
The $113M species-specific split among **E. coli / A. baumannii / K. pneumoniae** is in
Czaplewski appendix 1 p18 (Elsevier paywall; no PMC, no open supplement — three routes
tried and exhausted). It is bounded exactly (sum $113M, ordered E>A>K, each <$87M;
S. pneumoniae below the smallest shown bar) and **propagated as Monte-Carlo uncertainty**
in `gram_panel_alignment`, not hard-coded. The under-funded ranking is invariant to it.

### ⚠ Conflation guard (caught during this pull)
The figures **S. aureus 1.1M / E. coli 950k / S. pneumoniae 829k / K. pneumoniae 790k
/ P. aeruginosa 559k** are from a **different paper** — Ikuta KS, et al. "Global
mortality associated with 33 bacterial pathogens in 2019" (Lancet 2022, the "1 in 8
deaths" study) — and are **all-cause deaths associated with the pathogen, NOT AMR-
attributable/associated**. Do not use them as the AMR burden numerator. (Same class
of citation error as the MBIRA mix-up; logged so it isn't repeated.)

---

## 3. Closed-out index — computed results (`gram_panel_alignment`)

Monte-Carlo over the GRAM burden UIs and the one unfetched funding split, floor 0.02,
20,000 draws, seed `step_seed(4)`:

- **Cross-cutting: 57.8% non-pathogen-specific (flagged).**
- **Spearman ρ:** +0.14 (associated deaths) · +0.03 (attributable deaths) · −0.49
  (attributable DALYs); CI ≈ ±1 at n=6 → **no positive alignment** with burden.
- **Under-funded ranking (log2 mismatch median; robust to the split):**
  S. pneumoniae **+2.9** › K. pneumoniae **+2.0** › E. coli **+0.5** › A. baumannii **+0.1**
  › S. aureus **−1.0** › P. aeruginosa **−1.4**. S. pneumoniae + K. pneumoniae carry ≈99%
  of the "most under-funded" probability; S. aureus and P. aeruginosa are over-funded.

This is the descriptive Cross-Domain finding for report §3.4. Secure-env upgrade: drop in
the Czaplewski appendix p18 exact split (collapses the only remaining funding uncertainty;
ranking unchanged).

## Sources
- Czaplewski et al., Lancet Microbe 2026: https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(25)00216-2/fulltext · PDF: https://globalamrhub.org/wp-content/uploads/2026/01/Czaplewski-et-al-2026.pdf
- Murray et al., Lancet 2022 (GRAM): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8841637/
- GRAM project (Oxford): https://www.tropicalmedicine.ox.ac.uk/gram/news/global-burden-of-bacterial-antimicrobial-resistance
- Ikuta et al., Lancet 2022 (conflation guard): https://www.healthdata.org/news-events/newsroom/news-releases/lancet-one-eight-deaths-2019-linked-bacterial-infections-second
