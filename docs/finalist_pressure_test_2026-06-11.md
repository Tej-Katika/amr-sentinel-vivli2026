# Finalist-probability pressure-test (2026-06-11)

*An adversarially-verified reality-check on our odds of placing in the Vivli 2026 AMR
Surveillance Open Data Re-Use Data Challenge, validated against the actual 3-year history
of the challenge. Sources are primary Vivli pages unless noted; every base-rate figure
below survived 3-vote adversarial verification (deep-research run, 98 agents, 15 sources,
25 claims verified → 21 confirmed / 4 refuted).*

## Bottom line (and how it revises our prior)

Our standing memory said **"~50–60% finalist."** The verified history does **not** support a
number that high, for two specific reasons the pressure-test exposed:

1. **There is no published "finalist" denominator.** Vivli discloses *winners* (~6–7/yr) but
   **never publishes how many teams are shortlisted for the pitch round, nor any selection
   rate.** Two attempts to source a finalist count ("7 finalists in 2025"; "~a dozen in 2023")
   were **refuted 0–3** under verification. So any finalist percentage is an *inference on an
   unknown denominator*, not a fact. The old 18%→50% chain was built on that missing number.
2. **The "our archetype is favoured" multiplier is softer than we assumed.** Causal-inference /
   counterfactual methods *have won* (2025 grand prizes), but the stronger claim that **judges
   systematically favour causal inference was refuted 0–3.** Winning ≠ preferred. We should not
   double our odds on archetype-matching.

**Revised, defensible estimate:**

| Outcome | Verified base rate | Our quality-adjusted estimate | Confidence |
|---|---|---|---|
| Win **any** of the 6–7 awards | **~11%** (6.3/56 avg) | **~20–25%** | medium |
| Win the **Cross-Domain Award** specifically (our flagship path) | n/a (new 2026) | **~15–25%** | low-medium |
| Reach the **pitch/finalist shortlist** | **unpublished** | ~30–40% *if* shortlist ≈ 12–18 teams | low (denominator unverified) |
| Win a **Grand Prize** (Leadership/Visionary) | ~2/56 ≈ 3.5% | **~4–8%** | medium |

The honest headline: **~20–25% to place (win something), driven mostly by the Cross-Domain
niche; "finalist shortlist" is probably 30–40% but rests on an assumption Vivli won't confirm.**
This is a **downward revision** from 50–60%, and the revision is the *point* of the pressure-test.

## The verified reference class (2023–2025)

| Edition | Teams | Countries | Awards given | Winner/recognition rate |
|---|---|---|---|---|
| 2023 (1st annual) | 56 | 28 | 6 (Grand + 2 Impact + 2 Innovation + 1 Runner-up) | 10.7% |
| 2024 (2nd) | 55 | 27 | 6 + 1 Honorable Mention | 10.9–12.7% |
| 2025 (3rd) | 58 | 22 | 7 (2 Grand + Student×2 + Impact + Innovation + 1 HM) | 12.1% |
| **2026 (4th, current)** | TBD (~56 expected) | — | 6 named awards, >$50k total | — |

Participation is **remarkably flat (~56 teams)**, so the denominator is stable and predictable.
"Teams participated" is **not** cleanly split into EOI-only vs. full-submission counts — the
true scored-submission denominator could be smaller than 56, which would *raise* every rate
above. (Undisclosed; flagged as an open question.) No verified data exists before 2023.

## The judging reality

- **Process:** 300-word EOI gates data access → **5-page written submission is what gets scored**
  → an unstated number of **finalists are shortlisted for a live Zoom pitch + judge Q&A** → winners.
- **Criteria (published, unweighted):** **Methodology · Design · Innovation · Impact.** No public
  weighting; the panel may weight privately.
- **Panel:** independent, edition-specific (4 judges in 2023; ~8–9 in recent years), drawn from
  Pfizer / Paratek / academia, and for 2026 **includes a Global AMR R&D Hub expert (Dr Lesley
  Ogilvie)** — relevant to our Cross-Domain play.
- **What wins:** methodological novelty + policy relevance + LMIC focus recur. 2023 Grand Prize
  was an explicit **cross-dataset harmonization of six industry systems** vs. WHO GLASS; 2025
  grand prizes were causal-inference / counterfactual-policy work. LMIC teams (Kenya, Nigeria,
  Ghana, India) feature repeatedly — **independent / non-elite-institution teams do win.**

## Where WE sit against that reference class

**Above the median (push odds up):**
- **Reproducibility & rigor → Methodology criterion.** Pre-registered (OSF DOI), 135 unit tests,
  competing-risks estimators hand-built and tested. Few entries will be this disciplined.
- **Genuine multi-dataset integration → Innovation + Cross-Domain eligibility.** SPIDAAR + ATLAS +
  Global AMR R&D Hub. This *is* the 2023-Grand-Prize archetype (harmonize sources) and it
  satisfies the Cross-Domain requirement (Hub + ≥1 Register dataset) that most teams won't meet.
- **LMIC / sub-Saharan-Africa patient-outcome focus** — squarely in the winning pattern.
- **Novel reframing** (competing-risks excess-bed-days; the "burden paradox") → Innovation.
- **Policy-actionable deliverable** (stewardship what-if tool; R&D realignment) → Impact.

**Below the median (pull odds down — be honest):**
- **The headline is a null.** Resistance doesn't prolong stay; the mortality contrast is null.
  A judge skimming 56 submissions may read "null" as "underpowered," not "disciplined." **This is
  our single biggest liability** and no amount of method polish fully removes it.
- **n=21 susceptible arm.** Tiny. Direct hit on Methodology/Design even though we handle it openly.
- **Broken patient↔isolate linkage** → patient-level, not within-species. Structural, unfixable
  pre-secure-env.
- **Small independent team** vs. Oxford / St. Jude–Johns Hopkins / Vellore institutional entries.
  Not disqualifying (Kenya/Ibadan won) but the heavy hitters bring weight on Impact perception.

## The Cross-Domain Award: our highest-EV path (and its real size)

This is the one place the prior optimism *survives* scrutiny — with caveats:

- **Confirmed NEW for 2026.** $5,000 **travel grant, no cash**. Requires the **Global AMR R&D Hub
  dataset + ≥1 AMR Register dataset**; students and non-students judged **together**; **can be won
  alongside another award** (it's an extra shot, not mutually exclusive).
- **Why it favours us:** the eligibility filter (combine a *funding/R&D* dataset with *surveillance*)
  is non-obvious and new, so the **eligible pool is probably small** — most teams use a single
  surveillance register and never touch the Hub funding data. Our **Component 4 is purpose-built for
  exactly this** and a Hub expert sits on the panel.
- **Why not to over-count it:** (a) the eligible-pool size is *unverified* — if 15+ teams target it,
  base odds are ~7%, not high; (b) it's a $5k travel grant, **not** a grand prize; (c) our index is
  honestly *descriptive* (ρ≈0 at n=6) — rigorous but not flashy, which cuts against Innovation-seeking
  judges. Net: a realistic **15–25%** for this specific award — our best single line, but not a lock.

## What this implies for the next two weeks

1. **Lean into the Cross-Domain framing.** It is the highest expected-value, lowest-competition path
   and the one our prior optimism survives. Make the Hub + Register integration unmistakable in the
   5-page submission and ensure §3.4 reads as the *centerpiece-adjacent* contribution, not an appendix.
2. **Defuse the null, don't bury it.** The report already frames it as a "burden paradox" / systemic-
   leverage thesis — keep sharpening that so a skimming judge reads *discipline*, not *underpowered*.
   This is the highest-leverage editorial move on the general-award path.
3. **Bank the Methodology points.** Pre-registration + 135 tests + reproducible pipeline is a genuine
   differentiator; make it legible on page 1, not buried in §6.
4. **Don't bet the case on archetype-matching.** "We look like past winners" is weaker than assumed;
   sell the *substance* (rigor, multi-dataset, LMIC, actionable tool), not the resemblance.

## Caveats on this pressure-test itself

- The finalist (pitch-invited) **denominator is genuinely unpublished**; the 30–40% shortlist figure
  is an assumption (shortlist ≈ 12–18 of 56), not a sourced rate. If the shortlist is closer to the
  ~7 award winners, that figure drops toward 15–20%.
- Quality adjustments above are **judgement, not data** — there is no published rubric weighting to
  calibrate against, and no finalist-vs-non-finalist feature data exists.
- 2026 competition is unobserved; the new Cross-Domain pool size is the biggest single unknown.

*Full verified findings, votes, and source URLs: see the deep-research output for this run
(15 primary/secondary sources, 21 confirmed claims). Key primary pages: amr.vivli.org
data-challenge overview, FAQs, how-to-participate, and the 2023/2024/2025 winner pages.*
