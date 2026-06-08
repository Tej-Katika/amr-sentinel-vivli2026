# AMR Sentinel — Vivli 2026 Analysis

[![OSF Pre-registered](https://img.shields.io/badge/OSF-Pre--registered-337AB7)](https://doi.org/10.17605/OSF.IO/BFQDP)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue)](./LICENSE)
[![Python 3.12](https://img.shields.io/badge/python-3.12-3776AB)](./.python-version)

Pre-registered secondary analysis for the **2026 Vivli AMR Surveillance Open Data Re-Use Data Challenge**. It anchors on **SPIDAAR** — the only register in the Vivli AMR catalogue with patient length-of-stay and mortality outcomes (Ghana, Kenya, Malawi, Uganda) — and uses **Pfizer ATLAS** surveillance and the **Global AMR R&D Hub** public-funding database to ask what the resistance phenotype actually does to patients, and where the leverage to act really lies.

## The thesis — the AMR burden paradox

After data access, the headline **pivoted off a resistance→mortality bridge** (underpowered and contradicted by the authoritative SSA literature) **onto resistance-attributable excess length-of-stay** estimated with competing risks. The triangulated finding (SPIDAAR + MBIRA + Fiji) is that **resistance is *not* a clean per-patient killer** — so the actionable burden is **systemic**: whether patients receive *adequate empiric therapy*, not the phenotype itself. The empiric-adequacy g-formula is the centerpiece; the resistance-to-mortality angle is retained only as an honest, underpowered secondary.

> The pivot is logged as a post-data-access deviation in [`docs/deviation_log.md`](docs/deviation_log.md) and a supplementary [`docs/osf_addendum_2026-06-06.md`](docs/osf_addendum_2026-06-06.md), filed under the parent pre-registration's §11. SMART is excluded from the Challenge, so all external surveillance is **ATLAS-only**. Build-ready specs: [`docs/analysis_plan_2026.md`](docs/analysis_plan_2026.md); strategy: [`docs/strategy_2026.md`](docs/strategy_2026.md).

## Components

| # | Module | What it does |
|---|---|---|
| 1 (primary) | [`excess_los.py`](src/amr_sentinel_vivli/excess_los.py) | Resistance-attributable **excess bed-days**: competing-risks (Aalen-Johansen) restricted-mean occupancy to τ=28, death competing with discharge; crude + severity-standardized + CIF decomposition; stratified bootstrap CI |
| 1b (co-primary) | [`excess_los_sensitivity.py`](src/amr_sentinel_vivli/excess_los_sensitivity.py) | Power/precision simulation (what CI width & power the n_R≈135 / n_S≈21 design delivers, and how re-linkage would tighten it) + exposure-ascertainment selection sensitivity (IPW reweight, assignment envelope, MAR anchor) |
| 1 (figure) | [`excess_los_figures.py`](src/amr_sentinel_vivli/excess_los_figures.py) | State-occupation bed-occupancy curves + competing-outcome panel |
| 3 | [`bayesian_projection.py`](src/amr_sentinel_vivli/bayesian_projection.py) | **ATLAS catchment nowcast** (empirical-Bayes beta-binomial partial pooling, MC credible intervals) + data-availability matrix + borrowed regional trend + **SPIDAAR frame-contrast** (severe-HAI 3GC-R vs mixed-surveillance ceftazidime-R); figures in [`projection_figures.py`](src/amr_sentinel_vivli/projection_figures.py) |
| 4 (Cross-Domain) | [`rd_alignment.py`](src/amr_sentinel_vivli/rd_alignment.py) | Descriptive **R&D mismatch index** (global GRAM burden share vs Global AMR R&D Hub public+philanthropic funding share, log2 + Spearman ρ; cross-cutting magnitude first; low-denominator floor) |
| 5 (centerpiece) | [`stewardship_gformula.py`](src/amr_sentinel_vivli/stewardship_gformula.py) + [`app/`](app/) | Empiric-adequacy **g-formula**: counterfactual avertable bed-days/deaths, WHO-CHOICE costed, with a shipped **Streamlit what-if tool** (de-identified calibration artifact + synthetic demo) |

Cross-cutting: [`config.py`](src/amr_sentinel_vivli/config.py) (locked constants), [`data_loading.py`](src/amr_sentinel_vivli/data_loading.py) (SPIDAAR + ATLAS loaders), [`pipeline.py`](src/amr_sentinel_vivli/pipeline.py) (orchestration). `cox_mortality.py` and `mortality_bridge.py` are **superseded** by the excess-LOS pivot (mortality demoted to secondary) and retained for provenance.

## The Streamlit stewardship tool

```bash
streamlit run app/streamlit_app.py
```

Loads a de-identified calibration artifact (population-aggregate constants only — no patient records) and lets a catchment-region programme scale them to local inputs (patients treated, current vs target empiric adequacy, bed-day cost) → projected deaths averted, bed-days added, and cost. Ships a synthetic demo artifact so it runs publicly. See [`app/README.md`](app/README.md). Every output is labelled an *ecological-calibration scenario, not an individual effect*.

## Layout

```
src/amr_sentinel_vivli/   analysis modules
app/                      Streamlit stewardship what-if tool + synthetic demo artifact
tests/                    unit + config-lock tests (synthetic data only)
data/                     Vivli DUA-restricted inputs — gitignored, NEVER committed
figures/  results/        regenerated outputs — gitignored
docs/                     pre-registration, deviation log, OSF addendum, plans
notebooks/                exploratory notebooks (flagged exploratory per pre-reg §8.5)
```

## Reproducibility (pre-reg §12)

- Dependencies pinned in `pyproject.toml`; analysis runs inside the Vivli secure environment.
- Master seed `20260526`; per-step seeds via `numpy.random.SeedSequence` (`config.step_seed`).
- Estimators are hand-rolled in numpy/pandas (competing-risks RMST/CIF, empirical-Bayes pooling, g-formula, Monte-Carlo costing) so they run without lifelines/scipy/pymc; tests use synthetic frames only.
- Global AMR R&D Hub data is pulled once at a **locked snapshot date** (`config.RD_HUB_SNAPSHOT_DATE`, cited in every figure caption), as the Hub revises retrospectively.

## Data access

Vivli-derived data is **not** redistributed (Vivli Data Use Agreement). Inputs live under `data/` (gitignored); derived outputs in `results/`/`figures/` are also gitignored and regenerated from code. Only aggregated outputs (counts, rates, parameter estimates) leave the secure environment, subject to egress review. See [`data/README.md`](data/README.md).

## Positionality

Neither author is based in the SPIDAAR catchment region; the framework is built as parameterized infrastructure for catchment-region researchers to re-run with local parameters, not as a verdict for any country. No stewardship recommendation is final without local validation. Full statement in the pre-registration §16.

## Authors

- **Tejashwar Reddy Katika** — Independent Researcher (University of North Texas); Lead Applicant — ORCID [0009-0006-3015-5697](https://orcid.org/0009-0006-3015-5697)
- **Akhilesh Reddy Katika** — MS in Data Science, Flinders University — ORCID [0009-0009-4517-6126](https://orcid.org/0009-0009-4517-6126)

## License

[Apache License 2.0](./LICENSE).
