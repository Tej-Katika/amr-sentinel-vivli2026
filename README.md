# AMR Sentinel — Vivli 2026 Analysis

[![OSF Pre-registered](https://img.shields.io/badge/OSF-Pre--registered-337AB7)](https://doi.org/10.17605/OSF.IO/BFQDP)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue)](./LICENSE)
[![Python 3.12](https://img.shields.io/badge/python-3.12-3776AB)](./.python-version)

Pre-registered secondary analysis for the **2026 Vivli AMR Surveillance Open Data Re-Use Data Challenge**. The analysis links the **SPIDAAR** registry — the only register in the Vivli AMR catalogue with patient mortality outcomes (Ghana, Kenya, Malawi, Uganda) — with **ATLAS** and **SMART** resistance trajectories and the **Global AMR R&D Hub** public funding database, to estimate AMR-attributable mortality and characterize the alignment of public R&D investment with burden.

> ### ⚠️ Pre-data scaffold
> This repository is committed **before any Vivli data has been accessed**. The analysis plan is locked in the OSF pre-registration ([10.17605/OSF.IO/BFQDP](https://doi.org/10.17605/OSF.IO/BFQDP)), deposited **2026-05-29**, before the Vivli Data Use Agreement was signed or the secure environment entered. The module functions are intentionally stubs; implementations land only after data access, and any deviation from the pre-registered plan is logged per its §11.

## The five pre-registered steps

| Step | Module | What it does (pre-reg §7) |
|---|---|---|
| 1 | [`cox_mortality.py`](src/amr_sentinel_vivli/cox_mortality.py) | SPIDAAR Cox PH mortality model, IPTW, e-value (tests **H1**) |
| 2 | [`bayesian_projection.py`](src/amr_sentinel_vivli/bayesian_projection.py) | Bayesian hierarchical resistance projection 2025–2030 (tests **H2**) |
| 3 | [`mortality_bridge.py`](src/amr_sentinel_vivli/mortality_bridge.py) | PAF bridge: HR × projected resistance → attributable deaths |
| 4 | [`rd_alignment.py`](src/amr_sentinel_vivli/rd_alignment.py) | Public R&D landscape alignment vs burden (tests **H3**) |
| 5 | [`stewardship_gformula.py`](src/amr_sentinel_vivli/stewardship_gformula.py) | Empiric-adequacy g-formula: counterfactual avertable bed-days/deaths, WHO-CHOICE costed (exploratory) |

Cross-cutting: [`config.py`](src/amr_sentinel_vivli/config.py) (locked constants), [`sensitivity.py`](src/amr_sentinel_vivli/sensitivity.py) (pre-reg §8), [`pipeline.py`](src/amr_sentinel_vivli/pipeline.py) (orchestration).

## Layout

```
src/amr_sentinel_vivli/   analysis modules (one per pre-registered step)
tests/                    smoke + config-lock tests
data/                     Vivli DUA-restricted inputs — gitignored, NEVER committed
figures/  results/        regenerated outputs — gitignored
notebooks/                exploratory notebooks (flagged exploratory per pre-reg §8.5)
docs/                     pre-registration pointer + methods notes
```

## Reproducibility (pre-reg §12)

- Python 3.12, dependencies pinned in `pyproject.toml` (+ `uv.lock` once resolved)
- Master seed `20260526`; per-step seeds via `numpy.random.SeedSequence`
- DVC tracking on derived artifacts; MLflow run logging per fitted model
- Every figure/table in the final report cites the commit hash that produced it
- Global AMR R&D Hub data is pulled once at a **locked snapshot date** (cited in every figure caption), as the Hub revises retrospectively

## Data access

Vivli-derived data is **not** redistributed (Vivli Data Use Agreement). Inputs live under `data/` (gitignored). Only aggregated outputs (counts, rates, parameter estimates) and a synthetic-mirror dataset with matched marginals are released. See [`data/README.md`](data/README.md).

## Positionality

Neither author is based in the SPIDAAR catchment region; the framework is built as parameterized infrastructure for catchment-region researchers to re-run with local parameters, not as a verdict for any country. Full statement in the pre-registration §16.

## Authors

- **Tejashwar Reddy Katika** — Independent Researcher (University of North Texas); Lead Applicant — ORCID [0009-0006-3015-5697](https://orcid.org/0009-0006-3015-5697)
- **Akhilesh Reddy Katika** — MS in Data Science, Flinders University — ORCID [0009-0009-4517-6126](https://orcid.org/0009-0009-4517-6126)

## License

[Apache License 2.0](./LICENSE).
