# `data/` — Vivli DUA-restricted inputs (gitignored)

**Nothing in this directory except this file is committed.** Vivli-derived data
must not be redistributed under the Vivli Data Use Agreement.

At analysis time, the following live here (inside the Vivli secure environment):

- `spidaar/` — SPIDAAR cohort (mortality outcomes; Ghana, Kenya, Malawi, Uganda)
- `atlas/` — Pfizer ATLAS isolates
- `smart/` — Merck SMART isolates
- `rd_hub/` — locked Global AMR R&D Hub snapshot (public + philanthropic only).
  Record the snapshot date in `config.RD_HUB_SNAPSHOT_DATE` and never change it.
- `who_ghe/` — WHO Global Health Estimates population denominators

## What *is* released (per pre-reg §12)

- Aggregated outputs only: counts, rates, parameter estimates
- A synthetic-mirror dataset with matched marginals (no real patient records)

SHA-256 hashes of input files are recorded for reproducibility; the raw inputs
never leave the secure environment.
