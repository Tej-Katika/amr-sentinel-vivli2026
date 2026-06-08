# `data/` — Vivli DUA-restricted inputs (gitignored)

**Nothing in this directory except this file is committed.** Vivli-derived data
must not be redistributed under the Vivli Data Use Agreement.

At analysis time, the following live here (inside the Vivli secure environment):

- `spidaar/` — SPIDAAR cohort: patient file (length-of-stay + mortality outcomes;
  Ghana, Kenya, Malawi, Uganda) and the isolate file (per-mechanism resistance)
- `atlas/` — Pfizer ATLAS isolates (`atlas_vivli_2004_2024.csv`) — the **only**
  external surveillance source: **Merck SMART is excluded from the Challenge**, so
  there is no `smart/`.
- `rd_hub/` — locked Global AMR R&D Hub snapshot (public + philanthropic only),
  for the Cross-Domain R&D mismatch (Component 4). Record the snapshot date in
  `config.RD_HUB_SNAPSHOT_DATE` and never change it.

Country bed-day unit costs (WHO-CHOICE) and GRAM burden figures are external,
published constants documented/cited in `docs/reference_verified_2026-06-06.md` and
the analysis modules — not files staged here.

## What *is* released (per pre-reg §12)

- Aggregated outputs only: counts, rates, parameter estimates
- A synthetic-mirror dataset with matched marginals (no real patient records)

SHA-256 hashes of input files are recorded for reproducibility; the raw inputs
never leave the secure environment.
