# Secure-environment confirmatory-run runbook

Step-by-step for reproducing and egressing the AMR Sentinel headline numbers inside the
Vivli AMR secure environment. This is the **Step 6** referenced in the project plan.

## What this does

`scripts/confirmatory_run.py` runs the full `pipeline.run()` on the real delivered files,
flattens every result, and diffs each value against `scripts/headline_manifest.json` (the
report's headline numbers, transcribed with tolerances). It prints a **PASS / DRIFT** table
and exits **0** when everything reproduces, **1** on any *unexpected* drift. Because the
delivered data already produced the development numbers, this is a *reproduce-and-egress*
check, not a re-analysis.

## Prerequisites

1. **Code** present in the secure env (repo clone/copy or mounted).
2. **Python 3.10+** (3.13 is the tested interpreter).
3. **Packages — only these are needed at run time:** `numpy`, `pandas`, `xlrd`
   (+ `matplotlib` for figures). The heavy libs declared in `pyproject.toml`
   (scipy/statsmodels/scikit-learn/lifelines/pymc/arviz/preliz) are **not imported at
   run time** — every estimator is hand-rolled in NumPy precisely because the secure env
   lacks them. `pip install -e .` also works if the env allows it.
4. **The three delivered data files**, at these exact repo-relative paths:
   - `data/spidaar/spidaar_patientdata.xls`
   - `data/spidaar/spidaar_isolatedata.xls`
   - `data/atlas/atlas_vivli_2004_2024.csv`

   No R&D Hub file is required: the funding denominators are in code
   (`rd_alignment.RD_HUB_SNAPSHOT_2026` primary, `RD_HUB_GENUS_SNAPSHOT` robustness), so the
   run is self-contained on those three files.

## Steps

**1 — Verify the data is in place**

```bash
ls -la data/spidaar/spidaar_patientdata.xls \
       data/spidaar/spidaar_isolatedata.xls \
       data/atlas/atlas_vivli_2004_2024.csv
```

**2 — Set up Python**

```bash
python -m pip install numpy pandas xlrd matplotlib   # matplotlib only for figures
```

**3 — (Sanity) run the test suite** — should be green before trusting the numbers:

```bash
PYTHONPATH=src python -m pytest -q          # expect: 145 passed
```

**4 — Run the confirmatory harness**

```bash
PYTHONPATH=src python scripts/confirmatory_run.py
```

Prints the PASS/DRIFT table and writes `confirmatory_results.json`.

**5 — Read the result**

- **All `PASS`, or only `DRIFT*` (flagged) rows** → success (exit 0); the report's headline
  numbers reproduce on the real files.
- **Any unflagged `DRIFT`** (exit 1) → a genuine discrepancy. For that row decide: stale
  report figure (fix the report) vs. real change (fix/explain the code). Each row prints its
  report §-ref to localise it.

**6 — Reconcile the *expected* `DRIFT*` rows by hand** (they never fail the run; still check):

- **`bayesian_excess_los.*`** (mean −0.94, HDI, P(excess>0)=0.30, P(>1d)=0.14) — flagged
  because the local Bayesian step is a Laplace approximation. If you do not upgrade to NUTS
  they should simply `PASS`; if you do run NUTS, confirm they land within tolerance
  (±0.25 on the mean).
- **`rd_alignment.*` / `rd_alignment_catchment.*` Gram-negative rows** — flagged for the
  unfetched $113M Gram-negative funding split (Monte-Carlo'd locally). Largely moot now that
  the live in-scope Hub extract is wired into the genus robustness check, but the *primary*
  index still uses the MC split; dropping in the real split shifts these within tolerance.

**7 — Regenerate the egress-reviewed figures** (the 4 PNGs the report embeds, into `figures/`):

```bash
PYTHONPATH=src python scripts/make_figures.py
```

Writes `figures/excess_los_stateoccupation.png`, `evidence_forest.png`,
`spidaar_framecontrast.png`, `rd_mismatch_global_vs_catchment.png` (needs `matplotlib`).

**8 — Egress** — submit `confirmatory_results.json` and the four PNGs through Vivli egress
review. These are aggregate outputs only; no row-level data leaves the environment.

**9 — Finalize** — if step 5 is clean (PASS / DRIFT*-only), the report's "Final Submission"
status is validated. If an unflagged DRIFT appeared, correct it before submitting.

## Notes

- Run from the repo root with `PYTHONPATH=src` (src-layout package `amr_sentinel_vivli`).
- All randomness is seeded from `config.MASTER_SEED` (`20260526`); reruns are deterministic.
- The development numbers were produced on this same delivered data locally, so a clean
  confirmatory run is the expected outcome — the value is the auditable, egressed record.
