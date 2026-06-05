# Deviation Log

Records any departure from the pre-registered analysis plan
(OSF: https://doi.org/10.17605/OSF.IO/BFQDP, deposited 2026-05-29, before data access).
Per pre-reg §11, every deviation — methodological or operational — is logged here
with the date, what changed, and why. This file is the audit trail referenced in
the final report's deviation note.

| Date | Type | Change | Rationale |
|------|------|--------|-----------|
| 2026-05-31 | Tooling | Added `xlrd>=2.0.1` to `pyproject.toml` dependencies. | SPIDAAR was delivered as legacy `.xls` files (patient, isolate, definitions); pandas requires the `xlrd` engine to read that format. No change to any analysis method, parameter, or model — purely enables reading the format the data shipped in. |
