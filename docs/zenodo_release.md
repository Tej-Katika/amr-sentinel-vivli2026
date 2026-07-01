# Minting a Zenodo DOI for the code archive

This repo is **Zenodo-ready**: `.zenodo.json` (deposit metadata Zenodo reads on archival) and
`CITATION.cff` (GitHub "Cite this repository") are committed. Minting the DOI is a few clicks you
do yourself — it publishes an archive under your Zenodo identity, so it can't be automated here.

## Prerequisites
- The GitHub repo **must be public** (`https://github.com/Tej-Katika/amr-sentinel-vivli2026`).
  The challenge requires public outputs anyway. The Zenodo archive will contain the **code + docs
  only** — the licensed SPIDAAR/ATLAS data and the figures they generate are gitignored under the
  DUA and will **not** be included (correct and intended).
- The two metadata files are pushed to GitHub (see step 0).

## Steps
0. **Push the metadata** (and everything else) to GitHub:
   ```bash
   git push origin main
   ```
1. Go to **https://zenodo.org** → **Log in with GitHub** and authorize the Zenodo app
   (grant access to the repository).
2. Zenodo → **Account → GitHub** (https://zenodo.org/account/settings/github/). Find
   `Tej-Katika/amr-sentinel-vivli2026` in the list and flip its toggle **ON**. (If it isn't
   listed, click *Sync*; the repo must be public.) This tells Zenodo to archive future releases.
3. On GitHub, **create a release**: repo → *Releases* → *Draft a new release* → tag e.g.
   **`v1.0.0`**, title it, publish. (Only releases created *after* step 2 are archived.)
4. Zenodo automatically archives that tagged snapshot and **mints a DOI**, using `.zenodo.json`
   for the title/authors/license/keywords/related-DOI. It appears at
   https://zenodo.org/account/settings/github/ next to the repo within a minute or two.
5. Zenodo issues **two DOIs**: a *version* DOI (this release) and a *concept* DOI ("all versions").
   **Cite the concept DOI** in the report — it always resolves to the latest version.
6. Add the DOI to the submission:
   - Report header (`docs/final_report_2026.md`): replace the `Code:` line's Zenodo mention with
     the concept DOI.
   - `README.md`: add the Zenodo DOI badge (Zenodo shows the Markdown snippet on the deposit page).
   - Then rebuild the PDF: `scripts/build_report_pdf.sh`.

## Notes
- To update after minting: push changes, cut a new release (e.g. `v1.0.1`) — Zenodo archives it
  under the same concept DOI with a new version DOI.
- `.zenodo.json` already links the OSF pre-registration (`10.17605/OSF.IO/BFQDP`) as a related
  identifier, so the deposit records the pre-reg provenance.
