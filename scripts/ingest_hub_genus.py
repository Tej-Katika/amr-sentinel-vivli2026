"""Ingest a Global AMR R&D Hub "by Genus" Data Table export into the per-pathogen denominator.

The Hub Dynamic Dashboard's *project-level* "Export" button leaves the InfectiousAgent
column blank, so it cannot give per-pathogen funding. The **"See Data Table" -> Excel**
export of the *by Genus* panel does carry it. This script parses that export, maps each
genus onto the six GRAM panel species, writes a tidy CSV under ``data/rd_hub/``, and prints
a scope check + the resulting Component-4 mismatch (via
``rd_alignment.genus_robustness_alignment``) so you can see immediately whether the export
is correctly scoped.

    PYTHONPATH=src python scripts/ingest_hub_genus.py [path/to/by-Genus.xlsx] [extract_date]

If no path is given it looks for the newest ``*by Genus*.xlsx`` in the user's Downloads
``data`` folder. ``extract_date`` (YYYY-MM-DD) defaults to today-less form read from the
export's own "Last updated" filter text when present, else must be passed.

SCOPE: for the index to be comparable to the locked Czaplewski denominator (~$2.51bn),
the dashboard filters must be set BEFORE exporting: Type of Funder = Public-Government +
Public-Other + Private-Non Profit; Research Area = Therapeutics; Year = 2017-2023; Human.
The script prints the export's recorded "Applied filters" and flags an out-of-scope pull
(e.g. all funders / all years -> total nearer $17.9bn) rather than silently using it.
"""

from __future__ import annotations

import glob
import sys
from pathlib import Path

import pandas as pd

from amr_sentinel_vivli import config, rd_alignment

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

_DOWNLOADS_DATA = Path.home() / "Downloads" / "data"

# Genus label (lowercased substring) -> GRAM panel species. Genus-level aggregation for all
# six; Escherichia/Klebsiella/Acinetobacter/Pseudomonas ~= the panel species, while
# Staphylococcus/Streptococcus genera are broader than S. aureus / S. pneumoniae (documented
# caveat, mirrored in rd_alignment.RD_HUB_GENUS_SNAPSHOT).
_GENUS_TO_PANEL = {
    "escherichia": "Escherichia coli",
    "klebsiella": "Klebsiella pneumoniae",
    "staphylococcus": "Staphylococcus aureus",
    "streptococcus": "Streptococcus pneumoniae",
    "acinetobacter": "Acinetobacter baumannii",
    "pseudomonas": "Pseudomonas aeruginosa",
}

# A correctly-scoped public+philanthropic therapeutics 2017-2023 pull should total near this
# (Czaplewski locked snapshot); a full all-funder/all-year pull lands near the latter.
_SCOPE_TARGET_MUSD = 2510.0
_FULL_SCOPE_MUSD = 17900.0


def _find_default_input() -> Path | None:
    hits = glob.glob(str(_DOWNLOADS_DATA / "*by Genus*.xlsx"))
    hits.sort(key=lambda p: -Path(p).stat().st_mtime)
    return Path(hits[0]) if hits else None


def parse_by_genus(path: Path) -> tuple[pd.DataFrame, str]:
    """Parse a Hub by-Genus Data Table export -> (genus, investment_usd, priority) + filters text.

    The sheet has an "Applied filters:" preamble, then a header row whose first cell is
    "Type of Genus" (col0=genus, col1=USD amount, col3=priority). Rows past the header with a
    parseable amount are kept.
    """
    raw = pd.read_excel(path, header=None)
    applied = ""
    header_row = None
    for i, val in raw[0].items():
        text = str(val)
        if text.startswith("Applied filters"):
            applied = text.replace("\\n", " | ")
        if text.strip() == "Type of Genus":
            header_row = i
            break
    if header_row is None:
        raise ValueError(
            f"{path.name} does not look like a by-Genus Data Table (no 'Type of Genus' "
            "header). Export the 'See Data Table' of the by-Genus panel, not the project list."
        )

    body = raw.iloc[header_row + 1 :, [0, 1, 3]].copy()
    body.columns = ["genus", "investment_usd", "priority"]
    body["genus"] = body["genus"].astype("string").str.strip()
    body["investment_usd"] = pd.to_numeric(body["investment_usd"], errors="coerce")
    body = body[body["genus"].notna() & body["investment_usd"].notna()]
    return body.reset_index(drop=True), applied


def to_panel_funding(body: pd.DataFrame) -> dict[str, float]:
    """Map genus rows onto the six panel species (US$ millions)."""
    funding: dict[str, float] = {}
    for _, row in body.iterrows():
        low = row["genus"].lower()
        for key, species in _GENUS_TO_PANEL.items():
            if low.startswith(key):
                funding[species] = funding.get(species, 0.0) + float(row["investment_usd"]) / 1e6
                break
    return funding


def main(argv: list[str]) -> int:
    in_path = Path(argv[1]) if len(argv) > 1 else _find_default_input()
    if in_path is None or not in_path.exists():
        print("No by-Genus export found. Pass a path, or drop a '*by Genus*.xlsx' into "
              f"{_DOWNLOADS_DATA}. Export it via the dashboard's 'See Data Table' on the "
              "by-Genus panel (filters: Public+Philanthropic, Therapeutics, 2017-2023, Human).")
        return 2

    body, applied = parse_by_genus(in_path)
    funding = to_panel_funding(body)
    total_all = float(body["investment_usd"].sum()) / 1e6

    print(f"Parsed {in_path.name}: {len(body)} genus rows, total ${total_all:,.0f}M\n")
    print(f"Applied filters (from export): {applied or '<none recorded>'}\n")

    missing = [s for s in _GENUS_TO_PANEL.values() if s not in funding]
    print("Panel funding (US$ M), genus-mapped:")
    for species in rd_alignment.GRAM_BURDEN_2019:
        print(f"  {species:28} {funding.get(species, float('nan')):>10.1f}")
    if missing:
        print(f"  WARNING: panel species absent from export: {missing}")

    # Scope check.
    near = min((_SCOPE_TARGET_MUSD, "in-scope (~public+phil therapeutics 2017-23)"),
               (_FULL_SCOPE_MUSD, "OUT OF SCOPE (looks like all funders/areas/years)"),
               key=lambda t: abs(total_all - t[0]))
    print(f"\nScope check: total ${total_all:,.0f}M -> {near[1]}")
    if "OUT OF SCOPE" in near[1]:
        print("  -> Re-export with the dashboard filters applied before using as the "
              "primary denominator (this is fine as the robustness cross-check).")

    # Write tidy CSV under data/rd_hub/ (gitignored; Hub data is public, not DUA-bound).
    extract_date = argv[2] if len(argv) > 2 else (config.RD_HUB_SNAPSHOT_DATE or "undated")
    out_dir = config.DATA_DIR / "rd_hub"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"rd_hub_genus_{extract_date}.csv"
    out = pd.DataFrame(
        [{"pathogen": s, "investment_musd": round(v, 6)} for s, v in funding.items()]
    )
    out.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")

    # Show the resulting Component-4 mismatch using this parsed funding.
    print("\nComponent-4 mismatch under this export (genus_robustness_alignment):")
    res = rd_alignment.genus_robustness_alignment(funding_musd=funding, draws=3000, seed=4)
    for species in res["underfunded_ranking"]:
        m = res["per_pathogen"][species]["log2_mismatch_median"]
        tag = "under-funded" if m > 0 else "over-funded"
        print(f"  {species:28} log2={m:+.2f}  {tag}")
    print(f"  Spearman rho = {res['spearman']['spearman_rho']:+.3f} (n=6, descriptive)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
