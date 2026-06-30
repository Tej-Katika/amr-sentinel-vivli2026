"""Render the four report figures (PNG) for secure-environment egress.

The figure renderers are deliberately NOT wired into ``pipeline.run()`` (matplotlib is an
optional dependency and the PNGs are gitignored). This driver computes just the inputs each
figure needs from the delivered Vivli data and writes all four PNGs that
``docs/final_report_2026.md`` embeds, so the secure-env figure step is a single command:

    PYTHONPATH=src python scripts/make_figures.py [out_dir]

``out_dir`` defaults to ``figures/`` (repo root; the paths the report links to). Requires
``matplotlib`` in addition to the numpy/pandas/xlrd runtime. Filenames are kept EXACTLY as
the report references them — do not rename without updating the report.
"""

from __future__ import annotations

import sys
from pathlib import Path

from amr_sentinel_vivli import data_loading
from amr_sentinel_vivli.bayesian_projection import frame_contrast
from amr_sentinel_vivli.evidence_figures import plot_evidence_forest
from amr_sentinel_vivli.evidence_synthesis import run_evidence_synthesis
from amr_sentinel_vivli.excess_los_figures import plot_excess_los
from amr_sentinel_vivli.projection_figures import plot_frame_contrast
from amr_sentinel_vivli.rd_alignment import catchment_alignment, gram_panel_alignment
from amr_sentinel_vivli.rd_alignment_figures import plot_mismatch_global_vs_catchment

# Some figure labels carry non-ASCII (τ, →); force UTF-8 so stdout does not die on cp1252.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def main(argv: list[str]) -> int:
    out_dir = Path(argv[1]) if len(argv) > 1 else Path("figures")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading delivered data (SPIDAAR patient + isolates, ATLAS) ...")
    spidaar = data_loading.load_spidaar()
    spidaar_isolates = data_loading.load_spidaar_isolates()
    atlas = data_loading.load_atlas()

    print("Computing figure inputs ...")
    evidence = run_evidence_synthesis(spidaar)
    contrast = frame_contrast(atlas, spidaar_isolates)
    rd_global = gram_panel_alignment()
    rd_catchment = catchment_alignment(spidaar_isolates)

    # (report-embed filename -> renderer). Names MUST match docs/final_report_2026.md.
    figures = [
        ("excess_los_stateoccupation.png", lambda p: plot_excess_los(spidaar, p)),
        ("evidence_forest.png", lambda p: plot_evidence_forest(evidence, p)),
        ("spidaar_framecontrast.png", lambda p: plot_frame_contrast(contrast, p)),
        ("rd_mismatch_global_vs_catchment.png",
         lambda p: plot_mismatch_global_vs_catchment(rd_global, rd_catchment, p)),
    ]

    print(f"Rendering {len(figures)} figures to {out_dir}/ ...")
    for name, render in figures:
        render(out_dir / name)
        print(f"  wrote {out_dir / name}")

    print(f"\nDone: {len(figures)} figures in {out_dir}/. Submit these through egress review "
          "alongside confirmatory_results.json.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
