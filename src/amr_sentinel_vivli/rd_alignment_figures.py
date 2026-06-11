"""Figure for the Cross-Domain R&D mismatch (Component 4): global vs catchment.

Rendering only (matplotlib imported lazily); the numbers come from ``rd_alignment``
(``gram_panel_alignment`` and ``catchment_alignment``), unit-tested there.
"""

from __future__ import annotations

import numpy as np

# Short display labels for the six GRAM panel pathogens.
_SHORT = {
    "Escherichia coli": "E. coli",
    "Staphylococcus aureus": "S. aureus",
    "Klebsiella pneumoniae": "K. pneumoniae",
    "Streptococcus pneumoniae": "S. pneumoniae",
    "Acinetobacter baumannii": "A. baumannii",
    "Pseudomonas aeruginosa": "P. aeruginosa",
}


def plot_mismatch_global_vs_catchment(global_result: dict, catchment_result: dict, out_path,
                                      title: str | None = None):
    """Paired log2-mismatch bars (global GRAM vs SSA-catchment) to ``out_path`` (PNG).

    For each pathogen, the log2 burden/funding mismatch under the global GRAM burden and under
    the SSA severe-HAI catchment burden, ordered by the global value. Positive = under-funded
    relative to burden; the zero line is parity. Makes the catchment re-ordering visible —
    S. pneumoniae falls from top (global) to over-funded (catchment) while the Gram-negatives
    rise. Returns the path.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = global_result["per_pathogen"]
    c = catchment_result["per_pathogen"]
    pathogens = sorted(g, key=lambda p: g[p]["log2_mismatch_median"])  # ascending -> top is most
    gv = np.array([g[p]["log2_mismatch_median"] for p in pathogens])
    cv = np.array([c[p]["log2_mismatch_median"] for p in pathogens])

    y = np.arange(len(pathogens))
    h = 0.38
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    ax.axvline(0.0, color="#444444", lw=1, ls="--", zorder=0)
    ax.barh(y + h / 2, gv, h, color="#2471a3", label="Global (GRAM burden)")
    ax.barh(y - h / 2, cv, h, color="#c0392b", label="Catchment (SSA severe-HAI)")

    ax.set_yticks(y)
    ax.set_yticklabels([_SHORT.get(p, p) for p in pathogens], fontsize=9)
    ax.set_xlabel("log2 mismatch  (>0 = under-funded relative to burden)")
    ax.set_title(title or "Burden–funding mismatch: global vs catchment")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(y=0.04)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
