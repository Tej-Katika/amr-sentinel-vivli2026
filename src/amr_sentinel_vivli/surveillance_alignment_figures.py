"""Figure for Component 6: the surveillance blind-spot axis.

Rendering only (matplotlib imported lazily); numbers come from ``surveillance_alignment``
(``run_surveillance_alignment``), unit-tested there. Left: per-pathogen burden / funding /
surveillance shares (the pathogen neglected on all three stands out). Right: where the
surveillance physically is — the top single contributor vs all of sub-Saharan Africa vs the
study catchment — the geographic blind spot.
"""

from __future__ import annotations

import numpy as np

_SHORT = {
    "Escherichia coli": "E. coli",
    "Staphylococcus aureus": "S. aureus",
    "Klebsiella pneumoniae": "K. pneumoniae",
    "Streptococcus pneumoniae": "S. pneumoniae",
    "Acinetobacter baumannii": "A. baumannii",
    "Pseudomonas aeruginosa": "P. aeruginosa",
}


def plot_surveillance_blindspot(result: dict, out_path, title: str | None = None):
    """Render the two-panel surveillance blind-spot figure to ``out_path`` (PNG).

    ``result`` is a :func:`surveillance_alignment.run_surveillance_alignment` dict. Returns
    the path. matplotlib is imported here so the package does not require it.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ta = result["three_axis"]["per_pathogen"]
    geo = result["geographic"]

    # Order pathogens by burden share (descending) for the left panel.
    pathogens = sorted(ta, key=lambda p: ta[p]["burden_share"], reverse=True)
    burden = np.array([ta[p]["burden_share"] for p in pathogens])
    funding = np.array([ta[p]["funding_share"] for p in pathogens])
    surveil = np.array([ta[p]["surveillance_share"] for p in pathogens])
    labels = [_SHORT[p] for p in pathogens]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6),
                                   gridspec_kw={"width_ratios": [1.7, 1]})

    x = np.arange(len(pathogens))
    w = 0.27
    ax1.bar(x - w, burden, w, label="burden (GRAM deaths)", color="#7b241c")
    ax1.bar(x, funding, w, label="R&D funding (Hub)", color="#c0392b")
    ax1.bar(x + w, surveil, w, label="surveillance (ATLAS)", color="#2471a3")
    ax1.axhline(1 / len(pathogens), ls=":", color="grey", lw=1)  # equal-share line
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax1.set_ylabel("share of panel (each axis sums to 1)")
    ax1.set_title("Burden vs funding vs surveillance, by pathogen")
    ax1.legend(frameon=False, fontsize=8)
    for p in result["three_axis"]["neglected_on_all_axes"]:
        i = pathogens.index(p)
        ax1.annotate("neglected\non all 3", (i, max(burden[i], funding[i], surveil[i]) + 0.02),
                     ha="center", fontsize=7, color="#7b241c")

    # Right: geographic concentration of surveillance.
    reg_labels = [geo["top_country"]["country"], "Sub-Saharan\nAfrica", "Study\ncatchment"]
    reg_share = np.array([geo["top_country"]["share"], geo["ssa"]["share"],
                          geo["catchment"]["share"]]) * 100
    colors = ["#95a5a6", "#c0392b", "#7b241c"]
    yb = np.arange(len(reg_labels))[::-1]
    ax2.barh(yb, reg_share, color=colors)
    ax2.set_yticks(yb)
    ax2.set_yticklabels(reg_labels, fontsize=8)
    ax2.set_xlabel("% of all ATLAS isolates")
    ax2.set_title("Where surveillance is\n(SSA has the highest AMR death rate)")
    for yi, v in zip(yb, reg_share, strict=True):
        ax2.annotate(f"{v:.1f}%", (v, yi), xytext=(4, 0), textcoords="offset points",
                     va="center", fontsize=8)
    ratio = geo["top_country_to_ssa_ratio"]
    ax2.annotate(f"{geo['top_country']['country']} alone =\n{ratio:.0f}× all of SSA",
                 (0.97, 0.05), xycoords="axes fraction", ha="right", fontsize=7.5,
                 color="#7b241c")

    fig.suptitle(title or "The surveillance blind spot: the highest-burden region is the "
                 "least watched", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
