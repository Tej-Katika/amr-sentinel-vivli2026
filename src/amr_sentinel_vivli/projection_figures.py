"""Figures for the ATLAS catchment nowcast + SPIDAAR frame-contrast (Component 3).

Rendering only (matplotlib imported lazily); the underlying numbers come from
``bayesian_projection`` (``nowcast`` / ``frame_contrast``), which are unit-tested there.
"""

from __future__ import annotations

import numpy as np


def plot_frame_contrast(contrast: dict, out_path, title: str | None = None):
    """Render the SPIDAAR-vs-ATLAS 3GC-R frame-contrast to ``out_path`` (PNG).

    Grouped bars per country (plus "overall") of ATLAS mixed-surveillance ceftazidime-R
    against SPIDAAR severe-HAI 3GC-R, each with its Wilson interval — making the
    severe-HAI/severity frame-shift visible. Returns the path.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    countries = sorted(set(contrast["atlas"]["by_country"])
                       | set(contrast["spidaar"]["by_country"]))
    labels = countries + ["overall"]
    _empty = {"prevalence": np.nan, "ci_lower": np.nan, "ci_upper": np.nan}

    def _series(source):
        rows = []
        for c in countries:
            cell = contrast[source]["by_country"].get(c)
            rows.append(cell if cell else _empty)
        rows.append(contrast[source]["overall"])
        vals = np.array([r["prevalence"] for r in rows], dtype=float)
        lo = np.array([r["ci_lower"] for r in rows], dtype=float)
        hi = np.array([r["ci_upper"] for r in rows], dtype=float)
        return vals, np.vstack([vals - lo, hi - vals])

    atlas_v, atlas_err = _series("atlas")
    spid_v, spid_err = _series("spidaar")

    x = np.arange(len(labels))
    w = 0.38
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(x - w / 2, atlas_v, w, yerr=atlas_err, capsize=3, color="#2471a3",
           label="ATLAS surveillance (ceftazidime-R)")
    ax.bar(x + w / 2, spid_v, w, yerr=spid_err, capsize=3, color="#c0392b",
           label="SPIDAAR severe-HAI (3GC-R)")
    ax.set(xticks=x, ylabel="3GC resistance prevalence", ylim=(0, 1),
           title=title or f"Frame-contrast: {contrast['panel']}")
    ax.set_xticklabels(labels)
    ax.legend(frameon=False, fontsize=8)
    ax.margins(x=0.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
