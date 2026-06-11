"""Figure for the Bayesian evidence synthesis (Component 1c).

Rendering only (matplotlib imported lazily); the numbers come from ``evidence_synthesis``
(``run_evidence_synthesis``), which is unit-tested there.
"""

from __future__ import annotations

import numpy as np

from .evidence_synthesis import EXTERNAL_MORTALITY_STUDIES


def plot_evidence_forest(synthesis: dict, out_path, title: str | None = None):
    """Forest plot of the resistance→in-hospital-death synthesis to ``out_path`` (PNG).

    Rows: the adjusted external cohorts (MBIRA, Fiji), our crude cohort estimate, the
    random-effects pooled HR (diamond), and the prediction interval for a new setting — all
    on a log HR axis with the null (HR=1) marked. Returns the path.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cohort = synthesis["cohort_crude"]
    pool = synthesis["primary_pool_adjusted"]

    rows = []
    for name, s in EXTERNAL_MORTALITY_STUDIES.items():
        rows.append((name, s["hr"], s["ci_low"], s["ci_high"], "adjusted"))
    rows.append(("This cohort (crude)", cohort["hr"], cohort["ci_low"], cohort["ci_high"], "crude"))
    rows.append(("Pooled (random-effects)", pool["pooled_hr"], pool["hr_ci"][0],
                 pool["hr_ci"][1], "pool"))
    rows.append(("Prediction interval", pool["pooled_hr"], pool["prediction_hr_ci"][0],
                 pool["prediction_hr_ci"][1], "pred"))

    colors = {"adjusted": "#2471a3", "crude": "#7f8c8d", "pool": "#c0392b", "pred": "#c0392b"}
    y = np.arange(len(rows))[::-1]  # first row at the top

    fig, ax = plt.subplots(figsize=(8.5, 4.0))
    ax.axvline(1.0, color="#444444", lw=1, ls="--", zorder=0)
    for yi, (_name, hr, lo, hi, kind) in zip(y, rows, strict=True):
        c = colors[kind]
        if kind == "pred":
            ax.plot([lo, hi], [yi, yi], color=c, lw=1.4, ls=":", zorder=2)
        else:
            ax.plot([lo, hi], [yi, yi], color=c, lw=2, zorder=2)
            ax.scatter([hr], [yi], s=(120 if kind == "pool" else 55),
                       marker=("D" if kind == "pool" else "o"), color=c, zorder=3)
        ax.annotate(f"{hr:.2f} ({lo:.2f}–{hi:.2f})", xy=(1.0, yi), xytext=(0, 9),
                    textcoords="offset points", ha="center", va="bottom", fontsize=7.5,
                    color=c)

    ax.set_xscale("log")
    ax.set_xlim(0.25, 25)
    ax.set_xticks([0.5, 1, 2, 5, 10, 20])
    ax.set_xticklabels(["0.5", "1", "2", "5", "10", "20"])
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xlabel("In-hospital death hazard ratio (3GC-R vs 3GC-S), log scale")
    ax.set_title(title or "Resistance→mortality: random-effects evidence synthesis")
    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(y=0.12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
