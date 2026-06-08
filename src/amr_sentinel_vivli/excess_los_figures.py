"""Figures for the excess-LOS competing-risks analysis (Component 1).

Separates the *figure data* (pure, deterministic, unit-tested — per-arm state-occupation
and cumulative-incidence curves) from the *rendering* (matplotlib, imported lazily so the
package imports without it). The headline figure shows, per resistance arm, the
probability of still occupying a bed over time (its restricted-mean difference IS the
excess bed-days) alongside the competing discharge/death cumulative incidences at τ.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .excess_los import (
    DEFAULT_TAU,
    EXIT_DEATH,
    EXIT_DISCHARGE,
    build_exit_frame,
)


def state_occupation_curves(df: pd.DataFrame, tau: int = DEFAULT_TAU) -> dict:
    """Per-arm Aalen-Johansen state-occupation curves on integer days [0, τ].

    For each arm returns, at each day, the cumulative incidence of discharge and of death
    and the residual probability of still being admitted (the bed-day occupancy curve).
    Restricts to the ascertained resistant-vs-susceptible cohort. Pure/deterministic.
    """
    asc = df[pd.to_numeric(df["resistant"], errors="coerce").isin([0, 1])]
    ef = build_exit_frame(asc, tau)
    ef = ef.assign(resistant=pd.to_numeric(asc["resistant"], errors="coerce").to_numpy())
    days = np.arange(0, tau + 1)

    out: dict = {"days": days.tolist()}
    for arm, name in ((1, "resistant"), (0, "susceptible")):
        a = ef[ef["resistant"] == arm]
        t = a["time"].to_numpy(dtype=float)
        et = a["exit_type"].to_numpy(dtype=int)

        surv = 1.0
        inc: list = []  # (event_time, d_cif_discharge, d_cif_death)
        for ti in np.unique(t[et != 0]):
            if ti > tau:
                break
            at_risk = np.sum(t >= ti)
            d_disc = np.sum((t == ti) & (et == EXIT_DISCHARGE))
            d_death = np.sum((t == ti) & (et == EXIT_DEATH))
            inc.append((ti, surv * d_disc / at_risk, surv * d_death / at_risk))
            surv *= 1.0 - (d_disc + d_death) / at_risk

        cif_disc = np.array([sum(i[1] for i in inc if i[0] <= d) for d in days])
        cif_death = np.array([sum(i[2] for i in inc if i[0] <= d) for d in days])
        out[name] = {
            "n": int(len(a)),
            "cif_discharge": cif_disc.tolist(),
            "cif_death": cif_death.tolist(),
            "p_admitted": (1.0 - cif_disc - cif_death).tolist(),
        }
    return out


def plot_excess_los(df: pd.DataFrame, out_path, tau: int = DEFAULT_TAU):
    """Render the headline excess-LOS figure to ``out_path`` (PNG). Returns the path.

    Left panel: bed occupancy P(admitted) over time by arm — the area between the curves
    to τ is the excess bed-days. Right panel: competing discharge vs death cumulative
    incidence at τ by arm. matplotlib is imported here so the package does not require it.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    curves = state_occupation_curves(df, tau)
    days = np.asarray(curves["days"], dtype=float)
    colors = {"resistant": "#c0392b", "susceptible": "#2471a3"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    for name in ("resistant", "susceptible"):
        ax1.plot(days, curves[name]["p_admitted"], color=colors[name], lw=2,
                 label=f"{name} (n={curves[name]['n']})")
    ax1.fill_between(days, curves["resistant"]["p_admitted"],
                     curves["susceptible"]["p_admitted"], color="grey", alpha=0.18)
    ax1.set(xlabel="Days since enrolment", ylabel="P(still admitted)",
            title="Bed occupancy by resistance arm\n(area between curves = excess bed-days)",
            xlim=(0, tau), ylim=(0, 1))
    ax1.legend(frameon=False)

    arms = ["resistant", "susceptible"]
    disc = [curves[a]["cif_discharge"][-1] for a in arms]
    death = [curves[a]["cif_death"][-1] for a in arms]
    admit = [curves[a]["p_admitted"][-1] for a in arms]
    x = np.arange(len(arms))
    ax2.bar(x, disc, 0.6, label="discharged alive", color="#27ae60")
    ax2.bar(x, death, 0.6, bottom=disc, label="died", color="#7b241c")
    ax2.bar(x, admit, 0.6, bottom=np.add(disc, death), label="still admitted", color="#bdc3c7")
    ax2.set(xticks=x, ylabel=f"Cumulative incidence at day {tau}",
            title="Competing outcomes at τ", ylim=(0, 1))
    ax2.set_xticklabels(arms)
    ax2.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
