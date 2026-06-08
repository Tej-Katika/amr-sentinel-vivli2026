"""Empiric-adequacy stewardship what-if tool (Vivli 2026 AMR Data Challenge, Step 5).

A re-runnable local-calibration calculator for catchment-region stewardship programs.
It loads a small DE-IDENTIFIED calibration artifact (population-aggregate per-patient
constants from the SPIDAAR competing-risks g-formula — no patient records) and scales it
to local parameters: how many resistant-HAI patients you treat, your current vs target
empiric-therapy adequacy, your country's bed-day cost.

The honest, defensible message it surfaces (the "AMR burden paradox"): improving empiric
adequacy **averts deaths but adds bed-days**, because patients who would have died now
survive to discharge. Every output is an *ecological-calibration scenario*, never an
estimated individual treatment effect — the firewall is shown on screen at all times.

Run:  streamlit run app/streamlit_app.py
Status: EXPLORATORY; the txadp adequacy coding is Gate-A-unverified (see the sidebar).
Apache-2.0.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import streamlit as st

# Make the src/ package importable when run via `streamlit run app/streamlit_app.py`
# (no install step required for the demo).
_SRC = Path(__file__).resolve().parents[1] / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from amr_sentinel_vivli.stewardship_gformula import (  # noqa: E402
    adequacy_scenario,
    monte_carlo_cost,
)

_DEMO_ARTIFACT = Path(__file__).resolve().parent / "demo_calibration.json"
_ARM_LABELS = {
    "pooled": "All patients (pooled)",
    "resistant": "Resistant infections only",
    "susceptible": "Susceptible infections only",
}


def _load_artifact(uploaded) -> dict:
    """Load a calibration artifact: an uploaded file if given, else the bundled demo."""
    if uploaded is not None:
        return json.load(uploaded)
    with open(_DEMO_ARTIFACT) as f:
        return json.load(f)


def main() -> None:
    st.set_page_config(page_title="AMR Empiric-Adequacy Stewardship What-If", layout="wide")
    st.title("Empiric-adequacy stewardship what-if")
    st.caption(
        "Sub-Saharan Africa HAI cohort (SPIDAAR-calibrated) · competing-risks g-formula · "
        "Vivli 2026 AMR Data Challenge"
    )

    # --- The firewall: always on screen -------------------------------------
    st.warning(
        "**Ecological-calibration scenario — not an estimated individual effect.** "
        "These projections scale population-aggregate constants to your local inputs. "
        "They are exploratory and depend on a Gate-A-unverified treatment coding (see "
        "sidebar). Treat them as a planning calculator, not a causal estimate for any "
        "individual patient or a final stewardship recommendation.",
        icon="⚠️",
    )

    # --- Sidebar: artifact provenance + global selectors --------------------
    with st.sidebar:
        st.header("Calibration source")
        uploaded = st.file_uploader(
            "Calibration artifact (JSON)", type="json",
            help="Leave empty to use the bundled synthetic demo. In the Vivli secure "
                 "environment, generate a real artifact with "
                 "stewardship_gformula.build_calibration_artifact() after egress review.",
        )
        artifact = _load_artifact(uploaded)
        prov = artifact["provenance"]
        st.write(f"**Source:** {prov['source']}")
        st.write(
            f"Ascertained n = {prov['n_ascertained']} · on-support n = {prov['n_on_support']} · "
            f"strata = {prov['n_strata_supported']} "
            f"(suppressed < min n: {prov['n_cells_below_min_n_suppressed']})"
        )
        st.info(f"**Gate A:** {prov['gate_a_note']}", icon="🔬")

        st.header("Scenario settings")
        arm = st.radio(
            "Population", options=list(_ARM_LABELS), format_func=_ARM_LABELS.get,
            help="Resistance determines adequacy, so the principled view is the controlled "
                 "direct effect within a resistance stratum.",
        )
        currency = st.radio(
            "Currency", options=["usd2010", "ppp"],
            format_func={"usd2010": "2010 USD (market FX)", "ppp": "PPP int-$ (2010)"}.get,
        )

    if artifact["per_patient"].get(arm) is None:
        st.error(
            f"The loaded calibration has no usable '{arm}' arm "
            "(the contrast was not identified at that resolution). Pick another population."
        )
        st.stop()

    countries = sorted(artifact["unit_cost_usd2010"])

    # --- Inputs -------------------------------------------------------------
    st.subheader("Your local programme")
    c1, c2 = st.columns(2)
    with c1:
        country = st.selectbox("Country", countries)
        n_patients = st.number_input(
            "Resistant-HAI patients treated per year", min_value=0, value=200, step=10,
        )
    with c2:
        current_pct = st.slider("Current empiric adequacy (%)", 0, 100, 50)
        target_pct = st.slider("Target empiric adequacy (%)", 0, 100, 80)

    if target_pct < current_pct:
        st.info("Target adequacy is below current — no patients are upgraded in this scenario.")

    scenario = adequacy_scenario(
        artifact,
        n_patients=float(n_patients),
        current_adequacy=current_pct / 100.0,
        target_adequacy=target_pct / 100.0,
        country=country,
        currency=currency,
        arm=arm,
    )

    # --- Results ------------------------------------------------------------
    st.subheader("Projected effect of the adequacy improvement")
    m1, m2, m3 = st.columns(3)
    m1.metric("Patients upgraded / yr", f"{scenario['patients_upgraded']:.0f}")
    m2.metric(
        "Deaths averted / yr", f"{scenario['averted_deaths']:.1f}",
        help="Positive = lives saved by adequate empiric therapy.",
    )
    m3.metric(
        "Additional bed-days / yr", f"{scenario['added_bed_days']:+.0f}",
        help="Competing-risks paradox: averted deaths become survivors who occupy beds. "
             "Positive = more occupancy.",
    )

    unit = "USD (2010)" if currency == "usd2010" else "int-$ (2010 PPP)"
    st.metric(
        f"Additional bed-day cost / yr — {unit}",
        f"{scenario['added_cost']:+,.0f}",
        help=f"At {scenario['unit_cost_per_bed_day']:.2f} {unit} per bed-day (WHO-CHOICE).",
    )

    # Cost uncertainty from the WHO-CHOICE 95% UI (PPP int-$).
    if abs(scenario["added_bed_days"]) > 0:
        mc = monte_carlo_cost({country: scenario["added_bed_days"]}, draws=10000, seed=5)
        st.caption(
            f"Cost uncertainty (WHO-CHOICE PPP int-$ 95% UI): median "
            f"{mc['median']:+,.0f}, 95% CI [{mc['ci_lower']:+,.0f}, {mc['ci_upper']:+,.0f}] int-$."
        )

    st.markdown(
        "**How to read this.** Improving empiric adequacy is projected to *save lives* "
        "while *adding* bed-days and cost — not because resistance is harmless, but because "
        "effective therapy converts deaths into survivors who stay to be discharged. The "
        "leverage is systemic (getting the right empiric drug to the patient), and its bed-day "
        "cost is the price of averted mortality, not a failure."
    )

    with st.expander("Calibration constants & provenance"):
        st.json({"per_patient": artifact["per_patient"][arm], "provenance": prov})
        st.caption(scenario["firewall"])

    st.caption("Apache-2.0 · No stewardship recommendation is final without local validation.")


if __name__ == "__main__":
    main()
