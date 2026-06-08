# Empiric-adequacy stewardship what-if tool

A re-runnable local-calibration calculator (Step 5 of the AMR Sentinel pipeline). It
takes the SPIDAAR competing-risks g-formula's **population-aggregate** per-patient
constants and lets a catchment-region stewardship programme scale them to local
parameters — patients treated, current vs target empiric-therapy adequacy, and the
country bed-day cost — to project **deaths averted**, **bed-days added**, and **cost**.

The defensible message it surfaces (the *AMR burden paradox*): improving empiric adequacy
**saves lives but adds bed-days**, because patients who would have died now survive to
discharge. The leverage is systemic — getting the right empiric drug to the patient —
not the resistance phenotype itself.

## Run

```bash
streamlit run app/streamlit_app.py
```

By default it loads the bundled **synthetic demo** artifact (`app/demo_calibration.json`)
— no real patient data, safe to run and share publicly. Upload your own calibration JSON
in the sidebar to use a different cohort.

## The firewall (why this is safe to ship)

The tool never carries raw patient data. It reads only a small **de-identified calibration
artifact**: standardized scalars computed over ≥150 patients, with any stratum cell below
`min_cell_n` suppressed and never exported (Vivli k-anonymity / egress review). Every
on-screen output is labelled an **ecological-calibration scenario, not an estimated
individual treatment effect**.

## Generating a real calibration artifact (inside the Vivli secure environment)

```python
import json
from amr_sentinel_vivli.data_loading import load_spidaar
from amr_sentinel_vivli.stewardship_gformula import build_calibration_artifact

artifact = build_calibration_artifact(load_spidaar(), min_cell_n=5,
                                      source="SPIDAAR (secure env)")
with open("app/demo_calibration.json", "w") as f:   # subject to egress review
    json.dump(artifact, f, indent=2)
```

To regenerate the **synthetic demo** artifact (public):

```python
import json
from amr_sentinel_vivli.stewardship_gformula import synthetic_cohort, build_calibration_artifact
art = build_calibration_artifact(synthetic_cohort(seed=20260607, per_cell=10),
                                 source="synthetic-demo (no real patient data)")
json.dump(art, open("app/demo_calibration.json", "w"), indent=2)
```

## Status

**Exploratory.** The treatment node (`txadp` empiric-therapy adequacy) coding direction is
a **Gate-A** item — unverified against the official SPIDAAR codebook, and it flips the sign
of every projection. The sidebar shows the current assumed coding. No stewardship
recommendation is final without local validation. Apache-2.0.
