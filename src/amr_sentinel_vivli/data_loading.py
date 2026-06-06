"""Data loading for the Vivli secure research environment.

All inputs are Vivli DUA-restricted and live under ``data/`` (gitignored). These
loaders run *inside* the Vivli secure environment after the DUA is signed. The
raw inputs never leave that environment, so the raw→analysis mapping is factored
into pure, side-effect-free helpers (``_restrict_to_catchment``,
``_build_spidaar_analysis_frame``) that are unit-tested against small *synthetic*
frames — no real patient record is ever committed or used in a test.

SPIDAAR linkage deviation (logged in ``docs/deviation_log.md``)
--------------------------------------------------------------
The pre-registered Step-1 unit of analysis is the isolate (resistant vs.
susceptible *of the same species*; pre-reg §6-7). As delivered, SPIDAAR ships two
unlinked workbooks — ``spidaar_isolatedata.xls`` (244 isolates: organism, site,
country, sample month/year, resistance) and ``spidaar_patientdata.xls`` (336
patients: 30-day mortality, age, sex, patient-level resistance summary). They
share **no join key** (no ``pid`` in the isolate file, no ``iid`` in the patient
file, no crosswalk sheet), and the counts do not reconcile for a positional join
(244 isolate rows vs. ``sum(nisol)=235`` vs. 207 patients with ≥1 isolate). The
mortality outcome and age/sex exist only at the patient level. We therefore run
Step 1 at the **patient level**: each patient's resistance exposure is the
pre-derived patient-level summary ``amrp`` (≥1 resistant isolate vs. all
susceptible), with the patient's 30-day mortality outcome. This is a deviation
from the isolate-level plan; the calendar-``year`` propensity covariate is also
unavailable patient-side and is dropped. Both are logged per pre-reg §11.
"""

from __future__ import annotations

import pandas as pd

from . import config

_NOT_YET = (
    "Pre-data scaffold — no Vivli data has been accessed (pre-reg §15). "
    "Implement loaders only after the DUA is signed and the secure environment entered."
)

# 30-day horizon for the Step-1 survival outcome (config.OUTCOME == "mortality_30d").
MORTALITY_HORIZON_DAYS = 30

# --- SPIDAAR codebook maps (mirror the workbook ``codebook`` sheet) ---------
# Sex (patient ``sex``): 3 = "Missing" -> NaN.
_SEX_LABELS = {0: "Male", 1: "Female"}

# Age class (patient ``agegr``); ordinal 0..14. Kept as an ordinal code for the
# propensity model (true age in years is not released — only the class).
_AGE_GROUP_LABELS = {
    0: "<1", 1: "1-5", 2: "6-10", 3: "11-15", 4: "16-20", 5: "21-25",
    6: "26-30", 7: "31-35", 8: "36-40", 9: "41-45", 10: "46-50",
    11: "51-55", 12: "56-60", 13: "61-65", 14: "65+",
}

# Confirmed HAI category (patient ``chaicat``); 10 = "No HAI confirmed" -> NaN.
_INFECTION_SITE_LABELS = {
    1: "BSI", 2: "cUTI", 3: "cSSTI", 4: "HAP", 5: "cIAI",
    6: "BSI+cUTI", 7: "BSI+cSSTI", 8: "HAP+cSSTI", 9: "BSI+cUTI+cIAI",
}

# Patient-level resistance summary (``amrp``): exposure for the patient-level H1.
#   2  -> ≥1 resistant isolate to treatment AB class   -> resistant (1)
#   0  -> all isolates susceptible                      -> susceptible (0)
#   1  -> mixed S, untested (no resistant isolate seen) -> susceptible (0)
#  -1  -> no susceptibility result / no isolate         -> unascertainable (NaN)
_RESISTANT_FROM_AMRP = {2: 1.0, 0: 0.0, 1: 0.0}

# Isolate-file binary resistance flags (``c3r``, ``mdr``, ``mrsa``, ``rx``):
#   1 -> resistant / positive, 0 -> susceptible / negative. Any other code is a
# sentinel for "not tested / not applicable" (e.g. ``9`` for 3GC-R on a
# Gram-positive, ``99`` for MRSA on a non-staphylococcus) and maps to NaN.
_RESISTANT_BINARY = {1: 1.0, 0: 0.0}

# Mortality status for analysis (``dead``): 0 alive, 1 deceased, 9 deceased-censored.
_DEAD_DECEASED_CODES = (1, 9)
_DEAD_VALID_CODES = (0, 1, 9)

# Analysis-frame columns returned by load_spidaar (no ``year`` — see module docstring).
SPIDAAR_COLUMNS = (
    "pid", "country", "organism", "infection_site", "age", "age_group", "sex",
    "resistant", "amrp", "mortality_30d", "time_at_risk",
    "dead", "days_to_death", "days_observed",
)

# Analysis-frame columns returned by load_spidaar_isolates. Unlike the patient
# frame this carries NO outcome (mortality is patient-side only and unlinkable)
# but DOES carry the per-mechanism resistance detail the patient summary collapses
# away (``c3r`` / ``mdr`` / ``mrsa``) plus sample timing for the Step-2 anchor.
ISOLATE_COLUMNS = (
    "iid", "country", "organism", "organism_group", "specimen",
    "infection_site", "clinical_relevance", "sample_month", "sample_year",
    "resistant", "c3r", "mdr", "mrsa", "amrtx_resistant",
)


def _read_data_sheet(path, sheet: str = "data") -> pd.DataFrame:
    """Read one sheet of a legacy SPIDAAR ``.xls`` workbook (needs the xlrd engine)."""
    return pd.read_excel(path, sheet_name=sheet)


def _restrict_to_catchment(df: pd.DataFrame, country_col: str = "ctry") -> pd.DataFrame:
    """Keep only rows in the pre-registered SPIDAAR catchment (``config.CATCHMENT_COUNTRIES``)."""
    return df[df[country_col].isin(config.CATCHMENT_COUNTRIES)].copy()


def _build_spidaar_analysis_frame(raw: pd.DataFrame) -> pd.DataFrame:
    """Map a raw SPIDAAR ``patientdata`` 'data' sheet to the Step-1 analysis frame.

    Pure function (no I/O) so it can be unit-tested on synthetic rows. Expects the
    raw patient-level columns ``pid, ctry, agegr, sex, chaicat, isol, amrp, dead,
    dtpta, nobsd`` and returns the columns in ``SPIDAAR_COLUMNS``.

    Outcome (``mortality_30d`` / ``time_at_risk``): a 30-day-horizon survival pair.
    A death observed on/before day 30 (``dead`` in {1, 9} with days-to-death
    ``dtpta`` ≤ 30) is the event; later deaths and survivors are censored at
    ``min(follow-up, 30)`` days, where follow-up is ``dtpta`` for deaths and the
    observation length ``nobsd`` for survivors.

    Exposure (``resistant``): from ``amrp`` per ``_RESISTANT_FROM_AMRP``; NaN where
    resistance is unascertainable (``amrp == -1``) so the Step-1 model can restrict
    to the resistant-vs-susceptible comparison cohort.
    """
    raw = raw.reset_index(drop=True)
    dead = pd.to_numeric(raw["dead"], errors="coerce")
    dtpta = pd.to_numeric(raw["dtpta"], errors="coerce")    # days admission->death
    nobsd = pd.to_numeric(raw["nobsd"], errors="coerce")    # days under observation

    deceased = dead.isin(_DEAD_DECEASED_CODES) & dtpta.notna()
    event = deceased & (dtpta <= MORTALITY_HORIZON_DAYS)
    follow_up = dtpta.where(deceased, nobsd)                 # death time, else obs length
    time_at_risk = follow_up.clip(upper=MORTALITY_HORIZON_DAYS)

    agegr = pd.to_numeric(raw["agegr"], errors="coerce")
    chaicat = pd.to_numeric(raw["chaicat"], errors="coerce")
    sex_code = pd.to_numeric(raw["sex"], errors="coerce")
    amrp = pd.to_numeric(raw["amrp"], errors="coerce")

    out = pd.DataFrame({
        "pid": raw["pid"],
        "country": raw["ctry"],
        "organism": raw["isol"],
        "infection_site": chaicat.map(_INFECTION_SITE_LABELS),
        "age": agegr.astype("Int64"),
        "age_group": agegr.map(_AGE_GROUP_LABELS),
        "sex": sex_code.map(_SEX_LABELS),
        "resistant": amrp.map(_RESISTANT_FROM_AMRP).astype("Int64"),
        "amrp": amrp.astype("Int64"),
        "mortality_30d": event.astype(int),
        "time_at_risk": time_at_risk,
        "dead": dead.astype("Int64"),
        "days_to_death": dtpta,
        "days_observed": nobsd,
    })
    return out[list(SPIDAAR_COLUMNS)]


def _validate_spidaar(df: pd.DataFrame) -> None:
    """Fail fast on the invariants the Step-1 model relies on."""
    stray = set(df["country"].dropna()) - set(config.CATCHMENT_COUNTRIES)
    if stray:
        raise ValueError(f"SPIDAAR rows outside the catchment: {sorted(stray)}")
    bad_dead = set(df["dead"].dropna().astype(int)) - set(_DEAD_VALID_CODES)
    if bad_dead:
        raise ValueError(f"Unexpected `dead` codes: {sorted(bad_dead)}")
    if not set(df["mortality_30d"].unique()) <= {0, 1}:
        raise ValueError("mortality_30d must be binary 0/1.")
    if not set(df["resistant"].dropna().astype(int)) <= {0, 1}:
        raise ValueError("resistant must be binary 0/1 (NaN where unascertainable).")


def load_spidaar() -> pd.DataFrame:
    """Load the SPIDAAR patient-level cohort for the Step-1 mortality model (pre-reg §5).

    Reads ``data/spidaar/spidaar_patientdata.xls``, restricts to
    ``config.CATCHMENT_COUNTRIES``, and returns the analysis frame
    (``SPIDAAR_COLUMNS``) carrying the 30-day mortality outcome
    (``config.OUTCOME``), the patient-level resistance exposure, and the available
    propensity covariates in ``config.PROPENSITY_COVARIATES`` (all except
    ``year``, which is not recorded patient-side).

    The unit is the patient rather than the isolate, and ``year`` is omitted; both
    are deviations from the pre-registered isolate-level plan and are documented in
    the module docstring and ``docs/deviation_log.md`` (pre-reg §11). Rows where
    resistance is unascertainable (``amrp == -1``) keep a NaN ``resistant`` and are
    excluded by the Step-1 model, not dropped here.
    """
    path = config.DATA_DIR / "spidaar" / "spidaar_patientdata.xls"
    raw = _read_data_sheet(path)
    cohort = _restrict_to_catchment(raw)
    frame = _build_spidaar_analysis_frame(cohort)
    _validate_spidaar(frame)
    return frame.reset_index(drop=True)


def _build_spidaar_isolate_frame(raw: pd.DataFrame) -> pd.DataFrame:
    """Map a raw SPIDAAR ``isolatedata`` 'data' sheet to the isolate-level frame.

    Pure function (no I/O) so it can be unit-tested on synthetic rows. Expects the
    raw isolate columns ``iid, isolate, group, chaicat, clinrel, stype, sampdat,
    c3r, mdr, mrsa, amrtx, rx, ctry`` and returns the columns in ``ISOLATE_COLUMNS``.

    Resistance exposure (``resistant``): the workbook composite ``rx`` ("MDR
    and/or 3GC-R"), the isolate-level analog of the patient-level ``amrp`` (1 ->
    resistant, 0 -> not, NaN where not ascertained). The per-mechanism flags
    ``c3r`` / ``mdr`` / ``mrsa`` and the administered-class resistance
    ``amrtx_resistant`` are carried through (NaN where not tested / not applicable).

    Code semantics for the binary flags (``_RESISTANT_BINARY``) and ``amrtx``
    (``_RESISTANT_FROM_AMRP``) are inferred from the patient-side codebook and
    SPIDAAR conventions; the isolate workbook ships variable *labels* but not
    answer *codes*. VERIFY against the official codebook inside the secure
    environment and log any correction (pre-reg §11).
    """
    raw = raw.reset_index(drop=True)
    sampled = pd.to_datetime(raw["sampdat"], format="%m/%Y", errors="coerce")
    chaicat = pd.to_numeric(raw["chaicat"], errors="coerce")
    c3r = pd.to_numeric(raw["c3r"], errors="coerce")
    mdr = pd.to_numeric(raw["mdr"], errors="coerce")
    mrsa = pd.to_numeric(raw["mrsa"], errors="coerce")
    rx = pd.to_numeric(raw["rx"], errors="coerce")
    amrtx = pd.to_numeric(raw["amrtx"], errors="coerce")

    out = pd.DataFrame({
        "iid": raw["iid"],
        "country": raw["ctry"],
        "organism": raw["isolate"],
        "organism_group": pd.to_numeric(raw["group"], errors="coerce").astype("Int64"),
        "specimen": raw["stype"],
        "infection_site": chaicat.map(_INFECTION_SITE_LABELS),
        "clinical_relevance": raw["clinrel"],
        "sample_month": sampled.dt.month.astype("Int64"),
        "sample_year": sampled.dt.year.astype("Int64"),
        "resistant": rx.map(_RESISTANT_BINARY).astype("Int64"),
        "c3r": c3r.map(_RESISTANT_BINARY).astype("Int64"),
        "mdr": mdr.map(_RESISTANT_BINARY).astype("Int64"),
        "mrsa": mrsa.map(_RESISTANT_BINARY).astype("Int64"),
        "amrtx_resistant": amrtx.map(_RESISTANT_FROM_AMRP).astype("Int64"),
    })
    return out[list(ISOLATE_COLUMNS)]


def _validate_spidaar_isolates(df: pd.DataFrame) -> None:
    """Fail fast on the invariants the isolate-level analyses rely on."""
    stray = set(df["country"].dropna()) - set(config.CATCHMENT_COUNTRIES)
    if stray:
        raise ValueError(f"SPIDAAR isolates outside the catchment: {sorted(stray)}")
    for col in ("resistant", "c3r", "mdr", "mrsa", "amrtx_resistant"):
        vals = set(df[col].dropna().astype(int))
        if not vals <= {0, 1}:
            raise ValueError(
                f"{col} must be binary 0/1 (NaN where not ascertained); got {sorted(vals)}"
            )
    yrs = df["sample_year"].dropna().astype(int)
    if len(yrs) and (yrs.min() < 2000 or yrs.max() > 2100):
        raise ValueError(f"Implausible sample_year range: {yrs.min()}-{yrs.max()}")


def load_spidaar_isolates() -> pd.DataFrame:
    """Load the SPIDAAR isolate-level resistance file (pre-reg §5).

    Reads ``data/spidaar/spidaar_isolatedata.xls``, restricts to
    ``config.CATCHMENT_COUNTRIES``, and returns the analysis frame
    (``ISOLATE_COLUMNS``): 244 isolates (2021-2022) with organism, specimen, HAI
    site, clinical relevance, sample timing, and resistance — both the composite
    exposure ``resistant`` (``rx`` = MDR and/or 3GC-R) and the per-mechanism flags.

    This is the **isolate-level** SPIDAAR resource. It carries no patient outcome
    and shares no join key with the patient file (see the module docstring and
    ``docs/deviation_log.md``), so it powers the Step-2 local projection anchor and
    the standalone resistance-burden description — NOT the Step-1 mortality model,
    until a re-linked (``pid``↔``iid``) extract is obtained.
    """
    path = config.DATA_DIR / "spidaar" / "spidaar_isolatedata.xls"
    raw = _read_data_sheet(path)
    isolates = _restrict_to_catchment(raw)
    frame = _build_spidaar_isolate_frame(isolates)
    _validate_spidaar_isolates(frame)
    return frame.reset_index(drop=True)


def load_atlas() -> pd.DataFrame:
    """Load Pfizer ATLAS isolates (~917K; pre-reg §5) for the projection backbone."""
    raise NotImplementedError(_NOT_YET)


def load_smart() -> pd.DataFrame:
    """Load Merck SMART isolates (~300K+; pre-reg §5)."""
    raise NotImplementedError(_NOT_YET)


def load_rd_hub_snapshot() -> pd.DataFrame:
    """Load the locked Global AMR R&D Hub snapshot (public + philanthropic only).

    Requires ``config.RD_HUB_SNAPSHOT_DATE`` to be set; the window is
    ``config.RD_HUB_WINDOW`` (2017-2024). Private-sector R&D is out of scope.
    """
    if config.RD_HUB_SNAPSHOT_DATE is None:
        raise ValueError(
            "config.RD_HUB_SNAPSHOT_DATE is not set. Lock the snapshot date before "
            "loading R&D Hub data (pre-reg §6/§7; Hub data is retrospectively revised)."
        )
    raise NotImplementedError(_NOT_YET)
