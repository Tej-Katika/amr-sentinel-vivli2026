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
#   0  -> all isolates susceptible (confirmed)          -> susceptible (0)
#   1  -> mixed S, untested (no resistant isolate seen) -> NOT confirmed -> NaN
#  -1  -> no susceptibility result / no isolate         -> unascertainable -> NaN
# NOTE (deviation 2026-06-06): ``amrp == 1`` was previously mapped to susceptible
# (0); it is now excluded from the resistant-vs-susceptible contrast (-> NaN),
# because "mixed, S-untested" is not a confirmed-susceptible patient and lumping it
# into the comparison arm makes the causal contrast indefensible. Susceptible is now
# ``amrp == 0`` ONLY. See docs/deviation_log.md and docs/analysis_plan_2026.md.
# CODEBOOK-CONFIRMED (2026-06-08): the official patient codebook defines ``amrp``
# ("Drug susceptibility result to administered AB class, patient level") as
# -1 = No RX result, 0 = All isolates S, 1 = Mixed S-untested, 2 = >=1 R isolate —
# matching this mapping exactly.
_RESISTANT_FROM_AMRP = {2: 1.0, 0: 0.0}

# Isolate-file binary resistance flags (``c3r``, ``mdr``, ``mrsa``, ``rx``):
#   1 -> resistant / positive, 0 -> susceptible / negative. Any other code is a
# sentinel for "not tested / not applicable" (e.g. ``9`` for 3GC-R on a
# Gram-positive, ``99`` for MRSA on a non-staphylococcus) and maps to NaN.
_RESISTANT_BINARY = {1: 1.0, 0: 0.0}

# Mortality status for analysis (``dead``): 0 alive, 1 deceased, 9 deceased-censored.
_DEAD_DECEASED_CODES = (1, 9)
_DEAD_VALID_CODES = (0, 1, 9)

# Analysis-frame columns returned by load_spidaar (no ``year`` — see module docstring).
# Includes the excess-length-of-stay (Component 1) fields surfaced 2026-06-06:
# ``los`` (discharge LOS, the competing-risks event time), ``severity``/``ward``/
# ``days_to_enrolment`` (confounders / left-truncation), and ``treatment_adequacy``
# (raw ``txadp`` for the Step-5 g-formula; coding verified against the codebook there).
SPIDAAR_COLUMNS = (
    "pid", "country", "organism", "infection_site", "age", "age_group", "sex",
    "ward", "severity", "days_to_enrolment",
    "resistant", "amrp", "treatment_adequacy",
    "mortality_30d", "time_at_risk",
    "los", "dead", "days_to_death", "days_observed",
)

# --- ATLAS (Pfizer) panel for the Step-2 catchment nowcast -------------------
# The delivered ATLAS export stores, per antibiotic, a MIC column (e.g. ``Ceftazidime``)
# and an interpretation column (e.g. ``Ceftazidime_I`` -> Susceptible/Intermediate/
# Resistant). Verified in the delivered file (2026-06-07): the catchment is 1,519 of
# 1,011,168 isolates; Enterobacterales = E. coli + K. pneumoniae (665, all with a
# ceftazidime interpretation); **ceftriaxone interpretation is blank in the catchment**,
# so ceftazidime is the only viable 3GC panel cell. Catchment years: Ghana/Malawi/Uganda
# 2021-2023, Kenya 2013/2014 + 2021-2023, none after 2023 (deviation log 2026-06-07).
ATLAS_PANEL_SPECIES: tuple[str, ...] = ("Escherichia coli", "Klebsiella pneumoniae")
ATLAS_PANEL_DRUG = "Ceftazidime"
# Resistance = Resistant vs ascertained (S+I+R); Intermediate counts as non-resistant.
# Blank/other interpretations -> NaN (not ascertained, excluded from the denominator).
_ATLAS_INTERP_TO_RESISTANT = {"Resistant": 1.0, "Susceptible": 0.0, "Intermediate": 0.0}
ATLAS_COLUMNS = (
    "isolate_id", "species", "country", "year", "source", "speciality",
    "drug", "mic", "interpretation", "resistant",
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
    dtpta, nobsd, los, disev, ward, enrtpt, txadp`` and returns the columns in
    ``SPIDAAR_COLUMNS``.

    Outcome (``mortality_30d`` / ``time_at_risk``): a 30-day-horizon survival pair.
    A death observed on/before day 30 (``dead`` in {1, 9} with days-to-death
    ``dtpta`` ≤ 30) is the event; later deaths and survivors are censored at
    ``min(follow-up, 30)`` days, where follow-up is ``dtpta`` for deaths and the
    observation length ``nobsd`` for survivors.

    Excess-LOS fields (Component 1, surfaced 2026-06-06): ``los`` (discharge
    length-of-stay = the competing-risks event time; NaN for patients who died or
    were not discharged), ``severity`` (``disev``), ``ward``, ``days_to_enrolment``
    (``enrtpt``), and raw ``treatment_adequacy`` (``txadp``).

    Exposure (``resistant``): from ``amrp`` per ``_RESISTANT_FROM_AMRP`` — resistant
    is ``amrp == 2``, susceptible is ``amrp == 0`` ONLY. ``amrp`` in {1, -1} maps to
    NaN (mixed-untested / unascertainable) and is excluded from the
    resistant-vs-susceptible contrast (deviation 2026-06-06).
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

    # Excess-LOS (Component 1) fields: discharge length-of-stay, severity, ward, and
    # the admission->enrolment day-count (enables the left-truncation sensitivity);
    # plus raw treatment-adequacy ``txadp`` for the Step-5 g-formula.
    los = pd.to_numeric(raw["los"], errors="coerce")
    disev = pd.to_numeric(raw["disev"], errors="coerce")
    ward = pd.to_numeric(raw["ward"], errors="coerce")
    enrtpt = pd.to_numeric(raw["enrtpt"], errors="coerce")
    txadp = pd.to_numeric(raw["txadp"], errors="coerce")

    out = pd.DataFrame({
        "pid": raw["pid"],
        "country": raw["ctry"],
        "organism": raw["isol"],
        "infection_site": chaicat.map(_INFECTION_SITE_LABELS),
        "age": agegr.astype("Int64"),
        "age_group": agegr.map(_AGE_GROUP_LABELS),
        "sex": sex_code.map(_SEX_LABELS),
        "ward": ward.astype("Int64"),
        "severity": disev.astype("Int64"),
        "days_to_enrolment": enrtpt.astype("Int64"),
        "resistant": amrp.map(_RESISTANT_FROM_AMRP).astype("Int64"),
        "amrp": amrp.astype("Int64"),
        "treatment_adequacy": txadp.astype("Int64"),
        "mortality_30d": event.astype(int),
        "time_at_risk": time_at_risk,
        "los": los,
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


def _build_atlas_frame(raw: pd.DataFrame, drug: str = ATLAS_PANEL_DRUG) -> pd.DataFrame:
    """Map a raw ATLAS export to the tidy panel frame for one antibiotic.

    Pure function (no I/O) so it can be unit-tested on synthetic rows. Expects the raw
    ATLAS columns ``Isolate Id, Species, Country, Year, Source, Speciality`` plus the
    drug's MIC column (``<drug>``) and interpretation column (``<drug>_I``), and returns
    ``ATLAS_COLUMNS``. ``resistant`` is 1 where the interpretation is "Resistant", 0 for
    "Susceptible"/"Intermediate", and NaN where the interpretation is missing (excluded
    from the resistance denominator).

    Resistance is taken from the ATLAS-supplied interpretation (complete for the catchment
    Enterobacterales × ceftazidime cell). Re-deriving R from the raw MIC under EUCAST v15.0
    (``config.PRIMARY_BREAKPOINT_REGIME``) is a documented secure-environment refinement
    (Gate X — MIC string format + breakpoint table to be confirmed); the MIC is carried
    through (``mic``) so that step can run without re-reading the source.
    """
    raw = raw.reset_index(drop=True)
    interp = raw[f"{drug}_I"]
    out = pd.DataFrame({
        "isolate_id": raw["Isolate Id"],
        "species": raw["Species"],
        "country": raw["Country"],
        "year": pd.to_numeric(raw["Year"], errors="coerce").astype("Int64"),
        "source": raw["Source"],
        "speciality": raw["Speciality"],
        "drug": drug,
        "mic": raw[drug],
        "interpretation": interp,
        "resistant": interp.map(_ATLAS_INTERP_TO_RESISTANT).astype("Int64"),
    })
    return out[list(ATLAS_COLUMNS)]


def _validate_atlas(df: pd.DataFrame, catchment_only: bool) -> None:
    """Fail fast on the invariants the Step-2 nowcast relies on."""
    if catchment_only:
        stray = set(df["country"].dropna()) - set(config.CATCHMENT_COUNTRIES)
        if stray:
            raise ValueError(f"ATLAS rows outside the catchment: {sorted(stray)}")
    vals = set(df["resistant"].dropna().astype(int))
    if not vals <= {0, 1}:
        raise ValueError(
            f"resistant must be binary 0/1 (NaN where not ascertained); got {sorted(vals)}")
    yrs = df["year"].dropna().astype(int)
    if len(yrs) and (yrs.min() < 2000 or yrs.max() > 2100):
        raise ValueError(f"Implausible ATLAS year range: {yrs.min()}-{yrs.max()}")


def load_atlas(
    drug: str = ATLAS_PANEL_DRUG,
    species: tuple[str, ...] | None = ATLAS_PANEL_SPECIES,
    catchment_only: bool = True,
) -> pd.DataFrame:
    """Load Pfizer ATLAS isolates (pre-reg §5) for the Step-2 catchment nowcast.

    Reads ``data/atlas/atlas_vivli_2004_2024.csv`` (only the columns needed for ``drug``),
    optionally restricts to ``species`` (default the Enterobacterales panel carriers) and
    to ``config.CATCHMENT_COUNTRIES`` (default True), and returns the tidy ``ATLAS_COLUMNS``
    frame with a binary ceftazidime ``resistant`` flag. With the defaults this is the ~665
    catchment Enterobacterales × ceftazidime isolates the nowcast and frame-contrast use;
    pass ``catchment_only=False`` / ``species=None`` for the full ATLAS backbone.

    SMART is excluded from the Challenge (deviation log 2026-06-07); ATLAS is the only
    external isolate-surveillance source.
    """
    path = config.DATA_DIR / "atlas" / "atlas_vivli_2004_2024.csv"
    usecols = ["Isolate Id", "Species", "Country", "Year", "Source", "Speciality",
               drug, f"{drug}_I"]
    raw = pd.read_csv(path, usecols=usecols, low_memory=False)
    if species is not None:
        raw = raw[raw["Species"].isin(species)]
    if catchment_only:
        raw = _restrict_to_catchment(raw, country_col="Country")
    frame = _build_atlas_frame(raw, drug)
    _validate_atlas(frame, catchment_only)
    return frame.reset_index(drop=True)


def load_atlas_backbone() -> pd.DataFrame:
    """Full ATLAS Species/Country/Year backbone (every isolate) for Component 6.

    Unlike :func:`load_atlas` (one drug cell, catchment-only), this reads the whole register
    but only the three columns the surveillance blind-spot analysis needs — so it is cheap
    despite the ~1M rows. Returns a tidy frame with ``Species``, ``Country`` and integer
    ``Year`` (rows with an unparseable year are dropped). No resistance columns.
    """
    path = config.DATA_DIR / "atlas" / "atlas_vivli_2004_2024.csv"
    raw = pd.read_csv(path, usecols=["Species", "Country", "Year"], low_memory=False)
    raw = raw.assign(Year=pd.to_numeric(raw["Year"], errors="coerce"))
    raw = raw.dropna(subset=["Species", "Country"])
    if raw.empty:
        raise ValueError("ATLAS backbone is empty after dropping rows without species/country.")
    return raw.reset_index(drop=True)


def load_smart() -> pd.DataFrame:
    """Load Merck SMART isolates (~300K+; pre-reg §5)."""
    raise NotImplementedError(_NOT_YET)


# --- Global AMR R&D Hub investment export -----------------------------------
# A dated export pulled by hand from the Hub Dynamic Dashboard / Investment Gallery
# (globalamrhub.org/dynamic-dashboard/) inside the secure environment and saved to
# ``data/rd_hub/rd_hub_investment_<RD_HUB_SNAPSHOT_DATE>.csv``. The Hub User Guide
# (dashboard *Library* > Data Sources) sanctions "extract and use"; we freeze a dated
# copy because the live dashboard is retrospectively revised (pre-reg §6/§7). The
# published Czaplewski et al. extract (``rd_alignment.RD_HUB_SNAPSHOT_2026``) is kept
# as the reconciliation cross-check, NOT replaced — see
# ``rd_alignment.reconcile_hub_snapshot``.
#
# Expected tidy schema (one row per pathogen x year x funder slice). Column names are
# matched case-insensitively and stripped; only the four required columns must exist.
RD_HUB_EXPORT_COLUMNS = ("pathogen", "year", "funder_type", "investment_usd")
_RD_HUB_OPTIONAL_COLUMNS = ("sector", "funder")

# Funder types kept under config.RD_HUB_SCOPE ("public + philanthropic only"). Anything
# else (private, public-private partnership, unspecified) is dropped — documented, not
# silently mixed in.
_RD_HUB_KEPT_FUNDER_TYPES = frozenset({"public", "philanthropic"})

# Pathogen labels that denote NON-species-specific (cross-cutting) investment. Folded
# into the cross-cutting bucket rather than any single pathogen.
_RD_HUB_CROSS_CUTTING_LABELS = frozenset({
    "cross-cutting", "cross cutting", "crosscutting", "multiple", "pathogen-agnostic",
    "pathogen agnostic", "not pathogen specific", "not pathogen-specific", "unspecified",
    "none", "n/a", "na", "", "various", "broad-spectrum", "broad spectrum",
})


def _build_rd_hub_frame(raw: pd.DataFrame) -> pd.DataFrame:
    """Pure raw->tidy transform for a Hub investment export (unit-tested on synthetic rows).

    Normalises column names, restricts to the in-scope funder types and (if present) the
    human sector, clips to ``config.RD_HUB_WINDOW``, coerces the amount to float, drops
    non-positive/unparseable amounts, and flags cross-cutting rows. Returns a tidy frame
    with columns ``pathogen, year, funder_type, investment_usd, is_cross_cutting`` (+
    ``funder`` if the export carried it). No side effects; never touches disk.
    """
    rename = {c: c.strip().lower().replace(" ", "_") for c in raw.columns}
    frame = raw.rename(columns=rename)

    missing = [c for c in RD_HUB_EXPORT_COLUMNS if c not in frame.columns]
    if missing:
        raise ValueError(
            f"Hub export is missing required column(s) {missing}; expected at least "
            f"{list(RD_HUB_EXPORT_COLUMNS)} (got {list(frame.columns)})."
        )

    keep = list(RD_HUB_EXPORT_COLUMNS) + [
        c for c in _RD_HUB_OPTIONAL_COLUMNS if c in frame.columns
    ]
    frame = frame[keep].copy()

    frame["pathogen"] = frame["pathogen"].astype("string").str.strip()
    frame["funder_type"] = frame["funder_type"].astype("string").str.strip().str.lower()
    frame["year"] = pd.to_numeric(frame["year"], errors="coerce").astype("Int64")
    frame["investment_usd"] = pd.to_numeric(frame["investment_usd"], errors="coerce")

    # Scope filters (each documented in config / deviation log).
    frame = frame[frame["funder_type"].isin(_RD_HUB_KEPT_FUNDER_TYPES)]
    if "sector" in frame.columns:
        frame = frame[frame["sector"].astype("string").str.strip().str.lower() == "human"]
    lo, hi = config.RD_HUB_WINDOW
    frame = frame[(frame["year"] >= lo) & (frame["year"] <= hi)]

    # Drop rows that carry no usable amount.
    frame = frame[frame["investment_usd"].notna() & (frame["investment_usd"] > 0)]

    label = frame["pathogen"].fillna("").str.lower()
    frame["is_cross_cutting"] = label.isin(_RD_HUB_CROSS_CUTTING_LABELS)

    frame["year"] = frame["year"].astype(int)
    frame["investment_usd"] = frame["investment_usd"].astype(float)
    return frame.reset_index(drop=True)


def _validate_rd_hub(frame: pd.DataFrame) -> None:
    """Guard the loaded Hub frame: non-empty, in-window, in-scope, positive amounts."""
    if frame.empty:
        raise ValueError(
            "Hub export is empty after applying the scope/window filters. Check that the "
            "export covers public+philanthropic funders within "
            f"{config.RD_HUB_WINDOW} (human sector)."
        )
    lo, hi = config.RD_HUB_WINDOW
    if not frame["year"].between(lo, hi).all():
        raise ValueError(f"Hub export has rows outside the locked window {config.RD_HUB_WINDOW}.")
    if not frame["funder_type"].isin(_RD_HUB_KEPT_FUNDER_TYPES).all():
        raise ValueError("Hub export retained an out-of-scope funder type.")
    if not (frame["investment_usd"] > 0).all():
        raise ValueError("Hub export retained a non-positive investment amount.")


def load_rd_hub_snapshot() -> pd.DataFrame:
    """Load the dated Global AMR R&D Hub investment export (public + philanthropic only).

    Reads ``data/rd_hub/rd_hub_investment_<RD_HUB_SNAPSHOT_DATE>.csv`` — a dated export
    hand-pulled from the Hub Dynamic Dashboard in the secure environment (see the module
    note above for provenance). Requires ``config.RD_HUB_SNAPSHOT_DATE`` to be set; the
    window is ``config.RD_HUB_WINDOW`` and private-sector R&D is out of scope.

    Returns a tidy frame (one row per pathogen x year x funder slice). Collapse it into the
    index denominator with :func:`rd_alignment.snapshot_from_hub_export`, and cross-check it
    against the published Czaplewski snapshot with
    :func:`rd_alignment.reconcile_hub_snapshot`.
    """
    if config.RD_HUB_SNAPSHOT_DATE is None:
        raise ValueError(
            "config.RD_HUB_SNAPSHOT_DATE is not set. Lock the snapshot date before "
            "loading R&D Hub data (pre-reg §6/§7; Hub data is retrospectively revised)."
        )
    path = config.DATA_DIR / "rd_hub" / f"rd_hub_investment_{config.RD_HUB_SNAPSHOT_DATE}.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"Hub investment export not found at {path}. Pull a dated export from the Hub "
            "Dynamic Dashboard (filters: public+philanthropic funders, "
            f"{config.RD_HUB_WINDOW[0]}-{config.RD_HUB_WINDOW[1]}, human sector) inside the "
            "secure environment and save it there. Columns expected: "
            f"{list(RD_HUB_EXPORT_COLUMNS)} (+ optional {list(_RD_HUB_OPTIONAL_COLUMNS)})."
        )
    raw = pd.read_csv(path)
    frame = _build_rd_hub_frame(raw)
    _validate_rd_hub(frame)
    return frame
