"""Data loading for the Vivli secure research environment.

All inputs are Vivli DUA-restricted and live under ``data/`` (gitignored). These
loaders run *inside* the Vivli secure environment after the DUA is signed. No
data has been accessed at the time this module is committed (pre-reg §15).
"""

from __future__ import annotations

import pandas as pd

from . import config

_NOT_YET = (
    "Pre-data scaffold — no Vivli data has been accessed (pre-reg §15). "
    "Implement loaders only after the DUA is signed and the secure environment entered."
)


def load_spidaar() -> pd.DataFrame:
    """Load the SPIDAAR cohort (~244 isolates with mortality outcome; pre-reg §5).

    Expected columns include the outcome (``config.OUTCOME``), the resistance
    exposure, and the propensity covariates in ``config.PROPENSITY_COVARIATES``.
    Scope is restricted to ``config.CATCHMENT_COUNTRIES``.
    """
    raise NotImplementedError(_NOT_YET)


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
