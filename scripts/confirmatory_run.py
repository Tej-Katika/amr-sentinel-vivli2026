"""Self-checking confirmatory run of the pivoted AMR Sentinel pipeline.

This is the artifact re-executed in the Vivli secure environment. It:

  1. runs ``pipeline.run()`` on the delivered data,
  2. flattens the nested result dict to dotted-key leaves,
  3. diffs every headline value against ``scripts/headline_manifest.json``
     (the expected development-run values transcribed from the report), and
  4. prints a PASS / DRIFT table and saves a JSON record.

Exit status is 0 when every entry is within tolerance, 1 otherwise -- so the
secure-environment operator gets an immediate, auditable yes/no on whether the
report's headline numbers reproduce on the real files, with each drift localised
to a single report line.

Usage:

    PYTHONPATH=src python scripts/confirmatory_run.py [output.json]

Drift on entries flagged ``secure_env_upgrade`` in the manifest is expected
(e.g. the Bayesian Laplace approximation becomes full NUTS; the unfetched $113M
Gram-negative funding split is dropped in). Those rows are reported as
``DRIFT*`` and do NOT, on their own, set a non-zero exit status -- reconcile
them by hand against the report. Any *unflagged* drift is a genuine failure.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from amr_sentinel_vivli import pipeline

# Some leaf strings carry non-ASCII (e.g. the "prior<->posterior" arrow U+2194);
# force UTF-8 so the run does not die on Windows' default cp1252 console encoding.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

_MANIFEST = Path(__file__).with_name("headline_manifest.json")
_DEFAULT_OUTPUT = Path("confirmatory_results.json")


def flatten(obj: object, prefix: str = "") -> dict[str, object]:
    """Flatten nested dicts (and lists-of-dicts) to dotted-key scalar leaves.

    Dicts descend by key; a list/tuple whose elements are *all* dicts descends
    by integer index (so ``power_simulation.by_effect.3.power`` resolves). Every
    other value -- scalars, numeric lists (confidence intervals), string lists
    (rankings), ndarrays, DataFrames -- is kept as a leaf so it can be compared
    whole.
    """
    out: dict[str, object] = {}
    if isinstance(obj, dict):
        for key, val in obj.items():
            child = f"{prefix}.{key}" if prefix else str(key)
            out.update(flatten(val, child))
    elif isinstance(obj, list | tuple) and obj and all(isinstance(x, dict) for x in obj):
        for idx, val in enumerate(obj):
            child = f"{prefix}.{idx}" if prefix else str(idx)
            out.update(flatten(val, child))
    else:
        out[prefix] = obj
    return out


def _as_float_list(value: object) -> list[float] | None:
    if isinstance(value, np.ndarray | list | tuple):
        try:
            return [float(x) for x in value]
        except (TypeError, ValueError):
            return None
    return None


def compare(expected: object, actual: object, tol: float) -> tuple[bool, object]:
    """Return (within_tolerance, diff). ``diff`` is a max abs diff, 0/1, or a tag."""
    if actual is _MISSING:
        return False, "MISSING"

    # String (or list-of-strings) expectations -> exact match.
    if isinstance(expected, str):
        return (str(actual) == expected), (0.0 if str(actual) == expected else "str")
    if isinstance(expected, list) and all(isinstance(e, str) for e in expected):
        actual_list = list(actual) if isinstance(actual, list | tuple) else None
        ok = actual_list == list(expected)
        return ok, (0.0 if ok else "list")

    # Numeric list (confidence interval) -> element-wise abs tolerance.
    if isinstance(expected, list):
        actual_list = _as_float_list(actual)
        if actual_list is None or len(actual_list) != len(expected):
            return False, "shape"
        diffs = [abs(float(e) - a) for e, a in zip(expected, actual_list, strict=True)]
        worst = max(diffs)
        return (worst <= tol), worst

    # Scalar numeric.
    try:
        diff = abs(float(expected) - float(actual))
    except (TypeError, ValueError):
        return False, "type"
    return (diff <= tol), diff


_MISSING = object()


def _lookup(flat: dict[str, object], path: str) -> object:
    return flat.get(path, _MISSING)


def _fmt(value: object) -> str:
    if value is _MISSING:
        return "<missing>"
    if isinstance(value, float | np.floating):
        return f"{float(value):.4g}"
    if isinstance(value, np.ndarray | list | tuple):
        floats = _as_float_list(value)
        if floats is not None:
            return "[" + ", ".join(f"{x:.4g}" for x in floats) + "]"
        return repr(list(value))
    return repr(value)


def main(argv: list[str]) -> int:
    output_path = Path(argv[1]) if len(argv) > 1 else _DEFAULT_OUTPUT

    manifest = json.loads(_MANIFEST.read_text(encoding="utf-8"))
    entries = manifest["entries"]

    print(f"Running pipeline.run() and checking {len(entries)} headline values "
          f"against {_MANIFEST.name} ...\n")
    flat = flatten(pipeline.run())

    rows: list[dict[str, object]] = []
    hard_fail = 0
    soft_drift = 0
    for entry in entries:
        path = entry["path"]
        expected = entry["expected"]
        tol = float(entry.get("tol", 0.0))
        upgrade = bool(entry.get("secure_env_upgrade", False))
        actual = _lookup(flat, path)
        ok, diff = compare(expected, actual, tol)

        if ok:
            status = "PASS"
        elif upgrade:
            status = "DRIFT*"
            soft_drift += 1
        else:
            status = "DRIFT"
            hard_fail += 1

        rows.append({
            "path": path,
            "ref": entry.get("ref", ""),
            "report": entry.get("report", ""),
            "expected": expected,
            "actual": None if actual is _MISSING else _jsonable(actual),
            "tol": tol,
            "diff": diff if isinstance(diff, str) else round(float(diff), 6),
            "status": status,
            "secure_env_upgrade": upgrade,
            "note": entry.get("note", ""),
        })

    _print_table(rows)

    total = len(rows)
    n_pass = sum(1 for r in rows if r["status"] == "PASS")
    print(f"\n{n_pass}/{total} PASS · {soft_drift} expected-drift (DRIFT*) · "
          f"{hard_fail} unexpected DRIFT")

    record = {
        "manifest": _MANIFEST.name,
        "master_seed": manifest["_meta"].get("master_seed"),
        "summary": {"total": total, "pass": n_pass,
                    "expected_drift": soft_drift, "unexpected_drift": hard_fail},
        "rows": rows,
    }
    output_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(f"Saved {output_path}")

    if hard_fail:
        print("\nFAIL: unexpected drift in headline numbers -- reconcile the "
              "report (stale figure) or the code (genuine change) before lifting "
              "the DRAFT marker.")
        return 1
    if soft_drift:
        print("\nOK with expected drift: only secure-env-upgrade rows (DRIFT*) "
              "moved; reconcile those by hand against the report.")
    else:
        print("\nOK: every headline number reproduces within tolerance.")
    return 0


def _jsonable(value: object) -> object:
    if isinstance(value, float | np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    floats = _as_float_list(value)
    if floats is not None:
        return floats
    if isinstance(value, list | tuple):
        return [str(x) for x in value]
    return str(value)


def _print_table(rows: list[dict[str, object]]) -> None:
    status_w = 6
    path_w = min(max((len(str(r["path"])) for r in rows), default=4), 60)
    header = (f"{'STATUS':<{status_w}}  {'PATH':<{path_w}}  "
              f"{'EXPECTED':>12}  {'ACTUAL':>12}  REPORT / REF")
    print(header)
    print("-" * len(header))
    for r in rows:
        exp = _fmt(r["expected"])
        act = _fmt(_MISSING if r["actual"] is None else r["actual"])
        ref = f"{r['report']}  ({r['ref']})" if r["report"] else str(r["ref"])
        line = (f"{r['status']:<{status_w}}  {str(r['path']):<{path_w}}  "
                f"{exp:>12}  {act:>12}  {ref}")
        print(line)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
