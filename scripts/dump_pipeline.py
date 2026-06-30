"""Run the full pivoted pipeline on the delivered data and dump every leaf value.

A read-only inspection harness: it imports nothing it shouldn't and writes no
files. It exists so a human can see *exactly* what the pipeline produces today
(the "dev numbers") in a flat, labelled, diff-friendly form before wiring up the
self-checking confirmatory run (scripts/confirmatory_run.py).

    PYTHONPATH=src python scripts/dump_pipeline.py

Numbers are rounded for legibility only; the underlying values are unchanged.
"""

from __future__ import annotations

import sys

import numpy as np

from amr_sentinel_vivli import pipeline

# Some leaf strings carry non-ASCII (e.g. the "prior<->posterior" arrow U+2194); force
# UTF-8 so the dump does not die on Windows' default cp1252 console encoding.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def _fmt(value: object) -> str:
    if isinstance(value, float | np.floating):
        return f"{float(value):.4g}"
    if isinstance(value, np.ndarray | tuple | list):
        try:
            arr = np.asarray(value, dtype=float).ravel()
        except (ValueError, TypeError):
            return repr(value)
        if arr.size <= 8:
            return "[" + ", ".join(f"{x:.4g}" for x in arr) + "]"
        return f"<array shape={np.asarray(value).shape}>"
    return repr(value)


def walk(obj: object, prefix: str = "") -> None:
    if isinstance(obj, dict):
        for key, val in obj.items():
            walk(val, f"{prefix}.{key}" if prefix else str(key))
    else:
        print(f"{prefix:<55} = {_fmt(obj)}")


if __name__ == "__main__":
    result = pipeline.run()
    walk(result)
