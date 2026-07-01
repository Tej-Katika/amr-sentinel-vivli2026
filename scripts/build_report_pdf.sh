#!/usr/bin/env bash
# Build the 5-page submission PDF from docs/final_report_2026.md.
#
#   scripts/build_report_pdf.sh
#
# Pipeline: pandoc (markdown -> typst, standalone) then typst (typst -> PDF), with a compact
# layout tuned to the Vivli hard limit (max 5 pages excluding references; no appendices). The
# figures are the licensed-data PNGs under figures/ (run scripts/make_figures.py first); the
# resulting PDF therefore embeds DUA-restricted material and is written to build/ (gitignored).
#
# Requires pandoc >= 3 and typst >= 0.12 on PATH (override with PANDOC=/path TYPST=/path).
set -euo pipefail

PANDOC="${PANDOC:-pandoc}"
TYPST="${TYPST:-typst}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/docs/final_report_2026.md"
OUT_DIR="$ROOT/build"
BODY="$OUT_DIR/report.typ"
PDF="$OUT_DIR/final_report_2026.pdf"
TWEAKS="$OUT_DIR/layout.typ"

mkdir -p "$OUT_DIR"

# Page geometry must go THROUGH pandoc's typst template (its conf() wrapper sets page/text/par
# last, so header-includes can't override them). Tuned to the 5-page limit; still legible.
META="$OUT_DIR/meta.yaml"
cat > "$META" <<'YAML'
papersize: a4
fontsize: 10pt
linestretch: 0.92
margin:
  x: 1.5cm
  y: 1.4cm
YAML

# Show-rules DO cascade to later content, so heading sizes / paragraph spacing / figure width
# go here (injected into the preamble).
cat > "$TWEAKS" <<'TYP'
#show heading.where(level: 1): set text(size: 13pt)
#show heading.where(level: 2): set text(size: 11.5pt)
#show heading.where(level: 3): set text(size: 10.5pt)
#show heading: set block(above: 0.7em, below: 0.3em)
#set par(spacing: 0.62em)
#set image(width: 78%)
TYP

# markdown -> standalone typst (keeps the relative ../figures/*.png paths; no media extraction)
"$PANDOC" "$SRC" -t typst -s -o "$BODY" \
  --metadata-file="$META" \
  --include-in-header="$TWEAKS"

# typst -> PDF, with the project root as the file-access root so ../figures/ resolves
"$TYPST" compile "$BODY" "$PDF" --root "$ROOT"

echo "Wrote $PDF"
if command -v python >/dev/null 2>&1; then
  python - "$PDF" <<'PY' 2>/dev/null || true
import sys
try:
    import pypdf
    n = len(pypdf.PdfReader(sys.argv[1]).pages)
    print(f"Pages: {n}  ({'OK <=5' if n <= 5 else 'OVER the 5-page limit'})")
except Exception:
    pass
PY
fi
