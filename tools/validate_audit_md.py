#!/usr/bin/env python3
"""
tools/validate_audit_md.py -- schema validator for LEXICON-AUDIT-*.md.

Enforces the 9 rules in .planning/phases/06.1-*/06.1-RESEARCH.md
Section "Validation Architecture > Schema validator":

    1. Table has 9 columns matching the D-10 schema exactly.
    2. Row count within [80, 250] (widened from D-11's 100-150
       because total reachable is ~230 per RESEARCH Example 5).
    3. lexicon_path column matches regex `\\d+;(?:\\d+;)*`.
    4. frequency column values in {high, med, low, "--"}.
    5. qaction column values in {yes, no}.
    6. toolbar column values in {yes, no}.
    7. Exactly 5 rows have toolbar=yes (D-06 top-5 toolbar slice).
    8. Rows with qaction=yes have frequency in {high, med} (D-05/D-11).
    9. icon_source cells ending in .svg have a matching file under
       xvue/qt/resources/icons/mail/ (Plan 02 deliverable).
       WARN-only by default; set AUDIT_STRICT_ICONS=1 to make FAIL.

Usage:
    python3 tools/validate_audit_md.py <path-to-audit.md>

Exit codes:
    0 = all rules pass (WARN messages allowed for rule 9)
    1 = schema violation (first violation printed to stderr)
"""
import os
import re
import sys
from pathlib import Path

EXPECTED_COLS = [
    "lexicon_path", "description_fr", "description_en", "frequency",
    "qaction", "toolbar", "icon_source", "shortcut", "notes",
]
FREQ_OK = {"high", "med", "low", "\u2014"}  # em-dash
BOOL_OK = {"yes", "no"}
LEX_RX = re.compile(r"^\d+;(?:\d+;)*$")
ROW_BOUND = (80, 250)
TOOLBAR_YES_EXACT = 5
ICONS_ROOT = (
    Path(__file__).resolve().parent.parent
    / "xvue" / "qt" / "resources" / "icons"
)


def resolve_icon_path(audit_path, icon_filename):
    """Resolve an SVG filename against the icons/ tree.

    1. Deduce the owning module from the audit filename stem
       (e.g. "LEXICON-AUDIT-mail" -> "mail"). Try
       xvue/qt/resources/icons/<mod>/<filename> first.
    2. If not present (e.g. icon was shipped by a different module's
       phase and reused -- see 6.2 reusing mesh-draw.svg from mail/),
       scan every module subdirectory under xvue/qt/resources/icons/
       and return the first match.
    3. Return the owning-module path even if missing, so the error
       message points at where the icon SHOULD have been shipped.
    """
    stem = audit_path.stem  # e.g. "LEXICON-AUDIT-mail"
    m = re.match(r"LEXICON-AUDIT-(\w+)$", stem)
    owning_mod = m.group(1) if m else None
    if owning_mod:
        p = ICONS_ROOT / owning_mod / icon_filename
        if p.exists():
            return p
    if ICONS_ROOT.exists():
        for mod_dir in ICONS_ROOT.iterdir():
            if mod_dir.is_dir():
                p = mod_dir / icon_filename
                if p.exists():
                    return p
    return ICONS_ROOT / (owning_mod or "unknown") / icon_filename


def fail(msg):
    sys.stderr.write("validate_audit_md: %s\n" % msg)
    sys.exit(1)


def rows(path):
    """Yield (header, row) tuples for the first markdown table in the file."""
    in_table = False
    header = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("|") and not in_table:
            header = [c.strip() for c in line.strip("|").split("|")]
            in_table = True
            continue
        if in_table:
            if not line.startswith("|"):
                break
            if set(line.replace("|", "").strip()) <= set("-: "):
                continue
            yield header, [c.strip() for c in line.strip("|").split("|")]


def main():
    if len(sys.argv) != 2:
        fail("usage: validate_audit_md.py <audit.md>")
    p = Path(sys.argv[1])
    if not p.is_file():
        fail("not a file: %s" % p)

    header = None
    data = []
    for h, r in rows(p):
        header = h
        data.append(r)

    # Rule 1 -- 9 columns, exact schema
    if header != EXPECTED_COLS:
        fail("rule 1 (9-col schema) failed: got %r" % header)

    # Rule 2 -- row count bound
    lo, hi = ROW_BOUND
    if not (lo <= len(data) <= hi):
        fail("rule 2 (row count in [%d,%d]) failed: got %d" % (lo, hi, len(data)))

    # Rules 0/3/4/5/6 -- per-row enum and regex checks
    for i, row in enumerate(data, start=1):
        if len(row) != 9:
            fail("rule 0 (9 cols per row) failed at row %d: %r" % (i, row))
        lex, df, en, freq, qa, tb, icon, sc, notes = row
        if not LEX_RX.match(lex):
            fail("rule 3 (lexicon_path regex) failed at row %d: %r" % (i, lex))
        if freq not in FREQ_OK:
            fail("rule 4 (frequency enum) failed at row %d: %r" % (i, freq))
        if qa not in BOOL_OK:
            fail("rule 5 (qaction enum) failed at row %d: %r" % (i, qa))
        if tb not in BOOL_OK:
            fail("rule 6 (toolbar enum) failed at row %d: %r" % (i, tb))

    # Rule 7 -- exactly 5 toolbar=yes
    n_tb = sum(1 for r in data if r[5] == "yes")
    if n_tb != TOOLBAR_YES_EXACT:
        fail("rule 7 (toolbar=yes == 5) failed: got %d" % n_tb)

    # Rule 8 -- qaction=yes implies frequency in {high, med}
    for i, r in enumerate(data, start=1):
        if r[4] == "yes" and r[3] not in {"high", "med"}:
            fail(
                "rule 8 (qaction=yes => freq in high/med) failed at row %d: freq=%r"
                % (i, r[3])
            )

    # Rule 9 -- icon_source cells ending .svg have a matching file under
    # xvue/qt/resources/icons/mail/ (cross-check with Plan 02 deliverable).
    # Plan 01 runs BEFORE Plan 02 commits its icons; allow missing files
    # as WARN by default. Toggle via env var AUDIT_STRICT_ICONS=1.
    strict_icons = os.environ.get("AUDIT_STRICT_ICONS", "0") == "1"
    for i, r in enumerate(data, start=1):
        icon = r[6]
        if icon.endswith(".svg"):
            f = resolve_icon_path(p, icon)
            if not f.exists():
                msg = "rule 9 (svg exists) failed at row %d: %s" % (i, f)
                if strict_icons:
                    fail(msg)
                else:
                    sys.stderr.write("WARN: %s\n" % msg)

    print(
        "OK -- %d rows, %d toolbar=yes, schema valid" % (len(data), n_tb)
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
