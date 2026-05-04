#!/bin/bash
# bin/test_no_imagemagick_in_qt.sh — Phase 7 Plan 06 (EXPORT-06)
#
# CI-style assertion: the Qt path must NEVER reference ImageMagick.
#
# Per CONTEXT.md D-16/D-17 the grep is scoped to xvue/qt/ ONLY.
# bin/convertepsgif, xvue/video1.f, xvue/videofin.f, xvue/videonm.f
# (the LVIDEO + legacy-postprocess pipeline) stay alive until Phase 9
# RETIRE-03 — those files MUST NOT be matched here.
#
# Exit codes:
#   0  — gate passes (no ImageMagick references in xvue/qt/)
#   1  — gate fails (one or more references found; details printed)
#   2  — environment misconfigured ($MEFISTO unset)
#
# Allowlist rationale: the literal "convert" appears in many legitimate
# Qt API names (QString::convertTo*, QPageSize convertToOther, etc.)
# whose substring would otherwise trigger a false positive. The grep
# uses a negative-filter pipeline (\b word boundaries + a Qt-API
# exclusion list) to avoid false positives without weakening the gate
# against any genuine ImageMagick usage.

set -e

if [ -z "$MEFISTO" ]; then
    # Fall back to walking up from the script location to find the
    # repo root — useful when the script is invoked from CMake with
    # ${CMAKE_SOURCE_DIR}/../../bin/test_no_imagemagick_in_qt.sh and
    # MEFISTO is not in the environment.
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
    REPO_ROOT="$(dirname "$SCRIPT_DIR")"
    if [ -d "$REPO_ROOT/xvue/qt" ]; then
        MEFISTO="$REPO_ROOT"
    else
        echo "test_no_imagemagick_in_qt.sh: \$MEFISTO is unset and "
        echo "the script could not auto-detect the repo root from "
        echo "$SCRIPT_DIR. Export MEFISTO=/path/to/mefistosource."
        exit 2
    fi
fi

cd "$MEFISTO"

# The literal "convert" word, "ImageMagick" word, or the ".magick"
# domain — none of these may appear under xvue/qt/.
#
# Search scope:
#   - C/C++ sources, headers, CMake files, shell scripts, markdown,
#     and Qt resource XML (.qrc). build/ is excluded by --include
#     (only matching extensions are searched, not directories).
#
# Allowlist (legitimate Qt-API token usage):
#   - QPageSize          (PDF page-size enum tokens)
#   - convertToOther     (QImage::convertToFormat overloads etc.)
#   - convertTo(         (QString::convertTo*, QImage::convertTo*)
#   - QString::convertTo (defensive duplicate of the above)
#   - convertFromUtf     (QString::convertFromUtf*)
#   - Qt::ConvertibleTo  (Qt::ConvertibleToFloat / Bool / Int)
#
# The \b word-boundary anchor on the matcher itself + the explicit
# allowlist keeps the gate precise.
hits=$(grep -rEn --include='*.cpp' --include='*.h' --include='*.hpp' \
        --include='*.cmake' --include='CMakeLists.txt' \
        --include='*.sh' --include='*.md' --include='*.qrc' \
        '\b(convert|ImageMagick|magick)\b' xvue/qt/ \
        | grep -vE '(QPageSize|convertToOther|convertTo\(|QString::convertTo|convertFromUtf|Qt::ConvertibleTo)' \
        || true)

if [ -n "$hits" ]; then
    echo "FAIL: ImageMagick references found in xvue/qt/ (forbidden by EXPORT-06):"
    echo "$hits"
    exit 1
fi

echo "EXPORT-06 PASS: no ImageMagick references in xvue/qt/"
exit 0
