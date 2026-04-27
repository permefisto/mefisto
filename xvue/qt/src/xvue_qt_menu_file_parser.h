// xvue/qt/src/xvue_qt_menu_file_parser.h — Phase 6.1 Plan 02
//   updated 06.2 Plan 05 to honour the anglais flag at runtime.
// Bilingual menu-label parser for LIMTCL files under $MEFISTO/td/.
// Single source of truth for QAction labels is the Fortran menu file
// tree (CONTEXT D-12).
//
// Language resolution: at runtime we consult xvueIsEnglish() which
// probes $MEFISTO/td/m/anglais. When the flag is set AND
// $MEFISTO/td/ma/<name> is readable, that EN file is parsed.
// Otherwise we read $MEFISTO/td/m/<name> (legacy contract — also
// serves the FR install where td/m/ is the FR copy).
//
// Rationale (06.2-HUMAN-UAT.md gap 2): the legacy install copies
// ONE language into td/m/ at install time; users who toggle the
// anglais flag mid-life would need to re-run the install to refresh
// td/m/, which they often don't. Reading the canonical td/ma/ tree
// when the flag says EN closes the gap without requiring a re-install.
#pragma once

#include <QString>
#include <QHash>

class MenuFile {
public:
    const QString& title() const { return title_; }

    // Returns the label for numeric code `code`. Fallback is "N;" per
    // CONTEXT D-12 final bullet — e.g. label(80) on a file that never
    // defines code 80 returns "80;" so the QAction still shows SOMETHING
    // meaningful (the typed-lexicon shortcut). This is also the
    // RESEARCH Pitfall 5 log-and-fallback behaviour for the
    // td/m/modifm2d `80;` typo.
    QString label(int code) const;

    bool ok() const { return !labels_.isEmpty(); }

private:
    friend class XvueMenuFileParser;
    QString title_;
    QHash<int, QString> labels_;
};

class XvueMenuFileParser {
public:
    // Load the menu file for `name`. When xvueIsEnglish() is true and
    // $MEFISTO/td/ma/<name> is readable, that EN file is parsed; else
    // we fall back to $MEFISTO/td/m/<name> (legacy/FR-install path).
    // Cache is process-wide static (one XvueMenuFileParser session is
    // sufficient per RESEARCH Open Question #4). On parse failure,
    // returns a MenuFile whose .ok()==false and whose .label(N)
    // yields "N;".
    static const MenuFile& loadFor(const QString& name);

    // Test-only helper — clears the cache so per-test state is fresh.
    // Guarded by XVUE_QT_TESTING; not called from production paths.
    static void clearCacheForTesting();
};
