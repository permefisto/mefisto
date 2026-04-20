// xvue/qt/src/xvue_qt_menu_file_parser.h — Phase 6.1 Plan 02
// Bilingual menu-label parser for LIMTCL files under $MEFISTO/td/m/.
// Single source of truth for QAction labels is the Fortran menu file
// tree (CONTEXT D-12). Language resolution follows util/langue.f —
// $MEFISTO/td/m/anglais present => EN copy has been installed into
// td/m/; absent => FR copy. This parser ALWAYS reads td/m/<name>.
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
    // Load $MEFISTO/td/m/<name>. Cache is process-wide static (one
    // XvueMenuFileParser session is sufficient per RESEARCH Open
    // Question #4). On parse failure, returns a MenuFile whose
    // .ok()==false and whose .label(N) yields "N;".
    static const MenuFile& loadFor(const QString& name);

    // Test-only helper — clears the cache so per-test state is fresh.
    // Guarded by XVUE_QT_TESTING; not called from production paths.
    static void clearCacheForTesting();
};
