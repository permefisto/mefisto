// xvue/qt/src/xvue_qt_menu_file_parser.cpp — Phase 6.1 Plan 02
// Implementation of XvueMenuFileParser — see header for rationale.
//
// Grammar (empirically verified by reading td/m/debut + 7 other files):
//     variable <TAG> '<title>' entier
//     ( <NUM1> : '<label1>' , <NUM2> : '<label2>' , ... ) ;
// - Embedded single-quotes use Fortran doubling: 'd''un' == d'un.
// - Separator is strictly ':' — a rogue ';' (e.g. td/m/modifm2d line 17)
//   means the row is *not* matched and MenuFile::label() returns the
//   lexicon-path fallback "N;". RESEARCH Pitfall 5 log-and-fallback.
//
// Security: T-6.1-PATH-INJECT is mitigated because `name` is populated
// ONLY from hard-coded string literals in xvue_qt_mail_actions.cpp.
// No user-controlled path traversal is possible.
#include "xvue_qt_menu_file_parser.h"
#include "xvue_qt_i18n.h"   // 06.2 Plan 05: xvueIsEnglish() probe.

#include <QByteArray>
#include <QFile>
#include <QRegularExpression>
#include <QTextStream>

namespace {
QHash<QString, MenuFile>& cache() {
    static QHash<QString, MenuFile> c;
    return c;
}
} // namespace

QString MenuFile::label(int code) const {
    auto it = labels_.find(code);
    if (it != labels_.end()) return *it;
    return QStringLiteral("%1;").arg(code);
}

const MenuFile& XvueMenuFileParser::loadFor(const QString& name) {
    auto& c = cache();
    auto it = c.find(name);
    if (it != c.end()) return *it;

    MenuFile mf;

    const QByteArray mefisto = qgetenv("MEFISTO");
    if (mefisto.isEmpty()) {
        // Insert empty MenuFile so we don't retry for every call.
        return *c.insert(name, mf);
    }

    const QString home = QString::fromLocal8Bit(mefisto);

    // 06.2 Plan 05 gap-2 fix: when xvueIsEnglish() is true, prefer
    // $MEFISTO/td/ma/<name> (the canonical EN tree) over the
    // generic td/m/<name> (which often holds whichever language
    // the install copied at setup time, regardless of the flag's
    // current state).
    //
    // The probe is cheap because xvueIsEnglish() caches its
    // result via static-local initialisation in xvue_qt_i18n.cpp.
    QString chosen_path;
    if (xvueIsEnglish()) {
        const QString en_path = home + QStringLiteral("/td/ma/") + name;
        if (QFile::exists(en_path)) {
            chosen_path = en_path;
        }
    }
    if (chosen_path.isEmpty()) {
        chosen_path = home + QStringLiteral("/td/m/") + name;
    }

    QFile f(chosen_path);
    if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return *c.insert(name, mf);
    }

    QTextStream ts(&f);
    QString content = ts.readAll();

    // Title regex — capture the FIRST single-quoted string after
    // `variable <TAG>`. Per RESEARCH §Pattern 2 grammar.
    static const QRegularExpression titleRx(
        QStringLiteral(R"(variable\s+\S+\s+'((?:[^']|'')*)'\s+entier)"),
        QRegularExpression::CaseInsensitiveOption);
    auto tm = titleRx.match(content);
    if (tm.hasMatch()) {
        mf.title_ = tm.captured(1)
                     .replace(QStringLiteral("''"), QStringLiteral("'"));
    }

    // Pair regex — <digits> : '<label>' repeated. RESEARCH Pitfall 5
    // says td/m/modifm2d has a typo `80;` — we accept ONLY `:` as the
    // separator so the typo row falls back to "80;" via MenuFile::label().
    static const QRegularExpression pairRx(
        QStringLiteral(R"((\d+)\s*:\s*'((?:[^']|'')*)')"));
    auto mi = pairRx.globalMatch(content);
    while (mi.hasNext()) {
        auto m = mi.next();
        int code = m.captured(1).toInt();
        QString label = m.captured(2)
                         .replace(QStringLiteral("''"),
                                  QStringLiteral("'"));
        mf.labels_.insert(code, label);
    }

    return *c.insert(name, mf);
}

void XvueMenuFileParser::clearCacheForTesting() {
    cache().clear();
}
