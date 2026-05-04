// xvue/qt/src/xvue_qt_postscript.cpp
// Phase 7 Plan 02 (EXPORT-04, D-05/D-06/D-07): handleLasops() — verbatim
// port of xvuelc.c:1187-1304. Plan 03 fills per-primitive helpers
// (line/face/etc.) — those bodies stay stubbed in Plan 02 by design.
//
// EVERY branch, EVERY assignment, EVERY format string must match the
// original byte-for-byte. Anti-pattern: simplifying with QString /
// std::string. DO NOT.
//
// Three documented divergences from the byte-verbatim semantics:
//   1. T-07-09 fopen failure: xvuelc.c aborts via the C library exit
//      (status 1); we call QMessageBox::critical + qApp->quit() so
//      unsaved Fortran state has a chance to flush via the normal
//      teardown path.
//   2. Pitfall 3 chaine[-1] underflow: xvuelc.c writes through chaine[-1]
//      unconditionally (UB) when lasopsc==-3 (idx == -lasopsc-4 == -1).
//      We guard with an explicit idx>=0 && idx<MXRECT_ check. The bytes
//      inside chaine[idx] were never read back into the file output
//      stream, so PostScript byte output is unchanged. ASAN-clean.
//   3. Destructor + chaine free: in the lasops==0 branch, xvuelc.c reads
//      `chaine[i] = '\0'; free(chaine[i])` which (because chaine[i] is
//      char*) actually NULLs the pointer BEFORE freeing — yielding a
//      no-op free + memory leak. We port the *intended* semantics
//      (clobber-then-free), then add a defensive `chaine_[i] = nullptr`
//      to guard against double-free in our destructor (which also frees
//      chaine_[i]). The bytes inside chaine_[i] were never read after
//      the clobber in the legacy code, so PostScript output is unchanged.
#include "xvue_qt_postscript.h"
#include "xvue_qt_app.h"
#include "xvue_qt_api.h"   // XVUE_QT_ASSERT_MAIN_THREAD
#include "xvue_qt_i18n.h"  // xvueIsEnglish() — replaces xvuelc.c `langage`

#include <QApplication>
#include <QMessageBox>
#include <QString>
#include <QStringLiteral>
#include <algorithm>   // std::min, std::max
#include <cstdio>      // fopen, fclose, fprintf, remove
#include <cstdlib>     // calloc, free

const char* PsEmitter::kPostscriptFilename = "TEMPORAIRE.EPS";
const char* PsEmitter::kQualityFilename    = "TEMPORAIRE.QUA";

PsEmitter::PsEmitter() = default;

PsEmitter::~PsEmitter() {
    // Close any still-open TEMPORAIRE.EPS so on-disk PostScript is not
    // truncated. teardown_atexit() in xvue_qt_app.cpp resets ps_emitter_
    // BEFORE the QApplication is leaked, so this destructor runs in a
    // well-defined Qt state.
    if (fpo_ != nullptr) {
        std::fclose(fpo_);
        fpo_ = nullptr;
    }
    // Free any chaine_[i] buffers still allocated. xvuelc.c never reaches
    // this teardown path (program exit happens with chaine[] still live
    // and the OS reclaims) — but the Qt port may construct/destruct
    // PsEmitter multiple times across tests, so we must clean up.
    for (int i = 0; i < MXRECT_; ++i) {
        if (chaine_[i] != nullptr) {
            std::free(chaine_[i]);
            chaine_[i] = nullptr;
        }
    }
}

// ----------------------------------------------------------------------------
// handleLasops — VERBATIM port of xvuelc.c:1187-1304.
// EVERY branch, EVERY assignment, EVERY format string matches the original.
// Anti-pattern: simplifying with QString / std::string. DO NOT.
//
// Re-entrant: effacer in xvuelc.c:1414, 1435 calls back with lasops+100,
// hitting the mode-100 fall-through branch. Depth never exceeds 2 in
// legitimate use — we are the SAME instance so calls stack on lasopsc_.
// ----------------------------------------------------------------------------
void PsEmitter::handleLasops(int lasops) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    int i;

    /* ouverture du postscript */
    if (lasops > 0 && lasopsc_ == 0) {
        lasopsc_ = lasops;
        if (fpo_ != nullptr) {
            std::fclose(fpo_);
            fpo_ = nullptr;
        }
        if (modepsc_ != 0) {
            if ((fpo_ = std::fopen(kQualityFilename, "r")) != nullptr) {
                std::fclose(fpo_);
                fpo_ = nullptr;
                std::remove(kQualityFilename);
            }
        }
        if ((fpo_ = std::fopen(kPostscriptFilename, "w")) == nullptr) {
            // T-07-09: replace the legacy abort with QMessageBox + clean
            // shutdown so unsaved Fortran state has a chance to flush
            // through the normal teardown path. xvuelc.c:1234 used the
            // C-library exit; we use qApp->quit() instead.
            const QString title = xvueIsEnglish()
                ? QStringLiteral("Error")
                : QStringLiteral("Erreur");
            const QString body = xvueIsEnglish()
                ? QStringLiteral("NO creation of file TEMPORAIRE.EPS — quitting")
                : QStringLiteral("NON creation du fichier TEMPORAIRE.EPS — sortie");
            QMessageBox::critical(nullptr, title, body);
            if (auto* qa = qApp) qa->quit();
            return;   // Do NOT continue with a null fpo_.
        }
        if (modepsc_ != 0) {
            for (i = 0; i < 8; ++i) {
                chaine_[i] = (char*)std::calloc(longchaine_[i], sizeof(char));
            }
        }
        nbrcon_ = 0;
        concat_[0] = '\0';
        iTe_ = 0; iFa_ = 0; ity_ = 0; iep_ = 0; iPo_ = 0;
        ire_ = 0; iRe_ = 0; iel_ = 0; iEl_ = 0; iFP_ = 0;
    }
    else {
        if (lasops == 0) {
            lasopsc_ = lasops;
            if (fpo_ != nullptr) {
                std::fprintf(fpo_, "%s", concat_);
                std::fclose(fpo_);
                fpo_ = nullptr;
                if (modepsc_ != 0) {
                    // VERBATIM xvuelc.c:1252 semantics — clobber-then-free.
                    // xvuelc.c reads `chaine[i] = '\0' ; free(chaine[i])`
                    // which (because chaine[i] is `char*`) actually NULLs
                    // the pointer BEFORE freeing — yielding a no-op free
                    // + memory leak. Plan 02 ports the *intended*
                    // semantics (D-05): clobber the first byte through
                    // the still-valid pointer, free the buffer, THEN null.
                    // The trailing `chaine_[i] = nullptr` is the ONLY
                    // non-verbatim addition; it is defensive against
                    // double-free in the destructor (also frees chaine_[i])
                    // and does NOT alter PostScript byte output — the bytes
                    // inside chaine_[i] were never read after free.
                    for (i = 0; i < 8; ++i) {
                        if (chaine_[i] != nullptr) {
                            *chaine_[i] = '\0';
                            std::free(chaine_[i]);
                            chaine_[i] = nullptr;
                        }
                    }
                }
                menu_ = 0;
            }
        }
        else {
            if (lasops > 100) {
                /* remise a zero des controles de macros */
                iTe_ = 0; iFa_ = 0; ity_ = 0; iep_ = 0; iPo_ = 0;
                ire_ = 0; iRe_ = 0; iel_ = 0; iEl_ = 0; iFP_ = 0;
                if (fpo_ != nullptr) {
                    std::fclose(fpo_);
                    fpo_ = nullptr;
                }
                if (modepsc_ != 0) {
                    if ((fpo_ = std::fopen(kQualityFilename, "r")) != nullptr) {
                        std::fclose(fpo_);
                        fpo_ = nullptr;
                        std::remove(kQualityFilename);
                    }
                }
                if ((fpo_ = std::fopen(kPostscriptFilename, "w")) == nullptr) {
                    // T-07-09: same QMessageBox + qApp->quit() pattern as
                    // the lasops>0 open arm above (xvuelc.c:1276-1283).
                    const QString title = xvueIsEnglish()
                        ? QStringLiteral("Error")
                        : QStringLiteral("Erreur");
                    const QString body = xvueIsEnglish()
                        ? QStringLiteral("NO creation of file TEMPORAIRE.EPS — quitting")
                        : QStringLiteral("NON creation du fichier TEMPORAIRE.EPS — sortie");
                    QMessageBox::critical(nullptr, title, body);
                    if (auto* qa = qApp) qa->quit();
                    return;
                }
                if (modepsc_ != 0) {
                    for (i = 0; i < 8; ++i) {
                        if (chaine_[i] != nullptr) *chaine_[i] = '\0';
                    }
                }
                lasopsc_ = lasopsc_ - 100;
            }
            else {
                std::fprintf(fpo_, "%s", concat_);
                nbrcon_ = 0;
                concat_[0] = '\0';
                if (lasops < -1) {
                    /* effacement du menu correspondant */
                    lasopsc_ = std::max(lasops, -11);
                    if (modepsc_ != 0) {
                        // VERBATIM xvuelc.c: chaine[-lasopsc-4] indexing.
                        // -lasopsc_ is in [1..11]; -lasopsc_-4 is in [-3..7].
                        //
                        // DIVERGENCE FROM xvuelc.c: defensive idx>=0 && idx<MXRECT
                        // guard added. Original xvuelc.c writes through chaine[-1]
                        // (UB) but this is harmless because the byte was never
                        // read back into the file output stream. Guard added so
                        // ASAN does not flag the access. See SUMMARY.md.
                        // This is the ONLY divergence from byte-verbatim semantics
                        // in the dispatch logic (the destructor + fopen-failure
                        // divergences live elsewhere — see file header).
                        int idx = -lasopsc_ - 4;
                        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
                            *chaine_[idx] = '\0';
                        }
                    }
                }
                else {
                    lasopsc_ = std::min(lasops, 11);
                    if (lasopsc_ == 2) menu_ = 1;
                }
            }
        }
    }
}

// Per-primitive emit helpers — bodies land in Plan 03 GREEN. Empty stubs
// here for the TDD RED phase so the test target links while the new
// per-primitive byte-output tests fail.
void PsEmitter::line(int, int, int, int)                                {}
void PsEmitter::traitcouleur(int,int,int,int,float,float,float)         {}
void PsEmitter::face(const int*, int)                                   {}
void PsEmitter::face(const MefistoPoint*, int)                          {}
void PsEmitter::faceisocouleur(const int*, int, float, float, float)    {}
void PsEmitter::flpt(int, int, float)                                   {}
void PsEmitter::ellipse(int, int, int, int, int, int)                   {}
void PsEmitter::fond(float, float, float)                               {}
void PsEmitter::clear()                                                  {}
void PsEmitter::epaisseur(int)                                           {}
void PsEmitter::typetrait(int)                                           {}
void PsEmitter::chargefonte(const QString&, int, int, int, bool, bool)  {}
void PsEmitter::texte(const char*, int, int, int)                       {}
void PsEmitter::traits(const MefistoPoint*, int)                        {}
void PsEmitter::facetraits(const MefistoPoint*, int,
                            float, float, float, float)                 {}
void PsEmitter::bordrectangle(int, int, int, int)                       {}
void PsEmitter::rectangle(int, int, int, int)                           {}
void PsEmitter::bordarcellipse(int, int, int, int, float, float)        {}
void PsEmitter::arcellipse(int, int, int, int, float, float)            {}
