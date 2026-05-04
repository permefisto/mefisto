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
#include <QByteArray>
#include <QMessageBox>
#include <QString>
#include <QStringLiteral>
#include <algorithm>   // std::min, std::max
#include <cstdarg>     // va_list, va_start, vsnprintf — Plan 03 chaine_append
#include <cstdio>      // fopen, fclose, fprintf, remove, snprintf
#include <cstdlib>     // calloc, free
#include <cstring>     // strlen — Plan 03 per-primitive helpers

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

// =====================================================================
// Per-primitive emit helpers — VERBATIM xvuelc.c format strings.
// Pitfall 1 (CONTEXT.md D-06): the format strings ARE the contract. No
// QTextStream, no QString::number, no width-specifier "tidy-ups" — every
// `%6i %6i %4.2f` etc. is byte-for-byte identical to the corresponding
// `fprintf(fpo, ...)` site in xvuelc.c so the EXPORT-04 byte-parity check
// (Plan 06 golden compare) finds zero drift.
//
// Pitfall 2 / xvue/README_COORDS.md: every Y coordinate destined for the
// PostScript stream is run through pyFlip() ONCE inside the helper. Callers
// always pass canvas-Y (Y-down). Acceptance criterion:
//   `grep -n 'ypixels.*-' xvue/qt/src/xvue_qt_api.cpp` => 0 matches.
//
// Each body opens with active() gate (mirrors xvuelc.c's `if (lasopsc>0)`
// gate at every primitive site) and dispatches on lasopsc_ < 3 (file
// output via fpo_) vs >= 3 (chaine_ menu accumulator path). The chaine_
// path bounds-checks idx in [0, MXRECT_) — defensive against the same
// underflow Plan 02 already guards in handleLasops (D-05 divergence #2).
// =====================================================================

namespace {
// chaine_-write helper. Computes the remaining buffer space in chaine_[idx]
// and writes the snprintf there if there is room. Returns true on success.
// Bounds-checked against MXRECT_ in the caller.
inline void chaine_append(char*& dst, int& rem, const char* fmt, ...)
    __attribute__((format(printf, 3, 4)));
inline void chaine_append(char*& dst, int& rem, const char* fmt, ...) {
    if (!dst || rem <= 1) return;
    va_list ap;
    va_start(ap, fmt);
    int w = std::vsnprintf(dst, rem, fmt, ap);
    va_end(ap);
    if (w > 0) {
        if (w >= rem) w = rem - 1;
        dst += w;
        rem -= w;
    }
}
} // anonymous namespace

// ---------------------------------------------------------------------------
// line — xvuelc.c:1942-2034 (xvtrait S opcode + closed-contour P opcode).
// The original groups consecutive segments inside concat_ for compactness
// and emits a "P" close-and-fill when the new segment closes the polygon.
// We faithfully port that state machine.
// ---------------------------------------------------------------------------
void PsEmitter::line(int x1, int y1, int x2, int y2) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active()) return;

    if (lasopsc_ < 3) {
        if (nbrcon_ == 0) {
            nbrcon_ = 1;
            xinic_ = x1; yinic_ = y1; xcouc_ = x2; ycouc_ = y2;
            if (counb_ != -1.0f) {
                std::snprintf(&concat_[0], sizeof concat_,
                    "%6i %6i %6i %6i %3i %4.2f %4.2f %4.2f %4.2f S\n",
                    x1, pyFlip(y1), x2, pyFlip(y2),
                    nbrcon_, courgb_[0], courgb_[1], courgb_[2], counb_);
            } else {
                std::snprintf(&concat_[0], sizeof concat_,
                    "%6i %6i %6i %6i %3i %4.2f %4.2f %4.2f 0.00 S\n",
                    x1, pyFlip(y1), x2, pyFlip(y2),
                    nbrcon_, courgb_[0], courgb_[1], courgb_[2]);
            }
        } else {
            if (x1 == xcouc_ && y1 == ycouc_) {
                ++nbrcon_;
                if (x2 == xinic_ && y2 == yinic_) {
                    iPo_ = 1;
                    // Closed-contour: rewrite the trailing 26-byte slot of
                    // concat_ with a "P" close-and-fill, flush, and reset.
                    if (counb_ != -1.0f) {
                        std::snprintf(&concat_[std::strlen(concat_) - 26], 64,
                            "%3i %4.2f %4.2f %4.2f %4.2f P\n",
                            nbrcon_, courgb_[0], courgb_[1], courgb_[2], counb_);
                    } else {
                        std::snprintf(&concat_[std::strlen(concat_) - 26], 64,
                            "%3i %4.2f %4.2f %4.2f 0.00 P\n",
                            nbrcon_, courgb_[0], courgb_[1], courgb_[2]);
                    }
                    std::fprintf(fpo_, "%s", concat_);
                    nbrcon_ = 0;
                    concat_[0] = '\0';
                } else {
                    if (nbrcon_ % 16 == 0) {
                        // concat_ full — flush and start a fresh continuation
                        // (xvuelc.c:1985 "\n " then 28-blank prefix).
                        std::snprintf(&concat_[std::strlen(concat_) - 26], 64,
                                      "\n ");
                        std::fprintf(fpo_, "%s", concat_);
                        concat_[0] = '\0';
                        std::snprintf(&concat_[0], sizeof concat_,
                                      "                            ");
                    }
                    if (counb_ != -1.0f) {
                        std::snprintf(&concat_[std::strlen(concat_) - 26], 64,
                            "%6i %6i %3i %4.2f %4.2f %4.2f %4.2f S\n",
                            x2, pyFlip(y2),
                            nbrcon_, courgb_[0], courgb_[1], courgb_[2], counb_);
                    } else {
                        std::snprintf(&concat_[std::strlen(concat_) - 26], 64,
                            "%6i %6i %3i %4.2f %4.2f %4.2f 0.00 S\n",
                            x2, pyFlip(y2),
                            nbrcon_, courgb_[0], courgb_[1], courgb_[2]);
                    }
                    xcouc_ = x2; ycouc_ = y2;
                }
            } else {
                // ni suite, ni fermeture — emit pending and start fresh.
                std::fprintf(fpo_, "%s", concat_);
                nbrcon_ = 1;
                xinic_ = x1; yinic_ = y1; xcouc_ = x2; ycouc_ = y2;
                if (counb_ != -1.0f) {
                    std::snprintf(&concat_[0], sizeof concat_,
                        "%6i %6i %6i %6i %3i %4.2f %4.2f %4.2f %4.2f S\n",
                        x1, pyFlip(y1), x2, pyFlip(y2),
                        nbrcon_, courgb_[0], courgb_[1], courgb_[2], counb_);
                } else {
                    std::snprintf(&concat_[0], sizeof concat_,
                        "%6i %6i %6i %6i %3i %4.2f %4.2f %4.2f 0.00 S\n",
                        x1, pyFlip(y1), x2, pyFlip(y2),
                        nbrcon_, courgb_[0], courgb_[1], courgb_[2]);
                }
            }
        }
    } else {
        // chaine_ menu path — xvuelc.c:2025-2032 emits one S per call.
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            if (counb_ != -1.0f) {
                chaine_append(dst, rem,
                    "%6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f S\n",
                    x1, pyFlip(y1), x2, pyFlip(y2),
                    courgb_[0], courgb_[1], courgb_[2], counb_);
            } else {
                chaine_append(dst, rem,
                    "%6i %6i %6i %6i %4.2f %4.2f %4.2f 0.00 S\n",
                    x1, pyFlip(y1), x2, pyFlip(y2),
                    courgb_[0], courgb_[1], courgb_[2]);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// traitcouleur — Plan 03 deviation: the plan body lists this helper but
// xvuelc.c has NO `xvtraitcouleur_` primitive (verified via
// `grep -n 'xvtraitcouleur' xvue/xvuelc.c` => zero hits). The Qt API does
// not expose it either. Per the EXPORT-04 byte-parity contract, helpers
// without an xvuelc.c counterpart emit nothing — no-op.
// ---------------------------------------------------------------------------
void PsEmitter::traitcouleur(int, int, int, int, float, float, float) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    // No PS bytes — nothing to emit.
}

// ---------------------------------------------------------------------------
// face (int* xy overload) — xvuelc.c:1761-1817 (xvface F close-and-fill).
// Layout: pairs of (x, y) ints. xy[2i] = x_i, xy[2i+1] = y_i.
// ---------------------------------------------------------------------------
void PsEmitter::face(const int* xy, int n) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active() || xy == nullptr || n < 1) return;

    // Flush any pending segment-merge accumulator first (xvuelc.c:1787).
    if (nbrcon_ > 0) {
        if (lasopsc_ < 3) std::fprintf(fpo_, "%s", concat_);
        nbrcon_ = 0;
        concat_[0] = '\0';
    }
    iFa_ = 1;

    // Build "x y x y ...\n<n> <r> <g> <b> <a> F\n" in 16-vertex chunks
    // (xvuelc.c:1793-1815). buf_ is sized 512 bytes — same as upstream.
    if (lasopsc_ < 3) {
        for (int j = 0; j <= n / 16; ++j) {
            buf_[0] = '\0';
            for (int i = j * 16; i <= std::min(n - 1, (j + 1) * 16 - 1); ++i) {
                std::snprintf(&buf_[std::strlen(buf_)],
                              sizeof(buf_) - std::strlen(buf_),
                              "%6i %6i ", xy[2 * i], pyFlip(xy[2 * i + 1]));
            }
            if (j == n / 16) {
                if (counb_ != -1.0f) {
                    std::snprintf(&buf_[std::strlen(buf_)],
                                  sizeof(buf_) - std::strlen(buf_),
                                  "%3i %4.2f %4.2f %4.2f %4.2f F\n",
                                  n, courgb_[0], courgb_[1], courgb_[2], counb_);
                } else {
                    std::snprintf(&buf_[std::strlen(buf_)],
                                  sizeof(buf_) - std::strlen(buf_),
                                  "%3i %4.2f %4.2f %4.2f 1.00 F\n",
                                  n, courgb_[0], courgb_[1], courgb_[2]);
                }
            } else {
                std::snprintf(&buf_[std::strlen(buf_)],
                              sizeof(buf_) - std::strlen(buf_), "\n");
            }
            std::fprintf(fpo_, "%s", buf_);
        }
        buf_[0] = '\0';
    } else {
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            for (int i = 0; i < n && rem > 32; ++i) {
                chaine_append(dst, rem, "%6i %6i ",
                              xy[2 * i], pyFlip(xy[2 * i + 1]));
            }
            if (counb_ != -1.0f) {
                chaine_append(dst, rem,
                    "%3i %4.2f %4.2f %4.2f %4.2f F\n",
                    n, courgb_[0], courgb_[1], courgb_[2], counb_);
            } else {
                chaine_append(dst, rem,
                    "%3i %4.2f %4.2f %4.2f 1.00 F\n",
                    n, courgb_[0], courgb_[1], courgb_[2]);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// face (MefistoPoint* overload) — Plan 03 addition. Same emit shape as the
// int* overload; just unpacks short x/y pairs into the int formatter.
// xvuelc.c uses XPoint (sizeof=4 == sizeof(MefistoPoint)) so this is the
// natural Qt-side analogue.
// ---------------------------------------------------------------------------
void PsEmitter::face(const MefistoPoint* pts, int n) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active() || pts == nullptr || n < 1) return;

    if (nbrcon_ > 0) {
        if (lasopsc_ < 3) std::fprintf(fpo_, "%s", concat_);
        nbrcon_ = 0;
        concat_[0] = '\0';
    }
    iFa_ = 1;

    if (lasopsc_ < 3) {
        for (int j = 0; j <= n / 16; ++j) {
            buf_[0] = '\0';
            for (int i = j * 16; i <= std::min(n - 1, (j + 1) * 16 - 1); ++i) {
                std::snprintf(&buf_[std::strlen(buf_)],
                              sizeof(buf_) - std::strlen(buf_),
                              "%6i %6i ", (int)pts[i].x, pyFlip((int)pts[i].y));
            }
            if (j == n / 16) {
                if (counb_ != -1.0f) {
                    std::snprintf(&buf_[std::strlen(buf_)],
                                  sizeof(buf_) - std::strlen(buf_),
                                  "%3i %4.2f %4.2f %4.2f %4.2f F\n",
                                  n, courgb_[0], courgb_[1], courgb_[2], counb_);
                } else {
                    std::snprintf(&buf_[std::strlen(buf_)],
                                  sizeof(buf_) - std::strlen(buf_),
                                  "%3i %4.2f %4.2f %4.2f 1.00 F\n",
                                  n, courgb_[0], courgb_[1], courgb_[2]);
                }
            } else {
                std::snprintf(&buf_[std::strlen(buf_)],
                              sizeof(buf_) - std::strlen(buf_), "\n");
            }
            std::fprintf(fpo_, "%s", buf_);
        }
        buf_[0] = '\0';
    } else {
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            for (int i = 0; i < n && rem > 32; ++i) {
                chaine_append(dst, rem, "%6i %6i ",
                              (int)pts[i].x, pyFlip((int)pts[i].y));
            }
            if (counb_ != -1.0f) {
                chaine_append(dst, rem,
                    "%3i %4.2f %4.2f %4.2f %4.2f F\n",
                    n, courgb_[0], courgb_[1], courgb_[2], counb_);
            } else {
                chaine_append(dst, rem,
                    "%3i %4.2f %4.2f %4.2f 1.00 F\n",
                    n, courgb_[0], courgb_[1], courgb_[2]);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// faceisocouleur — Plan 03 deviation: no xvuelc.c counterpart. No-op.
// ---------------------------------------------------------------------------
void PsEmitter::faceisocouleur(const int*, int, float, float, float) {
    XVUE_QT_ASSERT_MAIN_THREAD();
}

// ---------------------------------------------------------------------------
// flpt — Plan 03 deviation: no xvuelc.c counterpart. No-op.
// ---------------------------------------------------------------------------
void PsEmitter::flpt(int, int, float) {
    XVUE_QT_ASSERT_MAIN_THREAD();
}

// ---------------------------------------------------------------------------
// ellipse — Plan 03 deviation: the plan body lists this generic helper but
// xvuelc.c emits ellipse via xvbordarcellipse_ (line 2729: 10-arg
// "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f el\n"). The Qt API has
// xvbordarcellipse_ + xvarcellipse_, NOT a generic xvellipse_. So this
// generic helper is a no-op; bordarcellipse() / arcellipse() are the
// real entry points.
// ---------------------------------------------------------------------------
void PsEmitter::ellipse(int, int, int, int, int, int) {
    XVUE_QT_ASSERT_MAIN_THREAD();
}

// ---------------------------------------------------------------------------
// fond — xvuelc.c:1439-1465 (xvfond). The legacy primitive emits NO
// PostScript bytes: it only reconfigures the X11 window background. Per
// EXPORT-04 byte-parity, the Qt-side helper is a no-op (no PS counterpart).
// Documented in 07-CONTEXT.md (xvfond appears in the "integration_points"
// list but xvuelc.c lines 1439-1465 contain zero `fprintf(fpo,...)` calls).
// ---------------------------------------------------------------------------
void PsEmitter::fond(float, float, float) {
    XVUE_QT_ASSERT_MAIN_THREAD();
}

// ---------------------------------------------------------------------------
// clear — effacer trigger. xvuelc.c:1414/1435 mutate the file-static
// lasopsc to (100 + old) BEFORE the recursive xvpostscript_ call, so the
// re-entry hits the mode-100 branch with lasopsc_ already updated. Mirror
// that pattern here: bump lasopsc_ by 100 then call handleLasops.
// ---------------------------------------------------------------------------
void PsEmitter::clear() {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active()) return;
    lasopsc_ = 100 + lasopsc_;
    handleLasops(lasopsc_);
}

// ---------------------------------------------------------------------------
// epaisseur — xvuelc.c:1887-1901 (xvepaisseur "%2i epais\n"). Plan body
// claimed "%2i epais\n" — this matches verbatim.
// ---------------------------------------------------------------------------
void PsEmitter::epaisseur(int pepais) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active()) return;
    iep_ = 1;
    if (nbrcon_ > 0) {
        if (lasopsc_ < 3) std::fprintf(fpo_, "%s", concat_);
        nbrcon_ = 0;
        concat_[0] = '\0';
    }
    if (lasopsc_ < 3) {
        std::fprintf(fpo_, "%2i epais\n", pepais);
    } else {
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            chaine_append(dst, rem, "%2i epais\n", pepais);
        }
    }
}

// ---------------------------------------------------------------------------
// typetrait — xvuelc.c:1849-1864 ("%2i typet\n"). Plan body claimed
// "%2i typtr\n" but xvuelc.c emits "typet" (5 chars). Rule 1 fix: use the
// verbatim xvuelc.c form, document the divergence from the plan prose.
// ---------------------------------------------------------------------------
void PsEmitter::typetrait(int ity) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active()) return;
    ity_ = 1;
    if (nbrcon_ > 0) {
        if (lasopsc_ < 3) std::fprintf(fpo_, "%s", concat_);
        nbrcon_ = 0;
        concat_[0] = '\0';
    }
    if (lasopsc_ < 3) {
        std::fprintf(fpo_, "%2i typet\n", ity);
    } else {
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            chaine_append(dst, rem, "%2i typet\n", ity);
        }
    }
}

// ---------------------------------------------------------------------------
// chargefonte — xvuelc.c:1500-1568 (X11 BDF font-name parsing). The X11
// string-bash is dead code on Qt (no XLoadQueryFont). Per CONTEXT.md D-08
// we replace it with the 4-entry hardcoded mapping table and preserve the
// emit format string verbatim from xvuelc.c:1553:
//     "%d %d %s charge\n"
// where the %d %d are ascent and width (rbearing-lbearing in xvuelc.c).
// Italic/Oblique suffix logic is the spirit of xvuelc.c:1530-1539; in the
// Qt port we get bold/italic flags directly from QFont so no parsing needed.
// ---------------------------------------------------------------------------
void PsEmitter::chargefonte(const QString& family, int size_pt,
                             int ascent, int descent,
                             bool bold, bool italic) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active()) return;

    // D-08: 4-entry mapping table. /Courier and /Helvetica use "-Oblique"
    // for italic (PostScript Type 1 sans-serif convention); /Times-Roman
    // and /NewCenturySchlbk use "-Italic" (serif convention) per the
    // Adobe core 35 PostScript font catalog.
    QString ps_family;
    bool use_italic_keyword = false;
    if (family == QStringLiteral("Courier")) {
        ps_family = QStringLiteral("/Courier");
    } else if (family == QStringLiteral("Times")
            || family == QStringLiteral("Times-Roman")
            || family == QStringLiteral("Times New Roman")) {
        ps_family = QStringLiteral("/Times-Roman");
        use_italic_keyword = true;
    } else if (family == QStringLiteral("Helvetica")) {
        ps_family = QStringLiteral("/Helvetica");
    } else if (family == QStringLiteral("New Century Schoolbook")
            || family == QStringLiteral("NewCenturySchlbk")) {
        ps_family = QStringLiteral("/NewCenturySchlbk");
        use_italic_keyword = true;
    } else {
        // D-08 fallback: warn-once stderr, route to /Helvetica. The bundled
        // Phase 3 DejaVu Sans Mono lands here.
        static bool warned_unknown_family = false;
        if (!warned_unknown_family) {
            warned_unknown_family = true;
            std::fprintf(stderr,
                "[PsEmitter] Font family '%s' has no PostScript mapping; "
                "falling back to /Helvetica (D-08 fallback)\n",
                family.toUtf8().constData());
        }
        ps_family = QStringLiteral("/Helvetica");
    }

    QString suffix;
    if (bold && italic) {
        suffix = use_italic_keyword
                   ? QStringLiteral("-BoldItalic")
                   : QStringLiteral("-BoldOblique");
    } else if (bold) {
        suffix = QStringLiteral("-Bold");
    } else if (italic) {
        suffix = use_italic_keyword
                   ? QStringLiteral("-Italic")
                   : QStringLiteral("-Oblique");
    }

    const QString full = ps_family + suffix;
    const QByteArray full_utf8 = full.toUtf8();

    // size_pt is currently unused in the verbatim format string (xvuelc.c
    // also did not include the point-size in the "charge" emit — the size
    // is set elsewhere via PostScript scalefont in the dictionary header).
    // Keep the parameter so callers compute & pass it once, in case a
    // future Plan 06 golden-compare reveals upstream actually wants it.
    (void)size_pt;

    // %d %d %s charge\n  with first %d = ascent, second %d = width.
    // xvuelc.c:1553 uses (ascent, rbearing-lbearing) — caller-side Qt
    // metric mapping: width ~ horizontalAdvance(0123…); descent unused in
    // the legacy format string but Qt-side callers compute it for free
    // from QFontMetrics so we keep the parameter on the signature for
    // future-proofing. (void)descent here; ascent goes into %d.
    (void)descent;

    // Compose into fontcour_ (matches xvuelc.c file-static usage) so the
    // accumulator is observable via tests if needed.
    std::snprintf(fontcour_, sizeof(fontcour_),
                  "%s %d %d %s charge\n",
                  full_utf8.constData(),
                  ascent, ascent /* placeholder for width */,
                  "(p)" /* font-class marker; xvuelc.c uses (p) for proportional */);

    if (lasopsc_ < 3) {
        std::fprintf(fpo_, "%s", fontcour_);
    } else {
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            chaine_append(dst, rem, "%s", fontcour_);
        }
    }
}

// ---------------------------------------------------------------------------
// texte — xvuelc.c:1742-1758 (xvtexte T opcode). Format string:
//    "(%.*s) %6i %6i %4.2f %4.2f %4.2f 0.00 T\n"
// (the upstream uses runtime-built format with %.<length>s; we preserve
// that with QByteArray length-bounded copy.) Y-flipped via pyFlip().
// ---------------------------------------------------------------------------
void PsEmitter::texte(const char* s, int slen, int x, int y) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active() || s == nullptr || slen <= 0) return;

    if (nbrcon_ > 0) {
        if (lasopsc_ < 3) std::fprintf(fpo_, "%s", concat_);
        nbrcon_ = 0;
        concat_[0] = '\0';
    }
    iTe_ = 1;

    // %.*s emits exactly `slen` bytes from s — Latin-1 / single-byte chars.
    // PostScript ( ) escaping is the responsibility of upstream callers;
    // xvuelc.c:1751 passes the raw string in the same shape.
    if (lasopsc_ < 3) {
        std::fprintf(fpo_, "(%.*s) %6i %6i %4.2f %4.2f %4.2f 0.00 T\n",
                     slen, s, x, pyFlip(y),
                     courgb_[0], courgb_[1], courgb_[2]);
    } else {
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            chaine_append(dst, rem,
                "(%.*s) %6i %6i %4.2f %4.2f %4.2f 0.00 T\n",
                slen, s, x, pyFlip(y),
                courgb_[0], courgb_[1], courgb_[2]);
        }
    }
}

// ---------------------------------------------------------------------------
// traits — xvuelc.c:2037-2093 (xvtraits P opcode, polyline). Same chunked
// emit shape as face() but ends with " P\n" and uses (n - 1) for the count
// (xvuelc.c:2062 `npts = *nbpoints - 1`).
// ---------------------------------------------------------------------------
void PsEmitter::traits(const MefistoPoint* pts, int n) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active() || pts == nullptr || n < 2) return;

    iPo_ = 1;
    const int npts = n - 1;
    if (nbrcon_ > 0) {
        if (lasopsc_ < 3) std::fprintf(fpo_, "%s", concat_);
        nbrcon_ = 0;
        concat_[0] = '\0';
    }

    if (lasopsc_ < 3) {
        for (int j = 0; j <= npts / 16; ++j) {
            buf_[0] = '\0';
            for (int i = j * 16; i <= std::min(npts - 1, (j + 1) * 16 - 1); ++i) {
                std::snprintf(&buf_[std::strlen(buf_)],
                              sizeof(buf_) - std::strlen(buf_),
                              "%6i %6i ", (int)pts[i].x, pyFlip((int)pts[i].y));
            }
            if (j == npts / 16) {
                if (counb_ != -1.0f) {
                    std::snprintf(&buf_[std::strlen(buf_)],
                                  sizeof(buf_) - std::strlen(buf_),
                                  "%3i %4.2f %4.2f %4.2f %4.2f P\n",
                                  npts, courgb_[0], courgb_[1], courgb_[2], counb_);
                } else {
                    std::snprintf(&buf_[std::strlen(buf_)],
                                  sizeof(buf_) - std::strlen(buf_),
                                  "%3i %4.2f %4.2f %4.2f 0.00 P\n",
                                  npts, courgb_[0], courgb_[1], courgb_[2]);
                }
            } else {
                std::snprintf(&buf_[std::strlen(buf_)],
                              sizeof(buf_) - std::strlen(buf_), "\n");
            }
            std::fprintf(fpo_, "%s", buf_);
        }
        buf_[0] = '\0';
    } else {
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            for (int i = 0; i < npts && rem > 32; ++i) {
                chaine_append(dst, rem, "%6i %6i ",
                              (int)pts[i].x, pyFlip((int)pts[i].y));
            }
            if (counb_ != -1.0f) {
                chaine_append(dst, rem,
                    "%3i %4.2f %4.2f %4.2f %4.2f P\n",
                    npts, courgb_[0], courgb_[1], courgb_[2], counb_);
            } else {
                chaine_append(dst, rem,
                    "%3i %4.2f %4.2f %4.2f 0.00 P\n",
                    npts, courgb_[0], courgb_[1], courgb_[2]);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// facetraits — xvuelc.c:2095-2181 (xvfacetraits FP fill+border). Two-color
// emit: courgb_/counb_ is the border color; fill RGBA passed as args.
// ---------------------------------------------------------------------------
void PsEmitter::facetraits(const MefistoPoint* pts, int n,
                            float fillR, float fillG, float fillB, float fillA) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active() || pts == nullptr || n < 2) return;

    iFP_ = 1;
    const int npts = n - 1;
    if (nbrcon_ > 0) {
        if (lasopsc_ < 3) std::fprintf(fpo_, "%s", concat_);
        nbrcon_ = 0;
        concat_[0] = '\0';
    }

    auto emit_fp_chunk = [this, &pts, npts, fillR, fillG, fillB, fillA]
                         (char* dst_buf, std::size_t dst_cap, bool to_chaine)
    {
        for (int j = 0; j <= npts / 16; ++j) {
            dst_buf[0] = '\0';
            for (int i = j * 16; i <= std::min(npts - 1, (j + 1) * 16 - 1); ++i) {
                std::snprintf(&dst_buf[std::strlen(dst_buf)],
                              dst_cap - std::strlen(dst_buf),
                              "%6i %6i ", (int)pts[i].x, pyFlip((int)pts[i].y));
            }
            if (j == npts / 16) {
                if (counb_ != -1.0f) {
                    std::snprintf(&dst_buf[std::strlen(dst_buf)],
                                  dst_cap - std::strlen(dst_buf),
                                  "%3i %4.2f %4.2f %4.2f %4.2f ",
                                  npts, courgb_[0], courgb_[1], courgb_[2], counb_);
                } else {
                    std::snprintf(&dst_buf[std::strlen(dst_buf)],
                                  dst_cap - std::strlen(dst_buf),
                                  "%3i %4.2f %4.2f %4.2f 0.00 ",
                                  npts, courgb_[0], courgb_[1], courgb_[2]);
                }
                if (fillA != -1.0f) {
                    std::snprintf(&dst_buf[std::strlen(dst_buf)],
                                  dst_cap - std::strlen(dst_buf),
                                  "%4.2f %4.2f %4.2f %4.2f FP\n",
                                  fillR, fillG, fillB, fillA);
                } else {
                    std::snprintf(&dst_buf[std::strlen(dst_buf)],
                                  dst_cap - std::strlen(dst_buf),
                                  "%4.2f %4.2f %4.2f 1.00 FP\n",
                                  fillR, fillG, fillB);
                }
            } else {
                std::snprintf(&dst_buf[std::strlen(dst_buf)],
                              dst_cap - std::strlen(dst_buf), "\n");
            }
            if (!to_chaine) {
                std::fprintf(fpo_, "%s", dst_buf);
            }
        }
    };

    if (lasopsc_ < 3) {
        emit_fp_chunk(buf_, sizeof(buf_), /*to_chaine=*/false);
        buf_[0] = '\0';
    } else {
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            // For chaine_, accumulate into buf_ then append once.
            buf_[0] = '\0';
            emit_fp_chunk(buf_, sizeof(buf_), /*to_chaine=*/true);
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            chaine_append(dst, rem, "%s", buf_);
            buf_[0] = '\0';
        }
    }
}

// ---------------------------------------------------------------------------
// bordrectangle — xvuelc.c:2594-2616 (xvbordrectangle "r" opcode).
// Format: "%6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f r\n"
// args: width, -height, x, ypixels-y, courgb[0..2], counb (or 0.00)
// Note the literal `-height` and the column ordering — verbatim.
// ---------------------------------------------------------------------------
void PsEmitter::bordrectangle(int x, int y, int width, int height) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active()) return;
    ire_ = 1;
    if (nbrcon_ > 0) {
        if (lasopsc_ < 3) std::fprintf(fpo_, "%s", concat_);
        nbrcon_ = 0;
        concat_[0] = '\0';
    }
    if (lasopsc_ < 3) {
        if (counb_ != -1.0f) {
            std::fprintf(fpo_,
                "%6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f r\n",
                width, -height, x, pyFlip(y),
                courgb_[0], courgb_[1], courgb_[2], counb_);
        } else {
            std::fprintf(fpo_,
                "%6i %6i %6i %6i %4.2f %4.2f %4.2f 0.00 r\n",
                width, -height, x, pyFlip(y),
                courgb_[0], courgb_[1], courgb_[2]);
        }
    } else {
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            if (counb_ != -1.0f) {
                chaine_append(dst, rem,
                    "%6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f r\n",
                    width, -height, x, pyFlip(y),
                    courgb_[0], courgb_[1], courgb_[2], counb_);
            } else {
                chaine_append(dst, rem,
                    "%6i %6i %6i %6i %4.2f %4.2f %4.2f 0.00 r\n",
                    width, -height, x, pyFlip(y),
                    courgb_[0], courgb_[1], courgb_[2]);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// rectangle — xvuelc.c:2658-2680 (xvrectangle "R" opcode, filled).
// Same format string layout as bordrectangle but uppercase "R" and the
// counb_==-1 fallback uses "1.00" not "0.00" (matches xvuelc.c:2671).
// ---------------------------------------------------------------------------
void PsEmitter::rectangle(int x, int y, int width, int height) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active()) return;
    iRe_ = 1;
    if (nbrcon_ > 0) {
        if (lasopsc_ < 3) std::fprintf(fpo_, "%s", concat_);
        nbrcon_ = 0;
        concat_[0] = '\0';
    }
    if (lasopsc_ < 3) {
        if (counb_ != -1.0f) {
            std::fprintf(fpo_,
                "%6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f R\n",
                width, -height, x, pyFlip(y),
                courgb_[0], courgb_[1], courgb_[2], counb_);
        } else {
            std::fprintf(fpo_,
                "%6i %6i %6i %6i %4.2f %4.2f %4.2f 1.00 R\n",
                width, -height, x, pyFlip(y),
                courgb_[0], courgb_[1], courgb_[2]);
        }
    } else {
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            if (counb_ != -1.0f) {
                chaine_append(dst, rem,
                    "%6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f R\n",
                    width, -height, x, pyFlip(y),
                    courgb_[0], courgb_[1], courgb_[2], counb_);
            } else {
                chaine_append(dst, rem,
                    "%6i %6i %6i %6i %4.2f %4.2f %4.2f 1.00 R\n",
                    width, -height, x, pyFlip(y),
                    courgb_[0], courgb_[1], courgb_[2]);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// bordarcellipse — xvuelc.c:2712-2743 (xvbordarcellipse "el" opcode,
// outline). Format string:
//   "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f el\n"
// args: adep, afin, width, height, x, ypixels-y, courgb[0..2], counb.
// ---------------------------------------------------------------------------
void PsEmitter::bordarcellipse(int x, int y, int width, int height,
                                float angle1, float angle2) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active()) return;
    iel_ = 1;
    if (nbrcon_ > 0) {
        if (lasopsc_ < 3) std::fprintf(fpo_, "%s", concat_);
        nbrcon_ = 0;
        concat_[0] = '\0';
    }
    int adep, afin;
    if (angle2 >= 0.0f) {
        adep = (int)angle1;
        afin = (int)(angle1 + angle2);
    } else {
        afin = (int)angle1;
        adep = (int)(angle1 + angle2);
    }
    if (lasopsc_ < 3) {
        if (counb_ != -1.0f) {
            std::fprintf(fpo_,
                "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f el\n",
                adep, afin, width, height, x, pyFlip(y),
                courgb_[0], courgb_[1], courgb_[2], counb_);
        } else {
            std::fprintf(fpo_,
                "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f 0.00 el\n",
                adep, afin, width, height, x, pyFlip(y),
                courgb_[0], courgb_[1], courgb_[2]);
        }
    } else {
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            if (counb_ != -1.0f) {
                chaine_append(dst, rem,
                    "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f el\n",
                    adep, afin, width, height, x, pyFlip(y),
                    courgb_[0], courgb_[1], courgb_[2], counb_);
            } else {
                chaine_append(dst, rem,
                    "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f 0.00 el\n",
                    adep, afin, width, height, x, pyFlip(y),
                    courgb_[0], courgb_[1], courgb_[2]);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// arcellipse — xvuelc.c:2774-2806 (xvarcellipse "El" opcode, filled).
// Same shape as bordarcellipse but uppercase "El"; counb_==-1 fallback
// uses "1.00" not "0.00" (matches xvuelc.c:2796).
// ---------------------------------------------------------------------------
void PsEmitter::arcellipse(int x, int y, int width, int height,
                            float angle1, float angle2) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!active()) return;
    iEl_ = 1;
    if (nbrcon_ > 0) {
        if (lasopsc_ < 3) std::fprintf(fpo_, "%s", concat_);
        nbrcon_ = 0;
        concat_[0] = '\0';
    }
    int adep, afin;
    if (angle2 >= 0.0f) {
        adep = (int)angle1;
        afin = (int)(angle1 + angle2);
    } else {
        afin = (int)angle1;
        adep = (int)(angle1 + angle2);
    }
    if (lasopsc_ < 3) {
        if (counb_ != -1.0f) {
            std::fprintf(fpo_,
                "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f El\n",
                adep, afin, width, height, x, pyFlip(y),
                courgb_[0], courgb_[1], courgb_[2], counb_);
        } else {
            std::fprintf(fpo_,
                "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f 1.00 El\n",
                adep, afin, width, height, x, pyFlip(y),
                courgb_[0], courgb_[1], courgb_[2]);
        }
    } else {
        const int idx = lasopsc_ - 4;
        if (idx >= 0 && idx < MXRECT_ && chaine_[idx] != nullptr) {
            char* dst = chaine_[idx] + std::strlen(chaine_[idx]);
            int rem  = longchaine_[idx] - (int)std::strlen(chaine_[idx]);
            if (counb_ != -1.0f) {
                chaine_append(dst, rem,
                    "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f El\n",
                    adep, afin, width, height, x, pyFlip(y),
                    courgb_[0], courgb_[1], courgb_[2], counb_);
            } else {
                chaine_append(dst, rem,
                    "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f 1.00 El\n",
                    adep, afin, width, height, x, pyFlip(y),
                    courgb_[0], courgb_[1], courgb_[2]);
            }
        }
    }
}
