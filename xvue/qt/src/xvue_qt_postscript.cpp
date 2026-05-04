// xvue/qt/src/xvue_qt_postscript.cpp
// Phase 7 Plan 02 (EXPORT-04): handleLasops() — TDD RED stage stub.
// Real verbatim port of xvuelc.c:1187-1304 lands in the GREEN commit.
//
// Plan 03 fills per-primitive helpers (line/face/etc.) — those bodies stay
// stubbed in Plan 02 by design.
#include "xvue_qt_postscript.h"
#include "xvue_qt_app.h"
#include "xvue_qt_api.h"   // XVUE_QT_ASSERT_MAIN_THREAD

const char* PsEmitter::kPostscriptFilename = "TEMPORAIRE.EPS";
const char* PsEmitter::kQualityFilename    = "TEMPORAIRE.QUA";

PsEmitter::PsEmitter() = default;

PsEmitter::~PsEmitter() {
    // Empty in RED. GREEN commit replaces this with the verbatim free pattern.
}

void PsEmitter::handleLasops(int /*lasops*/) {
    // RED stub: empty body — tests must FAIL.
}

// Per-primitive emit helpers — bodies land in Plan 03. Empty in Plan 02.
void PsEmitter::line(int, int, int, int)                                {}
void PsEmitter::traitcouleur(int,int,int,int,float,float,float)         {}
void PsEmitter::face(const int*, int)                                   {}
void PsEmitter::faceisocouleur(const int*, int, float, float, float)    {}
void PsEmitter::flpt(int, int, float)                                   {}
void PsEmitter::ellipse(int, int, int, int, int, int)                   {}
void PsEmitter::fond(float, float, float)                               {}
void PsEmitter::clear()                                                  {}
void PsEmitter::epaisseur(int)                                           {}
void PsEmitter::typetrait(int)                                           {}
void PsEmitter::chargefonte(const QString&, int, int, int, bool, bool)  {}
void PsEmitter::texte(const char*, int, int, int)                       {}
