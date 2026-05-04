// xvue/qt/src/xvue_qt_postscript.h
// Phase 7 Plan 02 (EXPORT-04, D-05, D-06, D-07): PsEmitter — verbatim Qt-side
// port of xvuelc.c:1187-1304 (xvpostscript_) plus the per-primitive emit
// helpers ported in Plan 03. Owned by XvueApp; single instance for the
// QApplication lifetime. Bridge-as-XvueApp-singleton pattern (Phase 5 D-02).
//
// All format strings inside this class MUST be byte-for-byte the same as
// the corresponding xvuelc.c sites. Reviewer: any "could simplify with
// QString::number" suggestion gets pushed back firmly (Pitfall 1).
//
// Y-flip rule (xvue/README_COORDS.md): pyFlip() is the ONLY place ypixels_-y
// arithmetic appears. Callers always pass canvas-Y (Y-down).
#pragma once
#include <cstdio>
#include <QString>

class XvueApp;

class PsEmitter {
public:
    PsEmitter();
    ~PsEmitter();

    // Verbatim port of xvuelc.c:1187-1304 (xvpostscript_ body).
    // Re-entrant (effacer recurses with lasops+100 — xvuelc.c:1414, 1435);
    // depth never exceeds 2 in legitimate use.
    void handleLasops(int lasops);

    // Active-emit gate: per-primitive helpers (Plan 03) early-return when this
    // is false. Mirrors the `if (lasopsc > 0)` gate scattered through xvuelc.c.
    bool active() const { return lasopsc_ > 0; }

    // ---- Per-primitive helpers — DECLARATIONS only here. Bodies in Plan 03.
    // The signatures are frozen by Plan 02 so Plan 03 can implement them
    // independently without churning the header.
    void line(int x1, int y1, int x2, int y2);                  // xvuelc.c:1942-1985 (xvtrait)
    void traitcouleur(int x1, int y1, int x2, int y2,
                      float r, float g, float b);               // xvuelc.c:~2000-2050 (xvtraitcouleur)
    void face(const int* xy, int n);                            // xvuelc.c:2050-2150 (xvface)
    void faceisocouleur(const int* xy, int n,
                        float r, float g, float b);             // xvuelc.c:2050-2150 (xvfaceisocouleur)
    void flpt(int x, int y, float v);                           // xvuelc.c:2580-2680 (xvflpt)
    void ellipse(int x, int y, int rx, int ry,
                 int a1, int a2);                               // xvuelc.c:2680-2780 (xvellipse)
    void fond(float r, float g, float b);                       // xvuelc.c:~2200-2260 (xvfond)
    void clear();                                                // effacer trigger
    void epaisseur(int pepais);                                  // xvuelc.c:1880-1900 (xvepaisseur)
    void typetrait(int ity);                                     // xvuelc.c:~1820-1860 (xvtypetrait)
    void chargefonte(const QString& family, int size_pt,
                     int ascent, int descent,
                     bool bold, bool italic);                   // xvuelc.c:1500-1568 (xvchargefonte) — D-08 mapping
    void texte(const char* s, int slen, int x, int y);          // xvuelc.c:~3100-3180 (xvtexte)

    // Public for unit tests + Plan 03 emit helpers. Y-flip happens here ONLY
    // (xvue/README_COORDS.md). Caller always passes canvas-Y (Y-down).
    int pyFlip(int y) const { return ypixels_ - y; }

    // Public read-only state (used by per-primitive helpers + tests).
    int  lasopsc() const { return lasopsc_; }
    int  modepsc() const { return modepsc_; }
    int  ypixels() const { return ypixels_; }
    void setCanvasDims(int xp, int yp) { xpixels_ = xp; ypixels_ = yp; }

    // ---- Test-only accessors (Plan 02 unit tests). Read-only window onto
    // the verbatim-ported file-statics so the GTest+QTest slots can assert
    // state-machine transitions without making the private members public.
    FILE*       fpoForTesting()       const { return fpo_; }
    const char* concatForTesting()    const { return concat_; }
    char**      chaineForTesting()          { return chaine_; }
    void        setModepscForTesting(int v) { modepsc_ = v; }

    // Test-only: snapshot of TEMPORAIRE.EPS path for fail-injection harness.
    static const char* kPostscriptFilename;   // = "TEMPORAIRE.EPS"
    static const char* kQualityFilename;      // = "TEMPORAIRE.QUA"

private:
    // ---- File-static-port (xvuelc.c:170-189). Names + types match. ----
    int   lasopsc_ = 0;
    int   modepsc_ = 0;
    FILE* fpo_     = nullptr;
    float counb_   = 0.0f;
    float courgb_[3] = {0.0f, 0.0f, 0.0f};
    int   palcourc_ = 0;
    int   xpixels_ = 0;
    int   ypixels_ = 0;
    static constexpr int MXRECT_ = 8;          // matches xvuelc.c MXRECT
    char* chaine_[MXRECT_]  = {nullptr};
    int   longchaine_[MXRECT_] = {1024,1024,1024,1024,1024,1024,1024,1024};
    char  format_[255]      = {0};
    int   menu_   = 0;
    int   nbrcon_ = 0;
    int   xinic_ = 0, yinic_ = 0, xcouc_ = 0, ycouc_ = 0;
    int   iTe_=0, iFa_=0, ity_=0, iep_=0, iPo_=0, ire_=0, iRe_=0, iel_=0, iEl_=0, iFP_=0;
    char  buf_[512]      = {0};
    char  concat_[1024]  = {0};   // sized at 1024 per xvuelc.c sprintf budget
    char  fontcour_[512] = {0};
};
