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

// Plan 03: pull in MefistoPoint (typedef struct {short x; short y;}) from
// xvue_qt_api.h for the face(const MefistoPoint*, int) overload that
// Qt-side callers in xvue_qt_api.cpp use to pass arrays without int*
// conversion. Both face() overloads emit byte-equal PostScript bytes.
#include "xvue_qt_api.h"

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

    // ---- Per-primitive helpers — Plan 02 froze these signatures, Plan 03
    // fills the bodies. Format strings inside the .cpp bodies are byte-for-
    // byte copies from the corresponding xvuelc.c sites.
    //
    // Plan 03 deviation note (Rule 3 blocking): the plan body of 07-03 lists
    // helpers (traitcouleur, faceisocouleur, flpt, ellipse) that have NO
    // counterpart in xvuelc.c — they are signatures inherited from the
    // Phase 7 plan-prose enumeration but no `if(lasopsc>0) fprintf(fpo,...)`
    // emit site exists upstream for them. Per the EXPORT-04 byte-parity
    // contract ("Each Qt drawing primitive that has an xvuelc.c counterpart
    // emits the SAME bytes"), helpers WITHOUT a counterpart are no-ops —
    // they early-return without touching fpo_/concat_. Documented inline.
    void line(int x1, int y1, int x2, int y2);                  // xvuelc.c:1942-2035 (xvtrait S/P)
    void traitcouleur(int x1, int y1, int x2, int y2,
                      float r, float g, float b);               // Plan 03: no xvuelc.c counterpart — no-op
    void face(const int* xy, int n);                            // xvuelc.c:1761-1817 (xvface F) — int* overload
    void face(const MefistoPoint* pts, int n);                  // Plan 03 overload — Qt-side wiring uses MefistoPoint
    void faceisocouleur(const int* xy, int n,
                        float r, float g, float b);             // Plan 03: no xvuelc.c counterpart — no-op
    void flpt(int x, int y, float v);                           // Plan 03: no xvuelc.c counterpart — no-op
    void ellipse(int x, int y, int rx, int ry,
                 int a1, int a2);                               // Plan 03: no xvuelc.c counterpart — no-op
    void fond(float r, float g, float b);                       // xvuelc.c:1439 — xvfond emits NO PS (no-op)
    void clear();                                                // effacer trigger — delegates to handleLasops(lasopsc_+100)
    void epaisseur(int pepais);                                  // xvuelc.c:1895 ("%2i epais\n")
    void typetrait(int ity);                                     // xvuelc.c:1856 ("%2i typet\n" verbatim)
    void chargefonte(const QString& family, int size_pt,
                     int ascent, int descent,
                     bool bold, bool italic);                   // xvuelc.c:1553 — D-08 mapping table replaces X11 BDF parse
    void texte(const char* s, int slen, int x, int y);          // xvuelc.c:1722-1758 (xvtexte T)

    // Plan 03 byte-parity additions (Rule 2 — required to wire the actual
    // Qt drawing primitives that DO have xvuelc.c counterparts but whose
    // helper signatures the plan body did not enumerate):
    void traits(const MefistoPoint* pts, int n);                // xvuelc.c:2037-2093 (xvtraits P)
    void facetraits(const MefistoPoint* pts, int n,
                    float fillR, float fillG, float fillB,
                    float fillA);                               // xvuelc.c:2095-2181 (xvfacetraits FP — fill+border)
    void bordrectangle(int x, int y, int width, int height);    // xvuelc.c:2573-2617 (xvbordrectangle r)
    void rectangle(int x, int y, int width, int height);        // xvuelc.c:2637-2681 (xvrectangle R)
    void bordarcellipse(int x, int y, int width, int height,
                        float angle1, float angle2);            // xvuelc.c:2684-2743 (xvbordarcellipse el)
    void arcellipse(int x, int y, int width, int height,
                    float angle1, float angle2);                // xvuelc.c:2746-2806 (xvarcellipse El)

    // Public for unit tests + Plan 03 emit helpers. Y-flip happens here ONLY
    // (xvue/README_COORDS.md). Caller always passes canvas-Y (Y-down).
    int pyFlip(int y) const { return ypixels_ - y; }

    // Public read-only state (used by per-primitive helpers + tests).
    int  lasopsc() const { return lasopsc_; }
    int  modepsc() const { return modepsc_; }
    int  ypixels() const { return ypixels_; }
    void setCanvasDims(int xp, int yp) { xpixels_ = xp; ypixels_ = yp; }

    // ---- Test-only accessors. Read-only window onto the verbatim-ported
    // file-statics so the GTest+QTest slots can assert state-machine
    // transitions without making the private members public.
    FILE*       fpoForTesting()       const { return fpo_; }
    const char* concatForTesting()    const { return concat_; }
    char**      chaineForTesting()          { return chaine_; }
    void        setModepscForTesting(int v) { modepsc_ = v; }
    // Plan 03: explicit test-only setters for counb_ / courgb_ so byte-output
    // tests can assert the counb!=-1 and counb==-1 emit branches independently.
    void        setCounbForTesting(float v) { counb_ = v; }
    void        setCourgbForTesting(float r, float g, float b) {
        courgb_[0] = r; courgb_[1] = g; courgb_[2] = b;
    }

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
