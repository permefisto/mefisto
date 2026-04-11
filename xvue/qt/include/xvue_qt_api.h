// xvue/qt/include/xvue_qt_api.h
// Phase 0: Fortran-facing extern "C" ABI for libxvueqt.a
// Copied byte-identically from xvue/xvuelc.c — DO NOT edit signatures without
// user approval. Every deviation is a link-time or crash-time bug.
//
// Scope resolution (Planner Alert, Option A per user):
//   - 57 Fortran-facing entries declared here
//   - 3 C-internal helpers NOT declared here (they move to xvue/qt/src/ as
//     static C++ functions in a later phase if needed):
//       (C-internal palette/colormap helpers at xvuelc.c:358, :463, :503)

#ifndef XVUE_QT_API_H
#define XVUE_QT_API_H

#include <cstddef>

// -------------------------------------------------------------------
// Trailing-underscore name mangling — copied verbatim from
// xvue/xvuelc.c lines 60-70. Do not change.
// -------------------------------------------------------------------
#ifdef __GNUC__
#  define proc(x) x##_
#else
#  define proc(x) x/**/_
#endif
#undef  proc
#define proc(x) x##_

// -------------------------------------------------------------------
// MefistoPoint — byte-identical to Xlib XPoint (D-07, Pitfall 4).
// Used by xvface_, xvtraits_, xvfacetraits_.
// -------------------------------------------------------------------
typedef struct { short x; short y; } MefistoPoint;
static_assert(sizeof(MefistoPoint) == 4,
              "MefistoPoint must match Xlib XPoint layout (4 bytes)");

// -------------------------------------------------------------------
// Debug thread-affinity macro — Phase 1 body (SHELL-07, D-12/D-13).
// Null-guarded: qApp may be null on the very first call before
// XvueApp::ensure() has run (e.g. xvpxecran_ callable standalone).
// Every real Phase 1 entry point therefore calls XvueApp::ensure()
// as its first statement, then this macro; the null-guard keeps the
// assertion safe for any caller that violates that order.
// -------------------------------------------------------------------
#ifdef QT_DEBUG
#  include <QThread>
#  include <QApplication>
#  define XVUE_QT_ASSERT_MAIN_THREAD()                                    \
      do {                                                                \
          if (qApp) {                                                     \
              Q_ASSERT(QThread::currentThread() == qApp->thread());       \
          }                                                               \
      } while (0)
#else
#  define XVUE_QT_ASSERT_MAIN_THREAD() ((void)0)
#endif

// -------------------------------------------------------------------
// Fortran-facing ABI — 57 entry points. Order matches xvue/xvuelc.c
// line order for easy side-by-side diff.
// -------------------------------------------------------------------
extern "C" {

// ---- 1. Init / environment ----
void  proc(languemefisto)(int *langue);                                     // xvuelc.c:227
void *proc(dctnmc)(int *nboctets);                                          // xvuelc.c:242 (returns void*)
void  proc(dstnmc)(void *mcoctets);                                         // xvuelc.c:254
void  proc(nomrepmefisto)(char *chaine, int *size);                         // xvuelc.c:266

// ---- 2. Graphics init / screen geometry ----
void  proc(xvinitgraphique)(void);                                          // xvuelc.c:286
void  proc(xtinit)(void);                                                   // xvuelc.c:306
void  proc(xvpxecran)(int *xp, int *yp);                                    // xvuelc.c:319
void  proc(xvmmecran)(int *xmm, int *ymm);                                  // xvuelc.c:337

// ---- 3. Accrochage (crosshair/snap) ----
void  proc(initaccrochage)(void);                                           // xvuelc.c:561

// ---- 4. xvinfo (14-arg mega-init — signature verbatim from xvuelc.c:612) ----
void  proc(xvinfo)( int *ix, int *iy, int *maxfonts,
                    int *n1coref, int *ndcoref, int *n1coelf,
                    int *ndcoelf, int *n1coulf, int *ndcoulf, int *nbcolo,
                    char namefonts[][256], int nbchar[], int *nbfonts, int *visuclass );

// ---- 5. Color palette (Fortran-facing helpers only) ----
void  proc(xvrecuprgbdec)(int *nbcolor, float *r, float *g, float *b);      // xvuelc.c:1044
// xvactivervb signature verbatim from xvuelc.c:1072
void  proc(xvactivervb)( int *palcour, int *nbcells,
                         float r[], float g[], float b[] );
void  proc(xvcouleur)(int *icolor);                                         // xvuelc.c:1119
void  proc(xvpostscript)(int *lasops);                                      // xvuelc.c:1187

// ---- 6. Pixmap save / restore ----
void  proc(fenetremempx)(void);                                             // xvuelc.c:1307
void  proc(mempxfenetre)(void);                                             // xvuelc.c:1321
void  proc(sauvefenetre)(void);                                             // xvuelc.c:1335
void  proc(restaurefenetre)(void);                                          // xvuelc.c:1350
void  proc(sauvemempx)(void);                                               // xvuelc.c:1365
void  proc(restauremempx)(void);                                            // xvuelc.c:1380
void  proc(effacemempx)(void);                                              // xvuelc.c:1395

// ---- 7. Clear / background / font ----
void  proc(effacer)(void);                                                  // xvuelc.c:1413
void  proc(xvfond)(int *icolor);                                            // xvuelc.c:1434
void  proc(xvchargefonte)(int *nofont0, int *nofont, int *largpx, int *hautpx); // xvuelc.c:1463
void  proc(xvnbpixeltexte)(char *texte, int *length, int *nbpxla, int *nbpxha);  // xvuelc.c:1576

// ---- 8. Window close / query ----
void  proc(xvfermer)(void);                                                 // xvuelc.c:1602
void  proc(xvpxfenetre)(int *x, int *y);                                    // xvuelc.c:1619

// ---- 9. Text rendering ----
void  proc(xvftexte)(char string[], int *length, int *x1, int *y1);         // xvuelc.c:1642
void  proc(xvtexte)(char string[], int *length, int *x1, int *y1);          // xvuelc.c:1662

// ---- 10. Filled polygons (use MefistoPoint shim) ----
void  proc(xvface)(int *n, MefistoPoint *pts);                              // xvuelc.c:1701

// ---- 11. Pen style / width ----
void  proc(xvtypetrait)(int *ptype);                                        // xvuelc.c:1760
void  proc(xvepaisseur)(int *pepais);                                       // xvuelc.c:1807

// ---- 12. Line drawing ----
void  proc(xvftrait)(int *x1, int *y1, int *x2, int *y2);                   // xvuelc.c:1845
void  proc(xvtrait)(int *x1, int *y1, int *x2, int *y2);                    // xvuelc.c:1862
void  proc(xvtraits)(int *nbpoints, MefistoPoint *points);                  // xvuelc.c:1977
void  proc(xvfacetraits)(int *ncf, int *nca, int *n, MefistoPoint *pts);    // xvuelc.c:2035

// ---- 13. Event / mouse ----
void  proc(xvsouris)(int *notypeevent, int *nbc, int *x1, int *y1);         // xvuelc.c:2123
void  proc(xvsouris2)(int *items, int *pmin0, int *notypeevent, int *ibutton, int *x1, int *y1); // xvuelc.c:2230
void  proc(deplsouris)(int *x, int *y);                                     // xvuelc.c:2374
void  proc(xvvoir)(void);                                                   // xvuelc.c:2384
void  proc(xvpause)(void);                                                  // xvuelc.c:2408

// ---- 14. Rectangles / ellipses ----
void  proc(xvfbordrectangle)(int *x, int *y, int *width, int *height);      // xvuelc.c:2426
void  proc(xvbordrectangle)(int *x, int *y, int *width, int *height);       // xvuelc.c:2443
void  proc(xvfrectangle)(int *x, int *y, int *width, int *height);          // xvuelc.c:2489
void  proc(xvrectangle)(int *x, int *y, int *width, int *height);           // xvuelc.c:2507
void  proc(xvbordarcellipse)(int *x, int *y, int *width, int *height, int *a1, int *a2); // xvuelc.c:2554
void  proc(xvarcellipse)(int *x, int *y, int *width, int *height, int *a1, int *a2);     // xvuelc.c:2616

// ---- 15. Time / date / host / env ----
void  proc(tempscpu)(double *tclock);                                       // xvuelc.c:2678
void  proc(secondes1970)(double *secondes);                                 // xvuelc.c:2694
void  proc(secondes1969)(double *secondes);                                 // xvuelc.c:2712
void  proc(nomordinateurhote)(char *host, int *nbcar);                      // xvuelc.c:2728
void  proc(ladate)(int *a, int *m, int *j);                                 // xvuelc.c:2753
void  proc(heureminuteseconde)(int *h, int *m, int *s, int *millis);        // xvuelc.c:2779
// valvarenv signature verbatim from xvuelc.c:2805
void  proc(valvarenv)( char *nom, int *lval_admis,
                       char *val, int *lval_trouve );

// ---- 16. PostScript file-mode ----
void  proc(xvinitierps)(int *modeps);                                       // xvuelc.c:2844
void  proc(xvimprimerps)(char nomfichier[], int *length);                   // xvuelc.c:2906
void  proc(xvsauverps)(char nomfichier[], int *length);                     // xvuelc.c:2954

} // extern "C"

#endif // XVUE_QT_API_H
