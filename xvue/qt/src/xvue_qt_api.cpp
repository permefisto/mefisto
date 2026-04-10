// xvue/qt/src/xvue_qt_api.cpp
// Phase 0: warn-once no-op stubs for every Fortran-facing extern "C" entry
// point in xvue_qt_api.h. Each stub prints a single stderr line on first
// call, then is silent for the rest of the process.
//
// Per D-18: do NOT abort, do NOT set error flags, do NOT touch any Qt object.
// Non-void returns use the simplest safe default (nullptr, 0.0).
//
// Per D-17: warn-once diagnostic pattern is a per-function `static bool warned`.

#include <cstdio>
#include "xvue_qt_api.h"

namespace {

inline void warn_once(bool &flag, const char *name) {
    if (!flag) {
        std::fprintf(stderr, "xvue-qt: stub %s not implemented yet\n", name);
        flag = true;
    }
}

} // anonymous namespace

extern "C" {

// ---- 1. languemefisto_ ----
void proc(languemefisto)(int *langue) {
    static bool warned = false;
    warn_once(warned, "languemefisto_");
    (void)langue;
}

// ---- 2. dctnmc_ (returns void*) ----
void *proc(dctnmc)(int *nboctets) {
    static bool warned = false;
    warn_once(warned, "dctnmc_");
    (void)nboctets;
    return nullptr;
}

// ---- 3. dstnmc_ ----
void proc(dstnmc)(void *mcoctets) {
    static bool warned = false;
    warn_once(warned, "dstnmc_");
    (void)mcoctets;
}

// ---- 4. nomrepmefisto_ ----
void proc(nomrepmefisto)(char *chaine, int *size) {
    static bool warned = false;
    warn_once(warned, "nomrepmefisto_");
    (void)chaine; (void)size;
}

// ---- 5. xvinitgraphique_ ----
void proc(xvinitgraphique)(void) {
    static bool warned = false;
    warn_once(warned, "xvinitgraphique_");
}

// ---- 6. xtinit_ ----
void proc(xtinit)(void) {
    static bool warned = false;
    warn_once(warned, "xtinit_");
}

// ---- 7. xvpxecran_ ----
void proc(xvpxecran)(int *xp, int *yp) {
    static bool warned = false;
    warn_once(warned, "xvpxecran_");
    (void)xp; (void)yp;
}

// ---- 8. xvmmecran_ ----
void proc(xvmmecran)(int *xmm, int *ymm) {
    static bool warned = false;
    warn_once(warned, "xvmmecran_");
    (void)xmm; (void)ymm;
}

// ---- 9. initaccrochage_ ----
void proc(initaccrochage)(void) {
    static bool warned = false;
    warn_once(warned, "initaccrochage_");
}

// ---- 10. xvinfo_ ----
void proc(xvinfo)( int *ix, int *iy, int *maxfonts,
                   int *n1coref, int *ndcoref, int *n1coelf,
                   int *ndcoelf, int *n1coulf, int *ndcoulf, int *nbcolo,
                   char namefonts[][256], int nbchar[], int *nbfonts, int *visuclass ) {
    static bool warned = false;
    warn_once(warned, "xvinfo_");
    (void)ix; (void)iy; (void)maxfonts;
    (void)n1coref; (void)ndcoref; (void)n1coelf;
    (void)ndcoelf; (void)n1coulf; (void)ndcoulf; (void)nbcolo;
    (void)namefonts; (void)nbchar; (void)nbfonts; (void)visuclass;
}

// ---- 11. xvrecuprgbdec_ ----
void proc(xvrecuprgbdec)(int *nbcolor, float *r, float *g, float *b) {
    static bool warned = false;
    warn_once(warned, "xvrecuprgbdec_");
    (void)nbcolor; (void)r; (void)g; (void)b;
}

// ---- 12. xvactivervb_ ----
void proc(xvactivervb)( int *palcour, int *nbcells,
                        float r[], float g[], float b[] ) {
    static bool warned = false;
    warn_once(warned, "xvactivervb_");
    (void)palcour; (void)nbcells; (void)r; (void)g; (void)b;
}

// ---- 13. xvcouleur_ ----
void proc(xvcouleur)(int *icolor) {
    static bool warned = false;
    warn_once(warned, "xvcouleur_");
    (void)icolor;
}

// ---- 14. xvpostscript_ ----
void proc(xvpostscript)(int *lasops) {
    static bool warned = false;
    warn_once(warned, "xvpostscript_");
    (void)lasops;
}

// ---- 15. fenetremempx_ ----
void proc(fenetremempx)(void) {
    static bool warned = false;
    warn_once(warned, "fenetremempx_");
}

// ---- 16. mempxfenetre_ ----
void proc(mempxfenetre)(void) {
    static bool warned = false;
    warn_once(warned, "mempxfenetre_");
}

// ---- 17. sauvefenetre_ ----
void proc(sauvefenetre)(void) {
    static bool warned = false;
    warn_once(warned, "sauvefenetre_");
}

// ---- 18. restaurefenetre_ ----
void proc(restaurefenetre)(void) {
    static bool warned = false;
    warn_once(warned, "restaurefenetre_");
}

// ---- 19. sauvemempx_ ----
void proc(sauvemempx)(void) {
    static bool warned = false;
    warn_once(warned, "sauvemempx_");
}

// ---- 20. restauremempx_ ----
void proc(restauremempx)(void) {
    static bool warned = false;
    warn_once(warned, "restauremempx_");
}

// ---- 21. effacemempx_ ----
void proc(effacemempx)(void) {
    static bool warned = false;
    warn_once(warned, "effacemempx_");
}

// ---- 22. effacer_ ----
void proc(effacer)(void) {
    static bool warned = false;
    warn_once(warned, "effacer_");
}

// ---- 23. xvfond_ ----
void proc(xvfond)(int *icolor) {
    static bool warned = false;
    warn_once(warned, "xvfond_");
    (void)icolor;
}

// ---- 24. xvchargefonte_ ----
void proc(xvchargefonte)(int *nofont0, int *nofont, int *largpx, int *hautpx) {
    static bool warned = false;
    warn_once(warned, "xvchargefonte_");
    (void)nofont0; (void)nofont; (void)largpx; (void)hautpx;
}

// ---- 25. xvnbpixeltexte_ ----
void proc(xvnbpixeltexte)(char *texte, int *length, int *nbpxla, int *nbpxha) {
    static bool warned = false;
    warn_once(warned, "xvnbpixeltexte_");
    (void)texte; (void)length; (void)nbpxla; (void)nbpxha;
}

// ---- 26. xvfermer_ ----
void proc(xvfermer)(void) {
    static bool warned = false;
    warn_once(warned, "xvfermer_");
}

// ---- 27. xvpxfenetre_ ----
void proc(xvpxfenetre)(int *x, int *y) {
    static bool warned = false;
    warn_once(warned, "xvpxfenetre_");
    (void)x; (void)y;
}

// ---- 28. xvftexte_ ----
void proc(xvftexte)(char string[], int *length, int *x1, int *y1) {
    static bool warned = false;
    warn_once(warned, "xvftexte_");
    (void)string; (void)length; (void)x1; (void)y1;
}

// ---- 29. xvtexte_ ----
void proc(xvtexte)(char string[], int *length, int *x1, int *y1) {
    static bool warned = false;
    warn_once(warned, "xvtexte_");
    (void)string; (void)length; (void)x1; (void)y1;
}

// ---- 30. xvface_ ----
void proc(xvface)(int *n, MefistoPoint *pts) {
    static bool warned = false;
    warn_once(warned, "xvface_");
    (void)n; (void)pts;
}

// ---- 31. xvtypetrait_ ----
void proc(xvtypetrait)(int *ptype) {
    static bool warned = false;
    warn_once(warned, "xvtypetrait_");
    (void)ptype;
}

// ---- 32. xvepaisseur_ ----
void proc(xvepaisseur)(int *pepais) {
    static bool warned = false;
    warn_once(warned, "xvepaisseur_");
    (void)pepais;
}

// ---- 33. xvftrait_ ----
void proc(xvftrait)(int *x1, int *y1, int *x2, int *y2) {
    static bool warned = false;
    warn_once(warned, "xvftrait_");
    (void)x1; (void)y1; (void)x2; (void)y2;
}

// ---- 34. xvtrait_ ----
void proc(xvtrait)(int *x1, int *y1, int *x2, int *y2) {
    static bool warned = false;
    warn_once(warned, "xvtrait_");
    (void)x1; (void)y1; (void)x2; (void)y2;
}

// ---- 35. xvtraits_ ----
void proc(xvtraits)(int *nbpoints, MefistoPoint *points) {
    static bool warned = false;
    warn_once(warned, "xvtraits_");
    (void)nbpoints; (void)points;
}

// ---- 36. xvfacetraits_ ----
void proc(xvfacetraits)(int *ncf, int *nca, int *n, MefistoPoint *pts) {
    static bool warned = false;
    warn_once(warned, "xvfacetraits_");
    (void)ncf; (void)nca; (void)n; (void)pts;
}

// ---- 37. xvsouris_ ----
void proc(xvsouris)(int *notypeevent, int *nbc, int *x1, int *y1) {
    static bool warned = false;
    warn_once(warned, "xvsouris_");
    (void)notypeevent; (void)nbc; (void)x1; (void)y1;
}

// ---- 38. xvsouris2_ ----
void proc(xvsouris2)(int *items, int *pmin0, int *notypeevent, int *ibutton, int *x1, int *y1) {
    static bool warned = false;
    warn_once(warned, "xvsouris2_");
    (void)items; (void)pmin0; (void)notypeevent; (void)ibutton; (void)x1; (void)y1;
}

// ---- 39. deplsouris_ ----
void proc(deplsouris)(int *x, int *y) {
    static bool warned = false;
    warn_once(warned, "deplsouris_");
    (void)x; (void)y;
}

// ---- 40. xvvoir_ ----
void proc(xvvoir)(void) {
    static bool warned = false;
    warn_once(warned, "xvvoir_");
}

// ---- 41. xvpause_ ----
void proc(xvpause)(void) {
    static bool warned = false;
    warn_once(warned, "xvpause_");
}

// ---- 42. xvfbordrectangle_ ----
void proc(xvfbordrectangle)(int *x, int *y, int *width, int *height) {
    static bool warned = false;
    warn_once(warned, "xvfbordrectangle_");
    (void)x; (void)y; (void)width; (void)height;
}

// ---- 43. xvbordrectangle_ ----
void proc(xvbordrectangle)(int *x, int *y, int *width, int *height) {
    static bool warned = false;
    warn_once(warned, "xvbordrectangle_");
    (void)x; (void)y; (void)width; (void)height;
}

// ---- 44. xvfrectangle_ ----
void proc(xvfrectangle)(int *x, int *y, int *width, int *height) {
    static bool warned = false;
    warn_once(warned, "xvfrectangle_");
    (void)x; (void)y; (void)width; (void)height;
}

// ---- 45. xvrectangle_ ----
void proc(xvrectangle)(int *x, int *y, int *width, int *height) {
    static bool warned = false;
    warn_once(warned, "xvrectangle_");
    (void)x; (void)y; (void)width; (void)height;
}

// ---- 46. xvbordarcellipse_ ----
void proc(xvbordarcellipse)(int *x, int *y, int *width, int *height, int *a1, int *a2) {
    static bool warned = false;
    warn_once(warned, "xvbordarcellipse_");
    (void)x; (void)y; (void)width; (void)height; (void)a1; (void)a2;
}

// ---- 47. xvarcellipse_ ----
void proc(xvarcellipse)(int *x, int *y, int *width, int *height, int *a1, int *a2) {
    static bool warned = false;
    warn_once(warned, "xvarcellipse_");
    (void)x; (void)y; (void)width; (void)height; (void)a1; (void)a2;
}

// ---- 48. tempscpu_ ----
void proc(tempscpu)(double *tclock) {
    static bool warned = false;
    warn_once(warned, "tempscpu_");
    (void)tclock;
}

// ---- 49. secondes1970_ ----
void proc(secondes1970)(double *secondes) {
    static bool warned = false;
    warn_once(warned, "secondes1970_");
    (void)secondes;
}

// ---- 50. secondes1969_ ----
void proc(secondes1969)(double *secondes) {
    static bool warned = false;
    warn_once(warned, "secondes1969_");
    (void)secondes;
}

// ---- 51. nomordinateurhote_ ----
void proc(nomordinateurhote)(char *host, int *nbcar) {
    static bool warned = false;
    warn_once(warned, "nomordinateurhote_");
    (void)host; (void)nbcar;
}

// ---- 52. ladate_ ----
void proc(ladate)(int *a, int *m, int *j) {
    static bool warned = false;
    warn_once(warned, "ladate_");
    (void)a; (void)m; (void)j;
}

// ---- 53. heureminuteseconde_ ----
void proc(heureminuteseconde)(int *h, int *m, int *s, int *millis) {
    static bool warned = false;
    warn_once(warned, "heureminuteseconde_");
    (void)h; (void)m; (void)s; (void)millis;
}

// ---- 54. valvarenv_ ----
void proc(valvarenv)( char *nom, int *lval_admis,
                      char *val, int *lval_trouve ) {
    static bool warned = false;
    warn_once(warned, "valvarenv_");
    (void)nom; (void)lval_admis; (void)val; (void)lval_trouve;
}

// ---- 55. xvinitierps_ ----
void proc(xvinitierps)(int *modeps) {
    static bool warned = false;
    warn_once(warned, "xvinitierps_");
    (void)modeps;
}

// ---- 56. xvimprimerps_ ----
void proc(xvimprimerps)(char nomfichier[], int *length) {
    static bool warned = false;
    warn_once(warned, "xvimprimerps_");
    (void)nomfichier; (void)length;
}

// ---- 57. xvsauverps_ ----
void proc(xvsauverps)(char nomfichier[], int *length) {
    static bool warned = false;
    warn_once(warned, "xvsauverps_");
    (void)nomfichier; (void)length;
}

} // extern "C"
