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
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <unistd.h>
#include "xvue_qt_api.h"
#include "xvue_qt_app.h"
#include "xvue_qt_window.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_state.h"
#include "xvue_qt_event.h"
#include "xvue_qt_menu_bridge.h"   // Phase 6.0 Plan 03: pre-AUTOEXIT pre-drain
#include "xvue_qt_console_dock.h"  // Phase 6.0 Plan 06: stdout redirect dispatch
#include "xvue_qt_prefs.h"         // Phase 6.0 Plan 06: per-module persistence init
#include "xvue_qt_postscript.h"    // Phase 7 Plan 02: PsEmitter dispatch
#include <QApplication>
#include <QCoreApplication>
#include <QCursor>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QFont>
#include <QFontMetrics>
#include <QGuiApplication>
#include <QLatin1Char>
#include <QPainter>
#include <QPixmap>
#include <QPoint>
#include <QPolygon>
#include <QRect>
#include <QScreen>
#include <QString>
#include <QWindow>
#include <vector>

// Forward declaration so the warn-once stubs in the anonymous namespace can
// reference XvueWindow/XvueMenuBridge by pointer. The full headers are
// already pulled in via the includes above.
class XvueWindow;
class XvueMenuBridge;

// Auto-batch on non-tty: ppmail/etc. read from stdin BEFORE xvinitgraphique_
// runs (project name, INTERA flag, ...). When launched without a controlling
// terminal (KDE launcher / kwin-mcp / pipeline / direct exec) those reads
// block forever and the Qt window appears alive but frozen. Detect non-tty
// at process start and set MEFISTO_BATCH_X11=1 + MEFISTO_XVSOURIS_AUTOEXIT=1
// (existing escape hatches from Phase 7/8) so Fortran takes the batch path.
// The user can override by setting MEFISTO_BATCH_X11=0 explicitly before
// launch — interactive terminal users keep stdin reads as today.
__attribute__((constructor))
static void xvue_qt_auto_batch_on_non_tty(void) {
    if (!isatty(fileno(stdin))) {
        if (getenv("MEFISTO_BATCH_X11") == nullptr) {
            setenv("MEFISTO_BATCH_X11", "1", 0);
        }
        if (getenv("MEFISTO_XVSOURIS_AUTOEXIT") == nullptr) {
            setenv("MEFISTO_XVSOURIS_AUTOEXIT", "1", 0);
        }
    }
}

namespace {

inline void warn_once(bool &flag, const char *name) {
    if (!flag) {
        std::fprintf(stderr, "xvue-qt: stub %s not implemented yet\n", name);
        flag = true;
    }
}

// Phase 6.0 Plan 06 + Phase 6.1 Plan 02: warn-once stubs for per-module
// action registration. 6.1..6.5 each ship a stronger symbol (a real
// registerXxxActions body in xvue_qt_<module>_actions.cpp) that displaces
// the stub here via the __attribute__((weak)) attribute below.
//
// Phase 6.1 Plan 02 deviation (Rule 3 — blocking fix): in 6.0 these bodies
// lived INSIDE the anonymous namespace, which gave them internal linkage
// and made the weak-override pattern non-functional (a strong symbol in
// another TU cannot displace an internal-linkage symbol). We move them
// OUT of the anonymous namespace, declare them `extern "C"` with
// `__attribute__((weak))`, and leave the signature unchanged. The
// dispatch call-site at xvue_module_init_ stays untouched — it resolves
// by unadorned name lookup from inside `extern "C" {`.
//
// Pure 6.0 builds (no module TU linked into libxvueqt.a) keep the
// warn-once no-op behaviour. 6.1+ builds displace the relevant symbol.
bool warn_register_mail_done_ = false;
bool warn_register_elas_done_ = false;
bool warn_register_flui_done_ = false;
bool warn_register_ther_done_ = false;
bool warn_register_nlse_done_ = false;

// Phase 6.0 hot-fix (Plan 05 + Plan 06 followup): suppress the empty-state
// hint as soon as Fortran issues any drawing primitive on the backing pixmap.
//
// Background: Plan 06's xvue_module_init_ flips state->has_user_content_,
// but pure 6.0 builds (mesher / solvers) don't yet CALL xvue_module_init_
// (that's added in 6.1..6.5). Meanwhile the legacy mesher paints chrome
// (POINTS & LINES menu, project header, button rasters) on the canvas
// during startup. Without this hook, the Plan 05 empty-state hint
// ("No project open / Choose File -> Open Project (Ctrl+O)") overlays
// the live mesher menu — visually misleading.
//
// Idempotent: first Fortran draw flips the flag and triggers a single
// repaint. Subsequent calls are a cheap branch + early return.
inline void xvue_qt_mark_user_content(XvueWindow* win) {
    if (!win) return;
    auto* st = win->state();
    if (!st || st->has_user_content_) return;
    st->has_user_content_ = true;
    if (win->canvas()) win->canvas()->update();
}

// File-local helper (D-13). Not part of the Fortran ABI.
// All four xv*rectangle_ entry points route through this one implementation.
//
// CR-01 (2026-04-11): legacy X11 distinguishes outline-only
// (XDrawRectangle: xvfbordrectangle_/xvbordrectangle_) from fill-only
// (XFillRectangle: xvfrectangle_/xvrectangle_). QPainter::drawRect strokes
// AND fills in one call, so the four symbols must dispatch on an explicit
// mode and use Qt::NoBrush / fillRect to keep the legacy semantics when
// pen and brush colors diverge (Phase 3).
enum class RectMode { Outline, Fill };

inline void xvue_qt_draw_rect_common(int x, int y, int w, int h, RectMode mode,
                                     bool emit_ps = true) {
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    const QRect r(x, y, w, h);
    if (mode == RectMode::Outline) {
        // XDrawRectangle: stroke only, no fill.
        QBrush saved = st->painter_->brush();
        st->painter_->setBrush(Qt::NoBrush);
        st->painter_->drawRect(r);
        st->painter_->setBrush(saved);
        // Phase 7 Plan 03: PS "r" opcode emit only for the persistent
        // (mempx-equivalent) variant. xvuelc.c xvfbordrectangle (line 2556)
        // draws to fenetre_mef and emits NO PS; xvbordrectangle (line 2573)
        // draws to mempx and DOES emit. Under Phase 2's single-backing
        // collapse both Qt entries share this helper, so the emit_ps flag
        // preserves the legacy PS-stream byte parity.
        if (emit_ps) XvueApp::psEmitter().bordrectangle(x, y, w, h);
    } else {
        // XFillRectangle: fill only, no outline.
        st->painter_->fillRect(r, st->painter_->brush());
        // Phase 7 Plan 03: PS "R" opcode emit for the persistent variant only.
        if (emit_ps) XvueApp::psEmitter().rectangle(x, y, w, h);
    }
    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) win->canvas()->update();
}

// Phase 4 (D-07, D-08, D-09): file-local helpers backing the four save/restore
// entry points. Anonymous-namespace visibility keeps them TU-local so the
// verify_abi ABI count stays at 57 (Pitfall 7). `inline` hints the linker.
//
// Qt 6 allows two QPainters on two DIFFERENT devices concurrently
// (doc.qt.io/qt-6/qpainter.html "A paint device can only be painted by
// one painter at a time" — per-device, not per-process). So the save
// direction uses a temporary scoped painter on saved_canvas_ while the
// long-lived painter_ stays bound to backing_ (Phase 2 D-05 invariant).
inline void xvue_qt_save_to_slot() {
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive() || !st->backing_) return;

    // D-02/D-14: lazy (re)allocation. If slot is null or wrong size,
    // drop it and allocate fresh. HiDPI via setDevicePixelRatio (02/D-06).
    if (!st->saved_canvas_ || st->saved_canvas_->size() != st->backing_->size()) {
        delete st->saved_canvas_;
        st->saved_canvas_ = new QPixmap(st->backing_->size());
    }
    // WR-02: Always refresh DPR, even on slot reuse. Qt 6 can change
    // backing_->devicePixelRatio() when the window moves between monitors
    // without changing backing_->size() (size is in device pixels), so the
    // reuse branch must re-sync DPR or restore would misscale. Cheap and
    // idempotent.
    st->saved_canvas_->setDevicePixelRatio(st->backing_->devicePixelRatio());

    // Scoped temporary painter on saved_canvas_; painter_ on backing_
    // is NOT touched (Phase 2 D-05; Pitfall 2).
    {
        QPainter tmp(st->saved_canvas_);
        tmp.drawPixmap(0, 0, *st->backing_);
    }
    // No canvas_->update() — save is invisible to the user.
}

inline void xvue_qt_restore_from_slot() {
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive() || !st->backing_) return;

    // D-12: size mismatch is auditable no-op with stderr warning.
    if (!st->saved_canvas_ || st->saved_canvas_->size() != st->backing_->size()) {
        std::fprintf(stderr, "xvue-qt: restore_from_slot: no slot or size mismatch\n");
        return;
    }

    // Restore direction uses the ACTIVE painter_ on backing_ (02/D-05).
    st->painter_->drawPixmap(0, 0, *st->saved_canvas_);
    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) win->canvas()->update();
}

// Phase 3 D-06 + Pitfall 7: xvtexte_ and xvftexte_ share this body under
// the Phase 2 single-backing model (02/D-05). In legacy xvuelc.c the two
// drew to DIFFERENT targets (mempx at :1658 vs fenetre_mef at :1678); the
// collapse here is legal ONLY because Phase 2 unified the draw target.
// If a future phase reintroduces multiple draw targets, revisit.
//
// Phase 7 Plan 03 (EXPORT-04): after the QPainter::drawText, dispatch into
// PsEmitter::texte for the PostScript "T" emit. Helper internally checks
// active() so this is a cheap no-op when no PS capture is in progress.
inline void xvue_qt_draw_text_common(const char* string, int length,
                                     int x1, int y1) {
    if (!string || length <= 0) return;              // WR-03 null-guard
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive() || !st->backing_) return;

    // Pitfall 4: TWO-ARG fromLatin1 — the legacy Fortran strings are NOT
    // NUL-terminated; explicit length is mandatory.
    QString qstr = QString::fromLatin1(string, length);
    st->painter_->setFont(st->current_font_);        // cheap no-op if unchanged
    st->painter_->drawText(x1, y1, qstr);            // baseline form (D-06)

    // Phase 7 Plan 03: PostScript emit after the QPainter call. Pass canvas-Y
    // (Y-down); the helper applies pyFlip internally per README_COORDS.md.
    XvueApp::psEmitter().texte(string, length, x1, y1);

    // Phase 2 D-01 epilogue: drawing primitives flush.
    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

} // anonymous namespace

// Phase 6.1 Plan 02 (Rule 3 blocking fix): the five per-module action
// registration stubs live at file scope with `extern "C"` linkage and
// `__attribute__((weak))` so a stronger symbol in
// xvue/qt/src/xvue_qt_<mod>_actions.cpp can displace them at link time.
//
// NOTE: declared `extern "C"` even though they are NOT part of the
// Fortran ABI (the call site is C++-side in xvue_module_init_ below,
// not a Fortran CALL). The C linkage is what lets the strong-symbol
// body in xvue_qt_mail_actions.cpp — also declared `extern "C"` —
// link against the same unmangled name. verify_abi.sh counts symbols
// whose name matches `[a-zA-Z_]\w*_$` AND is ` T ` (text / strong);
// these five remain ` W ` (weak) in 6.0 builds so they do NOT bump
// the frozen ABI count of 58. A 6.1 Qt build displaces `_mail` to
// ` T `, while `_elas / _flui / _ther / _nlse` stay ` W ` — so a grep
// for ` T ` still counts 58 (the pure Fortran ABI) + 1 (mail_actions
// strong) — see ACCEPTANCE grep in plan, which specifically targets
// the mail override only.
extern "C" void registerMailActions_stub_(XvueWindow*, XvueMenuBridge*)
    __attribute__((weak));
extern "C" void registerElasActions_stub_(XvueWindow*, XvueMenuBridge*)
    __attribute__((weak));
extern "C" void registerFluiActions_stub_(XvueWindow*, XvueMenuBridge*)
    __attribute__((weak));
extern "C" void registerTherActions_stub_(XvueWindow*, XvueMenuBridge*)
    __attribute__((weak));
extern "C" void registerNlseActions_stub_(XvueWindow*, XvueMenuBridge*)
    __attribute__((weak));

extern "C" void registerMailActions_stub_(XvueWindow*, XvueMenuBridge*) {
    if (!warn_register_mail_done_) {
        warn_register_mail_done_ = true;
        std::fprintf(stderr,
            "xvue-qt: registerMailActions stub (Phase 6.1 adds the real body).\n");
    }
}
extern "C" void registerElasActions_stub_(XvueWindow*, XvueMenuBridge*) {
    if (!warn_register_elas_done_) {
        warn_register_elas_done_ = true;
        std::fprintf(stderr,
            "xvue-qt: registerElasActions stub (Phase 6.2 adds the real body).\n");
    }
}
extern "C" void registerFluiActions_stub_(XvueWindow*, XvueMenuBridge*) {
    if (!warn_register_flui_done_) {
        warn_register_flui_done_ = true;
        std::fprintf(stderr,
            "xvue-qt: registerFluiActions stub (Phase 6.3 adds the real body).\n");
    }
}
extern "C" void registerTherActions_stub_(XvueWindow*, XvueMenuBridge*) {
    if (!warn_register_ther_done_) {
        warn_register_ther_done_ = true;
        std::fprintf(stderr,
            "xvue-qt: registerTherActions stub (Phase 6.4 adds the real body).\n");
    }
}
extern "C" void registerNlseActions_stub_(XvueWindow*, XvueMenuBridge*) {
    if (!warn_register_nlse_done_) {
        warn_register_nlse_done_ = true;
        std::fprintf(stderr,
            "xvue-qt: registerNlseActions stub (Phase 6.5 adds the real body).\n");
    }
}

extern "C" {

// ---- 1. languemefisto_ ----
// Returns 1 (English) if $MEFISTO/td/m/anglais exists, 0 (French) otherwise.
// Legacy xvuelc.c hardcoded /usr/local/mefisto; Qt port resolves $MEFISTO env.
void proc(languemefisto)(int *langue) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (langue == nullptr) return;
    const char *mefisto = getenv("MEFISTO");
    char path[1024];
    if (mefisto != nullptr) {
        std::snprintf(path, sizeof(path), "%s/td/m/anglais", mefisto);
    } else {
        std::strncpy(path, "/usr/local/mefisto/td/m/anglais", sizeof(path) - 1);
        path[sizeof(path) - 1] = '\0';
    }
    FILE *fp = std::fopen(path, "r");
    if (fp != nullptr) { *langue = 1; std::fclose(fp); }
    else               { *langue = 0; }
}

// ---- 2. dctnmc_ (returns void*) ----
// Allocate nboctets bytes; return pointer or NULL on failure.
void *proc(dctnmc)(int *nboctets) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (nboctets == nullptr || *nboctets <= 0) return nullptr;
    return std::malloc(static_cast<size_t>(*nboctets));
}

// ---- 3. dstnmc_ ----
// Free a buffer previously allocated by dctnmc_.
void proc(dstnmc)(void *mcoctets) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    std::free(mcoctets);
}

// ---- 4. nomrepmefisto_ ----
// Stash the MEFISTO root dir from Fortran so any future consumer can read it
// back via nomrepmefisto_buffer. Mirrors legacy xvuelc.c file-scope nom_homdir.
static char nomrepmefisto_buffer[256] = {0};
void proc(nomrepmefisto)(char *chaine, int *size) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (chaine == nullptr || size == nullptr || *size < 0) return;
    const int n = (*size < (int)(sizeof(nomrepmefisto_buffer) - 1))
                  ? *size : (int)(sizeof(nomrepmefisto_buffer) - 1);
    for (int i = 0; i < n; ++i) nomrepmefisto_buffer[i] = chaine[i];
    nomrepmefisto_buffer[n] = '\0';
}

// ---- 5. xvinitgraphique_ ----
void proc(xvinitgraphique)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();

    auto& win = XvueApp::window_slot();
    if (!win) {
        // D-07: lazy (re)allocation. QApplication is NOT recreated — call_once
        // guarantees exactly-one. D-02: title "MEFISTO", 800x600 set in ctor.
        win = std::make_unique<XvueWindow>();
    }
    win->show();
    win->raise();
    win->activateWindow();

    // Phase 9 Plan 9-06 (carry-forward #1): matched-dim Qt recapture.
    // Per RESEARCH §Pattern 4 + §Pitfall 5: honor MEFISTO_QT_WINDOW_SIZE=WIDTHxHEIGHT
    // env var, BUT ONLY in headless contexts (MEFISTO_BATCH_X11=1 OR
    // QT_QPA_PLATFORM=offscreen). Interactive sessions retain default sizing.
    // Use setMinimumSize+setMaximumSize on the canvas (NOT setFixedSize on the
    // window) so xvfermer_ can clear the constraints and a subsequent reopen
    // cycle does not inherit them. Constraints applied here cause the canvas
    // backing pixmap to be allocated at the env-specified dim during the
    // first XvueCanvas::resizeEvent triggered by win->adjustSize() / show(),
    // and clip any subsequent xvinfo_ -> win->resize(*ix, *iy) call. T-09-06
    // mitigation: implementation is env-var only (no new extern "C" entry —
    // ABI count remains 58). T-09-05-A mitigation: sscanf + 0 < {w,h} < 8192
    // bounds check; malformed silently ignored.
    {
        const char* batch = std::getenv("MEFISTO_BATCH_X11");
        const bool headless_batch = batch && std::strcmp(batch, "1") == 0;
        const char* qpa = std::getenv("QT_QPA_PLATFORM");
        const bool offscreen = qpa && std::strcmp(qpa, "offscreen") == 0;
        const char* qt_window_size = std::getenv("MEFISTO_QT_WINDOW_SIZE");
        if ((headless_batch || offscreen) && qt_window_size && qt_window_size[0]) {
            int w = 0, h = 0;
            if (std::sscanf(qt_window_size, "%dx%d", &w, &h) == 2 &&
                w > 0 && w < 8192 && h > 0 && h < 8192) {
                if (win && win->canvas()) {
                    win->canvas()->setMinimumSize(w, h);
                    win->canvas()->setMaximumSize(w, h);
                    win->adjustSize();
                }
            }
        }
    }

    // D-01 revised 2026-04-11 (debug phase-01-xvtest0-teardown-segfault):
    // the original "exactly one processEvents" rule is insufficient on X11
    // to actually realize a window (MapRequest, ConfigureNotify, Expose
    // each need their own trip through the event loop). We pump in a
    // bounded loop until QWindow::isExposed() is true or the 2 s timeout
    // fires (so a broken display cannot hang the Fortran main). The
    // ExcludeUserInputEvents flag remains mandatory (Pitfall 4).
    QElapsedTimer t;
    t.start();
    const int timeout_ms = 2000;
    QWindow* wh = win->windowHandle();
    while (t.elapsed() < timeout_ms) {
        QCoreApplication::processEvents(
            QEventLoop::ExcludeUserInputEvents, 20);
        wh = win->windowHandle();  // may become non-null after first pump
        if (wh && wh->isExposed()) break;
    }
}

// ---- 6. xtinit_ ----
// Phase 03.1 Pitfall 1 fix: promoted from warn-once stub to a real
// body so that XVOUVRIR -> XTINIT -> (XVINIT chain) opens a Qt window
// on the Qt backend the same way XVINITGRAPHIQUE does on the xvtest0
// code path. Idempotency is provided by window_slot()'s "if (!win)"
// lazy alloc — calling XVINITGRAPHIQUE after XTINIT (or in any order)
// does NOT double-create. XvueApp::ensure() is guarded by call_once
// per D-07. See .planning/phases/03.1-.../03.1-RESEARCH.md §Pitfall 1.
//
// The legacy X11 backend (xvue/xvuelc.c) retains its own independent
// xtinit_ body — this fix is Qt-side only, preserving BUILD-07
// bit-identity on the legacy path.
void proc(xtinit)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    // Forward to xvinitgraphique_'s body: lazy-alloc window_slot,
    // show/raise/activate, pump the bounded exposure loop. This is
    // idempotent — a subsequent xvinitgraphique_ call from any other
    // Fortran path will reuse the existing XvueWindow via the
    // "if (!win)" guard inside xvinitgraphique_.
    proc(xvinitgraphique)();
}

// ---- 7. xvpxecran_ ----
void proc(xvpxecran)(int *xp, int *yp) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    // D-16: logical pixels from the primary screen, callable before
    // xvinitgraphique_ (D-18). Multi-monitor awareness deferred.
    QScreen* s = QGuiApplication::primaryScreen();
    if (s && xp && yp) {
        *xp = s->size().width();
        *yp = s->size().height();
    }
}

// ---- 8. xvmmecran_ ----
void proc(xvmmecran)(int *xmm, int *ymm) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    // D-17: physical millimetres from QScreen::physicalSize(). Qt 6 guarantees
    // mm regardless of DPR.
    QScreen* s = QGuiApplication::primaryScreen();
    if (s && xmm && ymm) {
        QSizeF mm = s->physicalSize();
        *xmm = static_cast<int>(mm.width()  + 0.5);
        *ymm = static_cast<int>(mm.height() + 0.5);
    }
}

// ---- 9. initaccrochage_ ----
// Phase 5 Plan 05 Task 1 (EVENT-06, D-08). Allocate the 13x13 accrochage
// sprite held on XvueState::mempxaccro_. The sprite is a 3-pixel-thick black
// square border on a transparent background. xvsouris2_ blits this sprite
// onto the canvas at the nearest-item position during the accrochage path
// (Strategy B from 05-RESEARCH.md §6 — save the tile under the sprite with
// accroche_undo_tile_ on first draw, restore it before the next draw; reuses
// the Phase 4 saved_canvas_ ownership pattern instead of trying to emulate
// the X11 GXand/GXorInverted raster-op XOR trick).
//
// Pitfall 11: called before xvinitgraphique_ is tolerated silently — no
// window / no canvas / no state -> return without touching anything.
//
// Idempotency: if mempxaccro_ is already populated, delete and reallocate
// so a re-initialization does not leak the previous sprite. The sprite is
// a small fixed 13x13 pixmap so the cost of reallocation is negligible.
void proc(initaccrochage)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();

    auto& win = XvueApp::window_slot();
    if (!win) return;                 // Pitfall 11: no window yet
    auto* canvas = win->canvas();
    if (!canvas) return;              // Pitfall 11: no canvas yet
    auto* state = win->state();
    if (!state) return;               // defensive — cannot happen today

    // Idempotent: drop any previously-allocated sprite so we do not leak.
    if (state->mempxaccro_) {
        delete state->mempxaccro_;
        state->mempxaccro_ = nullptr;
    }

    // Allocate the 13x13 sprite. lmempxaccro/hmempxaccro in xvuelc.c:142-143.
    state->mempxaccro_ = new QPixmap(13, 13);
    state->mempxaccro_->fill(Qt::transparent);

    // Draw a 3-pixel-thick black square border. QPen width=3 with MiterJoin
    // produces square corners; drawRect(QRect(1,1,11,11)) strokes a
    // 11x11-edge rectangle centered on integer coordinates so the 3-pixel
    // stroke lands inside the 13x13 pixmap without clipping.
    {
        QPainter p(state->mempxaccro_);
        p.setRenderHint(QPainter::Antialiasing, false);  // crisp pixels
        QPen pen(Qt::black, 3);
        pen.setJoinStyle(Qt::MiterJoin);
        p.setPen(pen);
        p.setBrush(Qt::NoBrush);
        p.drawRect(QRect(1, 1, 11, 11));
    }
    // accroche_undo_tile_ is LEFT nullptr — it is allocated lazily by the
    // first xvsouris2_ motion when the filter actually needs to save a
    // tile under the sprite. See Plan 05 Task 2.
}

// ---- 10. xvinfo_ ----
void proc(xvinfo)( int *ix, int *iy, int *maxfonts,
                   int *n1coref, int *ndcoref, int *n1coelf,
                   int *ndcoelf, int *n1coulf, int *ndcoulf, int *nbcolo,
                   char namefonts[][256], int nbchar[], int *nbfonts, int *visuclass ) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();

    // D-03: window resize portion from Phase 1 stays.
    auto& win = XvueApp::window_slot();
    if (win && ix && iy) {
        win->resize(*ix, *iy);
    }

    // Phase 3 D-22/D-23: real font enumeration + palette ranges.
    // Pitfall 2 ordering: palette_init_once() already ran inside XvueState
    // ctor (via XvueWindow ctor via xvinitgraphique_), so any caller that
    // reads red/green/blue immediately after xvinfo_ returns (e.g.
    // xvue/xvinit.f:143 calling XVRECUPRGBDEC) sees populated data.
    //
    // Phase 9.1 maintainer fix 2026-05-06: Mefisto Fortran xvue/xvinit.f at
    // lines 222-289 recognizes ONLY X11 BDF font names of the form "WxH" to
    // populate NUFOHPX[]. If no name matches, NUFOHPX stays all-zero, NOFONT
    // stays 0, CHARGEFONTE(0) early-returns when current==target, NPHACA never
    // gets set, RECTTX line-height = 2*ECARLR ≈ 8px → all multi-line text
    // overlaps. Reported on screenshots after direct ppmail launch.
    //
    // Fix: report X11-BDF-compatible names whose H matches each kFontSizes[]
    // entry's approximate pixel height at 96dpi (Hpx = pt * 96/72 = pt * 4/3,
    // rounded). Qt port still renders the actual TrueType DejaVu Sans Mono;
    // these strings are purely for xvinit.f's pattern-detection. The width W
    // is informational only (RECTTX uses NPLACA from QFontMetrics, not W).
    // Each name advertises ALL BDF height tokens that should resolve to this
    // font index. xvue/xvinit.f INDEX(name, 'WxH') runs sequentially over its
    // pattern list and assigns NUFOHPX[H]=I on every match, so a single font
    // entry can claim multiple heights by listing multiple tokens. Crucial:
    // every height 7..24 the Fortran lahafo.f range can request must be
    // claimed by SOME entry, otherwise NUFOHPX[that_h]=0 + initial NOFONT=0
    // makes chargefonte.f short-circuit and NPHACA never populates -> luou.f
    // DYLGRC collapses to 2*ECARLR=8 -> RECTTX overlaps multi-line text.
    static const char* const kListFonts[XvueState::kNbFonts] = {
        "5x7 5x8 DejaVu Sans Mono 8pt",       // claims H=7,8
        "6x9 6x10 DejaVu Sans Mono 10pt",     // claims H=9,10
        "6x12 6x13 7x13 8x13 DejaVu Sans Mono 12pt", // claims H=12,13
        "7x14 DejaVu Sans Mono 14pt",         // claims H=14
        "9x15 8x16 DejaVu Sans Mono 16pt",    // claims H=15,16
        "10x20 DejaVu Sans Mono 18pt",        // claims H=20 (also covers 17,18,19 via xvinit fallthrough)
        "10x20 DejaVu Sans Mono 20pt",        // claims H=20 (last-write wins for 20)
        "12x24 DejaVu Sans Mono 24pt",        // claims H=24
        "12x24 DejaVu Sans Mono 28pt",
        "12x24 DejaVu Sans Mono 32pt",
    };

    int nfonts = XvueState::kNbFonts;
    if (maxfonts && *maxfonts > 0 && *maxfonts < nfonts) nfonts = *maxfonts;
    // Phase 03.1 Plan 03 — Rule 1 bug fix:
    //   Do NOT write back to *maxfonts. The sole live callers
    //   (xvue/xvinit.f:84, xvue/xvpxfe.f:65) pass the Fortran PARAMETER
    //   constant MAXFONTS (=512 from incl/xvfontes.inc). gfortran -O may
    //   place PARAMETER-sourced actual arguments in read-only memory, so
    //   writing through the pointer crashes with SIGSEGV at the first
    //   xvtest1..4 run. The legacy C backend (xvue/xvuelc.c:612-1042)
    //   never writes *maxfonts — it is a pure input/capacity parameter.
    //   Matching that contract is both correct and safe.
    if (nbfonts)  *nbfonts  = nfonts;
    if (namefonts && nbchar) {
        for (int k = 0; k < nfonts; ++k) {
            std::strncpy(namefonts[k], kListFonts[k], 255);
            namefonts[k][255] = '\0';
            nbchar[k] = static_cast<int>(std::strlen(kListFonts[k]));
        }
    }

    if (visuclass) *visuclass = 4;   // TrueColor — Qt is always 24-bit

    // Palette ranges (legacy xvinfo semantics — ref xvuelc.c:612-1042):
    //   [0..15]   = imposed defaults
    //   [16..255] = user-modifiable
    if (n1coref)  *n1coref  = 0;
    if (ndcoref)  *ndcoref  = 15;
    if (n1coelf)  *n1coelf  = 0;
    if (ndcoelf)  *ndcoelf  = 15;
    if (n1coulf)  *n1coulf  = 16;
    if (ndcoulf)  *ndcoulf  = 255;
    if (nbcolo)   *nbcolo   = 256;
}

// ---- 11. xvrecuprgbdec_ (Phase 3 D-18, TEXT-05) ----
// Read-only snapshot of the process-lifetime palette into Fortran arrays.
// Pitfall 2: palette_init_once() runs inside XvueState ctor, which is
// inside XvueWindow ctor, which is called from xvinitgraphique_ — strictly
// before xvue/xvinit.f:143's XVRECUPRGBDEC call (the only live caller per
// the call-site audit in 03-RESEARCH.md).
void proc(xvrecuprgbdec)(int *nbcolor, float *r, float *g, float *b) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!nbcolor || !r || !g || !b) return;            // WR-03

    int n = *nbcolor;
    if (n < 0) n = 0;
    if (n > XvueState::kMaxPalette) n = XvueState::kMaxPalette;

    for (int i = 0; i < n; ++i) {
        r[i] = XvueState::red[i];
        g[i] = XvueState::green[i];
        b[i] = XvueState::blue[i];
    }
}

// ---- 12. xvactivervb_ (Phase 3 D-17 AMENDED per A3, TEXT-05) ----
// Bulk palette load. Research correction over the original D-17 text:
// the single live caller xvue/palcde.f:619 passes NDCOUL+1 cells, so this
// is a full user-palette refresh, NOT a transient one-shot color write.
// Matches legacy xvuelc.c:1072-1116 bulk-copy semantics.
void proc(xvactivervb)( int *palcour, int *nbcells,
                        float r[], float g[], float b[] ) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    (void)palcour;  // palcourc was a PostScript state var, dropped (02/D-26)
    if (!nbcells || !r || !g || !b) return;            // WR-03

    int n = *nbcells;
    if (n < 0) n = 0;
    if (n > XvueState::kMaxPalette) n = XvueState::kMaxPalette;

    for (int i = 0; i < n; ++i) {
        XvueState::red[i]   = r[i];
        XvueState::green[i] = g[i];
        XvueState::blue[i]  = b[i];
        XvueState::palette_cache_dirty_[i] = true;
    }
    // No painter touch, no flush — xvcouleur_ picks up new values via the
    // dirty-flag rebuild on its next call.
}

// IN-02: shared palette-index clamp + dirty-cache rebuild used by
// xvcouleur_, xvfond_, and apply_palette_foreground. Returns the resolved
// QColor; callers decide whether to write foreground_/background_ and
// whether to call applyPen(). File-local (static) by design — all call
// sites are in this translation unit.
static const QColor& palette_resolve(int icolor, int fallback) {
    int i = icolor;
    if (i < 0 || i >= XvueState::kMaxPalette) i = fallback;
    if (XvueState::palette_cache_dirty_[i]) {
        XvueState::palette_cache_[i] = QColor::fromRgbF(
            XvueState::red[i],
            XvueState::green[i],
            XvueState::blue[i]);
        XvueState::palette_cache_dirty_[i] = false;
    }
    return XvueState::palette_cache_[i];
}

// ---- 13. xvcouleur_ (Phase 3 D-14, TEXT-04) ----
// State-change entry: install palette_cache_[i] as current foreground.
// NO flush — drawing primitives run their own D-01 epilogue.
void proc(xvcouleur)(int *icolor) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!icolor) return;                                // WR-03

    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st) return;

    // IN-02: legacy "fallback to red" (index 1) via shared helper.
    st->foreground_ = palette_resolve(*icolor, 1);
    st->applyPen();  // 02/D-20 syncs brush + pen from foreground_
    // D-14: no update(), no processEvents — state-change only.
}

// ---- 14. xvpostscript_ — Phase 7 Plan 02 (EXPORT-04, D-05, D-06): dispatch
//      one-liner; PsEmitter owns all real logic. ABI body changed only — the
//      proc(xvpostscript) extern "C" symbol is unchanged so the Fortran ABI
//      stays at 58 entry points. WR-03: null-guard *lasops dereference.
void proc(xvpostscript)(int *lasops) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!lasops) return;                            // WR-03 null-arg guard
    XvueApp::psEmitter().handleLasops(*lasops);
}

// ---- 15. fenetremempx_ ----
// D-04 (Phase 4): no-op — backing_ is the visible surface (Phase 2 D-04
// single-backing collapse). No fenetre/mempx distinction on Qt.
void proc(fenetremempx)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    // intentional no-op per Phase 4 D-04
}

// ---- 16. mempxfenetre_ ----
// D-04 (Phase 4): no-op — backing_ is the visible surface (Phase 2 D-04
// single-backing collapse). Phase 03-04 empirically validated that
// solver tracers (e.g. elas/trelas.f) produce correct displays under this.
void proc(mempxfenetre)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    // intentional no-op per Phase 4 D-04
}

// ---- 17. sauvefenetre_ ----
// D-07, D-08 (Phase 4): save backing_ -> saved_canvas_ via file-local helper.
void proc(sauvefenetre)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    xvue_qt_save_to_slot();
}

// ---- 18. restaurefenetre_ ----
// D-07, D-08 (Phase 4): restore saved_canvas_ -> backing_ via file-local helper.
void proc(restaurefenetre)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    xvue_qt_restore_from_slot();
}

// ---- 19. sauvemempx_ ----
// D-07, D-08 (Phase 4): bit-identical to sauvefenetre_ (single-slot model, D-01).
void proc(sauvemempx)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    xvue_qt_save_to_slot();
}

// ---- 20. restauremempx_ ----
// D-07, D-08 (Phase 4): bit-identical to restaurefenetre_ (single-slot model, D-01).
void proc(restauremempx)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    xvue_qt_restore_from_slot();
}

// ---- 21. effacemempx_ ----
// D-10, D-11 (Phase 4): same body as effacer_ — Qt has no separate
// mempx surface (Phase 2 D-04, D-15). Two symbols kept distinct for
// ABI preservation (D-08) but functionally identical.
void proc(effacemempx)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (st && st->painter_ && st->painter_->isActive() && st->backing_) {
        st->painter_->fillRect(st->backing_->rect(), st->background_);
    }
    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) win->canvas()->update();
}

// ---- 22. effacer_ (D-15) ----
void proc(effacer)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (st && st->painter_ && st->painter_->isActive() && st->backing_) {
        st->painter_->fillRect(st->backing_->rect(), st->background_);
    }
    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) {
        win->canvas()->update();
    }
    // WR-02: deferred flush -- xvvoir_/xvpause_ pump the event loop.
}

// ---- 23. xvfond_ (Phase 3: retired Phase 2 two-entry hack) ----
// Phase 2 used a hardcoded {black, white} 2-entry palette (D-14). Phase 3
// routes the color source through the real palette_cache_ the same way
// xvcouleur_ does, while preserving the Phase 2 D-24 fillRect mechanism.
// Phase 7 Plan 03: xvuelc.c xvfond_ (line 1439) emits NO PostScript bytes —
// it only mutates X11 window background. PsEmitter::fond is therefore a
// no-op; we still call it here to keep the wiring shape consistent so a
// future plan can add a PS-side equivalent without churning callers.
void proc(xvfond)(int *icolor) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!icolor) return;                               // WR-03

    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st) return;

    // IN-02: xvfond_ defaults to palette index 0 (black) via shared helper.
    // D-15: update XvueState::background_ through the live window, repaint.
    // Phase 2 D-24: re-fill the backing with the new background and flush.
    st->background_ = palette_resolve(*icolor, 0);
    if (st->painter_ && st->painter_->isActive() && st->backing_) {
        st->painter_->fillRect(st->backing_->rect(), st->background_);
    }
    // Phase 7 Plan 03: PsEmitter::fond is a no-op (legacy emits no PS).
    XvueApp::psEmitter().fond(
        st->background_.redF(), st->background_.greenF(),
        st->background_.blueF());
    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) {
        win->canvas()->update();
    }
    // WR-02: deferred flush -- xvvoir_/xvpause_ pump the event loop.
}

// ---- 24. xvchargefonte_ (Phase 3 D-04, TEXT-01) ----
// Bundled DejaVu Sans Mono at XvueState::kFontSizes[*nofont]. QFont is RAII
// (Pitfall 6) so *nofont0 (the legacy "previous font to free") is ignored.
// Phase 7 Plan 03 (EXPORT-04): emits PS "%s %d %d %s charge\n" via D-08
// hardcoded mapping table. DejaVu Sans Mono lands on /Helvetica fallback.
void proc(xvchargefonte)(int *nofont0, int *nofont, int *largpx, int *hautpx) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    (void)nofont0;  // Pitfall 6: QFont is RAII, no explicit free needed

    if (!nofont || !largpx || !hautpx) return;       // WR-03

    int idx = *nofont;
    if (idx < 0) idx = 0;
    if (idx >= XvueState::kNbFonts) idx = XvueState::kNbFonts - 1;

    auto& win = XvueApp::window_slot();
    if (!win) { *largpx = 0; *hautpx = 0; return; }
    auto* st = win->state();
    if (!st) { *largpx = 0; *hautpx = 0; return; }

    st->current_font_size_idx_ = idx;
    // Use setPixelSize so the BDF "WxH" prefix Fortran sees in xvinfo namefonts
    // matches the actual Qt-rendered pixel height. kFontPixelSizes[] mirrors
    // the BDF H values (8, 10, 13, 14, 15, 20, 20, 24, 24, 24) — Fortran's
    // xvue/xvinit.f populates NUFOHPX[H] from those, so a request for height H
    // pixels gets a font that actually renders at H pixels (avoids
    // Fortran's per-line Y-advance using a smaller value than the actual
    // glyph extent — which manifested as overlapping text on welcome /
    // impqua histogram / project info panels per maintainer screenshots
    // 2026-05-06).
    static constexpr int kFontPixelSizes[XvueState::kNbFonts] = {
        8, 10, 13, 14, 15, 20, 20, 24, 24, 24
    };
    QFont f(QStringLiteral("DejaVu Sans Mono"));
    f.setPixelSize(kFontPixelSizes[idx]);
    st->current_font_ = f;
    st->current_metrics_ = QFontMetrics(st->current_font_);

    if (st->painter_ && st->painter_->isActive()) {
        st->painter_->setFont(st->current_font_);
    }

    *largpx = st->current_metrics_.horizontalAdvance(QLatin1Char('0'));
    *hautpx = st->current_metrics_.height();
    // Floor: ensure at least 8px height so a tiny system-default rendering
    // doesn't give Fortran 0 (which would collapse DYLGRC line spacing).
    if (*hautpx < 8) *hautpx = 8;
    if (*largpx < 4) *largpx = 4;
    std::fprintf(stderr,
                 "xvchargefonte: idx=%d pixelSize=%d -> largpx=%d hautpx=%d\n",
                 idx, kFontPixelSizes[idx], *largpx, *hautpx);

    // Phase 7 Plan 03: PS font-load emit. Pass family + size (pt) + ascent/
    // descent / bold / italic; the D-08 mapping table inside chargefonte()
    // produces "/Courier" / "/Times-Roman" / "/Helvetica" / "/NewCenturySchlbk"
    // (with optional -Bold/-Italic/-Oblique suffix), and the format string
    // matches xvuelc.c:1553 byte-for-byte.
    XvueApp::psEmitter().chargefonte(
        st->current_font_.family(),
        XvueState::kFontSizes[idx],
        st->current_metrics_.ascent(),
        st->current_metrics_.descent(),
        st->current_font_.bold(),
        st->current_font_.italic());
}

// ---- 25. xvnbpixeltexte_ (Phase 3 D-05, TEXT-02) ----
void proc(xvnbpixeltexte)(char *texte, int *length, int *nbpxla, int *nbpxha) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!texte || !length || !nbpxla || !nbpxha) return;   // WR-03

    auto& win = XvueApp::window_slot();
    if (!win) { *nbpxla = 0; *nbpxha = 0; return; }
    auto* st = win->state();
    if (!st) { *nbpxla = 0; *nbpxha = 0; return; }

    // Pitfall 4: TWO-ARG fromLatin1 — explicit length from Fortran.
    QString qstr = QString::fromLatin1(texte, *length);
    *nbpxla = st->current_metrics_.horizontalAdvance(qstr);
    *nbpxha = st->current_metrics_.height();
}

// ---- 26. xvfermer_ ----
void proc(xvfermer)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    // Headless-test capture hook — matches the legacy xvuelc.c contract
    // (cf. xvue/xvuelc.c xvfermer_). Three optional env vars:
    //   MEFISTO_XVFERMER_READY_FILE  — path to touch right before
    //       destroying the window so an external screenshot tool can
    //       grab the final rendered state (root-window capture path).
    //   MEFISTO_XVFERMER_HOLD_MS     — milliseconds to sleep after
    //       flushing + touching the sentinel, before the window
    //       is destroyed. Default 0 (no hold).
    //   MEFISTO_QT_CAPTURE_PATH      — if set, call QWidget::grab()
    //       on the canvas and save directly to this PNG path before
    //       destroying the window. Works under QT_QPA_PLATFORM=offscreen
    //       with NO X11 dependency, so it does not need xcb-cursor0 or
    //       an Xvfb server — ideal for headless CI and for cases where
    //       the xcb platform plugin fails to initialise.
    // All three default to "no-op" when unset; interactive behavior is
    // unchanged. Qt-side counterpart of the legacy hook (commit e029b84)
    // with the extra in-process grab() path for full headless support.
    {
        auto& win = XvueApp::window_slot();
        if (win && win->canvas()) {
            win->canvas()->update();
        }
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);

        const char* qt_capture_path = std::getenv("MEFISTO_QT_CAPTURE_PATH");
        if (qt_capture_path && qt_capture_path[0] != '\0') {
            // In-process screenshot of the XvueState backing pixmap —
            // this is the authoritative surface every Phase 2 drawing
            // primitive paints on (see xvue_qt_canvas.cpp:35-36:
            // "paintEvent = drawPixmap blit ONLY" from state_->backing_).
            // Capturing it directly bypasses the widget paint pipeline,
            // which under QT_QPA_PLATFORM=offscreen does not always
            // trigger a full repaint between update() and grab(), and
            // guarantees the captured PNG matches exactly what the
            // drawing primitives produced up to this point. Also ends
            // any active painter session first so the pixmap is in a
            // consistent readable state.
            XvueState* st = (win ? win->state() : nullptr);
            if (st && st->backing_) {
                if (st->painter_ && st->painter_->isActive()) {
                    st->painter_->end();
                }
                if (!st->backing_->save(
                        QString::fromLatin1(qt_capture_path), "PNG")) {
                    std::fprintf(stderr,
                        "xvue-qt: failed to save backing to "
                        "MEFISTO_QT_CAPTURE_PATH=%s\n",
                        qt_capture_path);
                }
            } else if (win && win->canvas()) {
                // Fallback to widget grab if no backing is available
                // (e.g., the driver never called xvinfo_ to allocate it).
                QPixmap shot = win->canvas()->grab();
                if (!shot.save(
                        QString::fromLatin1(qt_capture_path), "PNG")) {
                    std::fprintf(stderr,
                        "xvue-qt: failed to save widget grab to "
                        "MEFISTO_QT_CAPTURE_PATH=%s\n",
                        qt_capture_path);
                }
            }
        }

        const char* ready_path = std::getenv("MEFISTO_XVFERMER_READY_FILE");
        const char* hold_env   = std::getenv("MEFISTO_XVFERMER_HOLD_MS");
        int hold_ms = 0;
        if (hold_env && hold_env[0] != '\0') {
            int hv = std::atoi(hold_env);
            if (hv >= 0 && hv <= 60000) hold_ms = hv;
        }
        if (ready_path && ready_path[0] != '\0') {
            FILE* f = std::fopen(ready_path, "w");
            if (f) { std::fputs("ready\n", f); std::fclose(f); }
        }
        if (hold_ms > 0) {
            // Pump the event loop during the hold so paint events
            // delivered to the canvas actually reach the X server.
            QElapsedTimer t; t.start();
            while (t.elapsed() < hold_ms) {
                QCoreApplication::processEvents(
                    QEventLoop::ExcludeUserInputEvents, 50);
                usleep(10 * 1000);
            }
        }
    }
    // Phase 9 Plan 9-06 (carry-forward #1, Pitfall 5 mitigation 2): clear
    // any MEFISTO_QT_WINDOW_SIZE-installed canvas size constraints before
    // window teardown so a subsequent xvinitgraphique_ reopen on the same
    // process does not inherit stale min/max from the previous open cycle.
    // Safe even when no env was set (constraints default to 0,0 / MAX,MAX).
    {
        auto& win_close = XvueApp::window_slot();
        if (win_close && win_close->canvas()) {
            win_close->canvas()->setMinimumSize(0, 0);
            win_close->canvas()->setMaximumSize(QWIDGETSIZE_MAX, QWIDGETSIZE_MAX);
        }
    }
    // D-06: destroy the window only. Do NOT touch qApp, do NOT call
    // qApp->quit(). The QApplication lives until the atexit handler.
    XvueApp::window_slot().reset();
    // D-06 addendum 2026-04-11 (debug phase-01-xvtest0-teardown-segfault):
    // drain DeferredDelete events queued by the widget-hierarchy teardown
    // immediately, while Qt is in a well-defined state. Otherwise these
    // events survive until the atexit handler and risk running against a
    // torn-down Qt. ExcludeUserInputEvents per Pitfall 4.
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

// ---- 27. xvpxfenetre_ (D-23) ----
void proc(xvpxfenetre)(int *x, int *y) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y) return;
    auto& win = XvueApp::window_slot();
    if (!win || !win->canvas()) { *x = 0; *y = 0; return; }
    *x = win->canvas()->width();   // logical pixels (SHELL-06)
    *y = win->canvas()->height();
}

// ---- 28. xvftexte_ (Phase 3 D-06, TEXT-03) ----
// Shares body with xvtexte_ under Phase 2 single-backing (02/D-05, Pitfall 7).
void proc(xvftexte)(char string[], int *length, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!string || !length || !x1 || !y1) return;      // WR-03
    xvue_qt_draw_text_common(string, *length, *x1, *y1);
}

// ---- 29. xvtexte_ (Phase 3 D-06, TEXT-03) ----
void proc(xvtexte)(char string[], int *length, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!string || !length || !x1 || !y1) return;      // WR-03
    xvue_qt_draw_text_common(string, *length, *x1, *y1);
}

// ---- 30. xvface_ (D-11, DRAW-03) ----
// Legacy: XFillPolygon -- fill only, no edge stroke (WR-01).
// Phase 7 Plan 03 (EXPORT-04): emits the PS "F" close-and-fill opcode after
// the QPainter call.
void proc(xvface)(int *n, MefistoPoint *pts) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    if (!n || *n < 3 || !pts) return;

    QPolygon poly;
    poly.reserve(*n);
    for (int i = 0; i < *n; ++i) {
        poly << QPoint(pts[i].x, pts[i].y);
    }
    // WR-01: drawPolygon strokes the outline with the current pen. Legacy
    // XFillPolygon only fills. Suppress the pen for the fill step so the
    // two colors stay independent once Phase 3 unlocks pen/brush divergence.
    QPen saved_pen = st->painter_->pen();
    st->painter_->setPen(Qt::NoPen);
    st->painter_->drawPolygon(poly, Qt::OddEvenFill);  // auto 1<->n close
    st->painter_->setPen(saved_pen);

    // Phase 7 Plan 03: PS "F" emit. MefistoPoint overload — helper Y-flips
    // each vertex internally per README_COORDS.md.
    XvueApp::psEmitter().face(pts, *n);

    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) win->canvas()->update();
}

// ---- 31. xvtypetrait_ (D-17, D-19, DRAW-06) ----
// Phase 7 Plan 03: emits PS "%2i typet\n" after the state update.
void proc(xvtypetrait)(int *ptype) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!ptype) return;   // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st) return;
    st->pen_style_ = *ptype;
    st->applyPen();   // applyPen() self-gates on painter_->isActive()
    XvueApp::psEmitter().typetrait(*ptype);
}

// ---- 32. xvepaisseur_ (D-18, DRAW-06) ----
// Phase 7 Plan 03: emits PS "%2i epais\n" after the state update.
void proc(xvepaisseur)(int *pepais) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!pepais) return;   // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st) return;
    st->pen_width_base_ = *pepais;
    st->applyPen();
    XvueApp::psEmitter().epaisseur(*pepais);
}

// ---- 34. xvtrait_ (D-09, DRAW-02) ----
// Phase 7 Plan 03: emits the PS "S" opcode after the QPainter::drawLine.
// Helper internally checks active() so this is a no-op when no PS capture
// is in progress. xvftrait_ does NOT emit (xvuelc.c:1905 is the
// fenetre_mef-only variant); see xvftrait body below.
void proc(xvtrait)(int *x1, int *y1, int *x2, int *y2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x1 || !y1 || !x2 || !y2) return;   // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    st->painter_->drawLine(*x1, *y1, *x2, *y2);
    // Phase 7 Plan 03 (EXPORT-04): PS "S" opcode emit. Pass canvas-Y;
    // the helper applies pyFlip internally per README_COORDS.md.
    XvueApp::psEmitter().line(*x1, *y1, *x2, *y2);
    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) win->canvas()->update();
    // WR-02: deferred flush -- xvvoir_/xvpause_ pump the event loop.
}

// ---- 33. xvftrait_ (D-09 -- semantically identical to xvtrait_ under
// the single-backing model; legacy window-vs-mempx split obsolete) ----
// Phase 7 Plan 03: xvuelc.c:1905 draws to fenetre_mef and emits NO
// PostScript. Inline the painter call instead of routing through xvtrait_
// so the PS-emit byte parity matches the legacy stream.
void proc(xvftrait)(int *x1, int *y1, int *x2, int *y2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x1 || !y1 || !x2 || !y2) return;   // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    st->painter_->drawLine(*x1, *y1, *x2, *y2);
    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) win->canvas()->update();
    // Note: NO psEmitter emit — legacy xvftrait_ never wrote PS bytes.
}

// ---- 35. xvtraits_ (D-10, DRAW-02) ----
// Phase 7 Plan 03: emits PS "P" polyline opcode after the QPainter call.
void proc(xvtraits)(int *nbpoints, MefistoPoint *points) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    if (!nbpoints || *nbpoints < 2 || !points) return;

    constexpr int STACK_LIMIT = 128;
    if (*nbpoints <= STACK_LIMIT) {
        QPoint qpts[STACK_LIMIT];
        for (int i = 0; i < *nbpoints; ++i) {
            qpts[i] = QPoint(points[i].x, points[i].y);
        }
        st->painter_->drawPolyline(qpts, *nbpoints);
    } else {
        std::vector<QPoint> qpts;
        qpts.reserve(*nbpoints);
        for (int i = 0; i < *nbpoints; ++i) {
            qpts.emplace_back(points[i].x, points[i].y);
        }
        st->painter_->drawPolyline(qpts.data(),
                                   static_cast<int>(qpts.size()));
    }
    // Phase 7 Plan 03: PS "P" emit (xvuelc.c:2037-2093). Helper handles
    // the (n-1)-segments count + Y-flip + chunked emit internally.
    XvueApp::psEmitter().traits(points, *nbpoints);
    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) win->canvas()->update();
    // WR-02: deferred flush -- xvvoir_/xvpause_ pump the event loop.
}

// Internal helper — inline palette install for callers that already hold
// a resolved XvueState*. Mirrors proc(xvcouleur) body but skips
// XvueApp::ensure()/thread-assert/window-slot re-lookup. IN-02: delegates
// clamp + dirty-cache rebuild to the shared palette_resolve() helper.
static void apply_palette_foreground(XvueState* st, int icolor) {
    st->foreground_ = palette_resolve(icolor, 1);  // legacy fallback to red
    st->applyPen();  // rebuilds pen from pen_style_ + foreground_, brush from foreground_
}

// ---- 36. xvfacetraits_ (D-12, DRAW-03; 02.1/D-1: honor ncf/nca) ----
// Phase 7 Plan 03 (EXPORT-04): emits PS "FP" fill+border opcode after the
// QPainter calls. Per xvuelc.c:2095-2181 the outline RGB is courgb_ (post-
// nca) and the fill RGB is courgbm (saved before the second xvcouleur
// call). We read the ncf-resolved QColor BEFORE the apply_palette_foreground
// pair so we have both colors at PS-emit time.
void proc(xvfacetraits)(int *ncf, int *nca, int *n, MefistoPoint *pts) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!ncf || !nca || !n || *n < 3 || !pts) return;   // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;

    QPolygon poly;
    poly.reserve(*n);
    for (int i = 0; i < *n; ++i) {
        poly << QPoint(pts[i].x, pts[i].y);
    }

    // Phase 7 Plan 03: capture the ncf-resolved fill color BEFORE we mutate
    // st->foreground_ via apply_palette_foreground, so the PS "FP" emit can
    // include both fill (ncf) and border (post-nca) RGB triplets.
    const QColor fill_qc = palette_resolve(*ncf, 1);

    // 02.1/D-1: mirror xvue/xvuelc.c:2055,2064 — install ncf, fill; install
    // nca, outline. Dashed style set by xvtypetrait_ survives the outline
    // step because pen_style_ is threaded through XvueState::applyPen()
    // (see xvue_qt_state.cpp:136-142), which rebuilds pen_ from
    // pen_style_ + foreground_ on every apply_palette_foreground() call.
    // WR-01: fill must not stroke (setPen(Qt::NoPen)); outline must not
    // re-fill (setBrush(Qt::NoBrush)).
    apply_palette_foreground(st, *ncf);
    st->painter_->setPen(Qt::NoPen);
    st->painter_->drawPolygon(poly, Qt::OddEvenFill);   // fill (D-12 order)
    // IN-01: no pen restore between fill and outline — the next
    // apply_palette_foreground(st, *nca) call invokes applyPen() which
    // unconditionally rebuilds pen_ from pen_style_ + foreground_, so any
    // intermediate setPen(Qt::NoPen) is immediately superseded.

    apply_palette_foreground(st, *nca);                 // rebuilds pen with nca color, preserves dash style
    st->painter_->setBrush(Qt::NoBrush);
    st->painter_->drawPolygon(poly);                    // outline (D-12 order) — uses current pen (incl dash)
    // IN-03 / WR-04: legacy leaves foreground_ at nca (we match that by
    // not touching foreground_). The outline step above left the painter
    // brush at Qt::NoBrush; while any subsequent xvcouleur_/applyPen()
    // would refresh it, WR-04 tightens the post-condition by restoring
    // the brush here so primitives that read the live painter brush
    // directly (e.g. xvface_) observe the canonical post-ncf->nca final
    // state without depending on an intervening state-change call.
    st->painter_->setBrush(QBrush(st->foreground_, Qt::SolidPattern));

    // Phase 7 Plan 03: PS "FP" emit. Pass fill RGB explicitly; the helper
    // pulls border RGB from PsEmitter::courgb_ (which would mirror the
    // current painter color set by xvcouleur_ — kept in sync once Plan 06
    // golden-compares end-to-end). For Plan 03 we pass fillA = -1.0f to
    // force the "1.00 FP" branch matching xvuelc.c:2166 default.
    XvueApp::psEmitter().facetraits(pts, *n,
        fill_qc.redF(), fill_qc.greenF(), fill_qc.blueF(),
        /*fillA=*/-1.0f);

    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) win->canvas()->update();
}

// ---- 37. xvsouris_ ----
// Headless-test short-circuit: when MEFISTO_XVSOURIS_AUTOEXIT is set,
// return a synthetic SPACE keypress so driver loops of the form
//   IF (NOTYEV .NE. 2) GOTO 100
// exit on first iteration. Pairs with the MEFISTO_XVFERMER_* capture
// hook above and matches the legacy xvuelc.c contract.
//
// When the env var is not set, the Qt-side behavior stays the stub it
// was before (returns no event — the driver's input loop is expected
// to be handled by the Qt event pump). Interactive paths are never
// regressed.
void proc(xvsouris)(int *notypeevent, int *nbc, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();

    // Phase 6.0 Plan 03 (UX-03): menu-queue pre-drain. If a QAction handler
    // queued a synthetic lexicon character via XvueMenuBridge::queueLexicon,
    // return it as a notypeevent=2 (KeyPress) BEFORE the AUTOEXIT short-
    // circuit OR the bridge->waitForEvent call. Pre-drain MUST run before
    // AUTOEXIT so menu-driven sequences win over the headless escape hatch
    // when both are armed (Research §1 "Interaction with AUTOEXIT" — a test
    // that queues chars wants to see those chars first; AUTOEXIT is a
    // headless safety valve, not a behavior gate).
    //
    // Mirrors the matching pre-drain in XvueEventBridge::waitForEvent (which
    // covers direct waitForEvent() callers that bypass this wrapper). Having
    // it in BOTH places is intentional: this site beats AUTOEXIT in the
    // production xvsouris_ path; the waitForEvent site beats nothing but
    // covers the test-only direct-bridge path. Both check
    // win->menuBridge() — Plan 06 wires the production bridge; Plan 03's
    // unit tests inject one via XvueWindow::setMenuBridgeForTesting.
    {
        auto& win = XvueApp::window_slot();
        if (win) {
            if (auto* mb = win->menuBridge()) {
                if (auto c = mb->popChar()) {
                    if (notypeevent) *notypeevent = 2;
                    if (nbc)         *nbc         = static_cast<unsigned char>(*c);
                    if (x1)          *x1          = 0;
                    if (y1)          *y1          = 0;
                    return;
                }
            }
        }
    }

    const char* autoexit = std::getenv("MEFISTO_XVSOURIS_AUTOEXIT");
    if (autoexit && autoexit[0] != '\0') {
        int delay_ms = 500;
        const char* delay_env = std::getenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");
        if (delay_env && delay_env[0] != '\0') {
            int d = std::atoi(delay_env);
            if (d >= 0 && d <= 60000) delay_ms = d;
        }
        // Paint the canvas and pump events so the window is up to date.
        auto& win = XvueApp::window_slot();
        if (win && win->canvas()) win->canvas()->update();
        QElapsedTimer t; t.start();
        while (t.elapsed() < delay_ms) {
            QCoreApplication::processEvents(
                QEventLoop::ExcludeUserInputEvents, 25);
            usleep(5 * 1000);
        }
        if (notypeevent) *notypeevent = 2;
        if (nbc)         *nbc         = ' ';
        if (x1)          *x1          = 0;
        if (y1)          *y1          = 0;
        return;
    }
    // Phase 5 Plan 04 (EVENT-02): real body. Route through the Plan 02 bridge.
    // No window open → silent abandon (Pitfall 11 analogue). Keeps
    // prpr/xvtest1-style drivers safe when called standalone.
    auto& win = XvueApp::window_slot();
    if (!win || !win->bridge()) {
        if (notypeevent) *notypeevent = 0;
        if (nbc)         *nbc         = 0;
        if (x1)          *x1          = 0;
        if (y1)          *y1          = 0;
        return;
    }
    auto r = win->bridge()->waitForEvent(XvueEventBridge::WaitMode::Souris);
    if (notypeevent) *notypeevent = r.notypeevent;
    if (nbc)         *nbc         = r.nbc;
    if (x1)          *x1          = r.x;
    if (y1)          *y1          = r.y;
}

// ---- 38. xvsouris2_ ----
// Phase 5 Plan 05 Task 2 (EVENT-03). Real body for the accrochage variant of
// xvsouris_. Routes through the Plan 02 bridge with WaitMode::Souris2 so the
// eventFilter runs the Strategy B save/restore dance: on each motion (and
// on button press) the filter finds the nearest items[] point, erases the
// previously-drawn sprite if any, blits the mempxaccro_ sprite at the new
// location, and returns notypeevent=5 (or 1 on release, or 0 on Esc/@
// abort). See xvue_qt_event.cpp for the full logic.
//
// AUTOEXIT: same MEFISTO_XVSOURIS_AUTOEXIT byte-for-byte contract as
// xvsouris_ / xvpause_ (D-10) so headless capture harnesses never hang.
// No-window guard mirrors Pitfall 11. Null-pointer guards on the Fortran
// out-params defend against pathological callers.
void proc(xvsouris2)(int *items, int *pmin0, int *notypeevent, int *ibutton, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    const char* autoexit = std::getenv("MEFISTO_XVSOURIS_AUTOEXIT");
    if (autoexit && autoexit[0] != '\0') {
        int delay_ms = 500;
        const char* delay_env = std::getenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");
        if (delay_env && delay_env[0] != '\0') {
            int d = std::atoi(delay_env);
            if (d >= 0 && d <= 60000) delay_ms = d;
        }
        auto& win = XvueApp::window_slot();
        if (win && win->canvas()) win->canvas()->update();
        QElapsedTimer t; t.start();
        while (t.elapsed() < delay_ms) {
            QCoreApplication::processEvents(
                QEventLoop::ExcludeUserInputEvents, 25);
            usleep(5 * 1000);
        }
        if (notypeevent) *notypeevent = 2;
        if (ibutton)     *ibutton     = ' ';
        if (x1)          *x1          = 0;
        if (y1)          *y1          = 0;
        (void)items; (void)pmin0;
        return;
    }

    // Pitfall 11 analogue: no window → silent abandon with zero out-params.
    auto& win = XvueApp::window_slot();
    if (!win || !win->bridge()) {
        if (notypeevent) *notypeevent = 0;
        if (ibutton)     *ibutton     = 0;
        if (x1)          *x1          = 0;
        if (y1)          *y1          = 0;
        return;
    }

    auto r = win->bridge()->waitForEvent(
        XvueEventBridge::WaitMode::Souris2, items, pmin0);
    if (notypeevent) *notypeevent = r.notypeevent;
    if (ibutton)     *ibutton     = r.nbc;
    if (x1)          *x1          = r.x;
    if (y1)          *y1          = r.y;
}

// ---- 39. deplsouris_ ----
// Phase 5 Plan 04 (EVENT-05, D-09). Warp the cursor to canvas-local (x, y)
// via QCursor::setPos(canvas->mapToGlobal(QPoint(x,y))). Non-blocking — does
// not touch the event loop.
//
// Wayland caveat (D-09, Pitfall 5): QCursor::setPos is a no-op on most pure
// Wayland compositors; X11 / XWayland is the supported session type.
//
// T-05-04-01 mitigation: bounds-check |x| / |y| to reject pathological
// Fortran callers. Values outside ±32768 silently no-op so a bad caller
// cannot warp the cursor to arbitrary screen coordinates.
void proc(deplsouris)(int *x, int *y) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y) return;
    const int xi = *x;
    const int yi = *y;
    if (xi < -32768 || xi > 32768 || yi < -32768 || yi > 32768) return;
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* canvas = win->canvas();
    if (!canvas) return;
    QCursor::setPos(canvas->mapToGlobal(QPoint(xi, yi)));
}

// ---- 40. xvvoir_ (D-02) ----
void proc(xvvoir)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (win && win->canvas()) {
        win->canvas()->update();
    }
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

// ---- 41. xvpause_ ----
// Phase 5 Plan 04 (EVENT-04). Block until any event arrives on the canvas,
// then return. Mirrors the xvuelc.c:2516-2531 XNextEvent/KeyPress semantics
// but through the Plan 02 XvueEventBridge so the Fortran call remains
// imperative while Qt runs a nested event loop underneath.
//
// AUTOEXIT extension: MEFISTO_XVSOURIS_AUTOEXIT short-circuits xvpause_ too
// so bin/xvtest-capture.sh and other headless harnesses never hang on a
// CALL XVPAUSE. Same env var as xvsouris_/xvsouris2_ (§8 of 05-RESEARCH.md,
// D-10 preserves the var name; Plan 04 extends its scope to xvpause_).
void proc(xvpause)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    const char* autoexit = std::getenv("MEFISTO_XVSOURIS_AUTOEXIT");
    if (autoexit && autoexit[0] != '\0') {
        int delay_ms = 500;
        const char* delay_env = std::getenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");
        if (delay_env && delay_env[0] != '\0') {
            int d = std::atoi(delay_env);
            if (d >= 0 && d <= 60000) delay_ms = d;
        }
        auto& win = XvueApp::window_slot();
        if (win && win->canvas()) win->canvas()->update();
        QElapsedTimer t; t.start();
        while (t.elapsed() < delay_ms) {
            QCoreApplication::processEvents(
                QEventLoop::ExcludeUserInputEvents, 25);
            usleep(5 * 1000);
        }
        return;
    }
    auto& win = XvueApp::window_slot();
    if (!win || !win->bridge()) return;  // no window — silent no-op (xvpause has no out-params)
    (void)win->bridge()->waitForEvent(XvueEventBridge::WaitMode::Pause);
}

// ---- 42. xvfbordrectangle_ (D-13, DRAW-04) ----
// Legacy: XDrawRectangle -- outline only on fenetre_mef.
// Phase 7 Plan 03 (EXPORT-04): xvuelc.c:2556 draws to fenetre_mef and emits
// NO PostScript. Under Phase 2 single-backing collapse this Qt entry shares
// xvue_qt_draw_rect_common with xvbordrectangle, so emit_ps=false here
// preserves legacy PS byte parity.
void proc(xvfbordrectangle)(int *x, int *y, int *width, int *height) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y || !width || !height) return;   // WR-03
    xvue_qt_draw_rect_common(*x, *y, *width, *height, RectMode::Outline,
                             /*emit_ps=*/false);
}

// ---- 43. xvbordrectangle_ (D-13, DRAW-04) ----
// Legacy: XDrawRectangle -- outline only.
void proc(xvbordrectangle)(int *x, int *y, int *width, int *height) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y || !width || !height) return;   // WR-03
    xvue_qt_draw_rect_common(*x, *y, *width, *height, RectMode::Outline);
}

// ---- 44. xvfrectangle_ (D-13, DRAW-04) ----
// Legacy: XFillRectangle -- fill only on fenetre_mef.
// Phase 7 Plan 03: xvuelc.c:2619 draws to fenetre_mef, emits NO PostScript.
// emit_ps=false preserves legacy byte parity (single-backing collapse note).
void proc(xvfrectangle)(int *x, int *y, int *width, int *height) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y || !width || !height) return;   // WR-03
    xvue_qt_draw_rect_common(*x, *y, *width, *height, RectMode::Fill,
                             /*emit_ps=*/false);
}

// ---- 45. xvrectangle_ (D-13, DRAW-04) ----
// Legacy: XFillRectangle -- fill only.
void proc(xvrectangle)(int *x, int *y, int *width, int *height) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y || !width || !height) return;   // WR-03
    xvue_qt_draw_rect_common(*x, *y, *width, *height, RectMode::Fill);
}

// ---- 46. xvbordarcellipse_ (D-14, RESEARCH Q1: drawArc, DRAW-05) ----
// Legacy xvuelc.c:2554 -- uses XDrawArc, outline only, float* angles in degrees.
// X11 x64 (1/64 deg) -> Qt x16 (1/16 deg). NOT x64.
// Phase 7 Plan 03: emits PS "el" arc-outline opcode (xvuelc.c:2729).
void proc(xvbordarcellipse)(int *x, int *y, int *width, int *height,
                            float *angle1, float *angle2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y || !width || !height || !angle1 || !angle2) return;  // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;

    const QRect bbox(*x - *width, *y - *height,
                     *width * 2,  *height * 2);
    const int start_16 = static_cast<int>(*angle1 * 16.0f);
    const int span_16  = static_cast<int>(*angle2 * 16.0f);

    st->painter_->drawArc(bbox, start_16, span_16);  // outline -- matches XDrawArc

    // Phase 7 Plan 03: PS "el" emit. Pass canvas-Y; helper Y-flips internally.
    XvueApp::psEmitter().bordarcellipse(*x, *y, *width, *height,
                                         *angle1, *angle2);

    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) win->canvas()->update();
    // WR-02: deferred flush -- xvvoir_/xvpause_ pump the event loop.
}

// ---- 47. xvarcellipse_ (D-14, RESEARCH Q1 CORRECTION: drawPie, DRAW-05) ----
// Legacy xvuelc.c:2616 -- uses XFillArc (filled pie slice to center).
// Qt's equivalent of XFillArc is drawPie, NOT drawArc. drawArc would only
// stroke the curve; XFillArc fills the pie wedge. See 02-RESEARCH.md Q1.
// Phase 7 Plan 03: emits PS "El" arc-fill opcode (xvuelc.c:2791).
void proc(xvarcellipse)(int *x, int *y, int *width, int *height,
                        float *angle1, float *angle2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y || !width || !height || !angle1 || !angle2) return;  // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;

    const QRect bbox(*x - *width, *y - *height,
                     *width * 2,  *height * 2);
    const int start_16 = static_cast<int>(*angle1 * 16.0f);
    const int span_16  = static_cast<int>(*angle2 * 16.0f);

    st->painter_->drawPie(bbox, start_16, span_16);  // filled wedge -- matches XFillArc

    // Phase 7 Plan 03: PS "El" emit. Pass canvas-Y; helper Y-flips internally.
    XvueApp::psEmitter().arcellipse(*x, *y, *width, *height,
                                     *angle1, *angle2);

    xvue_qt_mark_user_content(win.get());
    if (win->canvas()) win->canvas()->update();
    // WR-02: deferred flush -- xvvoir_/xvpause_ pump the event loop.
}

// ---- 48. tempscpu_ ----
// Returns CPU time used in seconds (clock() / CLOCKS_PER_SEC).
void proc(tempscpu)(double *tclock) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (tclock == nullptr) return;
    *tclock = static_cast<double>(clock()) / static_cast<double>(CLOCKS_PER_SEC);
}

// ---- 49. secondes1970_ ----
// Seconds since epoch + microsecond counter to disambiguate same-second calls
// (legacy parity per xvuelc.c:2837-2838).
static long nbs1970_counter = 0;
void proc(secondes1970)(double *secondes) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (secondes == nullptr) return;
    nbs1970_counter += 1;
    *secondes = static_cast<double>(time(nullptr))
              + static_cast<double>(nbs1970_counter) * 0.000001;
}

// ---- 50. secondes1969_ ----
// Plain seconds since epoch (no microsecond disambiguation).
void proc(secondes1969)(double *secondes) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (secondes == nullptr) return;
    *secondes = static_cast<double>(time(nullptr));
}

// ---- 51. nomordinateurhote_ ----
// Returns the hostname into host[*nbcar].
void proc(nomordinateurhote)(char *host, int *nbcar) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (host == nullptr || nbcar == nullptr) return;
    char name[80];
    if (gethostname(name, sizeof(name)) != 0) {
        std::fprintf(stderr, "nomordinateurhote: gethostname failed\n");
        *nbcar = 0;
        return;
    }
    name[sizeof(name) - 1] = '\0';
    const int len = static_cast<int>(std::strlen(name));
    *nbcar = len;
    std::memcpy(host, name, static_cast<size_t>(len));
}

// ---- 52. ladate_ ----
// Returns current date as year (since 1900 — matches legacy struct tm), month, day.
void proc(ladate)(int *a, int *m, int *j) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (a == nullptr || m == nullptr || j == nullptr) return;
    time_t t = time(nullptr);
    struct tm *pt = localtime(&t);
    if (pt == nullptr) { *a = 0; *m = 0; *j = 0; return; }
    *a = pt->tm_year;
    *m = pt->tm_mon + 1;
    *j = pt->tm_mday;
}

// ---- 53. heureminuteseconde_ ----
// Returns current hour, minute, second, milliseconds (millis always 0 — legacy parity).
void proc(heureminuteseconde)(int *h, int *m, int *s, int *millis) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (h == nullptr || m == nullptr || s == nullptr || millis == nullptr) return;
    time_t t = time(nullptr);
    struct tm *pt = localtime(&t);
    if (pt == nullptr) { *h = 0; *m = 0; *s = 0; *millis = 0; return; }
    *h = pt->tm_hour;
    *m = pt->tm_min;
    *s = pt->tm_sec;
    *millis = 0;
}

// ---- 54. valvarenv_ ----
// Reads the named environment variable into val[*lval_admis]. lval_trouve set to
// -1 if not found, else the value's full length (may exceed lval_admis if truncated).
// Note: nom is expected to be CHAR(0)-terminated by the Fortran caller.
void proc(valvarenv)( char *nom, int *lval_admis,
                      char *val, int *lval_trouve ) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (nom == nullptr || lval_admis == nullptr ||
        val == nullptr || lval_trouve == nullptr) return;
    const char *ptenv = getenv(nom);
    if (ptenv == nullptr) {
        *lval_trouve = -1;
        return;
    }
    int i = 0;
    while (ptenv[i] != '\0') {
        if (i < *lval_admis) val[i] = ptenv[i];
        ++i;
    }
    *lval_trouve = i;
}

// ---- 55. xvinitierps_ ----
void proc(xvinitierps)(int *modeps) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvinitierps_");
    (void)modeps;
}

// ---- 56. xvimprimerps_ ----
void proc(xvimprimerps)(char nomfichier[], int *length) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvimprimerps_");
    (void)nomfichier; (void)length;
}

// ---- 57. xvsauverps_ ----
void proc(xvsauverps)(char nomfichier[], int *length) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvsauverps_");
    (void)nomfichier; (void)length;
}

// ---- 58. xvue_module_init_ ----
// Phase 6.0 Plan 06: dispatch body replacing the Plan 01 scaffold stub.
// Called once per pp*_qt process from prpr/pp*.f AFTER xvinitgraphique_
// (6.1..6.5 add the CALL XVUE_MODULE_INIT('...') lines). Step 3 of this
// body touches consoleDock_ via installStdoutRedirect which requires
// XvueWindow to exist, and XvueWindow is constructed inside
// xvinitgraphique_. Empirically confirmed by
// test_xvue_qt_window_chrome.cpp:234 which invokes xvue_module_init_
// AFTER xvinitgraphique_.
//
// Pure 6.0 builds do not call this — the menuBridge stays with no
// registered module, and the shell still comes up with the
// {File, View, Help} shared menus only.
//
// T-06.0-01-01 mitigation: clamp attacker-controlled name length to [0, 32];
// treat name_len <= 0 as zero; never deref a null name buffer. The module
// names actually used (mail, elas, flui, ther, nlse) are all 4 chars, so 32
// is generously above the maximum.
//
// Steps:
//   1. Initialize XvuePrefs with the per-module key group.
//   2. Apply the persisted color-scheme preference (UX-13).
//   3. Wire the console-dock stdout redirect (UX-04). Idempotent — Plan 04's
//      installStdoutRedirect short-circuits when readFd_ >= 0.
//   4. Dispatch to the per-module registerXxxActions stub. 6.1..6.5 each
//      override one stub with a stronger symbol.
//   5. Mark the menu bridge as "module registered" sentinel.
//   6. Flip canvas has_user_content_ true so Plan 05's empty-state hint
//      stops rendering.
void proc(xvue_module_init)(char *name, int *name_len) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();

    // T-06.0-01-01: clamp attacker-controlled length.
    int n = (name_len ? *name_len : 0);
    if (n < 0)  n = 0;
    if (n > 32) n = 32;
    QString mod;
    if (n > 0 && name) {
        mod = QString::fromLatin1(name, n).trimmed();
    }

    // 1. Initialize persistence with per-module key group. Safe even for
    //    empty mod (XvuePrefs::initialize tolerates null-or-empty per Plan 02).
    XvuePrefs::initialize(mod.toLocal8Bit().constData());

    // 2. Apply the persisted color-scheme preference. Re-applies even if
    //    XvueWindow::ctor already did so in case xvue_module_init_ runs after
    //    a flag-flip via Preferences dialog.
    XvueApp::applyColorSchemePreference();

    // 3. Wire stdout redirect to the console dock. Must happen AFTER window
    //    is built (dock owns the QSocketNotifier). Idempotent — Plan 04's
    //    installStdoutRedirect short-circuits when readFd_ >= 0.
    auto& win = XvueApp::window_slot();
    if (win && win->consoleDock()) {
        win->consoleDock()->installStdoutRedirect();
    }

    // 4. Dispatch to the module-specific registerXxxActions. Each stub is a
    //    warn-once no-op in 6.0; 6.1..6.5 override via stronger symbol.
    XvueWindow*     wptr = win.get();
    XvueMenuBridge* mbptr = (win ? win->menuBridge() : nullptr);
    if      (mod == QStringLiteral("mail"))  registerMailActions_stub_(wptr, mbptr);
    else if (mod == QStringLiteral("elas"))  registerElasActions_stub_(wptr, mbptr);
    else if (mod == QStringLiteral("flui"))  registerFluiActions_stub_(wptr, mbptr);
    else if (mod == QStringLiteral("ther"))  registerTherActions_stub_(wptr, mbptr);
    else if (mod == QStringLiteral("nlse"))  registerNlseActions_stub_(wptr, mbptr);
    else {
        std::fprintf(stderr,
            "xvue-qt: xvue_module_init_('%s'): unknown module; shell comes up "
            "with shared menus only.\n",
            mod.toUtf8().constData());
    }

    // 5. Sentinel: mark the module as registered so any future Layer-2 assert
    //    in buildMenuBar (if added) does not fire. Safe to mark even on the
    //    "unknown module" branch — the shell still works without per-module
    //    QActions.
    if (win && win->menuBridge()) {
        win->menuBridge()->markModuleRegistered();
    }

    // 6. Flip canvas empty-state off: once a module is loaded, the canvas is
    //    "live" in the user's mental model even if no drawing has happened
    //    yet. Plan 05's paintEvent consults state->has_user_content_.
    if (win && win->state()) {
        win->state()->has_user_content_ = true;
    }
    if (win && win->canvas()) {
        win->canvas()->update();
    }
}

} // extern "C"
