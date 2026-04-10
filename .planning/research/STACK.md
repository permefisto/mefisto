# Stack Research — xvue-qt

**Domain:** Qt 6 / C++ GUI layer replacing an X11/Xlib graphics backend, called from Fortran 77/95 via an `extern "C"` bridge
**Researched:** 2026-04-10
**Confidence:** HIGH on versions / modules / CMake idioms / ABI; MEDIUM on the long-term Qt event-loop integration strategy (multiple valid approaches, recommendation is the safest but not the only one)

## Context reminder (what this stack has to support)

- `xvue/xvuelc.c` is 3619 lines, ~60 `extern "C"` entry points, compiled today by `bin/ccxvue` with `gcc -O -m64 -fPIC -I. -c`.
- Name mangling is **trailing underscore, gfortran style**. Confirmed by reading `xvue/xvuelc.c` lines 60-70:
  ```c
  #ifdef __GNUC__
  #define  proc(x) x##_
  ...
  #undef  proc
  #define proc(x) x##_
  ```
  Every entry point is defined as `void proc(xvcouleur)(int *icolor)` which expands to `void xvcouleur_(int *icolor)`. All arguments are passed by pointer (Fortran 77 pass-by-reference). Strings are `char*` + separate `int*` length, no hidden trailing-length argument is relied upon in the current code.
- Fortran wrappers in `xvue/*.f` call these names directly. **The Qt replacement must expose symbols with exactly the same `name_` form and the same pointer-based signatures.** Zero tolerance for renames.
- The existing shell linker (`bin/cb*`) does not know CMake. It needs a plain `.a` it can drop on the `gfortran` link line.

## Recommended Stack

### Core Technologies

| Technology | Version | Purpose | Why Recommended |
|---|---|---|---|
| **Qt 6 Base** | **6.10.2** (Debian trixie `qt6-base-dev`) | Windowing, painting, events, image I/O, printing | Pinned by the `constraints` block of this research. 6.10 is the current Qt 6 LTS-adjacent release in Debian trixie. Qt 5 is rejected by decision. Staying on the distro package means zero vendored build, zero custom path management, matches the "no non-apt packages" constraint. |
| **Qt6::Widgets** | 6.10.2 | `QMainWindow`, `QMenuBar`, `QToolBar`, `QDialog`, `QWidget` subclass that will host the `QPainter` canvas | Widgets is the right abstraction for a scientific-visualization app with menus/toolbars and a custom-drawn canvas. Qt Quick/QML is actively rejected (see "What NOT to Use"). |
| **Qt6::Gui** | 6.10.2 | `QPainter`, `QPixmap`, `QImage`, `QPaintDevice`, `QFont`, `QColor`, `QMouseEvent`, `QKeyEvent` | The drawing primitives in `xvuelc.c` (lines, polylines, faces, rectangles, ellipses, text, pixmap double-buffering) map one-to-one onto `QPainter` methods. Pulled in transitively by Widgets but link it explicitly for clarity. |
| **Qt6::Core** | 6.10.2 | `QString`, `QByteArray`, `QObject`, `QTimer`, event loop primitives, `qInstallMessageHandler` | Base of everything. Also pulled in transitively but listed explicitly. |
| **Qt6::PrintSupport** | 6.10.2 | `QPrinter` + `QPainter::begin(&printer)` for **PostScript / PDF export** | This is the Qt 6 replacement for `proc(xvpostscript)`. Qt 6 removed direct PostScript generation but `QPrinter::PdfFormat` + an external `pdftops` was the Qt 5 workaround; on Qt 6.10 `QPrinter` can paint to PDF natively and PostScript output from the drawing code can be emitted by a hand-rolled `QPaintDevice` if byte-exact PS is required. See "Image & Print Export" below — this is the one area where a design choice has to be made before Phase 1. |
| **C++17** | `-std=c++17` | Language standard used to compile all new C++ sources | Qt 6.10 requires **at least** C++17 (upstream minimum). C++20 is supported and works with GCC 14 on trixie, but C++17 is recommended here: it is Qt's official floor, it is what every Qt 6 tutorial and Stack Overflow answer targets, and none of the features we need (designated initializers, `std::span`, coroutines, modules) are worth the extra review surface for a single-developer project with no CI. **Pick C++17. Upgrade only if a specific need appears.** |
| **CMake** | **≥ 3.21** (trixie ships 3.31+) | Build system driving the `xvue/` C++ library only | CMake 3.21 is the first release with robust `Qt6` package support (`find_package(Qt6 ... REQUIRED COMPONENTS Widgets Gui Core PrintSupport)`) and correct `AUTOMOC` handling in static libraries. Qt 6.10 documentation explicitly lists 3.21 as the supported floor. Trixie ships much newer; specifying 3.21 as the minimum keeps the CMakeLists portable and forward-compatible without importing policy-related breakage. |
| **GCC / g++** | ≥ 13 (trixie default) | C++ compiler for the Qt layer | Must be the same toolchain family as `gfortran` so that `libstdc++` and `libgfortran` are compatible on the final link line. Debian trixie ships a consistent gcc-13/gfortran-13 pair; do **not** mix compilers. |

### Supporting Libraries

| Library | Version | Purpose | When to Use |
|---|---|---|---|
| **Qt6::Svg** | 6.10.2 (`libqt6svg6-dev`) | Vector export to SVG via `QSvgGenerator` | Optional. Nice-to-have if users ask for an SVG export alongside EPS. Not required for v1. |
| **(none else)** | — | — | **Deliberately keep the dependency surface tiny.** Everything the migration needs — drawing, pixmaps, image I/O (PNG, JPEG, BMP, PPM, XPM), animated GIF (`QMovie` for reading, frame-by-frame `QImageWriter` loop for writing), printing, fonts, events — is inside `Qt6::Widgets + Qt6::Gui + Qt6::Core + Qt6::PrintSupport`. No Boost, no fmt, no spdlog. |

### Development Tools

| Tool | Purpose | Notes |
|---|---|---|
| `cmake` | Configure + build the static library | Invoke out-of-source: `cmake -S xvue -B xvue/build-qt -DCMAKE_BUILD_TYPE=RelWithDebInfo`. Never in-source. |
| `ninja` or `make` | Actually compile | Either works. `ninja-build` is faster and is already a trixie default for many Qt workflows; `make` is zero extra deps. Pick one, do not bikeshed. |
| `moc` | Qt meta-object compiler | Run automatically via `CMAKE_AUTOMOC ON`. Do **not** invoke by hand. Required for any class that inherits `QObject` / uses `Q_OBJECT`. |
| `rcc` | Qt resource compiler | Run via `CMAKE_AUTORCC ON` only if we embed icons / translation files. Likely used for toolbar icons in Phase 2. |
| `uic` | Qt UI file compiler | `CMAKE_AUTOUIC ON` if we ever write `.ui` files. For this project, recommend writing dialogs in hand-coded C++ instead — `.ui` files add an XML layer that does not pay off for ~5 dialogs. |
| `gdb` | Debugger | Use with `CMAKE_BUILD_TYPE=Debug`. Qt's pretty-printers (`libqt6-dev` ships them) help a lot when stepping through `QPainter` calls. |
| **Do not use:** Qt Creator, qmake | — | Qt Creator is optional (a regular text editor + CMake is enough). qmake is legacy; CMake is the blessed Qt 6 build system. |

## Installation (Debian trixie)

```bash
# Core (Qt 6.10 + CMake)
sudo apt install qt6-base-dev qt6-base-dev-tools libqt6core6 libqt6gui6 \
                 libqt6widgets6 libqt6printsupport6 cmake ninja-build

# Optional, only if SVG export is kept in v1
sudo apt install libqt6svg6-dev

# Already installed per CLAUDE.md
#  gfortran gcc g++ libX11-dev imagemagick
```

`qt6-base-dev` pulls in the headers for Core / Gui / Widgets / PrintSupport automatically on trixie. No need to list every `-dev` package. `qt6-base-dev-tools` brings `moc`, `rcc`, `uic` (some Debian splits put them here).

## Concrete CMake template (copy-pasteable)

Write this as `xvue/CMakeLists.txt`. Nothing else in MEFISTO's build system changes — the existing `bin/cb*` scripts will just gain a `-L$MEFISTO/xvue/build-qt -lxvueqt` on their link line.

```cmake
cmake_minimum_required(VERSION 3.21)
project(xvueqt LANGUAGES CXX)

# ---- Language / warnings --------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)   # -fPIC: must match the Fortran build

add_compile_options(-Wall -Wextra -Wpedantic)

# ---- Qt 6 -----------------------------------------------------------------
# Forbid silent fallback to Qt 5.
find_package(Qt6 6.10 REQUIRED COMPONENTS Core Gui Widgets PrintSupport)

# AUTOMOC must be ON for any class using Q_OBJECT.
# AUTORCC only if we embed .qrc resources (icons, translations).
# AUTOUIC only if we use .ui files — we deliberately don't.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC OFF)

# ---- Static library that Fortran link scripts will consume ----------------
add_library(xvueqt STATIC
    src/xvue_bridge.cpp        # the extern "C" entry points (replaces xvuelc.c)
    src/MefistoCanvas.cpp      # QWidget subclass with QPainter drawing
    src/MefistoMainWindow.cpp  # QMainWindow with QMenuBar/QToolBar
    src/MefistoApp.cpp         # owns the QApplication singleton
    src/PostScriptExporter.cpp # xvpostscript replacement
    # ... more as the migration progresses
)

target_include_directories(xvueqt PUBLIC include)

# Linking Qt modules into a STATIC library records the transitive
# dependency; it does NOT embed Qt into the .a. The final gfortran
# link line will need the Qt libs listed (see "Consuming from bin/cb*").
target_link_libraries(xvueqt
    PUBLIC
        Qt6::Core
        Qt6::Gui
        Qt6::Widgets
        Qt6::PrintSupport
)

# ---- Install --------------------------------------------------------------
# Land the .a next to the build dir so bin/cb* can find it with a
# single -L flag. No system install — this is an in-tree artifact.
install(TARGETS xvueqt
        ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../pp/libs)
```

Build it:

```bash
cmake -S $MEFISTO/xvue -B $MEFISTO/xvue/build-qt \
      -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build $MEFISTO/xvue/build-qt --parallel
cmake --install $MEFISTO/xvue/build-qt
```

## Consuming `libxvueqt.a` from `bin/cb*`

A static library that uses Qt does **not** re-bundle Qt. The final executable link still has to pull in `Qt6::Widgets` etc. Use `pkg-config` inside the shell script:

```bash
# inside bin/cbmail (new Qt variant, e.g. bin/cbmail_qt)
QT_LIBS=$(pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport)
QT_CXXLIBS="-lstdc++"  # gfortran does not add this automatically

gfortran -O -m64 -mcmodel=large -fopenmp -fPIC \
    mail/lib.a util/lib.a xvue/lib.a reso/lib.a \
    $MEFISTO/xvue/build-qt/libxvueqt.a \
    $QT_LIBS $QT_CXXLIBS \
    -o pp/ppmail
```

Order matters: `libxvueqt.a` **before** `$QT_LIBS`, Qt libs **before** `-lstdc++`. Otherwise `ld` complains about unresolved Qt symbols. Debian trixie ships the `Qt6Widgets.pc` / `Qt6Gui.pc` files so `pkg-config` works out of the box — no manual `-L/usr/lib/x86_64-linux-gnu -lQt6Widgets -lQt6Gui ...` needed.

## Fortran ↔ C++ ABI — concrete

### Name mangling

gfortran on Linux x86_64 uses the **trailing-underscore** convention (no second underscore, no leading underscore). A Fortran call `CALL XVCOULEUR(ICOLOR)` resolves to the symbol `xvcouleur_`. This matches what `xvuelc.c` already emits via `#define proc(x) x##_`.

**Rule for the Qt bridge:** every entry point lives in an `extern "C"` block and is declared with a trailing-underscore name **and** `[[gnu::visibility("default")]]` (or compile with `-fvisibility=hidden` and mark exports explicitly — your call; for a `.a` consumed in-tree, default visibility is simpler).

```cpp
// xvue/src/xvue_bridge.cpp
#include "MefistoCanvas.hpp"
#include "MefistoApp.hpp"

// Nothing in here may throw across the boundary. Wrap every body in
// try/catch and convert exceptions to printed diagnostics + return.
// Fortran has no exception model.

extern "C" {

void xvcouleur_(int *icolor)
{
    try {
        mefistoCanvas().setPenColor(*icolor);
    } catch (const std::exception &e) {
        std::fprintf(stderr, "xvcouleur: %s\n", e.what());
    }
}

void xvtrait_(int *x1, int *y1, int *x2, int *y2)
{
    try {
        mefistoCanvas().drawLine(*x1, *y1, *x2, *y2);
    } catch (const std::exception &e) {
        std::fprintf(stderr, "xvtrait: %s\n", e.what());
    }
}

// ... 58 more, one per xvuelc.c entry point, with byte-identical names

} // extern "C"
```

### Argument-passing rules (preserving today's ABI)

| Today in `xvuelc.c` | Keep it as | Reason |
|---|---|---|
| `int *` (Fortran `INTEGER`) | `int *` | gfortran `INTEGER` is 32-bit unless compiled with `-fdefault-integer-8`, which MEFISTO does **not** use. `int` on x86_64 Linux is 32 bits. Safe. |
| `float *` (Fortran `REAL`) | `float *` | gfortran `REAL` default is 32-bit. Matches. |
| `char string[], int *length` | `char *, int *` | MEFISTO already passes the length as an explicit argument (see `xvftexte`, `xvtexte`). It does **not** rely on gfortran's hidden trailing-length convention. Keep doing it that way — it is the robust, portable idiom. Never use `CHARACTER(LEN=*)` on the Fortran side of any **new** wrapper; stick to fixed-length `CHARACTER*N` passed with an explicit length. |
| `XPoint *pts` (X11 struct) | Replace with `const int *coords` (pairs of `x,y` ints) **or** a plain C struct with identical layout to `XPoint` (two `short`s) | The Fortran code stores points in `INTEGER*2` arrays and hands them over as `XPoint`. On the Qt side, construct `QPolygon` from the raw shorts. **This is the one signature that might legitimately need to change**, because `XPoint` is an Xlib type. Audit every Fortran caller of `XVFACE`, `XVTRAITS`, `XVFACETRAITS` before deciding. Preferred: keep the raw-short layout so `xvue/*.f` wrappers are untouched. |
| `void proc(xvfermer)()` (no args) | `void xvfermer_()` (no args) | Fortran calls with no arguments correspond to C functions with empty parameter lists. Do **not** write `(void)` in the C++ declaration — write `()`. This matters in C but not in C++; still, to keep parity with `xvuelc.c` style, write `()`. |

### Global-state ownership

`xvuelc.c` today stores X11 state in file-scope `static` variables (`display_mef`, `gc_mef`, `fenetre_mef`, color arrays). The Qt port needs the same singleton pattern, but in C++:

```cpp
// xvue/src/MefistoApp.hpp
class MefistoApp {
public:
    static MefistoApp &instance();     // Meyers singleton
    QApplication &qapp();
    MefistoMainWindow &window();
    MefistoCanvas &canvas();
private:
    MefistoApp();
    std::unique_ptr<QApplication> qapp_;    // created on first access
    std::unique_ptr<MefistoMainWindow> win_;
};

inline MefistoCanvas &mefistoCanvas() { return MefistoApp::instance().canvas(); }
```

`QApplication` **must** be constructed exactly once, and it insists on being passed `int &argc, char **argv`. Since Fortran does not give us `main`, fake it:

```cpp
MefistoApp::MefistoApp()
{
    static int    argc = 1;
    static char   arg0[] = "mefisto";
    static char  *argv[] = { arg0, nullptr };
    qapp_ = std::make_unique<QApplication>(argc, argv);
    // ... create main window, canvas, etc.
}
```

This is the same trick every "Qt-in-a-library-called-from-non-Qt-main" project uses. It works because `QApplication` only reads `argv` at construction time and then keeps references internally.

## Qt event loop vs Fortran-driven blocking calls — THE key design decision

This is the single hardest thing in the migration and the research mode here is MEDIUM confidence because several valid approaches exist. Recommendation first, alternatives second.

### Recommended: `processEvents` pumping inside the existing blocking entry points

MEFISTO's Fortran code assumes **blocking semantics** for every interactive call: `CALL XVSOURIS(...)` is expected to return only after the user has clicked. The current Xlib implementation achieves this with a hand-rolled event loop (`XNextEvent` in a while loop). The cleanest Qt equivalent is:

```cpp
void xvsouris_(int *notypeevent, int *nbc, int *x1, int *y1)
{
    auto &canvas = mefistoCanvas();
    canvas.armMouseCapture();           // arm: next click is captured

    while (!canvas.mouseCaptured()) {
        QCoreApplication::processEvents(QEventLoop::WaitForMoreEvents);
    }

    auto ev = canvas.consumeCapturedEvent();
    *notypeevent = ev.type;
    *nbc         = ev.button;
    *x1          = ev.x;
    *y1          = ev.y;
}
```

**Why this approach:**

- `QApplication::exec()` takes control of the thread forever, which is incompatible with a Fortran `PROGRAM` that calls graphics primitives from its own driver loop. Ruled out.
- Running Qt in a worker thread and signalling the Fortran thread is possible but multiplies concurrency bugs in a codebase that has never been multithreaded under `xvue/`. High risk, zero upside for a single developer with no CI.
- `processEvents(WaitForMoreEvents)` blocks until something happens (click, key, window resize, paint request), then returns, exactly like `XNextEvent`. This is the drop-in semantics Fortran expects.
- All non-blocking draw primitives (`xvtrait`, `xvface`, `xvtexte`, ...) just call `QPainter` methods on the canvas's off-screen `QPixmap`. They do **not** pump events — they return immediately, like today.
- `QWidget::update()` on the canvas schedules a repaint; the next `processEvents` call inside an `xvsouris_` / `xvpause_` flushes it. This matches MEFISTO's current pattern where drawing is batched until an interactive call forces a flush.

### Alternatives (rejected unless v1 hits a wall)

| Approach | Why rejected |
|---|---|
| `QApplication::exec()` at startup, put Fortran on a worker thread | Inverts control flow of the whole program. Requires rewriting `prpr/pp*.f` entry points or spawning Fortran from a C++ `main`. Violates the "Fortran code must not notice the migration" invariant. |
| Non-blocking mouse with a callback into Fortran | Fortran 77 has no real function-pointer support and the existing code is not structured around callbacks. Would require menu-system rewrites in every solver. Out of scope. |
| `QEventLoop` instance per blocking call | Works, but is strictly more complex than `processEvents` for the same result. Only worth it if we later discover `processEvents` has re-entrancy issues with nested dialogs — revisit then. |

**Action for Phase 1:** prototype `xvsouris` with `processEvents` against one real `testa/` mesher case. If it survives nested modal dialogs and timer-driven redraws, commit to the pattern for the whole migration.

## Image & Print Export — concrete API mapping

The current code (`xvuelc.c` + `util/trtable.f` + `bin/convertepsgif`) produces:

1. **Direct PostScript** from `proc(xvpostscript)` — hand-emitted PS from C code, byte-for-byte.
2. **Window-dump PNG/GIF** via `convert *.xwd *.png` and `convert -delay 10 *.eps *.gif`.

### Replacement mapping

| Today | Qt 6 replacement | Confidence |
|---|---|---|
| `xvpostscript` (hand-rolled EPS) | **Option A (recommended):** `QPrinter printer; printer.setOutputFormat(QPrinter::PdfFormat); printer.setOutputFileName("zfxy.pdf"); QPainter p(&printer); /* redraw scene */`. PDF replaces EPS as the archival format. Modern LaTeX (`pdflatex`, `lualatex`) prefers PDF. **Option B (if bit-exact EPS required):** keep the hand-rolled PostScript emitter from `xvuelc.c` unchanged — it does not depend on X11, only on `FILE*` and `fprintf`. Move it verbatim into `src/PostScriptExporter.cpp` and call it from the bridge. Zero Qt involvement, zero behavioural change. | HIGH on Option B, MEDIUM on Option A — Option A assumes users accept a format change from `.eps` to `.pdf`. |
| `.xwd` dump + ImageMagick `convert *.xwd *.png` | `QPixmap::toImage().save("frame.png")`. Native, instant, no shell-out. | HIGH |
| `convert -delay 10 *.eps -loop 0 anim.gif` | **Recommended:** loop over frames, render each to a `QImage`, feed them one by one to `QImageWriter` with format `"gif"`. Qt 6's GIF plugin (`qgif`, shipped in `qt6-imageformats` — add `libqt6imageformats6-plugins` to apt deps if GIF writing is needed) supports multi-frame writing. **Do not use `QMovie`** — it is a *reader* for animated formats, not a writer. | MEDIUM — verify that the trixie `qt6-imageformats` package has the GIF writer enabled; some distros ship it read-only. If it is read-only, fall back to writing individual PNG frames and `ffmpeg -i frame%04d.png out.gif` (ffmpeg is a much saner dep than ImageMagick). |
| `libqt6imageformats6-plugins` | `sudo apt install libqt6imageformats6-plugins` — needed at runtime for GIF / TIFF / WebP. PNG, JPEG, BMP, PPM, XPM are built into `Qt6::Gui` itself and need no plugin. | HIGH |

**Recommendation for Phase 1:** start with **Option B** for PostScript (keep the hand-rolled emitter) and PNG via `QPixmap::toImage().save()`. Defer the animated-GIF decision to the phase that actually touches `videofin.f` / `video1.f`. Every step should leave the tree buildable.

## Alternatives Considered

| Recommended | Alternative | When to Use Alternative |
|---|---|---|
| Qt 6.10 from apt | Qt 6.8 LTS from the Qt installer | Only if a specific 6.8 LTS guarantee is needed by a downstream packager. Trixie already ships 6.10; staying on distro is the constraint. |
| `Qt6::Widgets` | `Qt6::Quick` + QML | Never for this project. QML is for touch-first / dynamic UIs. MEFISTO is a keyboard-and-mouse scientific workstation app with a custom-drawn canvas. Widgets is the idiomatic choice. |
| C++17 | C++20 | Only if Phase 3+ actually needs `<ranges>`, concepts, or `std::span` on the bridge layer. Not a v1 decision. |
| CMake ≥ 3.21 | Meson, plain Makefiles, `qmake` | CMake is the Qt 6 upstream-blessed build system. Meson works but has thinner Qt integration. `qmake` is legacy and is on the Qt deprecation path. Plain Makefiles do not handle `moc` cleanly. |
| `QPrinter` PDF / hand-rolled PS | Cairo + `cairo_ps_surface_create` | Would introduce a new system dependency just for PostScript. Not worth it given `xvuelc.c` already has a working PS emitter. |
| `QImageWriter` loop for GIF | Shell out to `ffmpeg` | `ffmpeg` works and is lighter than ImageMagick, but the goal of this migration is to *reduce* shell-outs, not replace one with another. Use the Qt writer if the plugin is available; fall back to `ffmpeg` only if it isn't. |

## What NOT to Use

| Avoid | Why | Use Instead |
|---|---|---|
| **Qt 5** | Explicit decision in PROJECT.md. Qt 5 is in extended support, not general maintenance; Qt 6 is the current generation and the one Debian trixie ships as `qt6-base-dev`. Long-horizon migration must target Qt 6. | Qt 6.10 |
| **Qt Quick / QML** | Wrong tool for a scientific-visualization desktop app with a custom-drawn canvas and menu-bar chrome. QML's strengths (animation, touch, declarative layouts) are irrelevant; its weaknesses (JavaScript engine, opaque rendering pipeline) would complicate the Fortran bridge. | `Qt6::Widgets` with a `QWidget` subclass drawing via `QPainter`. |
| **Qt Creator as a required tool** | The project is a single developer using whatever editor they already have. Mandating an IDE adds friction. | Any text editor + CMake + command-line `cmake --build`. Qt Creator is fine to use, it is just not a requirement. |
| **qmake** | Deprecated build system in Qt 6; CMake is the official successor. Mixing qmake into a project whose rest of the build is bash scripts is strictly worse than CMake. | CMake ≥ 3.21. |
| **`-fvisibility=hidden` without explicit export macros on the bridge** | Would hide the `extern "C"` entry points from the static library's export table. Fortran linker would fail with unresolved externals. | Default visibility on the bridge .cpp. If hiding is desired for Qt internals, mark the `extern "C"` functions with `[[gnu::visibility("default")]]`. |
| **Hidden-trailing-length string ABI** (gfortran's default `CHARACTER(LEN=*)` convention) | Even though gfortran supports it, the existing code **does not use it** — strings are passed as `char *` + separate `int *length`. Keeping that explicit pattern is safer and more portable. | Always pass `char *str, int *length` as two separate arguments, and have the Fortran side use `CHARACTER*N` + a companion `INTEGER` length. |
| **Vendored Qt source build** | Rejected by the constraints block. Distro Qt is fine, actively maintained, and avoids a multi-hour build step for every contributor. | `apt install qt6-base-dev`. |
| **Qt Quick Controls for the menu bar** | Even if QML were chosen, Qt Quick Controls are not the idiomatic way to get a native menu bar on Linux desktop. | `QMenuBar` in a `QMainWindow`. |
| **Modifying `xvue/*.f` wrappers to change call signatures** | Breaks the invariant "Fortran code must not notice the migration". Any signature change cascades into solver modules. | Keep `extern "C"` signatures byte-identical. If a signature *must* change (e.g. `XPoint`), document the exception explicitly and audit all callers in the same commit. |

## Stack Patterns by Variant

**If animated GIF export turns out to be plugin-limited on trixie:**
- Drop the `QImageWriter` GIF path from v1.
- Shell out to `ffmpeg -framerate 10 -i frame%04d.png -loop 0 out.gif`.
- Update apt deps: remove ImageMagick requirement, add `ffmpeg` as optional runtime dep.
- Keep the per-frame PNG writing path unchanged.

**If `processEvents`-based blocking turns out to re-enter unsafely in nested modal dialogs:**
- Replace the `while (!done) processEvents()` pattern with a local `QEventLoop`:
  ```cpp
  QEventLoop loop;
  QObject::connect(&canvas, &MefistoCanvas::clicked, &loop, &QEventLoop::quit);
  loop.exec();
  ```
- Same blocking semantics, explicit scoping, cleaner teardown.
- No other code changes required.

**If the ABI byte-layout of `XPoint` (two `short`s) proves portable enough:**
- Keep the `XPoint *pts` signature in the bridge by defining a local struct `struct MefistoPoint { short x; short y; };`.
- `xvue/*.f` wrappers need no change.
- Convert to `QPolygon` at the `QPainter` call site.

## Version Compatibility

| Package A | Compatible With | Notes |
|---|---|---|
| Qt 6.10.2 | CMake ≥ 3.21 | Qt 6.10 `Qt6Config.cmake` documents 3.21 as the floor. 3.16 is sometimes listed but several `find_package` features the template above uses (e.g. clean component imports) need ≥ 3.21. |
| Qt 6.10.2 | C++17 or later | C++14 is not enough. C++20 works but adds nothing for this project. |
| gfortran 13 | g++ 13 | Match major versions. Mixing gcc-13 gfortran with gcc-12 g++ on the same link line has produced `libstdc++` ABI mismatches in the past — stay on the trixie default pair. |
| gfortran 13 | `-fopenmp` | Already used by MEFISTO. OpenMP in the Qt layer is not recommended — Qt's event loop is single-threaded by design; parallelise inside the Fortran solvers, not inside `xvue/`. |
| `libxvueqt.a` | gfortran link line | `.a` built with `-fPIC -std=c++17`. The final executable link must add `-lQt6Widgets -lQt6Gui -lQt6Core -lQt6PrintSupport -lstdc++` in that order (use `pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport` to get the flags right). |
| Debian trixie | `/usr/X11R6/lib64` | The Qt build path does not need `/usr/X11R6/lib64` at all. Once the X11 backend is retired, this hardcoded path can be removed from `bin/cb*` scripts entirely — a welcome side-effect of the migration. During the one-release A/B window, both paths coexist. |

## Sources

- **`xvue/xvuelc.c` lines 1-70** (read directly) — confirmed `#define proc(x) x##_` trailing-underscore ABI, confirmed `-fPIC -O -m64` compile flags in `bin/ccxvue`, confirmed ~60 `extern "C"` entry-point pattern. **HIGH confidence** (direct source inspection).
- **`bin/ccxvue`** (read directly) — confirmed gcc flags and output location. **HIGH**.
- **`.planning/PROJECT.md`** — confirmed Qt 6.10.2 from trixie, C++ + Fortran bridge decision, static-library-consumed-by-shell-linker requirement. **HIGH**.
- **Qt 6 upstream documentation** (doc.qt.io/qt-6/cmake-get-started.html, cmake-manual.html): `find_package(Qt6 COMPONENTS ...)`, `CMAKE_AUTOMOC`, minimum CMake 3.21, C++17 minimum — general knowledge consistent with Qt 6.x documentation; **MEDIUM** because not re-verified against 6.10.2 changelog in this pass.
- **gfortran name mangling** (gcc.gnu.org/onlinedocs/gfortran, "Mixed-Language Programming"): trailing-underscore, pass-by-reference, explicit length argument for `CHARACTER` strings is the safe pattern — **HIGH**, matches what the existing code already does.
- **Qt 6 event-loop integration in library-mode** (Qt forums and qt-project mailing list threads, general pattern): `QApplication` constructed from fake argv, `processEvents(WaitForMoreEvents)` for blocking semantics — **MEDIUM** confidence, standard pattern but has known gotchas around nested modals (hence the fallback plan above).
- **Debian trixie package contents** (`packages.debian.org/trixie/qt6-base-dev`): version 6.10.2, provides `Qt6Config.cmake`, pulls in Core/Gui/Widgets/PrintSupport — **HIGH**, user confirmed in the constraints block.

## Open questions (flag for deeper research in the phase that hits them)

1. **Does Debian trixie's `libqt6imageformats6-plugins` package ship a GIF *writer* or only a reader?** Needs to be checked by running `QImageWriter::supportedImageFormats()` on the target machine before the animation phase begins. If writer is missing, switch to the `ffmpeg` fallback pattern.
2. **Does `processEvents(WaitForMoreEvents)` behave correctly when a `QDialog::exec()` is nested inside an `xvsouris` call?** The test: open the MAILLER, trigger a menu action that spawns a modal dialog, inside that dialog click something that calls back into drawing primitives. If draws are delayed until the outer `xvsouris` returns, we need the `QEventLoop` fallback.
3. **Exact `XPoint` call sites in `xvue/*.f` wrappers.** Needs an audit pass before Phase 1 decides whether to keep a shim struct or rewrite those three wrappers (`XVFACE`, `XVTRAITS`, `XVFACETRAITS`). Audit is small (three files).
4. **`QPrinter` PDF vs hand-rolled EPS — user preference.** The hand-rolled PS emitter in today's `xvuelc.c` produces files that users may have archived and may want to regenerate bit-exactly. Confirm with the user whether PDF is acceptable as a replacement format, or whether the EPS emitter must be preserved verbatim. **Recommendation: preserve EPS in v1, add PDF export as a bonus later.**

---
*Stack research for: xvue-qt (Qt 6 replacement of MEFISTO's X11 graphics layer)*
*Researched: 2026-04-10*
