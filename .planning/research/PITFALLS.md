# Pitfalls Research

**Domain:** Qt 6 migration of an X11/Xlib graphics library called from Fortran 77/95 via `extern "C"`
**Researched:** 2026-04-10
**Confidence:** HIGH on ABI and build-integration traps (grounded in direct reads of `xvue/xvuelc.c` and `bin/ccxvue`); MEDIUM on event-loop and HiDPI traps (depend on empirical behavior of `processEvents` inside nested `QDialog::exec`); MEDIUM on image-export parity (depends on trixie's `qt6-imageformats` plugin bundle).

## Context

This project replaces `xvue/xvuelc.c` (3619 lines of Xlib) with a Qt 6 C++ implementation behind a byte-identical `extern "C"` ABI, keeps the X11 backend alive for one A/B release cycle, and is driven by a single developer with no CI and no automated regression harness. Mistakes propagate: a subtle ABI mismatch or an event-loop deadlock will surface only under interactive use in a `testa/` case, often minutes into a mesher session. Every pitfall below is mapped to the roadmap phase where it should be addressed **before** it can bite, and to the warning signs that will appear if it isn't.

The phase numbering below refers to the phase structure proposed by STACK.md and ARCHITECTURE.md:

0. CMake skeleton + empty bridge
1. `XvueApp`/`XvueWindow`/`XvueCanvas` shell
2. Drawing primitives (points, lines, faces, rectangles, ellipses, text)
3. Colors, fonts, pixmap double-buffering
4. Event bridge (mouse, keyboard, blocking `xvsouris*`)
5. Real Qt UI chrome (QMenuBar/QToolBar/QDialog) — Level 3 UX
6. Image / GIF / PostScript export
7. A/B validation on `testa/` subset
8. X11 backend retirement

---

## Critical Pitfalls

### Pitfall 1: Fortran trailing-underscore name mangling mismatch

**What goes wrong:**
Linker error `undefined reference to xvtrait_` (or every single entry point at once) when linking `pp/ppmail_qt` against `libxvueqt.a`.

**Why it happens:**
gfortran on Linux x86_64 emits lowercase names with a **single trailing underscore** (`CALL XVTRAIT(...)` in Fortran → `xvtrait_` at the linker symbol level). The current `xvue/xvuelc.c` handles this with `#define proc(x) x##_` (lines 60–70). The Qt reimplementation must honor the same convention or every Fortran call site fails to resolve.

**How to avoid:**
Preserve the exact `#define proc(x) x##_` macro (or an equivalent `_` suffix) at the top of the new `extern "C"` bridge header. Emit **every** entry point through the macro. Add a build-time check that runs `nm libxvueqt.a | grep -c '_$'` and fails the build if any Fortran-facing symbol is missing the trailing underscore. Also forbid `-fno-underscoring` in `CMakeLists.txt` explicitly — a stray compiler flag can flip the convention silently.

**Warning signs:**
Undefined reference errors at link time on every Fortran-facing function at once; `nm libxvueqt.a | grep xvtrait` shows `xvtrait` without the trailing underscore.

**Phase to address:** Phase 0 (during the one-draw-primitive end-to-end test)
**Severity:** Critical

---

### Pitfall 2: Fortran passes everything by reference — C++ declarations must take pointers

**What goes wrong:**
Crash with a segmentation fault, or silent numeric garbage, when Qt calls `painter->drawLine(x1, y1, x2, y2)` inside `xvtrait_` because the four ints were declared as `int` in C++ but Fortran pushed four `int*`.

**Why it happens:**
Fortran 77/95 passes scalar arguments **by reference**, always. A Fortran `CALL XVTRAIT(X1,Y1,X2,Y2)` with `INTEGER X1,Y1,X2,Y2` hands four pointers to the C++ side, not four ints. The existing `xvuelc.c` signatures all look like `void proc(xvtrait)( int *x1, int *y1, int *x2, int *y2 )` — **every scalar argument is a pointer**. Dereference to use.

**How to avoid:**
Preserve the existing signatures literally. Copy them from `xvuelc.c` into the Qt bridge header without simplification. Add `-Wno-unused-parameter` is fine, but **do not** "clean up" signatures to pass ints by value. Code review checklist item: every int/float argument is `int*`/`float*`, every string is `char* text, int* length`.

**Warning signs:**
Segfault on the first interactive draw; `gdb` backtrace shows `drawLine` with nonsensical coordinate values (looks like memory addresses cast to int).

**Phase to address:** Phase 0 — enforced by copying `xvuelc.c` signatures verbatim
**Severity:** Critical

---

### Pitfall 3: Fortran character strings carry a hidden length argument — but only sometimes

**What goes wrong:**
`xvtexte_` or `xvftexte_` receives garbage text or crashes, because the C++ side ignored or misread the length argument.

**Why it happens:**
gfortran passes a `CHARACTER*N` argument as either (a) an explicit trailing `int` length appended to the end of the argument list, or (b) an explicit `int*` length passed **in order** as a normal argument. The existing `xvuelc.c` uses form (b) — every string entry point is `void proc(xvtexte)( char string[], int *length, int *x1, int *y1 )` with `length` as a regular pointer argument in the declared order, not a hidden trailing int. This is because the Fortran callers pass the length explicitly (e.g. `CALL XVTEXTE(STR, LEN, X, Y)`). **Do not** add a hidden trailing length; the Fortran side isn't using `-fdefault-character-kind` magic.

**How to avoid:**
Again, copy `xvuelc.c` signatures literally. Do not add or remove length arguments. Do not assume gfortran's hidden-trailing-length convention. If a future Fortran caller does pass a `CHARACTER*N` without an explicit length, handle it as a one-off in a dedicated bridge function, do not retrofit the whole API.

**Warning signs:**
Garbled text in the Qt window; `QString::fromLatin1(string, *length)` returning an empty or truncated string; valgrind reports an out-of-bounds read in `xvtexte_`.

**Phase to address:** Phase 2 (when drawing primitives land)
**Severity:** High

---

### Pitfall 4: `XPoint*` in the existing ABI is a short/short struct, not a Qt type

**What goes wrong:**
Entry points `xvface_`, `xvtraits_`, `xvfacetraits_` expect an `XPoint *pts` argument where `XPoint` is an Xlib `typedef struct { short x, y; }`. If the Qt bridge uses `QPoint*` (which is two ints) the caller and callee disagree on element size and every polyline draws the wrong shape (or crashes).

**Why it happens:**
Xlib `XPoint` is two 16-bit shorts for historical reasons. Fortran wrappers have been passing arrays with that memory layout for decades. The existing comment in `xvuelc.c` even flags this: `/*  ATTENTION XPoint p => short  p.x et p.y */` on the `xvface` declaration.

**How to avoid:**
Define a local shim struct in the Qt bridge header that matches the existing layout byte-for-byte:

```c
typedef struct { short x; short y; } MefistoPoint;
```

Use `MefistoPoint*` in the three affected entry points. Inside the C++ body, convert each element to `QPoint` when calling `QPainter::drawPolyline(const QPoint*, int)`. Audit every Fortran call site for these three functions once during Phase 2 to confirm no caller has drifted to a different layout.

**Warning signs:**
Polylines draw the wrong shape (every other point is offset); coordinates look bit-shifted; `sizeof(MefistoPoint) != 4`.

**Phase to address:** Phase 2 (drawing primitives)
**Severity:** Critical

---

### Pitfall 5: `QApplication` double-initialization / lifetime leak

**What goes wrong:**
Crash on the second call into `xvinitgraphique_` with a Qt assertion "QApplication instance already exists", or subtle event-delivery failures because two `QApplication` singletons fight for the same `argc/argv`.

**Why it happens:**
Fortran can call `xvinitgraphique_` multiple times in one process — the mesher reopens a display after saving a project, the solver opens a new one for post-processing, etc. The old X11 code handled this gracefully because Xlib `XOpenDisplay` is idempotent-ish. Qt's `QApplication` is a hard singleton: constructing a second one crashes.

**How to avoid:**
Guard `QApplication` construction with `std::call_once` (or a static pointer + `QCoreApplication::instance()` check) inside `XvueApp::ensure()`. Construct with fake static `argc`/`argv` arrays that outlive the application:

```cpp
static int    g_argc = 1;
static char   g_arg0[] = "mefisto";
static char*  g_argv[] = { g_arg0, nullptr };
```

**Never** destroy the `QApplication` inside `xvfermer_`. Only the process exit (via `std::atexit` registered at first init) tears it down. Document the rule in a comment on `XvueApp::ensure()`.

**Warning signs:**
"QApplication: there can only be one" assertion; mesher crashes on reopening its window; `qApp` returns a dangling pointer mid-session.

**Phase to address:** Phase 1 (XvueApp shell)
**Severity:** Critical

---

### Pitfall 6: `QApplication::exec()` top-level call inverts control flow

**What goes wrong:**
Either the Fortran program hangs forever on the first graphics call (because Qt owns the main loop and Fortran never gets control back), or the graphics calls return immediately without anything being drawn (because no event loop is pumping paint events).

**Why it happens:**
Qt's default paradigm is "`main()` calls `QApplication::exec()` which never returns until the app quits." MEFISTO's paradigm is the inverse: Fortran drives everything, graphics calls are **blocking synchronous subroutine calls** that must return control to the caller as soon as the operation is done. Calling `exec()` at the top level is incompatible.

**How to avoid:**
**Never** call `QApplication::exec()` anywhere in the Qt bridge. Instead:
- For pure draw calls (`xvtrait`, `xvface`, etc.): construct/reuse a `QPainter`, draw, return. Request a repaint via `canvas->update()` which schedules the actual paint for the next event-loop iteration, and pump once with `QCoreApplication::processEvents(QEventLoop::AllEvents)` before returning so the user sees the update.
- For blocking interactive calls (`xvsouris`, `xvsouris2`, `xvpause`): use a **nested local `QEventLoop`** inside the entry point. Quit the loop from a Qt event filter when the expected event (click, key, timeout) arrives. This lets Qt's event pump run while Fortran waits, without ever calling `exec()`.

Document both patterns in `xvue/README_ARCHITECTURE.md` and add a pre-commit grep that fails on any `QApplication::exec` in `xvue/*.cpp`.

**Warning signs:**
Mesher window appears blank; solver hangs on first mouse pick; `strace` shows no `poll()` return; `xvsouris` never returns.

**Phase to address:** Phase 1 (prevent) + Phase 4 (full nested-loop validation)
**Severity:** Critical

---

### Pitfall 7: `processEvents` starvation during long Fortran compute loops

**What goes wrong:**
The Qt window freezes (title bar grays out on GNOME, becomes "Not Responding" on KDE) during a long mesh refinement or solver iteration that does not call any `xvue/` entry point for tens of seconds.

**Why it happens:**
The nested-loop pattern (Pitfall 6) only pumps events when a Fortran call enters the bridge. A tight numerical loop that does no graphics for 30 seconds will freeze the UI for 30 seconds. The X11 backend had the same problem in principle but X11's "not responding" detection is looser than Wayland compositors'.

**How to avoid:**
Do **not** try to fix this by running Qt on a separate thread — that's a massive architectural change and Qt's thread-affinity rules make it dangerous. Instead, document the limitation in `README` and `LISEZMOI`, and add a low-effort mitigation: the Fortran main driver loops already call `xvpause_` or equivalent between iterations for progress reporting. Make those entry points pump `processEvents(QEventLoop::AllEvents, 10)` (with a 10ms budget) so the event loop gets at least a sliver of time per iteration.

**Warning signs:**
Compositor grays out the window during long runs; clicking the close button during a solver iteration does nothing until the iteration finishes; mouse cursor shows the "waiting" spinner over the window.

**Phase to address:** Phase 4 (event bridge), revisited in Phase 7 (A/B validation)
**Severity:** High

---

### Pitfall 8: Modal `QDialog::exec()` inside a nested `QEventLoop` causes re-entrancy confusion

**What goes wrong:**
Opening a modal file-save dialog via a `QAction::triggered` slot while the user is in the middle of an `xvsouris` mouse-pick call causes events to be delivered out of order, mouse clicks to be "stolen" by the dialog, or — in the worst case — a deadlock where the dialog's own nested loop can't finish because the outer `xvsouris` loop is blocking.

**Why it happens:**
`QDialog::exec()` runs its own nested `QEventLoop`. Stacking two nested loops on top of the main loop is supported by Qt in principle, but every re-entrant call site has to be audited, and scientific code is very bad at audits because the interaction patterns are hard to enumerate.

**How to avoid:**
Adopt the "menu actions queue, do not execute" pattern from ARCHITECTURE.md:
- `QAction::triggered` slots **do not** call Fortran or open modal dialogs directly. They push a command string onto a `XvueMenuBridge` queue.
- The next call to `xvsouris_` / `xvsouris2_` returns a synthetic keyboard event (`notypeevent = 2`) carrying the queued command. Fortran dispatches via its existing text-lexicon code path.
- Modal dialogs (file open, save, properties) are only opened at "safe" points: when no `xvsouris` nested loop is active. The easiest enforcement is a counter `XvueApp::blockingDepth()` — actions that need a modal dialog refuse to run (and show a status-bar message "Cannot open dialog while picking — press Esc first") if depth > 0.
- Document the rule as an invariant in `xvue/README_ARCHITECTURE.md`.

**Warning signs:**
Clicking a menu item while in mesher-pick mode causes the click to vanish; the save dialog opens but the mesher view is frozen behind it; Qt logs "QEventLoop::exec: instance X has already called exec" warnings.

**Phase to address:** Phase 5 (menu chrome)
**Severity:** High

---

### Pitfall 9: CMake `AUTOMOC` silently skips `.cpp` files that need moc

**What goes wrong:**
Linker error `undefined reference to vtable for XvueCanvas` or `undefined reference to XvueCanvas::staticMetaObject`. Means `moc` did not generate the meta-object for a `QObject` subclass.

**Why it happens:**
`AUTOMOC` scans for `Q_OBJECT` macros. If a `QObject` subclass is declared **in a `.cpp` file** (rather than a `.h` file) `AUTOMOC` historically missed it unless `CMAKE_AUTOMOC_MOC_OPTIONS` was tuned. Also, a stale CMake cache from before `AUTOMOC` was enabled will not re-scan old files until touched.

**How to avoid:**
- Declare every `QObject` subclass in a `.h` file under `xvue/include/` — never in a `.cpp` file.
- Set `set(CMAKE_AUTOMOC ON)` **before** `find_package(Qt6 ...)` in `CMakeLists.txt`.
- After any `xvue/CMakeLists.txt` edit, run `rm -rf xvue/build` and rebuild from scratch to clear the CMake cache.
- Add a `make clean && make` step to the release checklist.

**Warning signs:**
Undefined vtable errors; signals/slots silently not firing; `xvue/build/xvue_autogen/` directory missing expected `moc_*.cpp` files.

**Phase to address:** Phase 0 (CMake skeleton)
**Severity:** Medium

---

### Pitfall 10: Static library built without `-fPIC` fails to link into executables that need PIE

**What goes wrong:**
Link error `relocation R_X86_64_32 against .rodata can not be used when making a PIE object; recompile with -fPIC` when `bin/cbmail_qt` links `libxvueqt.a` into `pp/ppmail_qt`.

**Why it happens:**
Debian trixie's `gcc`/`g++` default to Position-Independent Executables (PIE). A static library consumed by a PIE binary must itself be compiled with `-fPIC`. Qt 6 itself is built with `-fPIC` but a hand-rolled `xvue/` CMakeLists.txt can easily forget the setting.

**How to avoid:**
Add `set(CMAKE_POSITION_INDEPENDENT_CODE ON)` at the top of `xvue/CMakeLists.txt`. Add a build-time sanity check: `readelf -r libxvueqt.a 2>&1 | grep R_X86_64_32` should return nothing.

**Warning signs:**
`relocation R_X86_64_32 ... can not be used` at the Fortran-link stage; works in isolated Qt test executables but fails when linked into the real `pp/pp*` binary.

**Phase to address:** Phase 0 (CMake skeleton)
**Severity:** High

---

### Pitfall 11: `-fopenmp` from the Fortran side conflicts with Qt's threading

**What goes wrong:**
The `pp/ppelas_OMP` (OpenMP-parallel) executable either refuses to link (`multiple definition of __kmpc_fork_call` if two OpenMP runtimes collide) or crashes mysteriously inside `QApplication::sendEvent` when called from an OpenMP worker thread.

**Why it happens:**
All `bin/cb*` shell scripts link with `-fopenmp` unconditionally (see the canonical flag set in CONVENTIONS.md). If `libxvueqt.a` is also compiled with `-fopenmp` by CMake, two OpenMP runtimes get linked. Separately, Qt has strict thread-affinity rules — `QWidget` subclasses can only be touched from the main thread. An OpenMP worker thread calling into `xvue/` would violate this.

**How to avoid:**
- Do **not** pass `-fopenmp` to CMake when building `libxvueqt.a`. The Qt layer is strictly single-threaded; it does not need OpenMP. Guard the flag in `xvue/CMakeLists.txt` with an explicit `target_compile_options(xvueqt PRIVATE -fno-openmp)`.
- Document in CONVENTIONS.md: **all graphics calls must happen on the Fortran main thread**, never inside an `!$OMP PARALLEL` region. The existing X11 code has the same rule implicitly; make it explicit.
- Add a debug-build assertion `Q_ASSERT(QThread::currentThread() == qApp->thread())` at the top of every `extern "C"` entry point.

**Warning signs:**
`multiple definition` link errors; random crashes deep in Qt event code; crashes only in `_OMP` executables, never in non-OMP ones.

**Phase to address:** Phase 0 (caught early) + Phase 7 (A/B validation on `ELASTICER_OMP`)
**Severity:** High

---

### Pitfall 12: Animated GIF writer plugin may not be present in Debian trixie's `qt6-imageformats`

**What goes wrong:**
`QImageWriter` succeeds for PNG/JPEG but `QImageWriter::supportedImageFormats()` returns no `"gif"` entry, so the animated-GIF export silently writes a single-frame PNG with `.gif` extension.

**Why it happens:**
Qt's image-format plugin bundle is split: the core `Qt6Gui` supports PNG/JPEG/BMP/PPM. Animated GIF writing is in the separate `qt6-imageformats` package and may or may not ship a GIF **writer** (it historically shipped a reader first). Trixie's package may or may not have it — needs verification at Phase 6 start.

**How to avoid:**
- At Phase 6 kickoff, run `QImageWriter::supportedImageFormats()` from a 10-line test program and check for `"gif"`. If absent, do not attempt a mid-project fight with Qt plugins; fall back to the documented alternative: write a sequence of PNG frames to a temp dir and shell out to `ffmpeg` (not ImageMagick — ffmpeg is leaner and more reliable for video output).
- Document the chosen path in STACK.md.
- Either way, **drop the ImageMagick dependency** — that was one of the xvue-qt goals and neither fallback needs it.

**Warning signs:**
`QImageWriter::write()` returns false; animated GIF has one frame; `QImageWriter::errorString()` mentions "unsupported format".

**Phase to address:** Phase 6 (export) — run the probe at phase start
**Severity:** Medium

---

### Pitfall 13: PostScript / EPS output fidelity drift vs the hand-rolled emitter

**What goes wrong:**
A scientific user comparing the Qt-backend EPS to the X11-backend EPS notices that text is in a different font, stroke widths are slightly different, or the bounding box is off. They assume the Qt migration broke their plots and revert.

**Why it happens:**
The existing `xvpostscript` entry point in `xvuelc.c` **does not use Xlib** — it's ~120 lines of pure `fprintf` writing PostScript strings directly. It has no X11 dependency. If the Qt reimplementation "helpfully" switches to `QPrinter` with PostScript output, the generated file is byte-different even when visually similar: Qt 6 reduced PostScript support compared to Qt 5, fonts are embedded differently, and the bounding box calculation differs.

**How to avoid:**
**Keep the existing hand-rolled PostScript emitter verbatim in the Qt backend.** Extract the `xvpostscript` body from `xvuelc.c` into a new `xvue/xvue_postscript.cpp` without touching a single `fprintf`. Link it into `libxvueqt.a` unchanged. This is the lowest-risk path and was called out by both STACK.md and ARCHITECTURE.md research.

If PDF output is ever wanted as a **bonus** (not a replacement), add `QPrinter::PdfFormat` via a new entry point, not by modifying `xvpostscript`.

**Warning signs:**
User reports "my plots look different"; `diff` of old-vs-new EPS output shows different PostScript operators; bounding box `%%BoundingBox:` line changes.

**Phase to address:** Phase 6 (export)
**Severity:** High

---

### Pitfall 14: X11-indexed colormap vs Qt 32-bit RGBA color drift

**What goes wrong:**
Color-coded stress / temperature / velocity plots in the Qt backend look slightly different from the X11 backend. Users trust the colors (they encode numerical data), so any drift is a correctness bug, not cosmetic.

**Why it happens:**
The existing `xvue/xvuelc.c` uses X11 indexed colormaps with explicit RGB→colormap-cell allocation (`xvStockeRGBtoColormap`, `xvColormapToRGB`). Qt 6 uses 32-bit RGBA everywhere. Naively porting `xvactivervb_` to `QColor::setRgb(r,g,b)` can introduce 1–2 bit drift because the X11 path quantized to the display's color depth (on old 8-bit displays: dramatically; on 24-bit: not at all, but still via a different path).

**How to avoid:**
- Audit every color-related entry point (`xvcouleur`, `xvactivervb`, `xvfond`, `xvCouleursImposees`, `xvStockeRGBtoColormap`, `xvColormapToRGB`) and document the expected output for each input RGB. Test cases live in `prpr/xvtest[1-4].f` and in the solver post-processing colormap code under `util/`.
- Implement an internal `std::array<QColor, MAX_PALETTE> g_palette` that mirrors the X11 colormap semantically. Every call to `xvactivervb_` populates the palette; every draw call uses `g_palette[index]`. No direct `setRgb` inside draw primitives.
- Freeze the scientific colormaps: they must NOT react to system dark-mode/light-mode. Only the chrome follows the palette.

**Warning signs:**
A/B validation on a colored solver plot shows visibly different shades; stress contour bands at the boundary between two colors move by a pixel.

**Phase to address:** Phase 3 (colors/fonts) + Phase 7 (A/B validation)
**Severity:** High

---

### Pitfall 15: Font metrics drift (X11 font server vs Qt text layout)

**What goes wrong:**
Text labels in the mesher overshoot their intended bounding box, overlap neighboring labels, or get clipped at the window edge — because Qt's `QFontMetrics` reports different widths than Xlib's `XTextExtents` for "the same" font.

**Why it happens:**
The existing `xvchargefonte` / `xvnbpixeltexte` entry points use Xlib's font server to measure text. Qt uses its own text layout (HarfBuzz under the hood in Qt 6). Even for "the same" font name, the reported width differs by a few pixels due to hinting, kerning, and subpixel positioning. Fortran code that places labels based on `xvnbpixeltexte_` return values will render slightly off.

**How to avoid:**
- Keep `xvnbpixeltexte_` as the single source of truth for text measurement. Every Fortran caller already reads it and positions labels accordingly; as long as the returned widths match what the Qt text drawer actually produces, labels will be correct.
- Implement `xvnbpixeltexte_` using `QFontMetrics::horizontalAdvance()` and `QFontMetrics::height()`. Do not cache across font changes.
- Use a **fixed, bundled font** (Qt's own `QFont("Monospace")` or a specific embedded font) rather than whatever the X11 font server picked. This makes the Qt output reproducible and independent of installed system fonts.
- Test on a mesher session with many text labels (`testa/pan2d`, `testa/nafems_le1`).

**Warning signs:**
Labels overlap; node-number labels clip at window edge; A/B comparison shows text in slightly different positions.

**Phase to address:** Phase 3 (fonts) + Phase 7 (A/B validation)
**Severity:** Medium

---

### Pitfall 16: Coordinate system origin mismatch (Y-up vs Y-down)

**What goes wrong:**
Meshes render upside-down in the Qt backend. Every triangle that was at the top of the X11 window is at the bottom of the Qt window.

**Why it happens:**
Both X11 and Qt use Y-down (0,0 at top-left) by default, so this **should not** happen — but the existing Fortran / `xvuelc.c` code sometimes inverts Y via `y = ymax - y` before calling Xlib, because scientific Fortran conventions are often Y-up. If the Qt bridge forgets the same inversion (or applies it twice), the display flips.

**How to avoid:**
Audit every entry point for Y-axis handling during Phase 0. Document in `xvue/README_COORDS.md` whether the bridge is Y-up or Y-down. Pick one convention (recommend: match the existing `xvuelc.c` exactly — if it inverts, the Qt bridge inverts identically; if it doesn't, the Qt bridge doesn't either). Add one `xvtest` case that draws a known asymmetric shape (an "F") to visually confirm orientation.

**Warning signs:**
First mesher test shows the mesh flipped vertically; solver plots show contours inverted; text appears mirrored.

**Phase to address:** Phase 0 — document during the one-primitive end-to-end test
**Severity:** Medium

---

### Pitfall 17: HiDPI scaling silently doubles coordinates

**What goes wrong:**
On a HiDPI display (4K monitor, retina laptop) Qt reports a `devicePixelRatio` of 2.0 and "helpfully" doubles logical coordinates to physical pixels. The Fortran code, which queries `xvpxecran_` for the window size and positions meshes in physical pixels, suddenly sees the mesh drawn at half the intended size (or full size but overflowing the window if it misreads `xvpxecran_`).

**Why it happens:**
Qt 6 distinguishes logical pixels (device-independent) from physical pixels. `QWidget::width()` returns logical; `QWidget::width() * devicePixelRatio()` returns physical. Fortran callers of `xvpxecran_` historically expected physical pixels (because X11 didn't have the distinction). A naive reimplementation returns logical and the Fortran layout math goes wrong.

**How to avoid:**
- Decide the convention **up front** in Phase 0 and document it. Recommendation: return **logical** pixels from `xvpxecran_` and let Qt handle physical-pixel scaling transparently. This gives the best HiDPI experience with zero Fortran changes, because the Fortran layout math is relative to the window size — as long as window size and draw primitives are consistent, the ratio works.
- **Never mix** physical and logical coordinates in the same entry point. If `xvpxecran_` returns logical, all `xvtrait_` coordinates are also logical.
- Test on a 4K display (or fake one with `QT_SCALE_FACTOR=2 ./pp/ppmail_qt`) during Phase 7.

**Warning signs:**
Mesh appears at half size on a 4K display; window resize events report doubled sizes; text is the wrong size relative to the mesh.

**Phase to address:** Phase 1 (XvueWindow shell) + Phase 7 (A/B validation on HiDPI)
**Severity:** Medium

---

### Pitfall 18: Stale `.o` files after a common-block edit — unique to this codebase

**What goes wrong:**
A Qt-layer bug is "fixed" but reappears on the next build, because a stale Fortran `.o` file from before a shared `incl/*.inc` edit is still linked into the executable. The bug isn't in the Qt code at all — it's a desync between the Fortran object files and their common blocks.

**Why it happens:**
CONCERNS.md already documents this: the shell build has no dependency tracker, and changing an `.inc` file silently leaves every `.o` that includes it stale. Adding Qt on top does not fix this; if anything, the extra CMake build step in `xvue/` makes it more likely that the developer cleans one build and not the other.

**How to avoid:**
- Document: **any change to `incl/*.inc` requires `bin/cbl_tout` full rebuild.** Qt build does not change this rule.
- Add a safety net: `bin/cbl_tout_qt` always runs `rm -rf xvue/build pp/` before building, so at least the Qt-side objects are never stale.
- During Phase 7 A/B validation, enforce a clean build before every test-case run.

**Warning signs:**
A bug fix "doesn't stick"; the same test case behaves differently on `make` vs `make clean && make`; `bin/cbl_tout` produces different binaries than a per-module `cb*` script.

**Phase to address:** Phase 0 (CMake skeleton) + Phase 7 (validation discipline)
**Severity:** Medium

---

### Pitfall 19: ImageMagick dependency removed before all call sites are audited

**What goes wrong:**
Dropping the ImageMagick dependency in a late phase breaks a corner of the codebase that still shells out to `convert` for something unrelated to the main EPS→GIF pipeline. The install-script / README updates turn `libmagickcore-dev` from required to unreferenced, but some script in `bin/` still expects it.

**Why it happens:**
`bin/convertepsgif` is the obvious ImageMagick caller, but shell scripts in `bin/` and possibly `td/` may invoke `convert` for other tasks (thumbnail generation, format conversion, print preparation). A grep-and-remove without auditing call sites can leave the install fragile.

**How to avoid:**
Before removing the dependency (end of Phase 6 or start of Phase 8), run `grep -rn 'convert' bin/ td/ testa/ testf/` and audit every match. Replace each one either with a Qt-side equivalent or with a `ffmpeg`/`ImageMagick`-if-present fallback. Update `README` and `LISEZMOI` only **after** the audit passes.

**Warning signs:**
Install on a fresh machine without ImageMagick works for `MAILLER` / `ELASTICER` but fails for an edge-case workflow; a tutorial in `td/` silently produces no output.

**Phase to address:** Phase 8 (X11 retirement)
**Severity:** Low

---

### Pitfall 20: Single-developer no-CI regression drift

**What goes wrong:**
A bug fix in Phase 4 regresses a feature from Phase 2. With no automated tests, the regression isn't noticed until weeks later when running a `testa/` case by hand.

**Why it happens:**
Already documented in CONCERNS.md and TESTING.md: no CI, no automated regression harness. The xvue-qt project multiplies the risk because it has more moving parts than a typical bug-fix change: CMake + Qt + event loop + ABI bridge + UI chrome.

**How to avoid:**
- Pick **5 canonical `testa/` cases** at Phase 0 (one per module: mesher, elas, flui, ther, nlse). Run all 5 through both backends (X11 and Qt) **at the end of every phase**. Log the run in a `.planning/phase-N/VALIDATION.md` file even if the log is just "5/5 passed visually".
- Before merging a phase, the X11 build must also still pass the same 5 cases. This protects against accidental damage to the legacy backend during the A/B window.
- Use `prpr/xvtest[1-4].f` as continuous smoke tests — they are fast (seconds) and exercise the graphics primitives in isolation.

**Warning signs:**
A feature that worked two phases ago no longer works; changes to one `.cpp` file break a seemingly unrelated entry point; A/B comparison suddenly differs where it matched before.

**Phase to address:** Every phase — this is a process rule, not a code change
**Severity:** High

---

## Confidence Assessment

| Pitfall area | Confidence | Basis |
|---|---|---|
| Fortran ↔ C++ ABI (Pitfalls 1–4) | HIGH | Grounded in direct reads of `xvue/xvuelc.c` lines 60–70 and sample entry points; gfortran trailing-underscore convention is documented and unambiguous on Linux x86_64 |
| QApplication lifecycle (Pitfall 5) | HIGH | Qt 6 documentation is explicit about the singleton rule |
| Event-loop strategy (Pitfalls 6–8) | MEDIUM | `processEvents` + nested `QEventLoop` is the documented library-mode pattern, but nested-modal re-entrancy needs empirical testing in a real mesher session |
| CMake / static library / OpenMP (Pitfalls 9–11) | HIGH | Standard CMake + Qt 6 pitfalls documented widely; Debian trixie PIE default is known |
| Image export (Pitfalls 12–13) | MEDIUM | Depends on trixie's `qt6-imageformats` bundle contents; PostScript verbatim-keeping is a design decision already validated by STACK.md |
| Color / font / coordinate drift (Pitfalls 14–16) | MEDIUM | Qt 6 and Xlib have well-documented differences but the exact magnitude depends on the current `xvuelc.c` behavior |
| HiDPI (Pitfall 17) | MEDIUM | Qt 6 HiDPI is mature but the Fortran-side assumption about pixel units needs audit |
| Build hygiene / regression (Pitfalls 18–20) | HIGH | CONCERNS.md and TESTING.md already document the no-CI / stale-object-file fragility; xvue-qt inherits those risks unchanged |

## Roadmap Implications

**Phase 0 is disproportionately important.** Pitfalls 1, 2, 4, 5, 6, 9, 10, 16, 18 are all caught or prevented during the initial CMake skeleton + one-primitive end-to-end test. The phase should include:
- ABI checklist (trailing underscore, pointer args, `MefistoPoint` shim)
- PIE/PIC, `-fno-openmp` for the Qt library, AUTOMOC set before find_package
- Y-up vs Y-down documentation
- Singleton `QApplication` with `call_once`
- Clean-build discipline baked into `bin/cbl_tout_qt`

**Phase 4 (event bridge) is the architectural keystone.** Pitfalls 6, 7, 8 all cluster here. Allocate extra time for empirical testing with real `testa/` cases. Be ready to iterate the `processEvents`/nested-`QEventLoop` design.

**Phase 6 (export) needs an early probe.** Pitfall 12 means the first task of Phase 6 must be running a 10-line test to check `QImageWriter::supportedImageFormats()` for GIF write support. Do not discover this mid-phase.

**Phase 7 (A/B validation) is where color/font/HiDPI drift (Pitfalls 14, 15, 17) surfaces.** Plan for 1–2 iterations of tuning. Do not treat A/B validation as "check and done" — treat it as "compare, fix, re-compare."

**Phase 8 (X11 retirement) must not move until the one-release-cycle A/B window has closed.** Pitfall 19 specifically warns against removing the ImageMagick dependency before auditing all call sites.

**Pitfall 20 (regression drift) is a process rule that applies to every phase.** It is the single biggest risk in a single-developer no-CI project. The 5-canonical-case end-of-phase rule must be enforced, not optional.

## Open Questions

- **Does `xvuelc.c` actually invert Y internally, or does it pass Fortran Y through as-is?** Needs a one-hour grep + read during Phase 0 to document the convention before the Qt bridge commits to either.
- **Does Debian trixie's `qt6-imageformats` package ship a GIF writer plugin?** Needs a 10-line test at Phase 6 kickoff.
- **How does `xvnbpixeltexte_` currently handle missing fonts (fallback behavior)?** Determines whether the Qt version must also carry a fallback or whether it can assume the bundled font is always available.
- **Are any graphics calls currently made from inside `!$OMP PARALLEL` regions in `_OMP` executables?** Needs a grep on `elas/`, `flui/`, `ther/` sources for `!$OMP` next to `CALL XV*`. Expected to be none, but must verify before making the "main-thread-only" rule normative.
