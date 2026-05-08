# xvue/qt/CMakeLists.txt
# Phase 0: Build skeleton for libxvueqt.a — Qt 6 replacement for xvuelc.o
# Owner: Alain Perronnet project / xvue-qt migration
# Per D-11: AUTOMOC before find_package, -fPIC, C++17, no -fopenmp.
# Per D-12: verify_abi target fails build on symbol-count drift.

cmake_minimum_required(VERSION 3.21)
project(xvueqt LANGUAGES CXX)

# D-11, Pitfall 9: MUST precede find_package(Qt6)
set(CMAKE_AUTOMOC ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # D-11, Pitfall 10
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets PrintSupport)

add_library(xvueqt STATIC
    src/xvue_qt_api.cpp
    src/xvue_qt_app.cpp
    src/xvue_qt_canvas.cpp
    src/xvue_qt_event.cpp
    src/xvue_qt_state.cpp
    src/xvue_qt_window.cpp
    # Phase 6.0 additions (scaffold: stubs; bodies fill in Plans 02-06):
    src/xvue_qt_menu_bridge.cpp
    src/xvue_qt_i18n.cpp
    src/xvue_qt_prefs.cpp
    src/xvue_qt_console_dock.cpp
    src/xvue_qt_error_batcher.cpp
    src/xvue_qt_about_dialog.cpp
    src/xvue_qt_preferences.cpp
    src/xvue_qt_recent_projects.cpp
    src/xvue_qt_shortcut_audit.cpp
    # Phase 7 Plan 02 (EXPORT-04): PsEmitter — verbatim port of xvuelc.c
    # xvpostscript_ state machine plus per-primitive emit helpers (Plan 03
    # fills bodies). Owned by XvueApp via XvueApp::psEmitter() singleton.
    src/xvue_qt_postscript.cpp
    # Phase 7 Plan 04 (EXPORT-02, EXPORT-05): XvueExport — PNG/JPEG/PDF
    # static helpers driven by File → Export → … menu actions. ABI stays
    # at 58 (no new extern "C" symbols). Plan 05 will grow the file with
    # GIF assembly + auto-snapshot beginAnimation/endAnimation.
    src/xvue_qt_export.cpp
)

# Phase 6.1 Plan 02 additions: bilingual LIMTCL menu-file parser and the
# strong-symbol override of registerMailActions_stub_ that wires ~31
# QActions + 4 toolbar entries per LEXICON-AUDIT-mail.md.
# Phase 6.2 Plan 02 addition: registerElasActions_stub_ strong body
# (alphabetical order: elas before mail) wiring ~27 QActions + 4
# toolbar entries per LEXICON-AUDIT-elas.md.
# Phase 6.3 Plan 02 addition: registerFluiActions_stub_ strong body
# (alphabetical order: elas/flui/mail) wiring all qaction=yes rows +
# 4 toolbar entries per LEXICON-AUDIT-flui.md.
# Phase 6.4 Plan 02 addition: registerTherActions_stub_ strong body
# (alphabetical order: elas/flui/mail/ther) wiring all 37 qaction=yes
# rows + 4 toolbar entries per LEXICON-AUDIT-ther.md.
# Phase 6.5 Plan 02 addition: registerNlseActions_stub_ strong body
# (alphabetical order: elas/flui/mail/nlse/ther — nlse goes BETWEEN
# mail and ther) wiring all 19 qaction=yes nlse-leaves + 4 toolbar
# entries per LEXICON-AUDIT-nlse.md. The per-module set is now COMPLETE.
target_sources(xvueqt PRIVATE
    src/xvue_qt_menu_file_parser.cpp
    src/xvue_qt_elas_actions.cpp
    src/xvue_qt_flui_actions.cpp
    src/xvue_qt_mail_actions.cpp
    src/xvue_qt_nlse_actions.cpp
    src/xvue_qt_ther_actions.cpp
)

# Phase 3 D-01/D-02: bundle DejaVuSansMono.ttf into libxvueqt.a via Qt
# resource system. Static-lib case — see Pitfall 1 in 03-RESEARCH.md:
# Q_INIT_RESOURCE(xvue_fonts) MUST be called at runtime from
# XvueApp::ensure() (see Task 2) or addApplicationFont returns -1.
qt_add_resources(xvueqt xvue_fonts
    PREFIX "/xvue/qt"
    FILES fonts/DejaVuSansMono.ttf
)

# Phase 6.0: empty icon resource — Plans 6.1..6.5 append SVG entries via the
# qrc XML. AUTORCC picks the .qrc up because we add it to target_sources.
# The empty-FILES form of qt_add_resources is unreliable across Qt versions;
# the .qrc-file route is the Qt-idiomatic alternative.
set(XVUE_ICONS_QRC ${CMAKE_CURRENT_SOURCE_DIR}/resources/xvue_icons.qrc)
target_sources(xvueqt PRIVATE ${XVUE_ICONS_QRC})
set_target_properties(xvueqt PROPERTIES AUTORCC ON)

target_include_directories(xvueqt PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(xvueqt
    PUBLIC
        Qt6::Core
        Qt6::Gui
        Qt6::Widgets
        Qt6::PrintSupport
)

# D-11, Pitfall 11: graphics layer MUST NOT propagate -fopenmp
# Phase 6.0 Pitfall 2.4: -finput-charset=UTF-8 so accented FR i18n strings
# embedded in Plan 02's tr() table compile cleanly across host locales.
target_compile_options(xvueqt PRIVATE -fno-openmp -Wall -Wextra -finput-charset=UTF-8)

# 2026-05-08: silence false-positive -Wdangling-reference from gcc 13/14 on
# `const MenuFile& foo = XvueMenuFileParser::loadFor(...)` bindings in the
# *_actions.cpp files. loadFor returns a reference into the parser's static
# cache; gcc's heuristic does not see the cache and assumes a temporary.
# Only apply if the compiler accepts the flag (older gcc / clang lack it).
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-dangling-reference" XVUE_HAS_NO_DANGLING_REF)
if(XVUE_HAS_NO_DANGLING_REF)
    target_compile_options(xvueqt PRIVATE -Wno-dangling-reference)
endif()

# D-12: verify_abi — fail build if symbol count in libxvueqt.a drifts from the
# extern "C" declaration count in the public header. This is the whole-project
# ABI drift guard; it runs on every build.
# Expected count after Option A resolution of Planner Alert: 58.
add_custom_target(verify_abi ALL
    COMMAND ${CMAKE_COMMAND} -E echo "verify_abi: counting ABI symbols in libxvueqt.a..."
    COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake/verify_abi.sh
        $<TARGET_FILE:xvueqt>
        ${CMAKE_CURRENT_SOURCE_DIR}/include/xvue_qt_api.h
    DEPENDS xvueqt
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
)

# D-10: SHELL-03 enforcement — fail the build if QApplication::exec or
# qApp->exec() appears in xvue/qt/src or xvue/qt/include. Mirrors the Phase 0
# verify_abi pattern: a POST_BUILD grep that aborts if it finds a match.
add_custom_target(verify_no_exec ALL
    COMMAND ${CMAKE_COMMAND} -E echo "verify_no_exec: scanning for QApplication::exec..."
    COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake/verify_no_exec.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    DEPENDS xvueqt
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
)

# Phase 6.0 §10 layer 1: D-04 modifier-rule lint. Fails the build if any
# bareword-character QKeySequence appears under xvue/qt/src.
add_custom_target(verify_shortcut_modifiers ALL
    COMMAND ${CMAKE_COMMAND} -E echo "verify_shortcut_modifiers: scanning..."
    COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake/verify_shortcut_modifiers.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    DEPENDS xvueqt
    VERBATIM
)

# Phase 6.0: SVG-only enforcement for custom icons (D-10 + UI-SPEC §Iconography).
# Plans 6.1..6.5 will populate resources/icons/ with .svg entries — this guard
# fails the build if a raster file ever sneaks in.
add_custom_target(verify_icon_source ALL
    COMMAND ${CMAKE_COMMAND} -E echo "verify_icon_source: scanning..."
    COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake/verify_icon_source.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/resources/icons
    DEPENDS xvueqt
    VERBATIM
)

# Phase 7 Plan 06 (EXPORT-06): grep gate — fail the build if any
# legacy raster-tool reference appears under xvue/qt/. Per CONTEXT.md
# D-16/D-17 the gate is scoped to xvue/qt/ ONLY (NOT xvue/) — the
# LVIDEO pipeline (xvue/video1.f, xvue/videofin.f, xvue/videonm.f) and
# the legacy GIF post-processor under bin/ legitimately stay alive in
# the legacy tree until Phase 9 RETIRE-03.
#
# The script lives in $MEFISTO/bin/ rather than under xvue/qt/cmake/ so
# it can also be invoked standalone (CI hook, manual sanity check) by
# the same path used by ctest. ${CMAKE_CURRENT_SOURCE_DIR}/../../
# resolves the repo root from xvue/qt/.
add_custom_target(verify_no_imagemagick_in_qt ALL
    COMMAND ${CMAKE_COMMAND} -E echo "verify_no_imagemagick_in_qt: scanning..."
    COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/../../bin/test_no_imagemagick_in_qt.sh
    DEPENDS xvueqt
    VERBATIM
)

# Phase 9 Plan 9-03 (RETIRE-02): post-retirement gate. Mirrors
# verify_no_imagemagick_in_qt structure. Fails the build if any X11
# linker line or path resurfaces under bin/cb* + xvue/qt/.
add_custom_target(verify_no_x11_in_build ALL
    COMMAND ${CMAKE_COMMAND} -E echo "verify_no_x11_in_build: scanning..."
    COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/../../bin/test_no_x11_in_build.sh
    DEPENDS xvueqt
    VERBATIM
)

# Phase 9 Plan 9-04 (RETIRE-03): post-retirement gate. Mirrors verify_no_*.
# Delegates the actual scan to bin/test_no_lvideo.sh, which fails the
# build if any legacy LVIDEO callback or shell-out helper resurfaces in
# active Fortran sources (xvue/, flui/, ther/, util/, prpr/).
add_custom_target(verify_no_lvideo ALL
    COMMAND ${CMAKE_COMMAND} -E echo "verify_no_lvideo: scanning..."
    COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/../../bin/test_no_lvideo.sh
    DEPENDS xvueqt
    VERBATIM
)

# Phase 9 Plan 9-09 (Phase 8 D-09 carry-forward): optional pp/*
# freshness check. Per RESEARCH §Open Question 4: NOT in the ALL chain
# (would race build-order between libxvueqt.a and pp/* link). Primary
# enforcement is the end-section invocation in bin/cbl_tout. This CMake
# target is opt-in for developers who prefer build-time fast feedback —
# invoke explicitly via `cmake --build xvue/qt/build --target verify_pp_freshness`.
option(MEFISTO_VERIFY_PP_FRESHNESS_ALL
       "Add verify_pp_freshness to ALL target (off by default per RESEARCH §OQ4)"
       OFF)

add_custom_target(verify_pp_freshness
    COMMAND ${CMAKE_COMMAND} -E echo "verify_pp_freshness: scanning..."
    COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake/verify_pp_freshness.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/build/libxvueqt.a
        ${CMAKE_CURRENT_SOURCE_DIR}/../../pp
    DEPENDS xvueqt
    VERBATIM
)

if(MEFISTO_VERIFY_PP_FRESHNESS_ALL)
    # Force inclusion in default-build target. Off-by-default per OQ4.
    set_property(TARGET verify_pp_freshness PROPERTY EXCLUDE_FROM_ALL FALSE)
    # And add a depends-on chain so default ALL invokes it:
    add_custom_target(_verify_pp_freshness_all_hook ALL
        DEPENDS verify_pp_freshness
    )
else()
    set_property(TARGET verify_pp_freshness PROPERTY EXCLUDE_FROM_ALL TRUE)
endif()

# Phase 7 Plan 01 (EXPORT-01): probes/ subdirectory holds standalone probe
# binaries (e.g., qimagewriter_probe). Gated so the canonical bin/cbl_tout
# build (Qt-only since Phase 9 RETIRE-02 — formerly the parallel cbl_tout
# vs cbl_tout_qt pair) does not depend on it; bin/cb_probe_qt enables the
# option explicitly.
option(XVUE_QT_BUILD_PROBES "Build xvue/qt/probes/* standalone binaries" OFF)
if(XVUE_QT_BUILD_PROBES)
    add_subdirectory(probes)
endif()

# D-13: install target so the shell scripts can locate libxvueqt.a by a stable path
install(TARGETS xvueqt
    ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
)

# Phase 5 Wave 0 (05-01): optional test subdirectory. Gated by a flag so the
# normal build stays unchanged; bin/cbl_tout turns it ON to build the
# xvue_qt_event_tests target as part of the full Qt build.
# (Pre-Phase-9 this comment named bin/cbl_tout_qt — RETIRE-02 collapsed
# that to bin/cbl_tout.)
option(XVUE_QT_BUILD_TESTS "Build all xvue_qt QTest binaries" OFF)
if(XVUE_QT_BUILD_TESTS)
    # WR-01/WR-02 gap-closure: enable test-only API in the library itself
    # (clearCacheForTesting, xvueClearLanguageCacheForTesting) so that test
    # binaries can link against the symbols they declare via XVUE_QT_TESTING
    # guards in xvue_qt_menu_file_parser.h and xvue_qt_i18n.h.
    # The define is scoped to test builds only (XVUE_QT_BUILD_TESTS=ON);
    # production builds of libxvueqt.a never include this define.
    target_compile_definitions(xvueqt PRIVATE XVUE_QT_TESTING)
    # Phase 7 Plan 02 (EXPORT-04): enable_testing() must live in the parent
    # CMakeLists.txt so `ctest -R xvue_qt_postscript_tests` from the build
    # root finds tests registered via add_test() inside the tests/
    # subdirectory. Calling enable_testing() only in tests/CMakeLists.txt
    # registers tests in tests/CTestTestfile.cmake which the build-root
    # CTest config never references — VALIDATION.md row 07-02-01's
    # `ctest -R '^xvue_qt_postscript_tests$'` from xvue/qt/build/ would
    # then report "No tests were found".
    enable_testing()
    add_subdirectory(tests)
endif()
