# xvue/qt/tests/CMakeLists.txt
# Phase 5 Wave 0 (05-01 Task 1): QTest unit-test target for the xvue_qt event
# bridge. Gated by -DXVUE_QT_BUILD_TESTS=ON from the parent CMakeLists.
#
# Target: xvue_qt_event_tests (runs under xvfb-run, links Qt6::Test).
find_package(Qt6 REQUIRED COMPONENTS Test Widgets)

add_executable(xvue_qt_event_tests
    test_xvue_qt_event.cpp
    test_helpers.cpp
)

target_link_libraries(xvue_qt_event_tests
    PRIVATE
        xvueqt
        Qt6::Test
        Qt6::Widgets
        Qt6::Gui
        Qt6::Core
)

target_include_directories(xvue_qt_event_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# Phase 6.0 Plan 01 Task 1: scaffold sanity tests for the 9 new src/xvue_qt_*
# classes + the xvue_module_init_ ABI hook (entry 58). Mirrors the
# xvue_qt_event_tests setup so xvfb-run can exec it identically.
add_executable(xvue_qt_menu_tests
    test_xvue_qt_menu_scaffold.cpp
)

target_link_libraries(xvue_qt_menu_tests
    PRIVATE
        xvueqt
        Qt6::Test
        Qt6::Widgets
        Qt6::Gui
        Qt6::Core
)

target_include_directories(xvue_qt_menu_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# ---------------------------------------------------------------------------
# Phase 6.0 Wave 2: i18n / menu-bridge / prefs / dialogs / canvas gesture
# test binaries. Plan 02 owns ALL three add_executable() blocks below to
# avoid parallel-worktree merge conflicts on this CMakeLists.txt — Plans
# 04 and 05 only write their .cpp files, never edit this file.
# ---------------------------------------------------------------------------

# ---- Plan 02 test binary ----
add_executable(xvue_qt_i18n_menu_prefs_tests
    test_xvue_qt_i18n_menu_prefs.cpp
    test_helpers.cpp
)
target_link_libraries(xvue_qt_i18n_menu_prefs_tests
    PRIVATE
        xvueqt
        Qt6::Test
        Qt6::Widgets
        Qt6::Gui
        Qt6::Core
)
target_include_directories(xvue_qt_i18n_menu_prefs_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# ---- Plan 04 test binary (console dock + error batcher + dialogs +
#                          recent-projects).
add_executable(xvue_qt_dialogs_console_tests
    test_xvue_qt_dialogs_console.cpp
    test_helpers.cpp
)
target_link_libraries(xvue_qt_dialogs_console_tests
    PRIVATE
        xvueqt
        Qt6::Test
        Qt6::Widgets
        Qt6::Gui
        Qt6::Core
)
target_include_directories(xvue_qt_dialogs_console_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# ---- Plan 05 test binary (canvas gestures + view transform + empty-state).
add_executable(xvue_qt_canvas_gestures_tests
    test_xvue_qt_canvas_gestures.cpp
    test_helpers.cpp
)
target_link_libraries(xvue_qt_canvas_gestures_tests
    PRIVATE
        xvueqt
        Qt6::Test
        Qt6::Widgets
        Qt6::Gui
        Qt6::Core
)
target_include_directories(xvue_qt_canvas_gestures_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# ---------------------------------------------------------------------------
# Phase 6.0 Wave 3 (Plan 03): event-bridge menu-queue pre-drain tests.
# Distinct from xvue_qt_event_tests so a Plan 03 regression cannot
# masquerade as a Phase 5 regression — they fail (or skip) independently.
# ---------------------------------------------------------------------------
add_executable(xvue_qt_event_menu_predrain_tests
    test_xvue_qt_event_menu_predrain.cpp
    test_helpers.cpp
)
target_link_libraries(xvue_qt_event_menu_predrain_tests
    PRIVATE
        xvueqt
        Qt6::Test
        Qt6::Widgets
        Qt6::Gui
        Qt6::Core
)
target_include_directories(xvue_qt_event_menu_predrain_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# ---------------------------------------------------------------------------
# Phase 6.0 Wave 4 (Plan 06): integration tests for the composed shell —
# XvueApp menuBridge() + applyColorSchemePreference(), XvueWindow chrome
# (menu bar + toolbar + status bar + console dock + 13 shared QActions),
# and the xvue_module_init_ dispatch body.
# ---------------------------------------------------------------------------
add_executable(xvue_qt_window_chrome_tests
    test_xvue_qt_window_chrome.cpp
    test_helpers.cpp
)
target_link_libraries(xvue_qt_window_chrome_tests
    PRIVATE
        xvueqt
        Qt6::Test
        Qt6::Widgets
        Qt6::Gui
        Qt6::Core
)
target_include_directories(xvue_qt_window_chrome_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# ---------------------------------------------------------------------------
# Phase 6.1 Plan 03 (D-13): mesher menu round-trip tests.
# 4 required cases + 1 bonus D-09 closeEvent case.
# Pattern: QAction::trigger() (NOT postEvent) drives the production bridge;
# xvsouris_ drains the chars via the Plan 03 pre-drain path.
# ---------------------------------------------------------------------------
add_executable(xvue_qt_mail_menu_tests
    test_xvue_qt_mail_menu.cpp
    test_helpers.cpp
)
target_link_libraries(xvue_qt_mail_menu_tests
    PRIVATE
        xvueqt
        Qt6::Test
        Qt6::Widgets
        Qt6::Gui
        Qt6::Core
)
target_include_directories(xvue_qt_mail_menu_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)
# Enable test-only API in XvueMenuFileParser and xvue_qt_i18n.
target_compile_definitions(xvue_qt_mail_menu_tests PRIVATE XVUE_QT_TESTING)

# ---------------------------------------------------------------------------
# Phase 6.2 Plan 03 (D-13 elas port): elasticity menu round-trip tests.
# 4 required cases + 1 bonus D-09 closeEvent case. Mirrors the 6.1
# xvue_qt_mail_menu_tests target verbatim, substituting elas for mail.
# QAction::trigger() drives the production bridge; xvsouris_ drains the
# chars via the shared Plan 03 pre-drain path.
# ---------------------------------------------------------------------------
add_executable(xvue_qt_elas_menu_tests
    test_xvue_qt_elas_menu.cpp
    test_helpers.cpp
)
target_link_libraries(xvue_qt_elas_menu_tests
    PRIVATE
        xvueqt
        Qt6::Test
        Qt6::Widgets
        Qt6::Gui
        Qt6::Core
)
target_include_directories(xvue_qt_elas_menu_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)
# Enable test-only API in XvueMenuFileParser and xvue_qt_i18n.
target_compile_definitions(xvue_qt_elas_menu_tests PRIVATE XVUE_QT_TESTING)

# ---------------------------------------------------------------------------
# Phase 6.3 Plan 03 (D-13 flui port): fluider menu round-trip tests.
# 4 required cases + 1 bonus closeEvent + testMenuOrder + testBilingualLabelsEnglish.
# Mirrors the 6.2 xvue_qt_elas_menu_tests target verbatim, substituting flui
# for elas. testMenuOrder codifies the 6.2 Plan-04 {File, Fluid, View, Help}
# gate as an automated regression guard. testBilingualLabelsEnglish codifies
# the 6.2 Plan-05 EN-parser gate. Together these two slots replace the 6.2
# manual A/B human-verify checkpoint, making Phase 6.3 fully autonomous.
# ---------------------------------------------------------------------------
add_executable(xvue_qt_flui_menu_tests
    test_xvue_qt_flui_menu.cpp
    test_helpers.cpp
)
target_link_libraries(xvue_qt_flui_menu_tests
    PRIVATE
        xvueqt
        Qt6::Test
        Qt6::Widgets
        Qt6::Gui
        Qt6::Core
)
target_include_directories(xvue_qt_flui_menu_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)
# Enable test-only API in XvueMenuFileParser (clearCacheForTesting)
# and xvue_qt_i18n (xvueClearLanguageCacheForTesting).
# Required for testBilingualLabelsEnglish. Mirrors elas target.
target_compile_definitions(xvue_qt_flui_menu_tests PRIVATE XVUE_QT_TESTING)

# ---------------------------------------------------------------------------
# Phase 6.4 Plan 03 (D-13 ther port): thermal menu round-trip tests.
# 4 required cases + 1 bonus closeEvent + testMenuOrder + testBilingualLabelsEnglish.
# Mirrors the 6.3 xvue_qt_flui_menu_tests target verbatim, substituting ther
# for flui. testMenuOrder codifies the {File, Thermal, View, Help} gate.
# testBilingualLabelsEnglish codifies the EN-parser gate against td/ma/debuther.
# testHelpNoQueue allowlist: {98;} per LEXICON-AUDIT-ther row 98 (NOT {97;} --
# that's flui's value; codifies the 6.3 Plan 03 Auto-fix Rule 1 lesson).
# ---------------------------------------------------------------------------
add_executable(xvue_qt_ther_menu_tests
    test_xvue_qt_ther_menu.cpp
    test_helpers.cpp
)
target_link_libraries(xvue_qt_ther_menu_tests
    PRIVATE
        xvueqt
        Qt6::Test
        Qt6::Widgets
        Qt6::Gui
        Qt6::Core
)
target_include_directories(xvue_qt_ther_menu_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)
# Enable test-only API in XvueMenuFileParser (clearCacheForTesting)
# and xvue_qt_i18n (xvueClearLanguageCacheForTesting).
# Required for testBilingualLabelsEnglish. Mirrors elas/flui targets.
target_compile_definitions(xvue_qt_ther_menu_tests PRIVATE XVUE_QT_TESTING)

# No install(); tests run in-place from xvue/build/tests/.
