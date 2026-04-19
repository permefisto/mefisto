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

# No install(); tests run in-place from xvue/build/tests/.
