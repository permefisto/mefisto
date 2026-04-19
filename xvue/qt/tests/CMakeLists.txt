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

# Phase 6.0 Plan 05 Task 2: 10 QTest cases covering wheel zoom, middle-drag
# pan, context menu suppression, resetView, and empty-state paint branch.
# Note: this target is provisionally added here in the Plan 05 worktree to
# enable verification under the orchestrator's parallel-worktree run. Plan
# 02's worktree is the canonical owner of tests/CMakeLists.txt edits in
# Wave 2 — the orchestrator merge will resolve any conflict (the Plan 02
# version already includes this same target per the Plan 05 frontmatter
# note "tests/CMakeLists.txt additions for this plan's test binary are
# owned by Plan 02 to avoid Wave 2 parallel-worktree merge conflicts").
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
