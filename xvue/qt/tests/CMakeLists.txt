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

# No install(); tests run in-place from xvue/build/tests/.
