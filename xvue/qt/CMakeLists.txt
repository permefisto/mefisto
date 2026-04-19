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

# D-12: verify_abi — fail build if symbol count in libxvueqt.a drifts from the
# extern "C" declaration count in the public header. This is the whole-project
# ABI drift guard; it runs on every build.
# Expected count after Option A resolution of Planner Alert: 57.
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

# D-13: install target so the shell scripts can locate libxvueqt.a by a stable path
install(TARGETS xvueqt
    ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
)

# Phase 5 Wave 0 (05-01): optional test subdirectory. Gated by a flag so the
# normal build stays unchanged; bin/cbl_tout_qt turns it ON to build the
# xvue_qt_event_tests target as part of the full Qt build.
option(XVUE_QT_BUILD_TESTS "Build xvue_qt QTest binaries (xvue_qt_event_tests)" OFF)
if(XVUE_QT_BUILD_TESTS)
    add_subdirectory(tests)
endif()
