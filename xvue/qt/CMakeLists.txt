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
    src/xvue_qt_state.cpp
    src/xvue_qt_window.cpp
)

# Phase 3 D-01/D-02: bundle DejaVuSansMono.ttf into libxvueqt.a via Qt
# resource system. Static-lib case — see Pitfall 1 in 03-RESEARCH.md:
# Q_INIT_RESOURCE(xvue_fonts) MUST be called at runtime from
# XvueApp::ensure() (see Task 2) or addApplicationFont returns -1.
qt_add_resources(xvueqt xvue_fonts
    PREFIX "/xvue/qt"
    FILES fonts/DejaVuSansMono.ttf
)

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
target_compile_options(xvueqt PRIVATE -fno-openmp -Wall -Wextra)

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

# D-13: install target so the shell scripts can locate libxvueqt.a by a stable path
install(TARGETS xvueqt
    ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
)
