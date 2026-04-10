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
    COMMAND sh -c "nm $<TARGET_FILE:xvueqt> | grep -c ' T [a-zA-Z_][a-zA-Z0-9_]*_$$' > nm_count.txt"
    COMMAND sh -c "grep -c '^[[:space:]]*\\(void\\|int\\|float\\|double\\|long\\|short\\|unsigned\\|void[[:space:]]*\\*\\)[[:space:]]*proc(' ${CMAKE_CURRENT_SOURCE_DIR}/include/xvue_qt_api.h > hdr_count.txt"
    COMMAND sh -c "NM_COUNT=$$(cat nm_count.txt); HDR_COUNT=$$(cat hdr_count.txt); echo \"nm count: $$NM_COUNT  header count: $$HDR_COUNT\"; if [ \"$$NM_COUNT\" != \"$$HDR_COUNT\" ]; then echo 'ERROR: ABI symbol count drift between libxvueqt.a and xvue_qt_api.h' >&2; exit 1; fi; if [ \"$$NM_COUNT\" != \"57\" ]; then echo \"WARNING: expected 57 Fortran-facing entries (Option A of Planner Alert), got $$NM_COUNT\" >&2; fi"
    DEPENDS xvueqt
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
)

# D-13: install target so the shell scripts can locate libxvueqt.a by a stable path
install(TARGETS xvueqt
    ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
)
