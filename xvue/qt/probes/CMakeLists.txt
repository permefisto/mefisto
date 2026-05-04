# xvue/qt/probes/CMakeLists.txt
# Phase 7 Plan 01 (EXPORT-01): standalone probe binary. NOT a library
# dependency; built only on demand by bin/cb_probe_qt.

cmake_minimum_required(VERSION 3.21)
find_package(Qt6 REQUIRED COMPONENTS Core Gui)

add_executable(qimagewriter_probe qimagewriter_probe.cpp)

target_link_libraries(qimagewriter_probe
    PRIVATE
        Qt6::Core
        Qt6::Gui
)

target_compile_options(qimagewriter_probe PRIVATE -Wall -Wextra)
