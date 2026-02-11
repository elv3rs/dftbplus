# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

set(HSD_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(HSD_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(HSD_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)

add_subdirectory(${HSD_SOURCE_DIR} ${HSD_BINARY_DIR})

add_library(hsd::hsd INTERFACE IMPORTED GLOBAL)
target_link_libraries(hsd::hsd INTERFACE hsd)
