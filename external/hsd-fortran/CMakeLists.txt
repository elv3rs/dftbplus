# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# Disable hsd-fortran tests, examples, and benchmarks when built as a subproject
set(HSD_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(HSD_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(HSD_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)

add_subdirectory(${HSD_SOURCE_DIR} ${HSD_BINARY_DIR})
