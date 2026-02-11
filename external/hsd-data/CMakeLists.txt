# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

set(HSD_DATA_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(HSD_DATA_BUILD_APP OFF CACHE BOOL "" FORCE)

add_subdirectory(${HSD-DATA_SOURCE_DIR} ${HSD-DATA_BINARY_DIR})

add_library(hsd-data::hsd-data INTERFACE IMPORTED GLOBAL)
target_link_libraries(hsd-data::hsd-data INTERFACE hsd-data)
