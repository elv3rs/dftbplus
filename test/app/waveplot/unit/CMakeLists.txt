set(fypp_flags ${FYPP_BUILD_FLAGS})
get_target_property(_fortuno_incdir Fortuno::fortuno_include_dir INTERFACE_INCLUDE_DIRECTORIES)
list(APPEND fypp_flags "-I${_fortuno_incdir}" "-I${CMAKE_SOURCE_DIR}/src/dftbp/include")

set(unit-test-prefix "unit")

add_subdirectory(molorb)

