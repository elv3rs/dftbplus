set(projectdir ${PROJECT_SOURCE_DIR})

set(sources_lib_fpp
  gridcache.F90
  initwaveplot.F90
  slater.F90
)
set(sources_main_fpp
    waveplot.F90
)
add_subdirectory(molorb)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/src/dftbp/include)

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources_lib_fpp}" sources_lib_f90_preproc)
dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources_main_fpp}" sources_main_f90_preproc)

add_library(waveplot_lib STATIC ${sources_lib_f90_preproc})

# Library targets public such that we can link against this from elsewhere
# Unveil the module and source directories as well
target_include_directories(waveplot_lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/molorb>
)

target_link_libraries(waveplot_lib PUBLIC dftbplus)
if(WITH_MPI)
  target_link_libraries(waveplot_lib PUBLIC MPI::MPI_Fortran)
endif()

# Executable links against library
add_executable(waveplot ${sources_main_f90_preproc})
target_link_libraries(waveplot PRIVATE waveplot_lib)

install(TARGETS waveplot DESTINATION "${CMAKE_INSTALL_BINDIR}")


