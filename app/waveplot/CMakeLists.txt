set(projectdir ${PROJECT_SOURCE_DIR})

set(sources-lib-fpp
  gridcache.F90
  initwaveplot.F90
  molorb.F90
  molorb2.F90
  slater.F90
  )
set(sources-main-fpp
    waveplot.F90
)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/src/dftbp/include)

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources-lib-fpp}" sources-lib-f90-preproc)
dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources-main-fpp}" sources-main-f90-preproc)

add_library(waveplot-lib STATIC ${sources-lib-f90-preproc})
# Library targets public such that we can link against this from elsewhere
# Unveil the module and source directories as well
target_include_directories(waveplot-lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_link_libraries(waveplot-lib PUBLIC dftbplus)
if(WITH_MPI)
  target_link_libraries(waveplot-lib PUBLIC MPI::MPI_Fortran)
endif()

# Executable links against library
add_executable(waveplot ${sources-main-f90-preproc})
target_link_libraries(waveplot PRIVATE waveplot-lib)

install(TARGETS waveplot DESTINATION "${CMAKE_INSTALL_BINDIR}")


