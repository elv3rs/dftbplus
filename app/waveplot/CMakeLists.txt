set(projectdir ${PROJECT_SOURCE_DIR})

set(sources_main_fpp
  initwaveplot.F90
    waveplot.F90
)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/src/dftbp/include)

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources_main_fpp}" sources_main_f90_preproc)

# Link waveplot executable against libwavegrid
add_executable(waveplot ${sources_main_f90_preproc})
target_link_libraries(waveplot PRIVATE libwavegrid)

install(TARGETS waveplot DESTINATION "${CMAKE_INSTALL_BINDIR}")
