set(projectdir ${PROJECT_SOURCE_DIR})

set(sources_lib_fpp
  gridcache.F90
  slater.F90
)

add_subdirectory(molorb)
if(WITH_CUDA)
    set(sources_lib_cuda
    molorb/kernel.cu
    )
endif()

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/src/dftbp/include)

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources_lib_fpp}" sources_lib_f90_preproc)

add_library(libwavegrid STATIC
  ${sources_lib_f90_preproc}
  ${sources_lib_cuda}
)

# Library targets public such that we can link against this from elsewhere
# Unveil the module and source directories as well
target_include_directories(libwavegrid PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/molorb>
)

target_link_libraries(libwavegrid PUBLIC dftbplus)
if(WITH_MPI)
  target_link_libraries(libwavegrid PUBLIC MPI::MPI_Fortran)
endif()






