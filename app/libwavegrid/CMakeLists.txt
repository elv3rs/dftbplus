set(projectdir ${PROJECT_SOURCE_DIR})

set(sources_lib_fpp
  slater.F90
  libwavegrid.F90
)

add_subdirectory(molorb)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/src/dftbp/include)

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources_lib_fpp}" sources_lib_f90_preproc)

add_library(libwavegrid STATIC
  ${sources_lib_f90_preproc}
  ${sources_lib_cuda}
)


if(WITH_CUDA)
    add_library(wavegrid_cuda STATIC molorb/kernel.cu)

    if(WITH_OMP)
        target_compile_options(wavegrid_cuda PRIVATE "-Xcompiler=-fopenmp")
        target_link_libraries(wavegrid_cuda PUBLIC OpenMP::OpenMP_CXX)
    endif()

    target_link_libraries(libwavegrid PUBLIC wavegrid_cuda)
endif()


# Library targets public such that we can link against this from elsewhere
# Unveil the module and source directories as well
target_include_directories(libwavegrid PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/molorb>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/molorb>
)

target_link_libraries(libwavegrid PUBLIC dftbplus)

if(WITH_MPI) 
  target_link_libraries(libwavegrid PUBLIC MPI::MPI_Fortran)
endif()






