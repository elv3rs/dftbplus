set(projectdir ${PROJECT_SOURCE_DIR})

set(sources_lib_fpp
  gridcache.F90
  slater.F90
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
    find_package(CUDAToolkit REQUIRED)
    
    # Enable CXX language to ensure OpenMP for C++ is handled correctly.
    enable_language(CXX)
    
    # 1. Create a dedicated library for the CUDA code.
    add_library(wavegrid_cuda STATIC molorb/kernel.cu)

    # 2. THE FIX: Directly add the required compile flag to the CUDA target.
    #    This is the most reliable way to inject the flag.
    target_compile_options(wavegrid_cuda PRIVATE "-Xcompiler=-fopenmp")
    
    # 3. Find OpenMP and link against it. This is still necessary to tell
    #    the FINAL LINKER to add the OpenMP runtime library (e.g., -lgomp).
    find_package(OpenMP REQUIRED)
    target_link_libraries(wavegrid_cuda PUBLIC OpenMP::OpenMP_CXX)

    # Link against other CUDA libs as before.
    target_link_libraries(wavegrid_cuda PUBLIC CUDA::cublas)

    # Finally, link the main libwavegrid library against our new CUDA library.
    # This pulls in the compiled kernel.cu.o and its link dependencies.
    target_link_libraries(libwavegrid PUBLIC wavegrid_cuda)
endif()






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






