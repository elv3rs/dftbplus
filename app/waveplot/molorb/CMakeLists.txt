set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources_local_fpp
  ${curdir}/molorb.F90
  ${curdir}/pointwise.F90
  ${curdir}/parallel.F90
)

set(sources_lib_fpp "${sources_lib_fpp};${sources_local_fpp}" PARENT_SCOPE)

