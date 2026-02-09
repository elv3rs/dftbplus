set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/apicallback.F90
  ${curdir}/eigenvects.F90
  ${curdir}/elstattypes.F90
  ${curdir}/forcetypes.F90
  ${curdir}/hsdhelpers.F90
  ${curdir}/initprogram.F90
  ${curdir}/inputconversion.F90
  ${curdir}/inputdata.F90
  ${curdir}/main.F90
  ${curdir}/mainio.F90
  ${curdir}/outputfiles.F90
  ${curdir}/oldcompat.F90
  ${curdir}/parser.F90
  ${curdir}/parser_hamiltonian.F90
  ${curdir}/parser_customisation.F90
  ${curdir}/parser_analysis.F90
  ${curdir}/parser_electrostatics.F90
  ${curdir}/parser_dispersion.F90
  ${curdir}/parser_driver.F90
  ${curdir}/parser_excited.F90
  ${curdir}/parser_external.F90
  ${curdir}/parser_filling.F90
  ${curdir}/parser_general.F90
  ${curdir}/parser_geometry.F90
  ${curdir}/parser_hybrid.F90
  ${curdir}/parser_kpoints.F90
  ${curdir}/parser_parallel.F90
  ${curdir}/parser_parseroptions.F90
  ${curdir}/parser_reks.F90
  ${curdir}/parser_sccoptions.F90
  ${curdir}/parser_shared_utils.F90
  ${curdir}/parser_spin.F90
  ${curdir}/parser_skfiles.F90
  ${curdir}/parser_solver.F90
  ${curdir}/parser_transport.F90
  ${curdir}/qdepextpotgen.F90
  ${curdir}/qdepextpotproxy.F90
  ${curdir}/specieslist.F90
  ${curdir}/transportio.F90
)

if(WITH_API)
  list(APPEND sources-fpp ${curdir}/mainapi.F90)
endif()

add_subdirectory(input)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
