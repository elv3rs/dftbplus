set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/commontypes.F90
  ${curdir}/densedescr.F90
  ${curdir}/dynneighlist.F90
  ${curdir}/latpointiter.F90
  ${curdir}/eleccutoffs.F90
  ${curdir}/integral.F90
  ${curdir}/multipole.F90
  ${curdir}/oldskdata.F90
  ${curdir}/orbitals.F90
  ${curdir}/parallelks.F90
  ${curdir}/typegeometry.F90
  ${curdir}/typegeometryhsd.F90
  ${curdir}/wrappedintr.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
