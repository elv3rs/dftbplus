set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/blockpothelper.F90
  ${curdir}/bondpopulations.F90
  ${curdir}/boundarycond.F90
  ${curdir}/chargepenalty.F90
  ${curdir}/charges.F90
  ${curdir}/coordnumber.F90
  ${curdir}/coulomb.F90
  ${curdir}/dense.F90
  ${curdir}/densitymatrix.F90
  ${curdir}/determinants.F90
  ${curdir}/dftbplusu.F90
  ${curdir}/dftd3param.F90
  ${curdir}/dftd4param.F90
  ${curdir}/dispcommon.F90
  ${curdir}/dispdftd4.F90
  ${curdir}/dispersions.F90
  ${curdir}/dispiface.F90
  ${curdir}/dispslaterkirkw.F90
  ${curdir}/dispuff.F90
  ${curdir}/dispuffdata.F90
  ${curdir}/elecconstraints.F90
  ${curdir}/elstatpot.F90
  ${curdir}/energytypes.F90
  ${curdir}/encharges.F90
  ${curdir}/etemp.F90
  ${curdir}/extcharges.F90
  ${curdir}/extfields.F90
  ${curdir}/forces.F90
  ${curdir}/getenergies.F90
  ${curdir}/h5correction.F90
  ${curdir}/halogenx.F90
  ${curdir}/hamiltonian.F90
  ${curdir}/hybridxc.F90
  ${curdir}/mdftb.F90
  ${curdir}/nonscc.F90
  ${curdir}/onscorrection.F90
  ${curdir}/orbitalequiv.F90
  ${curdir}/periodic.F90
  ${curdir}/pmlocalisation.F90
  ${curdir}/populations.F90
  ${curdir}/potentials.F90
  ${curdir}/rshgamma.F90
  ${curdir}/scc.F90
  ${curdir}/sccinit.F90
  ${curdir}/shift.F90
  ${curdir}/shortgamma.F90
  ${curdir}/shortgammafuncs.F90
  ${curdir}/simpledftd3.F90
  ${curdir}/sk.F90
  ${curdir}/slakocont.F90
  ${curdir}/slakoeqgrid.F90
  ${curdir}/sparse2dense.F90
  ${curdir}/spin.F90
  ${curdir}/spinorbit.F90
  ${curdir}/stress.F90
  ${curdir}/thirdorder.F90
  ${curdir}/uniquehubbard.F90
  ${curdir}/vdwdata.F90)

if(WITH_MBD)
  list(APPEND sources-fpp ${curdir}/dispmbd.F90)
endif(WITH_MBD)

add_subdirectory(repulsive)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
