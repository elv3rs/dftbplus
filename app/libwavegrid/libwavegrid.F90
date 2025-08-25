module libwavegrid
  use libwavegrid_molorb, only: TSpeciesBasis, TMolecularOrbital, TMolecularOrbital_init
  use libwavegrid_molorb, only: getValue, getAtomicDensities, getTotalChrg
  implicit none

  public :: TSpeciesBasis, TMolecularOrbital, TMolecularOrbital_init
  public :: getValue, getAtomicDensities, getTotalChrg

end module libwavegrid

