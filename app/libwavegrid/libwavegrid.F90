module libwavegrid
  use libwavegrid_molorb, only: TSpeciesBasis, TMolecularOrbital, TMolecularOrbital_init
  use libwavegrid_molorb, only: getValue, getAtomicDensities, getTotalChrg
  use libwavegrid_slater, only: TSlaterOrbital
  implicit none

  public :: TSpeciesBasis, TMolecularOrbital, TSlaterOrbital
  public :: TMolecularOrbital_init
  public :: getValue, getAtomicDensities, getTotalChrg

end module libwavegrid

