!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module dftbp_wavegrid_basis
  use dftbp_common_accuracy, only : dp
  use dftbp_wavegrid_basis_spharmonics, only: realTessY
  use dftbp_wavegrid_basis_orbital, only: TOrbital, maxCutoff
  use dftbp_wavegrid_basis_slater, only: TSlaterOrbital, TSlaterOrbital_init
  use dftbp_wavegrid_basis_lut, only: TRadialTableOrbital
  use dftbp_wavegrid_basis_gaussian, only: TGaussianOrbital
  implicit none

  private

  public :: realTessY
  public :: TOrbital, maxCutoff
  public :: TSlaterOrbital, TGaussianOrbital, TRadialTableOrbital
  public :: TSlaterOrbital_init

end module dftbp_wavegrid_basis

