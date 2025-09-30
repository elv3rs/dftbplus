!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set pure = "" if defined('WITH_ASSERT') else "pure"
!> Holds TSlaterOrbital, a concrete implementation of TOrbital.
module dftbp_wavegrid_basis_slater
  use dftbp_common_accuracy, only : dp
  use dftbp_wavegrid_basis_orbital, only : TOrbital
  implicit none

  private

  public :: TSlaterOrbital
  
  !> Concrete class for a Slater-Type Orbital (STO).
  type, extends(TOrbital) :: TSlaterOrbital
    !> Maximum power of the radial distance
    integer :: nPow

    !> Number of exponential coefficients
    integer :: nAlpha

    !> Summation coefficients. Shape: [nPow, nAlpha]
    real(dp), allocatable :: aa(:,:)

    !> Exponential coefficients (stored negated)
    real(dp), allocatable :: alpha(:)
  contains
    procedure :: getRadial => TSlaterOrbital_getRadial
    procedure :: init => TSlaterOrbital_init
  end type TSlaterOrbital

contains

  !> Initialises using STO parameters.
  subroutine TSlaterOrbital_init(this, aa, alpha, ll, resolution, cutoff)

    !> TSlaterOrbital instance to initialise
    class(TSlaterOrbital), intent(out) :: this

    !> Summation coefficients (nCoeffPerAlpha, nAlpha)
    real(dp), intent(in) :: aa(:,:)

    !> Exponential coefficients
    real(dp), intent(in) :: alpha(:)

    !> Angular momentum of the orbital
    integer, intent(in) :: ll

    !> Grid distance for the orbital
    real(dp), intent(in) :: resolution

    !> Cutoff, after which orbital is assumed to be zero
    real(dp), intent(in) :: cutoff

    integer :: iGrid
    real(dp) :: r

    @:ASSERT(cutoff > 0.0_dp)

    this%angMom = ll
    this%cutoffSq = cutoff ** 2

    this%nAlpha = size(alpha)
    @:ASSERT(size(aa, dim=2) == this%nAlpha)
    this%nPow = size(aa, dim=1)

    allocate(this%aa(this%nPow, this%nAlpha), source=aa)

    ! The STO formula uses exp(-alpha * r).
    ! Directly store -alpha to avoid repeated negation when calculating.
    allocate(this%alpha(this%nAlpha))
    this%alpha(:) = -1.0_dp * alpha(:)

  end subroutine TSlaterOrbital_init


  !> Calculates the value of an STO analytically.
  ${pure}$ function TSlaterOrbital_getRadial(this, r) result(sto)

    !> SlaterOrbital instance
    class(TSlaterOrbital), intent(in) :: this

    !> Distance, where the STO should be calculated
    real(dp), intent(in) :: r

    !> Value of the STO on return
    real(dp) :: sto

    real(dp) :: pows(this%nPow)
    real(dp) :: rTmp
    integer :: ii, jj

    ! Avoid 0.0**0 as it may lead to arithmetic exception on pre-2008 compilers
    if (this%angMom == 0 .and. r < epsilon(1.0_dp)) then
      rTmp = 1.0_dp
    else
      rTmp = r**this%angMom
    end if

    do ii = 1, this%nPow
      pows(ii) = rTmp
      rTmp = rTmp * r
    end do
    sto = 0.0_dp
    do ii = 1, this%nAlpha
      rTmp = 0.0_dp
      do jj = 1, this%nPow
        rTmp = rTmp + this%aa(jj, ii) * pows(jj)
      end do
      sto = sto + rTmp * exp(this%alpha(ii) * r)
    end do

  end function TSlaterOrbital_getRadial


end module dftbp_wavegrid_basis_slater
