!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
!> This module defines TOrbital, which, aside from holding the angular momentum,
!! may be subclassed by one of three concrete radial function implementations:
!! Slater Type Orbitals (STO), Gaussian Type Orbitals (GTO) or a radial lookup table (LUT),
!! allowing arbitrary radial functions to be used. Any orbital may be resampled onto a LUT
!! with a new chosen resolution and cutoff.

#:include 'common.fypp'
#:set pure = "" if defined('WITH_ASSERT') else "pure"

module dftbp_wavegrid_basis_orbital
  use dftbp_common_accuracy, only : dp
  implicit none


  public :: TOrbital
  
  type, abstract :: TOrbital
    !> Angular momentum (l)
    integer :: angMom = -1
    !> Square of the Cutoff, after which the orbital is assumed to be zero
    real(dp) :: cutoffSq
  contains
    procedure(IGetRadial), deferred :: getRadial
    procedure(IAssign), deferred, pass(lhs) :: assign
    generic, public :: assignment(=) => assign
    procedure :: getNorm => TOrbital_getNorm
    procedure :: getThresholdCutoff => TOrbital_getThresholdCutoff
    procedure :: getDensityCutoff => TOrbital_getDensityCutoff
  end type TOrbital

  abstract interface
    ${pure}$ function IGetRadial(this, r) result(val)
      import :: TOrbital, dp
      class(TOrbital), intent(in) :: this
      real(dp), intent(in) :: r
      real(dp) :: val
    end function IGetRadial

    subroutine IAssign(lhs, rhs)
      import :: TOrbital
      class(TOrbital), intent(out) :: lhs
      class(TOrbital), intent(in) :: rhs
    end subroutine IAssign
  end interface
contains

  !> Returns the norm of the orbital (i.e. sqrt(integral of R(r)^2 * r^2 dr))
  ${pure}$ function TOrbital_getNorm(this, nGrid) result(norm)
    class(TOrbital), intent(in) :: this
    real(dp) :: norm
    integer, intent(in), optional :: nGrid
    integer :: iGrid, segments
    real(dp) :: dr, r, integral, fPrev, fCurr
    
    segments = 1000
    if (present(nGrid)) then
      segments = nGrid
    end if
   
    dr = sqrt(this%cutoffSq) / real(segments, dp)
    integral = 0.0_dp
    fPrev = 0.0_dp
    do iGrid = 1, segments
      r = iGrid * dr
      fCurr = this%getRadial(r)**2 * r**2
      ! Trapezoidal segments
      integral = integral + 0.5_dp * (fPrev + fCurr) * dr
      fPrev = fCurr
    end do
    norm = sqrt(integral)

  end function TOrbital_getNorm

  !> Returns the distance at which the orbital has dropped below <threshold>.
  !! Works backwards from the original cutoff, or <maxCutoff>, if passed.
  ${pure}$ function TOrbital_getThresholdCutoff(this, thresholdIn, maxCutoff, nGrid) result(cutoff)
    class(TOrbital), intent(in) :: this
    !> Threshold value used (default: 1.0e-5)
    real(dp), intent(in), optional :: thresholdIn
    !> Maximum cutoff considered (default: sqrt(this%cutoffSq))
    real(dp), intent(in), optional :: maxCutoff
    !> Number of grid points to use for searching (default: 1000)
    integer, intent(in), optional :: nGrid
    real(dp) :: cutoff
    integer :: segments, iGrid
    real(dp) :: dr, threshold
    
    segments = 1000
    if (present(nGrid)) segments = nGrid

    threshold = 1.0e-5_dp
    if (present(thresholdIn)) threshold = thresholdIn

    if (present(maxCutoff)) then
      cutoff = min(sqrt(this%cutoffSq), maxCutoff)
    else
      cutoff = sqrt(this%cutoffSq)
    end if

    dr = cutoff / real(segments, dp)
    ! Move inwards from original cutoff until threshold is exceeded
    do iGrid = segments, 1, -1
      if (abs(this%getRadial(iGrid * dr)) > threshold) exit
      cutoff = iGrid * dr
    end do

  end function TOrbital_getThresholdCutoff

  !> Returns the distance enclosing a certain fraction (default: 99.99%) of the total density.
  ${pure}$ function TOrbital_getDensityCutoff(this, densityFrac, nGrid) result(cutoff)
    class(TOrbital), intent(in) :: this
    !> Fraction of the total density to enclose (default: 0.9999)
    real(dp), intent(in), optional :: densityFrac
    !> Number of grid points to use for numerical integration (default: 1000)
    integer, intent(in), optional :: nGrid
    real(dp) :: cutoff
    integer :: segments
    real(dp) :: dr, integral, targetIntegral, frac, fPrev, fCurr

    frac = 0.9999_dp
    if (present(densityFrac)) frac = densityFrac
    segments = 1000
    if (present(nGrid)) segments = nGrid

    ! Do not require a normalized orbital to begin with
    targetIntegral = frac * this%getNorm(nGrid=segments)**2

    dr = sqrt(this%cutoffSq) / real(segments, dp)
    integral = 0.0_dp
    cutoff = 0.0_dp
    fPrev = 0.0_dp
    do while (cutoff < sqrt(this%cutoffSq))
      cutoff = cutoff + dr
      fCurr = this%getRadial(cutoff)**2 * cutoff**2
      integral = integral + 0.5_dp * (fPrev + fCurr) * dr
      fPrev = fCurr
      if (integral >= targetIntegral) exit
    end do

  end function TOrbital_getDensityCutoff

end module dftbp_wavegrid_basis_orbital

