
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set pure = "" if defined('WITH_ASSERT') else "pure"
!> Holds TRadialTableOrbital, a concrete implementation of TOrbital utilising a 1d interpolated lookup table.
!> This allows for arbitrary radial functions.
module dftbp_wavegrid_basis_lut
  use dftbp_common_accuracy, only : dp
  use dftbp_wavegrid_basis_orbital, only : TOrbital
  implicit none

  private

  public :: TRadialTableOrbital
  
  !> Concrete class for an orbital represented by a linearly
  !! interpolated lookup table.
  type, extends(TOrbital) :: TRadialTableOrbital
    !> Grid spacing (resolution)
    real(dp) :: gridDist
    !> Inverse of the grid spacing
    real(dp) :: invLutStep
    !> Orbital values on the grid
    real(dp), allocatable :: gridValue(:)
  contains
    procedure :: getRadial => TRadialTable_getRadial
    procedure :: initFromArray => TRadialTableOrbital_initFromArray
    procedure :: initFromOrbital => 
  end type TRadialTableOrbital

contains

  !> Initialises using a verbatim array of values.
  subroutine TRadialTable_initFromArray(this, gridValue, gridDist, angMom)
    class(TRadialTable), intent(out) :: this
    real(dp), intent(in) :: gridValue(:)
    real(dp), intent(in) :: gridDist
    integer, intent(in) :: angMom
    real(dp) :: cutoff

    this%useRadialLut = .true.
    this%angMom = angMom
    this%gridDist = gridDist
    this%invLutStep = 1.0_dp / gridDist

    cutoff = this%gridDist * real(size(gridValue) - 1, dp)
    this%cutoffSq = cutoff**2

    allocate(this%gridValue, source=gridValue)

  end subroutine TRadialTable_initFromArray

  !> Resamples another orbital onto a LUT with given resolution.
  subroutine TRadialTable_initFromOrbital(this, other, resolution)
    class(TRadialTable), intent(out) :: this
    class(TOrbital), intent(in) :: other
    real(dp), intent(in) :: resolution
    integer :: iGrid
    real(dp) :: r, cutoff

    @:ASSERT(resolution > 0.0_dp)

    ! Set parameters
    this%angMom = other%angMom
    this%cutoffSq = other%cutoffSq
    this%gridDist = resolution
    this%invLutStep = 1.0_dp / resolution
    
    ! Allocate LUT grid
    cutoff = sqrt(this%cutoffSq)
    allocate(this%gridValue(floor(cutoff / resolution) + 2))

    ! Populate LUT by sampling the other orbital
    do iGrid = 1, size(this%gridValue)
      r = real(iGrid - 1, dp) * resolution
      this%gridValue(iGrid) = other%getRadial(r)
    end do

  end subroutine TOrbital_resampleToLut


  !> Returns the value of the RadialFunction at a given point.
  !! Builds a 1d cache grid across which the result is interpolated
  !! in order to speed up evaluation for subsequent calls.
  ${pure}$ function TRadialTable_getRadial(this, r) result(sto)

    !> RadialTable instance
    class(TRadialTable), intent(in) :: this

    !> Distance, where STO should be calculated
    real(dp), intent(in) :: r

    !> Contains the value of the function on return
    real(dp) :: sto

    integer :: ind
    real(dp) :: frac, posOnGrid

    @:ASSERT(r >= 0.0_dp)

    ! ind = 1 means zero distance as r = (ind - 1) * gridDist
    posOnGrid = r * this%invLutStep
    ind = floor(posOnGrid) + 1
    if (ind < this%nGrid) then
      frac = posOnGrid - real(ind - 1, dp)
      sto = (1.0_dp - frac) * this%gridValue(ind) + frac * this%gridValue(ind+1)
    else
      sto = 0.0_dp
    end if

  end function TRadialTable_getRadial
end module dftbp_wavegrid_basis_lut
