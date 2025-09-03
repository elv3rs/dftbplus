!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines to calculate a Slater type orbital (STO).
module libwavegrid_slater
  use dftbp_common_accuracy, only : dp
  implicit none

  private

  public :: TSlaterOrbital, realTessY


  !> Data type for STOs.
  type TSlaterOrbital
    !> Angular momentum of the orbital
    integer :: angMom

    !> Square of the Cutoff, after which the orbital is assumed to be zero
    real(dp) :: cutoffSq

    !> Whether to use cached values instead of direct calculation.
    logical :: useRadialLut = .false.

    ! ##  LUT parameters ##
    !> Grid distance (resolution)
    real(dp) :: gridDist

    !> Inverse of the grid distance
    real(dp) :: invLutStep

    !> Number of grid points
    integer :: nGrid

    !> STO values on the distance grid
    real(dp), allocatable :: gridValue(:)

    ! ## Direct calculation STO parameters ##
    !> Maximum power of the radial distance
    integer :: nPow

    !> Number of exponential coefficients
    integer :: nAlpha

    !> Summation coefficients. Shape: [nPow, nAlpha]
    real(dp), allocatable :: aa(:,:)

    !> Exponential coefficients
    real(dp), allocatable :: alpha(:)


    !! TODO: remove this from struct and store elsewhere.
    !> Occupation of the orbital (for atomic density)
    real(dp) :: occupation


  contains

    !> Initialises a SlaterOrbital.
    procedure :: init => TSlaterOrbital_init

    !> Initialises using a LUT for radial values.
    procedure :: initFromLut => TSlaterOrbital_initFromLut

    !> Returns the value of the Lut in a given point.
    procedure :: getRadialCached => TSlaterOrbital_getRadialValueCached
  
    !> Non-interpolated version using direct STO calculation.
    procedure :: getRadialDirect => TSlaterOrbital_getRadialValueDirect

    !> Dispatch to cached or direct version based on this%useRadialLut
    procedure :: getRadial
  
  end type TSlaterOrbital



contains


  !> Returns the real tesseral spherical harmonics in a given point.
  !! This function only work for angular momenta between 0 and 3 (s-f).
  function realTessY(ll, mm, coord, rrOpt) result(rty)
    !$omp declare target

    !> Angular momentum of the spherical harmonics (0 <= ll <= 3)
    integer, intent(in) :: ll

    !> Magnetic quantum number
    integer, intent(in) :: mm

    !> Coordinate where the value should be calculated
    real(dp), intent(in) :: coord(:)

    !> Length of the coordinate vector, if known in advance
    real(dp), intent(in), optional :: rrOpt

    real(dp) :: rty

    real(dp) :: rr, xx, yy, zz

    if (present(rrOpt)) then
      rr = rrOpt
    else
      rr = norm2(coord)
    end if

    @:ASSERT(ll >= 0 .and. ll <= 3)
    @:ASSERT(abs(mm) <= ll)
    @:ASSERT(size(coord) == 3)
    @:ASSERT(rr >= 0.0_dp)

    xx = coord(1)
    yy = coord(2)
    zz = coord(3)

    if (rr < epsilon(1.0_dp) .and. ll /= 0) then
      rty = 0.0_dp
      return
    end if

    select case (ll)
    case(0)
      rty = 0.2820947917738782_dp
    case(1)
      select case(mm)
      case(-1)
        ! y
        rty = 0.4886025119029198_dp * yy / rr
      case(0)
        ! z
        rty = 0.4886025119029198_dp * zz / rr
      case(1)
        ! x
        rty = 0.4886025119029198_dp * xx / rr
      end select
    case(2)
      select case(mm)
      case(-2)
        ! xy
        rty = 1.092548430592079_dp * xx * yy / rr**2
      case(-1)
        ! yz
        rty = 1.092548430592079_dp * yy * zz / rr**2
      case(0)
        ! z**2
        rty = -0.3153915652525200_dp * (-2.0_dp * zz**2 + xx**2 + yy**2) / rr**2
      case(1)
        ! xz
        rty = 1.092548430592079_dp * xx * zz / rr**2
      case(2)
        ! x**2-y**2
        rty = 0.5462742152960395_dp * (xx**2 - yy**2) / rr**2
      end select
    case(3)

      !> general set for f orbitals (not cubic), see
      !> http://winter.group.shef.ac.uk/orbitron/AOs/4f/equations.html
      select case (mm)
      case(-3)
        ! y(3x**2-y**2)
        rty = 0.5900435899266435_dp * yy * (3.0_dp * xx**2 - yy**2) / rr**3
      case(-2)
        ! x**2+y**2+z**2
        rty = 2.890611442640554_dp * xx * yy *zz / rr**3
      case(-1)
        ! yz**2
        rty = -0.4570457994644658_dp * (-4.0_dp * zz**2 + xx**2 + yy**2) * yy / rr**3
      case(0)
        ! z**3
        rty = -0.3731763325901155_dp * zz * (-2.0_dp * zz**2 + 3.0_dp * xx**2 + 3.0_dp * yy**2)&
            & / rr**3
      case(1)
        ! xz**2
        rty = -0.4570457994644658_dp * (-4.0_dp * zz**2 + xx**2 + yy**2) * xx / rr**3
      case(2)
        ! z(x**2-y**2)
        rty = 1.445305721320277_dp * zz * (xx**2 - yy**2) / rr**3
      case(3)
        ! x(x**2-3y**2)
        rty = 0.5900435899266435_dp * xx * (xx**2 - 3.0_dp * yy**2) / rr**3
      end select
    end select

  end function realTessY

  !> Initialises using a LUT for radial values.
  subroutine TSlaterOrbital_initFromLut(this, gridValue, gridDist, angMom)
    class(TSlaterOrbital), intent(inout) :: this
    real(dp), intent(in) :: gridValue(:)
    real(dp), intent(in) :: gridDist
    integer, intent(in) :: angMom
    real(dp) :: cutoff

    this%useRadialLut = .true.
    this%angMom = angMom
    this%gridDist = gridDist
    this%invLutStep = 1.0_dp / gridDist
    this%nGrid = size(gridValue)

    cutoff = this%gridDist * real(this%nGrid - 1, dp)
    this%cutoffSq = cutoff**2

    allocate(this%gridValue(this%nGrid))
    this%gridValue(:) = gridValue(:)

  end subroutine TSlaterOrbital_initFromLut

  !> Initialises a SlaterOrbital.
  subroutine TSlaterOrbital_init(this, aa, alpha, ll, resolution, cutoff, useRadialLut)

    !> SlaterOrbital instance to initialise
    class(TSlaterOrbital), intent(inout) :: this

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

    !> Whether to use the cached grid for evaluation
    logical, intent(in), optional :: useRadialLut

    integer :: iGrid, ii
    real(dp) :: rr

    @:ASSERT(cutoff > 0.0_dp)

    this%angMom = ll
    this%cutoffSq = cutoff ** 2

    ! Store parameters for direct calculation
    @:ASSERT(size(aa, dim=2) == this%nAlpha)
    this%nAlpha = size(alpha)
    this%nPow = size(aa, dim=1)
    allocate(this%aa(this%nPow, this%nAlpha))
    allocate(this%alpha(this%nAlpha))
    this%aa(:,:) = aa
    this%alpha(:) = -1.0_dp * alpha

    if (present(useRadialLut)) then
      this%useRadialLut = useRadialLut
    else
      this%useRadialLut = .false.
    end if

    if (this%useRadialLut) then
      ! Obtain STO on a grid
      @:ASSERT(resolution > 0.0_dp)
      this%nGrid = floor(cutoff / resolution) + 2
      this%gridDist = resolution
      this%invLutStep = 1.0_dp / resolution

      allocate(this%gridValue(this%nGrid))
      do iGrid = 1, this%nGrid
        rr = real(iGrid - 1, dp) * resolution
        this%gridValue(iGrid) = this%getRadialDirect(rr)
      end do
    end if

  end subroutine TSlaterOrbital_init


  function getRadial(this, rr) result(sto)
    class(TSlaterOrbital), intent(in) :: this
    real(dp), intent(in) :: rr
    real(dp) :: sto

    if (this%useRadialLut) then
      sto = this%getRadialCached(rr)
    else
      sto = this%getRadialDirect(rr)
    end if

  end function getRadial

  !> Returns the value of the SlaterOrbital in a given point.
  !! Builds a 1d cache grid across which the result is interpolated
  !! in order to speed up evaluation for subsequent calls.
  function TSlaterOrbital_getRadialValueCached(this, rr) result(sto)

    !> SlaterOrbital instance
    class(TSlaterOrbital), intent(in) :: this

    !> Distance, where STO should be calculated
    real(dp), intent(in) :: rr

    !> Contains the value of the function on return
    real(dp) :: sto

    integer :: ind
    real(dp) :: frac

    @:ASSERT(rr >= 0.0_dp)

    ! ind = 1 means zero distance as rr = (ind - 1) * gridDist
    ind = floor(rr * this%invLutStep) + 1
    if (ind < this%nGrid) then
      frac = mod(rr, this%gridDist) * this%invLutStep
      sto = (1.0_dp - frac) * this%gridValue(ind) + frac * this%gridValue(ind+1)
    else
      sto = 0.0_dp
    end if

  end function TSlaterOrbital_getRadialValueCached


  !> Calculates the value of an STO analytically.
  function TSlaterOrbital_getRadialValueDirect(this, rr) result(sto)

    !> SlaterOrbital instance
    class(TSlaterOrbital), intent(in) :: this

    !> Distance, where the STO should be calculated
    real(dp), intent(in) :: rr

    !> Value of the STO on return
    real(dp) :: sto

    real(dp) :: pows(this%nPow)
    real(dp) :: rTmp
    integer :: ii, jj

    ! Avoid 0.0**0 as it may lead to arithmetic exception
    if (this%angMom == 0 .and. rr < epsilon(1.0_dp)) then
      rTmp = 1.0_dp
    else
      rTmp = rr**this%angMom
    end if
    do ii = 1, this%nPow
      pows(ii) = rTmp
      rTmp = rTmp * rr
    end do
    sto = 0.0_dp
    do ii = 1, this%nAlpha
      rTmp = 0.0_dp
      do jj = 1, this%nPow
        rTmp = rTmp + this%aa(jj, ii) * pows(jj)
      end do
      sto = sto + rTmp * exp(this%alpha(ii) * rr)
    end do

  end function TSlaterOrbital_getRadialValueDirect


end module libwavegrid_slater
