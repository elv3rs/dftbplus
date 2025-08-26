!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module libwavegrid_molorb_offloaded
  use dftbp_common_accuracy, only : dp
  use libwavegrid_molorb_types, only : TSystemParams, TPeriodicParams, TCalculationContext
  use libwavegrid_slater, only : TSlaterOrbital
  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_double_complex, c_bool, c_ptr, c_loc, &
      & c_null_ptr
  implicit none
  private

#:if WITH_CUDA
  public :: evaluateCuda
#:endif
  !> Data for the basis set in SoA format
  type TBasisParams
    integer :: nStos
    logical :: useRadialLut

    integer :: nLutPoints
    real(dp) :: invLutStep
    real(dp), allocatable :: lutGridValues(:, :)

    !! SoA
    integer :: maxNPows
    integer :: maxNAlphas
    integer, allocatable :: angMoms(:)
    real(dp), allocatable :: cutoffsSq(:)
    integer, allocatable :: nPows(:)
    integer, allocatable :: nAlphas(:)
    real(dp), allocatable :: coeffs(:,:,:)
    real(dp), allocatable :: alphas(:,:)

    logical :: isInitialized = .false.
  end type TBasisParams
contains
#:if WITH_CUDA
  subroutine evaluateCuda(system, stos, periodic, kIndexes, phases, ctx, &
      & eigVecsReal, eigVecsCmpl, valueReal, valueCmpl)

    !> System
    type(TSystemParams), intent(in), target :: system
    !> Basis set
    type(TSlaterOrbital), intent(in), target :: stos(:)
    !> Periodic boundary conditions
    type(TPeriodicParams), intent(in), target :: periodic
    integer, intent(in), target :: kIndexes(:)
    complex(dp), intent(in), target :: phases(:, :)
    !> Calculation flags
    type(TCalculationContext), intent(in) :: ctx
    !> Eigenvectors
    real(dp), intent(in), target :: eigVecsReal(:, :)
    complex(dp), intent(in), target :: eigVecsCmpl(:, :)
    !> Output grids
    real(dp), intent(out), target :: valueReal(:, :, :, :)
    complex(dp), intent(out), target :: valueCmpl(:, :, :, :)

    type, bind(c) :: TGridParamsC
      integer(c_int) :: nPointsX, nPointsY, nPointsZ
      type(c_ptr) :: origin, gridVecs
    end type

    type, bind(c) :: TSystemParamsC
      integer(c_int) :: nAtom, nCell, nSpecies, nOrb
      type(c_ptr) :: coords, species, iStos
    end type

    type, bind(c) :: TPeriodicParamsC
      integer(c_int) :: isPeriodic
      type(c_ptr) :: latVecs, recVecs2pi, kIndexes, phases
    end type

    type, bind(c) :: TSlaterOrbitalC
      integer(c_int) :: nStos, useRadialLut, nLutPoints
      real(c_double) :: inverseLutStep
      type(c_ptr) :: lutGridValues

      integer(c_int) :: maxNPows, maxNAlphas
      type(c_ptr) :: sto_angMoms, sto_nPows, sto_nAlphas
      type(c_ptr) :: sto_cutoffsSq, sto_coeffs, sto_alphas
    end type

    type, bind(c) :: TCalculationParamsC
      integer(c_int) :: nEigIn, nEigOut, isRealInput, calcAtomicDensity, calcTotalChrg
      type(c_ptr) :: eigVecsReal, eigVecsCmpl
      type(c_ptr) :: valueReal_out, valueCmpl_out
    end type


    interface
      subroutine evaluate_on_device_c(grid, system, periodic, basis, calc) bind(C, name ='evaluate_on_device_c')
        import
        type(TGridParamsC), intent(in) :: grid
        type(TSystemParamsC), intent(in) :: system
        type(TPeriodicParamsC), intent(in) :: periodic
        type(TSlaterOrbitalC), intent(in) :: basis
        type(TCalculationParamsC), intent(in) :: calc
      end subroutine evaluate_on_device_c
    end interface

    type(TBasisParams), target :: basis

    type(TGridParamsC) :: grid_p
    type(TSystemParamsC) :: system_p
    type(TPeriodicParamsC) :: periodic_p
    type(TSlaterOrbitalC) :: basis_p
    type(TCalculationParamsC) :: calc_p

    ! Populate the structs
    if (ctx%isRealOutput) then
      grid_p%nPointsX = size(valueReal, dim=1)
      grid_p%nPointsY = size(valueReal, dim=2)
      grid_p%nPointsZ = size(valueReal, dim=3)
    else
      grid_p%nPointsX = size(valueCmpl, dim=1)
      grid_p%nPointsY = size(valueCmpl, dim=2)
      grid_p%nPointsZ = size(valueCmpl, dim=3)
    end if
    grid_p%origin = c_loc(system%origin)
    grid_p%gridVecs = c_loc(system%gridVecs)

    system_p%nAtom = system%nAtom
    system_p%nCell = size(system%coords, dim=3)
    system_p%nSpecies = system%nSpecies
    system_p%nOrb = system%nOrb
    system_p%coords = c_loc(system%coords)
    system_p%species = c_loc(system%species)
    system_p%iStos = c_loc(system%iStos)

    periodic_p%isPeriodic = merge(1, 0, periodic%isPeriodic)
    periodic_p%latVecs = c_loc(periodic%latVecs)
    periodic_p%recVecs2pi = c_loc(periodic%recVecs2pi)
    periodic_p%kIndexes = c_loc(kIndexes)
    periodic_p%phases = c_loc(phases)


    call prepareBasisSet(basis, stos)
    basis_p%useRadialLut = merge(1, 0, basis%useRadialLut)
    if (basis%useRadialLut) then
      basis_p%nLutPoints = basis%nLutPoints
      basis_p%inverseLutStep = basis%invLutStep
      basis_p%lutGridValues = c_loc(basis%lutGridValues)
    else
      basis_p%maxNPows = basis%maxNPows
      basis_p%maxNAlphas = basis%maxNAlphas
      basis_p%sto_angMoms = c_loc(basis%angMoms)
      basis_p%sto_nPows = c_loc(basis%nPows)
      basis_p%sto_nAlphas = c_loc(basis%nAlphas)
      basis_p%sto_cutoffsSq = c_loc(basis%cutoffsSq)
      basis_p%sto_coeffs = c_loc(basis%coeffs)
      basis_p%sto_alphas = c_loc(basis%alphas)
    end if


    if (ctx%isRealInput) then
      calc_p%nEigIn = size(eigVecsReal, dim=2)
    else
      calc_p%nEigIn = size(eigVecsCmpl, dim=2)
    end if
    if (ctx%isRealOutput) then
      calc_p%nEigOut = size(valueReal, dim=4)
    else
      calc_p%nEigOut = size(valueCmpl, dim=4)
    end if
    if (ctx%calcTotalChrg) then
      @:ASSERT(calc_p%nEigOut == 1)
    end if
    calc_p%isRealInput = merge(1, 0, ctx%isRealInput)
    calc_p%calcAtomicDensity = merge(1, 0, ctx%calcAtomicDensity)
    calc_p%calcTotalChrg = merge(1, 0, ctx%calcTotalChrg)
    calc_p%eigVecsReal = c_loc(eigVecsReal)
    calc_p%eigVecsCmpl = c_loc(eigVecsCmpl)
    calc_p%valueReal_out = c_loc(valueReal)
    calc_p%valueCmpl_out = c_loc(valueCmpl)

    call evaluate_on_device_c(grid_p, system_p, periodic_p, basis_p, calc_p)

  end subroutine evaluateCuda
 
  ! Convert the basis set to SoA format or unified Lut table
  ! Currently, mixed lut/direct calculation is not supported.
  ! Additionally, all orbitals must use identical LUT settings.
  subroutine prepareBasisSet(this, stos)
    type(TBasisParams), intent(out) :: this
    type(TSlaterOrbital), intent(in) :: stos(:)
    integer :: iOrb

    this%nStos = size(stos)
    ! We will assert that all orbitals use identical LUT settings
    this%useRadialLut = stos(1)%useRadialLut
    this%nLutPoints = stos(1)%nGrid
    this%invLutStep = stos(1)%invLutStep

    if (this%useRadialLut) then
      print *, "Using radial LUT for STO evaluation with ", this%nLutPoints, " points."
      allocate(this%lutGridValues(this%nLutPoints, this%nStos))
      do iOrb = 1, this%nStos
        @:ASSERT(stos(iOrb)%useRadialLut)
        @:ASSERT(stos(iOrb)%nGrid == this%nLutPoints)
        @:ASSERT(abs(stos(iOrb)%invLutStep - this%invLutStep) < 1.0e-12_dp)

        this%lutGridValues(:, iOrb) = stos(iOrb)%gridValue
      end do
    else ! Direct evaluation
      ! Allocate SoA arrays
      allocate(this%angMoms(this%nStos))
      allocate(this%nPows(this%nStos))
      allocate(this%nAlphas(this%nStos))
      allocate(this%cutoffsSq(this%nStos))

      ! Populate SoA arrays
      do iOrb = 1, this%nStos
        @:ASSERT(.not. stos(iOrb)%useRadialLut)
        this%angMoms(iOrb) = stos(iOrb)%angMom
        this%cutoffsSq(iOrb) = stos(iOrb)%cutoffSq
        this%nPows(iOrb) = stos(iOrb)%nPow
        this%nAlphas(iOrb) = stos(iOrb)%nAlpha
      end do
      this%maxNPows = maxval(this%nPows)
      this%maxNAlphas = maxval(this%nAlphas)

      ! Allocate and populate coefficient/alpha matrices
      allocate(this%coeffs(this%maxNPows, this%maxNAlphas, this%nStos))
      allocate(this%alphas(this%maxNAlphas, this%nStos))
      do iOrb = 1, this%nStos
        this%coeffs(1:stos(iOrb)%nPow, 1:stos(iOrb)%nAlpha, iOrb) = stos(iOrb)%aa
        this%alphas(1:stos(iOrb)%nAlpha, iOrb) = stos(iOrb)%alpha
      end do
    end if

    this%isInitialized = .true.
  end subroutine prepareBasisSet
#:endif

end module libwavegrid_molorb_offloaded
