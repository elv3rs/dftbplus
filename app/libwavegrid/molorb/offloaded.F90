!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module libwavegrid_molorb_offloaded
  use dftbp_common_accuracy, only : dp
  use libwavegrid_molorb_types, only : TSystemParams, TPeriodicParams, TBasisParams, TCalculationContext
  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_double_complex, c_bool, c_ptr, c_loc, &
      & c_null_ptr
  implicit none
  private

#:if WITH_CUDA
  public :: prepareGPUCoefficients, evaluateCuda
#:endif

contains
#:if WITH_CUDA
  !> Prepare coefficient vectors for GPU calculation by, if required due to total charge calculation,
  !> scaling the eigenvectors by sqrt(occupationVec).
  subroutine prepareGPUCoefficients(ctx, eigVecsReal, eigVecsCmpl, occupationVec, coeffVecReal, coeffVecCmpl)
    type(TCalculationContext), intent(in) :: ctx
    real(dp), intent(in) :: eigVecsReal(:,:)
    complex(dp), intent(in) :: eigVecsCmpl(:,:)
    real(dp), intent(in), optional :: occupationVec(:)
    real(dp), allocatable, intent(out) :: coeffVecReal(:,:)
    complex(dp), allocatable, intent(out) :: coeffVecCmpl(:,:)

    integer :: iEig

    allocate(coeffVecReal(size(eigVecsReal, dim=1), size(eigVecsReal, dim=2)))
    allocate(coeffVecCmpl(size(eigVecsCmpl, dim=1), size(eigVecsCmpl, dim=2)))

    if (ctx%calcTotalChrg) then
      print *, "Baking occupationVec into GPU Eigenvector coefficients"
      @:ASSERT(size(occupationVec) == size(eigVecsReal, dim=2))
      do iEig = 1, size(eigVecsReal, dim=2)
        coeffVecReal(:, iEig) = eigVecsReal(:, iEig) * sqrt(occupationVec(iEig))
      end do
      do iEig = 1, size(eigVecsCmpl, dim=2)
        coeffVecCmpl(:, iEig) = eigVecsCmpl(:, iEig) * sqrt(occupationVec(iEig))
      end do
    else
      coeffVecReal = eigVecsReal
      coeffVecCmpl = eigVecsCmpl
    end if

  end subroutine prepareGPUCoefficients

  subroutine evaluateCuda(origin, gridVecs, &
      & system, basis, periodic, kIndexes, phases, &
      & ctx, &
      & eigVecsReal, eigVecsCmpl, valueReal, valueCmpl)

    !> Grid
    real(dp), intent(in), target :: origin(3)
    real(dp), intent(in), target :: gridVecs(3, 3)
    !> System
    type(TSystemParams), intent(in), target :: system
    !> Basis set
    type(TBasisParams), intent(in), target :: basis
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

    type, bind(c) :: TBasisParamsC
      integer(c_int) :: nStos, maxNPows, maxNAlphas
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
        type(TBasisParamsC), intent(in) :: basis
        type(TCalculationParamsC), intent(in) :: calc
      end subroutine evaluate_on_device_c
    end interface

    type(TGridParamsC) :: grid_p
    type(TSystemParamsC) :: system_p
    type(TPeriodicParamsC) :: periodic_p
    type(TBasisParamsC) :: basis_p
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
    grid_p%origin = c_loc(origin)
    grid_p%gridVecs = c_loc(gridVecs)
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
    basis_p%nStos = basis%nStos
    basis_p%maxNPows = basis%maxNPows
    basis_p%maxNAlphas = basis%maxNAlphas
    basis_p%sto_angMoms = c_loc(basis%angMoms)
    basis_p%sto_nPows = c_loc(basis%nPows)
    basis_p%sto_nAlphas = c_loc(basis%nAlphas)
    basis_p%sto_cutoffsSq = c_loc(basis%cutoffsSq)
    basis_p%sto_coeffs = c_loc(basis%coeffs)
    basis_p%sto_alphas = c_loc(basis%alphas)
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
#:endif

end module libwavegrid_molorb_offloaded
