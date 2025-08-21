!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module libwavegrid_molorb_parallel
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
  use libwavegrid_molorb_types, only : TSystemParams, TPeriodicParams, TBasisParams, TCalculationContext
  use libwavegrid_slater, only : realTessY, getRadial
  use omp_lib, only : omp_is_initial_device, omp_get_num_devices
  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_double_complex, c_bool, c_ptr, c_loc, &
      & c_null_ptr
  implicit none
  private
  !> Max powers expected in STO basis.
  integer, parameter :: MAX_STO_POWS = 16

  public :: evaluateParallel

contains


  !> Returns the values of several molecular orbitals on grids.
  !> This dispatches to either CPU / GPU implementation, and handles total Charge calculation using occupationVec if present.
  subroutine evaluateParallel(origin, gridVecs, system, periodic, kIndexes, phases, basis, &
      & ctx, eigVecsReal, eigVecsCmpl,  valueReal, valueCmpl, occupationVec)

    !> Origin of the grid
    real(dp), intent(in) :: origin(:)
    !> Grid vectors
    real(dp), intent(in) :: gridVecs(:,:)

    !> System geometry and composition
    type(TSystemParams), intent(in) :: system

    !> Periodic boundary conditions data
    type(TPeriodicParams), intent(in) :: periodic
    !> Index of the k-points for each orbital in kPoints
    integer, intent(in) :: kIndexes(:)
    !> Phase factors for periodic images
    complex(dp), intent(in) :: phases(:,:)

    !> Basis set data in SoA format
    type(TBasisParams), intent(in) :: basis

    !> Calculation control flags
    type(TCalculationContext), intent(in) :: ctx

    !> Real eigenvectors, or null-array
    real(dp), intent(in) :: eigVecsReal(:,:)
    !> Complex eigenvectors, or null-array
    complex(dp), intent(in) :: eigVecsCmpl(:,:)

    !> Contains the real grid on exit
    real(dp), intent(out) :: valueReal(:,:,:,:)
    !> Contains the complex grid on exit
    complex(dp), intent(out) :: valueCmpl(:,:,:,:)

    !> If this is present, calculate total charge. I.e. sum the squared states weighted by occupationVec(nEig).
    !> Output valueReal and valueCmpl will collapse to one slice in the last dimension, (x,y,z,1).
    real(dp), intent(in), optional :: occupationVec(:)
    !! Variables for total charge calculation
    integer :: iEig, iStart, iEnd, nChunk, iEigInChunk, nEigs, nEigsPerChunk
    real(dp), allocatable :: bufferReal(:,:,:,:), coeffVecReal(:,:)
    complex(dp), allocatable :: bufferCmpl(:,:,:,:), coeffVecCmpl(:,:)

    if (ctx%calcTotalChrg) then
      @:ASSERT(size(valueReal, dim=4) <= 1)
      @:ASSERT(size(valueCmpl, dim=4) <= 1)
    end if


    if (ctx%runOnGPU) then
      #:if WITH_CUDA
        ! GPU implementation passes occupation information by baking their sqrt into the eigenvectors
        call prepareGPUCoefficients(ctx, eigVecsReal, eigVecsCmpl, occupationVec, coeffVecReal, coeffVecCmpl)

        call evaluateCuda(origin, gridVecs, &
            & system, basis, periodic, kIndexes, phases, &
            & ctx, &
            & coeffVecReal, coeffVecCmpl, valueReal, valueCmpl)
      #:endif
    else ! CPU implementation
      #:if WITH_OMP
        print *, "Using CPU OMP molorb calculation."
      #:else
        print *, "Using CPU serial molorb calculation. OMP not available. (will run slower)"
      #:endif
      if (.not. ctx%calcTotalChrg) then
        call evaluateOMP(origin, gridVecs, &
            & system, basis, periodic, kIndexes, phases, &
            & ctx, &
            & eigVecsReal, eigVecsCmpl, valueReal, valueCmpl)
      else
        ! Number of eigenvectors to calculate at once in a chunk.
        ! We need a function to query free RAM to dynamically size this.
        ! Additionally, expose to user.
        nEigsPerChunk = -1 

        nEigs = size(occupationVec)
        if (nEigsPerChunk <= 0 .or. nEigsPerChunk > nEigs) then
          nEigsPerChunk = nEigs
        end if
        print *, "Using CPU molorb calculation with ", nEigsPerChunk, " eigenvectors per chunk."
        if (ctx%isRealInput) then
          allocate(bufferReal(size(valueReal, 1), size(valueReal, 2), size(valueReal, 3), nEigsPerChunk))
        else ! Complex input
          allocate(bufferCmpl(size(valueCmpl, 1), size(valueCmpl, 2), size(valueCmpl, 3), nEigsPerChunk))
        end if

        ! Zero accumulator
        valueReal(:, :, :, 1) = 0.0_dp

        lpChunk: do iStart = 1, nEigs, nEigsPerChunk
          iEnd = min(iStart + nEigsPerChunk - 1, nEigs)
          nChunk = iEnd - iStart + 1

          if (ctx%isRealInput) then
            call evaluateOMP(origin, gridVecs, &
                & system, basis, periodic, kIndexes, phases, &
                & ctx, &
                & eigVecsReal(:, iStart:iEnd), eigVecsCmpl, bufferReal, valueCmpl)

            do iEigInChunk = 1, nChunk
              iEig = iStart + iEigInChunk - 1
              valueReal(:,:,:,1) = valueReal(:,:,:,1) + bufferReal(:,:,:,iEigInChunk)**2 * occupationVec(iEig)
            end do
          else ! Complex input
            call evaluateOMP(origin, gridVecs, &
                & system, basis, periodic, kIndexes(iStart:iEnd), phases, &
                & ctx, &
                & eigVecsReal, eigVecsCmpl(:, iStart:iEnd), valueReal, bufferCmpl)

            do iEigInChunk = 1, nChunk
              iEig = iStart + iEigInChunk - 1
              valueReal(:,:,:,1) = valueReal(:,:,:,1) + abs(bufferCmpl(:,:,:,iEigInChunk))**2 * occupationVec(iEig)
            end do
          end if
        end do lpChunk
        if (allocated(bufferReal)) then
          deallocate(bufferReal)
        end if
        if (allocated(bufferCmpl)) then
          deallocate(bufferCmpl)
        end if
      end if
    end if

  end subroutine evaluateParallel



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


#:if WITH_CUDA
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
      integer(c_int) :: nEigIn, nEigOut, isRealInput, isDensityCalc, calcTotalChrg
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
    grid_p%nPointsX = size(valueReal, dim=1)
    grid_p%nPointsY = size(valueReal, dim=2)
    grid_p%nPointsZ = size(valueReal, dim=3)
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
    calc_p%isDensityCalc = merge(1, 0, ctx%isDensityCalc)
    calc_p%calcTotalChrg = merge(1, 0, ctx%calcTotalChrg)
    calc_p%eigVecsReal = c_loc(eigVecsReal)
    calc_p%eigVecsCmpl = c_loc(eigVecsCmpl)
    calc_p%valueReal_out = c_loc(valueReal)
    calc_p%valueCmpl_out = c_loc(valueCmpl)

    call evaluate_on_device_c(grid_p, system_p, periodic_p, basis_p, calc_p)

  end subroutine evaluateCuda
#:endif




  subroutine evaluateOMP(origin, gridVecs, &
      & system, basis, periodic, kIndexes, phases, &
      & ctx, &
      & eigVecsReal, eigVecsCmpl, valueReal, valueCmpl)

    !> Grid
    real(dp), intent(in) :: origin(3)
    real(dp), intent(in) :: gridVecs(3, 3)
    !> System
    type(TSystemParams), intent(in) :: system
    !> Basis set
    type(TBasisParams), intent(in) :: basis
    !> Periodic boundary conditions
    type(TPeriodicParams), intent(in) :: periodic
    integer, intent(in) :: kIndexes(:)
    complex(dp), intent(in) :: phases(:, :)
    !> Calculation flags
    type(TCalculationContext), intent(in) :: ctx
    !> Eigenvectors
    real(dp), intent(in) :: eigVecsReal(:, :)
    complex(dp), intent(in) :: eigVecsCmpl(:, :)
    !> Output grids
    real(dp), intent(out) :: valueReal(:, :, :, :)
    complex(dp), intent(out) :: valueCmpl(:, :, :, :)

    !! Thread private variables
    integer ::  ind, iSpecies
    real(dp) :: tmp_pows(MAX_STO_POWS), xyz(3), diff(3)
    real(dp) :: rSq, r, val, radialVal, tmp_rexp, frac(3)
    !! Loop Variables
    integer :: i1, i2, i3, iEig, iAtom, iOrb, iM, iL, iCell

    @:ASSERT(basis%maxNPows <= MAX_STO_POWS)

    !$omp parallel do collapse(3) &
    !$omp&    private(i1, i2, i3, iCell, iAtom, iOrb, iEig, iL, iM, xyz, diff, &
    !$omp&              r, val, radialVal, tmp_pows, tmp_rexp, ind, iSpecies, rSq) &
    !$omp&    shared(gridVecs, origin, system, basis, periodic, &
    !$omp&              eigVecsReal, eigVecsCmpl, phases, ctx, valueReal, valueCmpl)
    lpI3: do i3 = 1, size(valueReal, dim=3)
      lpI2: do i2 = 1, size(valueReal, dim=2)
        lpI1: do i1 = 1, size(valueReal, dim=1)
            valueReal(i1, i2, i3, :) = 0.0_dp
            xyz(:) = origin(:) + real(i1 - 1, dp) * gridVecs(:, 1) &
                             & + real(i2 - 1, dp) * gridVecs(:, 2) &
                             & + real(i3 - 1, dp) * gridVecs(:, 3)

            ! Map grid coordinates into unit cell
            if (periodic%isPeriodic) then
              frac(:) = matmul(xyz, periodic%recVecs2pi)
              xyz(:) = matmul(periodic%latVecs, frac - real(floor(frac), dp))
            end if

            ! Get contribution from every atom in every cell for current point
            lpCell: do iCell = 1, size(system%coords, dim=3)
              ind = 0
              lpAtom: do iAtom = 1, size(system%coords, dim=2)
                iSpecies = system%species(iAtom)
                diff(:) = xyz - system%coords(:, iAtom, iCell)
                rSq = dot_product(diff, diff)

                lpOrb: do iOrb = system%iStos(iSpecies), system%iStos(iSpecies + 1) - 1
                  iL = basis%angMoms(iOrb)
                  ! Calculate wave function only if atom is inside the cutoff
                  if (rSq > basis%cutoffsSq(iOrb)) then
                    ind = ind + 2*iL + 1
                    cycle lpOrb
                  end if
                  r = sqrt(rSq)

                  radialVal = getRadial(iL, &
                              & basis%nPows(iOrb), &
                              & basis%nAlphas(iOrb), &
                              & basis%coeffs(1:basis%nPows(iOrb), 1:basis%nAlphas(iOrb), iOrb), &
                              & basis%alphas(1:basis%nAlphas(iOrb), iOrb), &
                              & r)

                  lpM : do iM = -iL, iL
                    ind = ind + 1
                    val = radialVal * realTessY(iL, iM, diff, r)
                    if (ctx%isDensityCalc) val = val * val

                    if (ctx%isRealInput) then
                      do iEig = 1, size(eigVecsReal, dim=2)
                        valueReal(i1, i2, i3, iEig) = valueReal(i1, i2, i3, iEig) + val * eigVecsReal(ind, iEig)
                      end do
                    else ! Complex
                      do iEig = 1, size(eigVecsCmpl, dim=2)
                        valueCmpl(i1, i2, i3, iEig) = valueCmpl(i1, i2, i3, iEig) + val &
                            & * phases(iCell, kIndexes(iEig)) *  eigVecsCmpl(ind, iEig)
                      end do
                    end if
                  end do lpM
                end do lpOrb
              end do lpAtom
            end do lpCell
        end do lpI1
      end do lpI2
    end do lpI3
    !$omp end parallel do
  end subroutine evaluateOMP

end module libwavegrid_molorb_parallel
