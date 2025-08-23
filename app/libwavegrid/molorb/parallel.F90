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
#:if WITH_CUDA
  use libwavegrid_molorb_offloaded, only : prepareGPUCoefficients, evaluateCuda
#:endif
  use libwavegrid_slater, only : realTessY, getRadial
  use omp_lib, only : omp_is_initial_device, omp_get_num_devices
  implicit none
  private
  !> Max powers expected in STO basis.
  integer, parameter :: MAX_STO_POWS = 16

  public :: evaluateParallel

contains


  !> Returns the values of several molecular orbitals on grids.
  !> This dispatches to either CPU / GPU implementation, and handles total Charge calculation using occupationVec if present.
  subroutine evaluateParallel(system, periodic, kIndexes, phases, basis, &
      & ctx, eigVecsReal, eigVecsCmpl,  valueReal, valueCmpl, occupationVec)

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
        print *, "GPU offloaded molorb."
        ! GPU implementation passes occupation information by baking their sqrt into the eigenvectors
        call prepareGPUCoefficients(ctx, eigVecsReal, eigVecsCmpl, occupationVec, coeffVecReal, coeffVecCmpl)

        call evaluateCuda(system, basis, periodic, kIndexes, phases, ctx, &
            & coeffVecReal, coeffVecCmpl, valueReal, valueCmpl)
      #:else
        call error("Libwavegrid: GPU offloaded molorb requested, but compiled without CUDA support.")
      #:endif
    else ! CPU implementation
      #:if WITH_OMP
        print *, "OMP parallel CPU molorb."
      #:else
        print *, "Serial CPU molorb."
      #:endif
      if (.not. ctx%calcTotalChrg) then
        call evaluateOMP(system, basis, periodic, kIndexes, phases, ctx, &
            & eigVecsReal, eigVecsCmpl, valueReal, valueCmpl)
      else
        ! Number of eigenvectors to calculate at once in a chunk.
        ! We need a function to query free RAM to dynamically size this.
        ! TODO: Additionally, expose to user.
        nEigsPerChunk = -1

        nEigs = size(occupationVec)
        if (nEigsPerChunk <= 0 .or. nEigsPerChunk > nEigs) then
          nEigsPerChunk = nEigs
        end if
        if (ctx%isRealInput) then
          allocate(bufferReal(size(valueReal, 1), size(valueReal, 2), size(valueReal, 3), nEigsPerChunk))
        else ! Complex input
          allocate(bufferCmpl(size(valueReal, 1), size(valueReal, 2), size(valueReal, 3), nEigsPerChunk))
        end if

        ! Zero accumulator
        valueReal(:, :, :, 1) = 0.0_dp

        lpChunk: do iStart = 1, nEigs, nEigsPerChunk
          iEnd = min(iStart + nEigsPerChunk - 1, nEigs)
          nChunk = iEnd - iStart + 1
          print *, "Processing", nEigsPerChunk, "eigenvectors from", iStart, "to", iEnd, "of", nEigs

          if (ctx%isRealInput) then
            call evaluateOMP(system, basis, periodic, kIndexes, phases, ctx, &
                & eigVecsReal(:, iStart:iEnd), eigVecsCmpl, bufferReal, valueCmpl)

            do iEigInChunk = 1, nChunk
              iEig = iStart + iEigInChunk - 1
              valueReal(:,:,:,1) = valueReal(:,:,:,1) + bufferReal(:,:,:,iEigInChunk)**2 * occupationVec(iEig)
            end do
          else ! Complex input
            call evaluateOMP(system, basis, periodic, kIndexes(iStart:iEnd), phases, ctx, &
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


  subroutine evaluateOMP(system, basis, periodic, kIndexes, phases, ctx, &
      & eigVecsReal, eigVecsCmpl, valueReal, valueCmpl)

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
    integer :: nPoints(4)

    if(ctx%isRealInput) then
      valueReal = 0.0_dp
      nPoints = shape(valueReal)
    else 
      valueCmpl = 0.0_dp
      nPoints = shape(valueCmpl)
    end if

    @:ASSERT(basis%maxNPows <= MAX_STO_POWS)

    !$omp parallel do collapse(3) &
    !$omp&    private(i1, i2, i3, iCell, iAtom, iOrb, iEig, iL, iM, xyz, frac, diff, &
    !$omp&              r, val, radialVal, tmp_pows, tmp_rexp, ind, iSpecies, rSq) &
    !$omp&    shared(ctx, system, basis, periodic, phases, kIndexes, &
    !$omp&              eigVecsReal, eigVecsCmpl, valueReal, valueCmpl)
    lpI3: do i3 = 1, nPoints(3)
      lpI2: do i2 = 1, nPoints(2)
        lpI1: do i1 = 1, nPoints(1)
            xyz(:) = system%origin(:) + real(i1 - 1, dp) * system%gridVecs(:, 1) &
                                    & + real(i2 - 1, dp) * system%gridVecs(:, 2) &
                                    & + real(i3 - 1, dp) * system%gridVecs(:, 3)

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
                    if (ctx%calcAtomicDensity) val = val * val

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
