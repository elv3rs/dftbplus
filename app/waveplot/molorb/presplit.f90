!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------! 


module waveplot_molorb_parallel
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : imag
  use waveplot_slater, only : TSlaterOrbital, realTessY, getVerboseRadial
  use omp_lib, only : omp_is_initial_device, omp_get_num_devices
  implicit none  
  
  public :: evaluateParallel

contains

  !> Returns the values of several molecular orbitals on grids.
  !! Caveat: The flag tPeriodic decides if the complex or the real version is read/written for the
  !! various parameters.
  subroutine evaluateParallel(origin, gridVecs, eigVecsReal, eigVecsCmpl, nAtom, nOrb, coords,&
      & species, iStos, stos, tPeriodic, tReal, latVecs, recVecs2p, kPoints,&
      & kIndexes, nCell, cellVec, tAddDensities, valueReal, valueCmpl)

    !> Origin of the grid
    real(dp), intent(in) :: origin(:)

    !> Grid vectors
    real(dp), intent(in) :: gridVecs(:,:)

    !> Real eigenvectors, or null-array
    real(dp), intent(in) :: eigVecsReal(:,:)

    !> Complex eigenvectors, or null-array
    complex(dp), intent(in) :: eigVecsCmpl(:,:)

    !> Nr. of atoms
    integer, intent(in) :: nAtom

    !> Nr. of orbitals
    integer, intent(in) :: nOrb

    !> Coordinates of the atoms
    real(dp), intent(in) :: coords(:,:,:)

    !> Species for each atom
    integer, intent(in) :: species(:)

    !> Starting position of the STOs for each species
    integer, intent(in) :: iStos(:)

    !> Array containing the STOs
    type(TSlaterOrbital), intent(in) :: stos(:)

    !> If the system is periodic
    logical, intent(in) :: tPeriodic

    !> If the system is real
    logical, intent(in) :: tReal

    !> Lattice vectors or null-array
    real(dp), intent(in) :: latVecs(:,:)

    !> Reciprocal vectors divided by 2pi (periodic) or a null-array (molecular)
    real(dp), intent(in) :: recVecs2p(:,:)

    !> Kpoints or null-array
    real(dp), intent(in) :: kPoints(:,:)

    !> Index of the k-points for each orbital in KPoints
    integer, intent(in) :: kIndexes(:)

    !> Nr. of cells to consider
    integer, intent(in) :: nCell

    !> Translation vector of the considered cells
    real(dp), intent(in) :: cellVec(:,:)

    !> If densities should be added instead of wave funcs
    logical, intent(in) :: tAddDensities

    !> Contains the real grid on exit
    real(dp), intent(out) :: valueReal(:,:,:,:)

    !> Contains the complex grid on exit
    complex(dp), intent(out) :: valueCmpl(:,:,:,:)

    !!Complex
    complex(dp), allocatable :: atomOrbValCmpl(:)
    complex(dp) :: phases(nCell, size(kPoints, dim=2))
    integer :: nPoints(4)



    real(dp), allocatable :: atomOrbValReal(:)
    integer :: i1, i2, i3

    !! SOA for stos contents
    integer, allocatable :: angMoms(:)
    real(dp), allocatable :: cutoffs(:)

    integer, allocatable :: sto_nPows(:)
    integer, allocatable :: sto_nAlphas(:)
    real(dp), allocatable :: sto_coeffs(:,:,:)
    real(dp), allocatable :: sto_alphas(:,:)
    integer :: maxNPows, maxNAlphas

    !! Create arrays out of sto object data to simplify OMP parallelization
    block
      integer ::  iOrb = 0

      allocate(angMoms(size(stos)), source=0)
      allocate(sto_nPows(size(stos)), source=0)
      allocate(sto_nAlphas(size(stos)), source=0)
      allocate(cutoffs(size(stos)), source=0.0_dp)


      do iOrb = 1, size(stos)
        angMoms(iOrb) = stos(iOrb)%angMom
        cutoffs(iOrb) = stos(iOrb)%cutoff
        sto_nPows(iOrb) = stos(iOrb)%nPow
        sto_nAlphas(iOrb) = stos(iOrb)%nAlpha
      end do


      !! Determine max power/ coefficient array size
      maxNPows = maxval(sto_nPows)
      maxNAlphas = maxval(sto_nAlphas)

      allocate(sto_coeffs(maxNPows, maxNAlphas, size(stos)))
      allocate(sto_alphas(maxNAlphas, size(stos)))

      ! Copy data from the sto objects to the arrays
      do iOrb = 1, size(stos)
        sto_coeffs(1:sto_nPows(iOrb), 1:sto_nAlphas(iOrb), iOrb) = stos(iOrb)%aa
        sto_alphas(1:sto_nAlphas(iOrb), iOrb) = stos(iOrb)%alpha
      end do


      do iOrb = 1, size(stos)
        print *, 'Orbital:', iOrb, 'Shape(aa):', shape(stos(iOrb)%aa), &
          'Expected Shape:', sto_nPows(iOrb), sto_nAlphas(iOrb)
        print *, 'Orbital:', iOrb, 'Size(alpha):', size(stos(iOrb)%alpha), &
          'Expected Size:', sto_nAlphas(iOrb)
      end do

    end block


    ! Array for the contribution of each orbital (and its periodic images)
    allocate(atomOrbValReal(nOrb))
    nPoints = shape(valueReal)
    valueReal(:,:,:,:) = 0.0_dp


    print *, "Devices:", omp_get_num_devices()
    print *, "maxNPows:", maxNPows, "maxNAlphas:", maxNAlphas
    print *, "Max sto_nPows in data:", maxval(sto_nPows)
    print *, "iStos", iStos
    print *, "Stos:", size(stos), "nOrb:", nOrb
    print *, "species", species


    print *, "sto_alphas:", sto_alphas
    print *, "sto_nalphas:", sto_nAlphas

    !$omp target teams distribute parallel do collapse(3) &
    !$omp&    map(to: gridVecs, origin,  &
    !$omp&                  nCell, nAtom, species, coords, cutoffs, angMoms, iStos, &
    !$omp&                  sto_nPows, sto_nAlphas, sto_coeffs, sto_alphas, &
    !$omp&                  maxNAlphas, eigVecsReal) &
    !$omp&    map(from: valueReal) 
    ! Loop over all grid points
    lpI3: do i3 = 1, nPoints(3)
      lpI2: do i2 = 1, nPoints(2)
        lpI1: do i1 = 1, nPoints(1)
          block
            real(dp) :: sto_tmp_pows(16)
            real(dp) :: sto_tmp_rexp
            real(dp) :: curCoords(3,3), xyz(3), diff(3)
            real(dp) :: xx, val, radialVal
            integer :: ind, iEig, iAtom, iOrb, iM, iSpecies, iL, iCell, ii, jj


            curCoords(:, 3) = real(i3 - 1, dp) * gridVecs(:, 3)
            curCoords(:, 2) = real(i2 - 1, dp) * gridVecs(:, 2)
            curCoords(:, 1) = real(i1 - 1, dp) * gridVecs(:, 1)
            xyz(:) = sum(curCoords, dim=2) + origin

            ! Get contribution from every atom in every cell for current point
            lpCell: do iCell = 1, nCell
              ind = 0
              lpAtom: do iAtom = 1, nAtom
                iSpecies = species(iAtom)
                diff(:) = xyz - coords(:, iAtom, iCell)
                xx = norm2(diff)
                lpOrb: do iOrb = iStos(iSpecies), iStos(iSpecies + 1) - 1
                  iL = angMoms(iOrb)
                  ! Calculate wave function only if atom is inside the cutoff
                  if (xx > cutoffs(iOrb)) then
                    ! Skip this orbital
                    ind = ind + 2*iL + 1
                    cycle lpOrb
                  end if



                  ! Avoid 0.0**0 as it may lead to arithmetic exception
                  if (iL == 0 .and. xx < epsilon(1.0_dp)) then
                     sto_tmp_rexp = 1.0_dp
                  else
                    sto_tmp_rexp = xx**iL
                  end if

                  ! Compute radial powers
                  do ii = 1, sto_nPows(iOrb)
                    sto_tmp_pows(ii) = sto_tmp_rexp
                    sto_tmp_rexp = sto_tmp_rexp * xx
                  end do

                  radialVal = 0.0_dp
                  do ii = 1, sto_nAlphas(iOrb)
                    sto_tmp_rexp = 0.0_dp
                    do jj = 1, sto_nPows(iOrb)
                      sto_tmp_rexp = sto_tmp_rexp + sto_coeffs(jj, ii, iOrb) * sto_tmp_pows(jj)
                    end do
                    radialVal = radialVal + sto_tmp_rexp * exp((-0.5_dp)*xx) !*  exp(sto_alphas(ii, iOrb) * xx) 
                  end do
                  !---- End inlined radial calculation ---!

                  lpM : do iM = -iL, iL
                    ind = ind + 1
                    val =  radialVal !* realTessY(iL, iM, diff, xx)
                    do iEig = 1, nPoints(4)
                      valueReal(i1, i2, i3, iEig) = valueReal(i1, i2, i3, iEig) + val * eigVecsReal(ind, iEig)
                    end do
                  end do lpM
                end do lpOrb
              end do lpAtom
            end do lpCell
          end block
        end do lpI1
      end do lpI2
    end do lpI3
    !$omp end target teams distribute parallel do

  end subroutine evaluateParallel
end module waveplot_molorb_parallel
