!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------! 


module waveplot_molorb_pointwise
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : imag
  use dftbp_dftb_boundarycond, only : TBoundaryConditions
  use dftbp_dftb_periodic, only : getCellTranslations
  use dftbp_common_status, only : TStatus
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_typegeometry, only : TGeometry
  use waveplot_slater, only : TSlaterOrbital, realTessY
  use dftbp_math_lapackroutines, only: gesv
  implicit none  
  
  public :: evaluatePointwise

contains

  !> Returns the values of several molecular orbitals on grids.
  !! Caveat: The flag tPeriodic decides if the complex or the real version is read/written for the
  !! various parameters.
  subroutine evaluatePointwise(origin, gridVecs, eigVecsReal, eigVecsCmpl, nAtom, nOrb, coords,&
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

    integer, allocatable :: angMoms(:)
    real(dp), allocatable :: cutoffs(:)
    real(dp) :: curCoords(3,3), xyz(3), diff(3), frac(3)
    real(dp) :: atomAllOrbVal(nOrb, nCell)
    logical :: nonZeroMask(nOrb), allZero
    integer :: nNonZero
    integer, target :: nonZeroIndContainer(nOrb)
    integer, pointer :: nonZeroIndices(:)
    real(dp), allocatable :: atomOrbValReal(:)
    complex(dp), allocatable :: atomOrbValCmpl(:)
    complex(dp) :: phases(nCell, size(kPoints, dim=2))
    real(dp) :: xx, val
    integer :: nPoints(4)
    integer :: ind, i1, i2, i3, iEig, iAtom, iOrb, iM, iSpecies, iL, iCell


    !! Create arrays out of sto object data for better cache locality
    allocate(angMoms(size(iStos)))
    allocate(cutoffs(size(iStos)))

    do iOrb = 1, size(stos)
      angMoms(iOrb) = stos(iOrb)%angMom
      cutoffs(iOrb) = stos(iOrb)%cutoff
    end do



    ! Array for the contribution of each orbital (and its periodic images)
    if (tReal) then
      allocate(atomOrbValReal(nOrb))
      nPoints = shape(valueReal)
    else
      allocate(atomOrbValCmpl(nOrb))
      nPoints = shape(valueCmpl)
    end if

    ! Phase factors for the periodic image cell. Note: This will be conjugated in the scalar product
    ! below. This is fine as, in contrast to what was published, DFTB+ uses implicitly exp(-ikr) as
    ! a phase factor, as the unpack routines assemble the lower triangular matrix with exp(ikr) as
    ! factor.
    phases(:,:) = exp(imag * matmul(transpose(cellVec), kPoints))

    ! Loop over all grid points
    lpI3: do i3 = 1, nPoints(3)
      curCoords(:, 3) = real(i3 - 1, dp) * gridVecs(:, 3)
      lpI2: do i2 = 1, nPoints(2)
        curCoords(:, 2) =  real(i2 - 1, dp) * gridVecs(:, 2)
        lpI1: do i1 = 1, nPoints(1)
          curCoords(:, 1) = real(i1 - 1, dp) * gridVecs(:, 1)
          xyz(:) = sum(curCoords, dim=2) + origin
          if (tPeriodic) then
            frac(:) = matmul(xyz, recVecs2p)
            xyz(:) = matmul(latVecs, frac - real(floor(frac), dp))
          end if
          ! Get contribution from every atom in every cell for current point
          allZero = .true.
          lpCell: do iCell = 1, nCell
            ind = 1
            lpAtom: do iAtom = 1, nAtom
              iSpecies = species(iAtom)
              diff(:) = xyz - coords(:, iAtom, iCell)
              xx = norm2(diff)
              lpOrb: do iOrb = iStos(iSpecies), iStos(iSpecies + 1) - 1
                iL = angMoms(iOrb)
                ! Calculate wave function only if atom is inside the cutoff
                if (xx <= cutoffs(iOrb)) then
                  allZero = .false.
                  call stos(iOrb)%getRadial(xx, val)
                  do iM = -iL, iL
                    atomAllOrbVal(ind, iCell) = val * realTessY(iL, iM, diff, xx)
                    ind = ind + 1
                  end do
                else
                  atomAllOrbVal(ind:ind+2*iL, iCell) = 0.0_dp
                  ind = ind + 2 * iL + 1
                end if
              end do lpOrb
            end do lpAtom
          end do lpCell

          if (allZero) then
            if (tReal) then
              valueReal(i1, i2, i3, :) = 0.0_dp
            else
              valueCmpl(i1, i2, i3, :) = 0.0_dp
            end if
            cycle lpI1
          end if

          ! Establish mask and index of nonzero elements
          nonZeroMask = any(atomAllOrbVal /= 0.0_dp, dim=2)
          nNonZero = 0
          do iOrb = 1, nOrb
            if (nonZeroMask(iOrb)) then
              nNonZero = nNonZero + 1
              nonZeroIndContainer(nNonZero) = iOrb
            end if
          end do
          nonZeroIndices => nonZeroIndContainer(1:nNonZero)

          ! Sum the contribution from all cells and multiply by the provided coefficients (usually
          ! the eigenvector)
          if (tReal) then
            if (tAddDensities) then
              atomAllOrbVal(:,:) = atomAllOrbVal**2
            end if
            atomOrbValReal(:) = sum(atomAllOrbVal, dim=2)
            do iEig = 1, nPoints(4)
              valueReal(i1, i2, i3, iEig) = dot_product(atomOrbValReal(nonZeroIndices),&
                  & eigVecsReal(nonZeroIndices, iEig))
            end do
          else
            ind = 0
            do iEig = 1, nPoints(4)
              if (kIndexes(iEig) /= ind) then
                ind = kIndexes(iEig)
                atomOrbValCmpl(nonZeroIndices) = (0.0_dp, 0.0_dp)
                do iCell = 1, nCell
                  atomOrbValCmpl(nonZeroIndices) = atomOrbValCmpl(nonZeroIndices)&
                      & + atomAllOrbVal(nonZeroIndices, iCell) * phases(iCell, ind)
                end do
              end if
              valueCmpl(i1, i2, i3, iEig) = dot_product(atomOrbValCmpl(nonZeroIndices),&
                  & eigVecsCmpl(nonZeroIndices, iEig))
            end do
          end if
        end do lpI1
      end do lpI2
    end do lpI3

  end subroutine evaluatePointwise
end module waveplot_molorb_pointwise
