!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to calculate the value of one or more molecular orbitals composed from STOs on
!! an equidistant grid.
module libwavegrid_molorb
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_boundarycond, only : TBoundaryConds
  use dftbp_dftb_periodic, only : getCellTranslations
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_typegeometry, only : TGeometry
  use libwavegrid_slater, only : TSlaterOrbital, realTessY
  use libwavegrid_molorb_parallel, only: evaluateParallel

  implicit none

  private


  !> Data type containing information about the basis for a species.
  type TSpeciesBasis

    !> Atomic number of the species
    integer :: atomicNumber

    !> Nr. of orbitals
    integer :: nOrb

    !> STO for each orbital
    type(TSlaterOrbital), allocatable :: stos(:)

  end type TSpeciesBasis


  !> Data type containing information for molecular orbital calculator.
  type TMolecularOrbital

    !> Nr. of atoms
    integer :: nAtom

    !> Nr. of species
    integer :: nSpecies

    !> Species of each atom
    integer, allocatable :: species(:)

    !> Index array for STOs
    integer, allocatable :: iStos(:)

    !> All STOs sequentially
    type(TSlaterOrbital), allocatable :: stos(:)

    !> Nr. of orbitals in the system
    integer :: nOrb

    !> If system is periodic
    logical :: tPeriodic

    !> Lattice vectors
    real(dp), allocatable :: latVecs(:,:)

    !> Reciprocal vectors divided by 2pi
    real(dp), allocatable :: recVecs2p(:,:)

    !> Cell shift vectors
    real(dp), allocatable :: cellVec(:,:)

    !> Nr. of cell shift vectors
    integer :: nCell

    !> Coordinates in all cells
    real(dp), allocatable :: coords(:,:,:)

    !> If it is initialised
    logical :: tInitialised = .false.

  end type TMolecularOrbital


  !> Returns the value of one or more molecular orbitals on a grid
  interface getValue
    module procedure TMolecularOrbital_getValue_real
    module procedure TMolecularOrbital_getValue_cmpl
  end interface

  public :: TSpeciesBasis
  public :: TMolecularOrbital, TMolecularOrbital_init, getValue

contains


  !> Initialises MolecularOrbital instance.
  subroutine TMolecularOrbital_init(this, geometry, boundaryCond, basis)

    !> Molecular Orbital
    type(TMolecularOrbital), intent(out) :: this

    !> Geometrical information.
    type(TGeometry), intent(in) :: geometry

    !> Boundary condition
    type(TBoundaryConds), intent(in) :: boundaryCond

    !> Basis for each species.
    type(TSpeciesBasis), intent(in) :: basis(:)

    integer :: nOrb, ii, jj, ind, iSp
    real(dp) :: mCutoff
    real(dp), allocatable :: rCellVec(:,:)

    @:ASSERT(.not. this%tInitialised)
    @:ASSERT(geometry%nSpecies == size(basis))

    this%nAtom = geometry%nAtom
    this%nSpecies = geometry%nSpecies
    allocate(this%species(this%nAtom))
    this%species(:) = geometry%species

    ! Create sequential list of STOs
    nOrb = 0
    do ii = 1, this%nSpecies
      nOrb = nOrb + (basis(ii)%nOrb)
    end do
    allocate(this%iStos(this%nSpecies+1))
    allocate(this%stos(nOrb))

    mCutoff = -1.0_dp

    ind = 1
    do ii = 1, this%nSpecies
      this%iStos(ii) = ind
      nOrb = basis(ii)%nOrb
      this%stos(ind:ind+nOrb-1) = basis(ii)%stos(1:nOrb)
      do jj = 1, basis(ii)%nOrb
        mCutoff = max(mCutoff, basis(ii)%stos(jj)%cutoff)
      end do
      ind = ind + nOrb
    end do
    this%iStos(ii) = ind

    ! Count all orbitals (including m-dependence)
    nOrb = 0
    do ii = 1, this%nAtom
      iSp = this%species(ii)
      do jj = 1, basis(iSp)%nOrb
        nOrb = nOrb + 1 + 2 * basis(iSp)%stos(jj)%angMom
      end do
    end do
    this%nOrb = nOrb

    ! Get cells to look for when adding STOs from periodic images
    this%tPeriodic = geometry%tPeriodic
    if (this%tPeriodic) then
      allocate(this%latVecs(3, 3))
      allocate(this%recVecs2p(3, 3))
      this%latVecs(:,:) = geometry%latVecs
      call invert33(this%recVecs2p, this%latVecs)
      this%recVecs2p(:,:) = reshape(this%recVecs2p, [3, 3], order=[2, 1])
      call getCellTranslations(this%cellVec, rCellVec, this%latVecs, this%recVecs2p, mCutoff)
      this%nCell = size(this%cellVec,dim=2)
    else
      allocate(this%latVecs(3, 0))
      allocate(this%recVecs2p(3, 0))
      allocate(this%cellVec(3, 1))
      this%cellVec(:,:) = 0.0_dp
      allocate(rCellVec(3, 1))
      rCellVec(:,:) = 0.0_dp
      this%nCell = 1
    end if

    ! Create coordinates for central cell and periodic images
    allocate(this%coords(3, this%nAtom, this%nCell))
    this%coords(:,:,1) = geometry%coords
    call boundaryCond%foldCoordsToCell(this%coords(:,:,1), this%latVecs)
    if (this%tPeriodic) then
      do ii = 2, this%nCell
        do jj = 1, this%nAtom
          this%coords(:, jj, ii) = this%coords(:, jj, 1) + rCellVec(:, ii)
        end do
      end do
    end if

    this%tInitialised = .true.

  end subroutine TMolecularOrbital_init


  !> Returns molecular orbitals on a grid.
  subroutine TMolecularOrbital_getValue_real(this, origin, gridVecs, eigVecsReal, &
      & valueOnGrid, addDensities)

    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this

    !> Origin of the grid
    real(dp), intent(in) :: origin(:)

    !> Grid vectors
    real(dp), intent(in) :: gridVecs(:,:)

    !> Summation coefficients for the STOs
    real(dp), intent(in) :: eigVecsReal(:,:)

    !> Molecular orbitals on a grid
    real(dp), intent(out) :: valueOnGrid(:,:,:,:)

    !> Add densities instead of wave functions
    logical, intent(in), optional :: addDensities

    !> timing variables
    integer(dp) :: startTime, endTime, clockRate

    real(dp) :: kPoints(3, 0)
    integer :: kIndexes(0)
    complex(dp) :: valueCmpl(0, 0, 0, 0)
    complex(dp) :: eigVecsCmpl(0, 0)
    logical :: tAddDensities

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(origin) == 3)
    @:ASSERT(all(shape(gridVecs) == [3, 3]))
    @:ASSERT(size(eigVecsReal, dim=1) == this%nOrb)
    @:ASSERT(all(shape(valueOnGrid) > [1, 1, 1, 0]))
    @:ASSERT(size(eigVecsReal, dim=2) == size(valueOnGrid, dim=4))

    if (present(addDensities)) then
      tAddDensities = addDensities
    else
      tAddDensities = .false.
    end if

    call system_clock(count_rate=clockRate, count=startTime)

    call evaluateParallel(origin, gridVecs, eigVecsReal, eigVecsCmpl, this%nAtom, this%nOrb,&
        & this%coords, this%species, this%iStos, this%stos,&
        & this%tPeriodic, .true., this%latVecs, this%recVecs2p, kPoints, kIndexes, this%nCell,&
        & this%cellVec, tAddDensities, valueOnGrid, valueCmpl)

  
    call system_clock(count=endTime)

    print *, "MolorbTime:", real(endTime - startTime, dp) / clockRate
  end subroutine TMolecularOrbital_getValue_real


  !> Returns molecular orbitals on a grid.
  subroutine TMolecularOrbital_getValue_cmpl(this, origin, gridVecs, eigVecsCmpl, kPoints,&
      & kIndexes, valueOnGrid)

    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this

    !> Origin of the grid
    real(dp), intent(in) :: origin(:)

    !> Grid vectors
    real(dp), intent(in) :: gridVecs(:,:)

    !> Summation coefficients for the STOs
    complex(dp), intent(in) :: eigVecsCmpl(:,:)

    !> Array of k-points
    real(dp), intent(in) :: kPoints(:,:)

    !> Index of the k-points in kPoints for every mol.orbital
    integer, intent(in) :: kIndexes(:)

    !> Molecular orbitals on grid on exit.
    complex(dp), intent(out) :: valueOnGrid(:,:,:,:)

    real(dp) :: valueReal(0,0,0,0)
    real(dp) :: eigVecsReal(0,0)
    logical :: tAddDensities = .false.

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(origin) == 3)
    @:ASSERT(all(shape(gridVecs) == [3, 3]))
    @:ASSERT(size(eigVecsCmpl, dim=1) == this%nOrb)
    @:ASSERT(all(shape(valueOnGrid) > [0, 0, 0, 0]))
    @:ASSERT(size(eigVecsCmpl, dim=2) == size(valueOnGrid, dim=4))
    @:ASSERT(size(kPoints, dim=1) == 3)
    @:ASSERT(size(kPoints, dim=2) > 0)
    @:ASSERT(size(kIndexes) == size(eigVecsCmpl, dim=2))
    @:ASSERT(maxval(kIndexes) <= size(kPoints, dim=2))
    @:ASSERT(minval(kIndexes) > 0)
  
    call evaluateParallel(origin, gridVecs, eigVecsReal, eigVecsCmpl, this%nAtom, this%nOrb,&
        & this%coords, this%species, this%iStos, this%stos,&
        & this%tPeriodic, .false., this%latVecs, this%recVecs2p, kPoints, kIndexes, this%nCell,&
        & this%cellVec, tAddDensities, valueReal, valueOnGrid)

  end subroutine TMolecularOrbital_getValue_cmpl




end module libwavegrid_molorb
