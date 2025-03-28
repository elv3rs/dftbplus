!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to calculate the value of one or more molecular orbitals composed from STOs on
!! an equidistant grid.
module waveplot_molorb
  use waveplot_slater, only : getValue, realTessY, TSlaterOrbital
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : imag
  use dftbp_dftb_boundarycond, only : TBoundaryConds
  use dftbp_dftb_periodic, only : getCellTranslations
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_typegeometry, only : TGeometry
  implicit none

  private
  save


  !> Data type containing information about the basis for a species.
  type TSpeciesBasis

    !> Atomic number of the species
    integer :: atomicNumber

    !> Nr. of orbitals
    integer :: nOrb

    !> Angular momentum for each orbital
    integer, allocatable :: angMoms(:)

    !> Cutoff for each orbital
    real(dp), allocatable :: cutoffs(:)

    !> STO for each orbital
    type(TSlaterOrbital), allocatable :: stos(:)

    !> Occupation for each orbital
    real(dp), allocatable :: occupations(:)

  end type TSpeciesBasis


  !> Data type containing information for molecular orbital calculator.
  type TMolecularOrbital
    private

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

    !> Cutoff for each STO
    real(dp), allocatable :: cutoffs(:)

    !> Angular mometum for each STO
    integer, allocatable :: angMoms(:)

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

    integer :: nOrb
    integer :: ii, jj, ind, iSp
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
    allocate(this%cutoffs(nOrb))
    allocate(this%angMoms(nOrb))
    ind = 1
    do ii = 1, this%nSpecies
      this%iStos(ii) = ind
      nOrb = basis(ii)%nOrb
      this%stos(ind:ind+nOrb-1) = basis(ii)%stos(1:nOrb)
      this%cutoffs(ind:ind+nOrb-1) = basis(ii)%cutoffs(1:nOrb)
      this%angMoms(ind:ind+nOrb-1) = basis(ii)%angMoms(1:nOrb)
      ind = ind + nOrb
    end do
    this%iStos(ii) = ind

    ! Count all orbitals (including m-dependence)
    nOrb = 0
    do ii = 1, this%nAtom
      iSp = this%species(ii)
      nOrb = nOrb + sum(2*this%angMoms(this%iStos(iSp):this%iStos(iSp+1)-1)+1)
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
      mCutoff = maxval(this%cutoffs)
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
  subroutine TMolecularOrbital_getValue_real(this, origin, gridVecs, eigVecsReal, valueOnGrid,&
      & addDensities)

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

    real(dp), save :: kPoints(3, 0)
    integer, save :: kIndexes(0)
    complex(dp), save :: valueCmpl(0, 0, 0, 0)
    complex(dp), save :: eigVecsCmpl(0, 0)
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

    call localGetValue(origin, gridVecs, eigVecsReal, eigVecsCmpl, this%nAtom, this%nOrb,&
        & this%coords, this%species, this%cutoffs, this%iStos, this%angMoms, this%stos,&
        & this%tPeriodic, .true., this%latVecs, this%recVecs2p, kPoints, kIndexes, this%nCell,&
        & this%cellVec, tAddDensities, valueOnGrid, valueCmpl)

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

    real(dp), save :: valueReal(0,0,0,0)
    real(dp), save :: eigVecsReal(0,0)
    logical, save :: tAddDensities = .false.

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

    call localGetValue(origin, gridVecs, eigVecsReal, eigVecsCmpl, this%nAtom, this%nOrb,&
        & this%coords, this%species, this%cutoffs, this%iStos, this%angMoms, this%stos,&
        & this%tPeriodic, .false., this%latVecs, this%recVecs2p, kPoints, kIndexes, this%nCell,&
        & this%cellVec, tAddDensities, valueReal, valueOnGrid)

  end subroutine TMolecularOrbital_getValue_cmpl


  !> Returns the values of several molecular orbitals on grids.
  !! Caveat: The flag tPeriodic decides if the complex or the real version is read/written for the
  !! various parameters.
  subroutine localGetValue(origin, gridVecs, eigVecsReal, eigVecsCmpl, nAtom, nOrb, coords,&
      & species, cutoffs, iStos, angMoms, stos, tPeriodic, tReal, latVecs, recVecs2p, kPoints,&
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

    !> total Nr. of orbitals, i.e. ∑ atoms  ∑ orbitals
    integer, intent(in) :: nOrb

    !> Coordinates of the atoms
    real(dp), intent(in) :: coords(:,:,:)

    !> Species for each atom
    integer, intent(in) :: species(:)

    !> Cutoff for each STO
    real(dp), intent(in) :: cutoffs(:)

    !> Starting position of the STOs for each species
    integer, intent(in) :: iStos(:)

    !> Angular moment for each STO
    integer, intent(in) :: angMoms(:)

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

    real(dp) :: curCoords(3,3), xyz(3), diff(3), frac(3)
    complex(dp) :: phases(nCell, size(kPoints, dim=2))
    real(dp) :: xx, val
    integer :: nPoints(4)
    integer :: i1, i2, i3, iEig, iAtom, iOrb, iM, iSpecies, iL, iCell

    ! Cache variables

    ! Todo: This overrides targetResolution and resolutionFactor.
    integer, parameter :: resolutionFactorLead = 5

    ! Target subdivision resolution.
    ! Todo: Figure out what unit we are using here. A?
    real(dp) :: targetResolution(3)
    
    ! Resolution factors calculated to meet the targetResolution requirements
    integer :: resolutionFactor(3)

    ! Wavefunction cache overview:
    ! Wavefunctions are evaluated across a finer grid (resolutionFactor as many points as the main one).
    ! When aligned, points can simply be read using a stride of resolutionFactor.
    ! Main feature:
    ! The X-coordinate is subdivided into <resolutionFactor> chunked parts to ensure a contiguous memory layout.
    ! Y and Z have been subdivided as well, though only for convenience.
    ! Indices : [X, cacheInd, XChunked, Y, YChunked, Z, ZChunked]
    real(dp), allocatable, save :: wavefunctionCache(:,:,:,:,:,:,:)

    integer, allocatable, save, target :: chunkedIndices(:,:,:)
    integer, allocatable, save, target :: sliceIndicesMain(:,:,:,:)
    integer, allocatable, save, target :: sliceIndicesCache(:,:,:,:)

    logical, save :: tCacheInitialised = .false.

    ! cacheInd is the index of the orbital in the cache as a running sum of (iSpecies, iOrb, iM)
    integer :: cacheInd, cacheSize

    ! Cache grid vectors (x_fine, y_fine, z_fine)
    real(dp) :: cacheGridVecs(3,3)

    ! Used for decomposing Atom position into cacheGrid Basis vecs
    real(dp) :: cacheBasis(3,3)
    real(dp) :: pos(3, 1)
    
    ! Aligning the cache to the main grid for each atom position
    integer :: nUniqueOrb

    ! Mapping array describing cache layout: (iM, iOrb) -> cacheInd.
    integer, allocatable :: cacheIndexMap(:,:)


    integer :: i1Chunked, i2Chunked, i3Chunked
    integer :: coeffInd
    integer :: nPointsHalved(3)



    integer, pointer :: iChunk(:)
    integer, pointer :: iMain(:,:)
    integer, pointer :: iCache(:,:)

    logical, parameter :: tSparseCache = .true.
    real(dp) :: expectedSizeMB

    !---------------------------------
    
    nUniqueOrb = SIZE(angMoms)
    allocate(cacheIndexMap(-MAXVAL(angMoms):MAXVAL(angMoms), nUniqueOrb))

    cacheIndexMap(:,:) = -1 ! Ensure we dont go oob

    ! Calculate nr. of cache entries <cacheSize> and populate mapping array <cacheIndexMap>
    cacheInd = 1
    print *, "nUniqueOrb", nUniqueOrb
    do iOrb = 1, nUniqueOrb
      iL = angMoms(iOrb)
      do iM = -iL, iL
        cacheInd = cacheInd + 1
        cacheIndexMap(iM, iOrb) = cacheInd
      end do
    end do
    cacheSize = cacheInd
    
    ! Main grid size
    if (tReal) then
      nPoints = shape(valueReal)
    else
      nPoints = shape(valueCmpl)
      stop "Complex not implemented yet"
    end if

    targetResolution = norm2(gridVecs, dim=1) / resolutionFactorLead
    resolutionFactor = norm2(gridVecs, dim=1) / targetResolution
    
    ! Fine Grid Basis Vectors
    cacheGridVecs(:,1) = gridVecs(:,1) / resolutionFactor(1)
    cacheGridVecs(:,2) = gridVecs(:,2) / resolutionFactor(2)
    cacheGridVecs(:,3) = gridVecs(:,3) / resolutionFactor(3)

    ! Print cutoff size. We assume that maxval == minval.
    if (maxval(cutoffs) /= minval(cutoffs)) then
      print *, "Warn: Different cutoffs (max/min):", maxval(cutoffs), minval(cutoffs)
    end if


    ! Choose Cache size as min(cutoff, nPoints). This assumes that the basis is orthogonal.
    ! We might have to invert the Gram matrix to get the correct cutoffs.
    nPointsHalved(1) = MIN(int(MAXVAL(cutoffs) / norm2(cacheGridVecs(:,1))), nPoints(1))
    nPointsHalved(2) = MIN(int(MAXVAL(cutoffs) / norm2(cacheGridVecs(:,2))), nPoints(2))
    nPointsHalved(3) = MIN(int(MAXVAL(cutoffs) / norm2(cacheGridVecs(:,3))), nPoints(3))
    
    ! Todo: Calculate resolutionFactor to achieve a subdivision of targetResolution 1/A.
    !       -> Figure out a sensible targetResolution for testing
    !       -> Allow Different resolution factors per dimension
    ! Todo: Extend cutoff estimation to non-orthogonal basis by building Gram matrix
    ! Todo: Implement Complex Version
    ! Todo: Compare with non-chunked version, decide if chunking is worth the added complexity
    ! Todo: Benchmark dense/sparse variant. Can we remove the sparse version or does it provide a significant speedup for tiny systems?






    print "(*(G0, 1X))", "Main Grid Dimensions", nPoints
    print "(*(G0, 1X))", " ->", size(valueReal), "elements,",  sizeof(valueReal) / 1000.0 / 1000.0, "MB"
    print "(*(G0, 1X))", "tSparseCache:", tSparseCache



    if (tCacheInitialised) then
      print "(*(G0, 1X))", "Reusing saved wavefunctions"
      goto 591 ! Temporary quick fix
    end if


    expectedSizeMB = 0.0_dp
    expectedSizeMB = storage_size(1.0_dp) * product(nPointsHalved)
    expectedSizeMB  = expectedSizeMB * 2**3 * resolutionFactor**3 * cacheSize
    expectedSizeMB  = expectedSizeMB / 8.0 / 1024.0 / 1024.0

    print "(*(G0, 1X))", "Allocating Cache Grid of dimensions", nPointsHalved * 2
    print "(*(G0, 1X))", "expected Cache Allocation Size", expectedSizeMB, "MB"
    if (expectedSizeMB > 8000) then
      stop "Expected cache array size exceeds 8GB"
    end if

    ! We define (0,0,0) to be the origin using Fortrans fancy arbitrary indices feature.
    allocate(wavefunctionCache(   -nPointsHalved(1):nPointsHalved(1), &
                                &  1:cacheSize, &
                                &  0:resolutionFactor-1, &
                                & -nPointsHalved(2):nPointsHalved(2),&
                                &  0:resolutionFactor-1, &
                                & -nPointsHalved(3):nPointsHalved(3), &
                                &  0:resolutionFactor-1))
    ! zero out the cache
    wavefunctionCache(:,:,:,:,:,:,:) = 0.0_dp
    print "(*(G0, 1X))", " ->", size(wavefunctionCache), "elements,",  sizeof(wavefunctionCache) / 1000.0 / 1000.0, "MB"


    ! Allocate indice arrays
    allocate(chunkedIndices(3, nAtom, nCell))
    allocate(sliceIndicesMain(3,2, nAtom, nCell))
    allocate(sliceIndicesCache(3,2, nAtom, nCell))

    print *, "Aligning Arrays..."
    ! Factored this out of the main loop to get chunked information for sparse caching
    do iCell = 1, nCell
      coeffInd = 1
      do iAtom = 1, nAtom
        iSpecies = species(iAtom)
        ! Determine Array Offsets by aligning the wavefunction cache, then clamping to array bounds.
        pos(:,1) = coords(:, iAtom, iCell) + origin
        cacheBasis(:,:) = cacheGridVecs(:,:)
        ! Decompose Atom position onto basis
        call gesv(cacheBasis, pos)
        ! Atom Offset in terms of cacheGridVecs now stored in pos.
        ! Choose closest chunk
        chunkedIndices(1, iAtom, iCell) = ABS (int(MOD(pos(1,1), 1.0_dp) * real(resolutionFactor, dp)))
        chunkedIndices(2, iAtom, iCell) = ABS (int(MOD(pos(2,1), 1.0_dp) * real(resolutionFactor, dp)))
        chunkedIndices(3, iAtom, iCell) = ABS (int(MOD(pos(3,1), 1.0_dp) * real(resolutionFactor, dp)))
        ! Align to main grid, clamp to bounds
        sliceIndicesMain(3,1, iAtom, iCell) = MAX(1, int(pos(3,1)) - nPointsHalved(3))
        sliceIndicesMain(3,2, iAtom, iCell) = MIN(nPoints(3), int(pos(3,1)) + nPointsHalved(3))

        sliceIndicesMain(2,1, iAtom, iCell) = MAX(1, int(pos(2,1)) - nPointsHalved(2))
        sliceIndicesMain(2,2, iAtom, iCell) = MIN(nPoints(2), int(pos(2,1)) + nPointsHalved(2))

        sliceIndicesMain(1,1, iAtom, iCell) = MAX(1, int(pos(1,1)) - nPointsHalved(1))
        sliceIndicesMain(1,2, iAtom, iCell) = MIN(nPoints(1), int(pos(1,1)) + nPointsHalved(1))

        sliceIndicesCache(:,1, iAtom, iCell) = sliceIndicesMain(:,1, iAtom, iCell) - int(pos(:,1))
        sliceIndicesCache(:,2, iAtom, iCell) = sliceIndicesMain(:,2, iAtom, iCell) - int(pos(:,1))
      end do
    end do



  


    print "(*(G0, 1X))",  "Caching", nUniqueOrb,  "wavefunctions:" 
    ! Loop over all wavefunctions and cache them on the fine grid
    lpOrb: do iOrb = 1, nUniqueOrb
      print "(*(G0, 1X))", " -> Caching orbital ", iOrb
      iL = angMoms(iOrb)
      ! For every Point on the fine grid:
      lpI3Chunked : do i3Chunked = 0, resolutionFactor-2
        if (tSparseCache) then
          if (.not. any(i3Chunked == chunkedIndices(3, :, :))) cycle
        end if
        lpI3 : do i3 = -nPointsHalved(3), nPointsHalved(3)
          curCoords(:, 3) = real(i3*resolutionFactor + i3Chunked, dp) * cacheGridVecs(:, 3)
          lpI2Chunked : do i2Chunked = 0, resolutionFactor-2
            if (tSparseCache) then
              if (.not. any(i2Chunked == chunkedIndices(2, :, :))) cycle
            end if
            lpI2 : do i2 = -nPointsHalved(2), nPointsHalved(2)
              curCoords(:, 2) = real(i2*resolutionFactor + i2Chunked, dp) * cacheGridVecs(:, 2)
              lpI1Chunked : do i1Chunked = 0, resolutionFactor-2
                if (tSparseCache) then
                  if (.not. any(i1Chunked == chunkedIndices(1, :, :))) cycle
                end if
                lpI1 : do i1 = -nPointsHalved(1), nPointsHalved(1)
                  curCoords(:, 1) = real(i1*resolutionFactor + i1Chunked, dp) * cacheGridVecs(:, 1)
                  
                  ! Calculate position
                  xyz(:) = sum(curCoords, dim=2)
                  if (tPeriodic) then
                    frac(:) = matmul(xyz, recVecs2p)
                    xyz(:) = matmul(latVecs, frac - real(floor(frac), dp))
                  end if
                  xx = norm2(xyz)
                  if (xx <= cutoffs(iOrb)) then
                    ! Get radial dependence
                    call getValue(stos(iOrb), xx, val)
                    lpIM : do iM = -iL, iL
                      cacheInd = cacheIndexMap(iM, iOrb)
                      ! Combine with angular dependence and add to cache
                      wavefunctionCache(i1, cacheInd, i1Chunked, i2, i2Chunked, i3, i3Chunked) = val * realTessY(iL, iM, diff, xx)
                    end do lpIM
                  end if
                end do lpI1
              end do lpI1Chunked
            end do lpI2
          end do lpI2Chunked
        end do lpI3
      end do lpI3Chunked
    end do lpOrb
    tCacheInitialised = .true.


591 print *, "Applying wavefunctions"
    ! Apply wavefunctions. For each atom, determine the offsets, then align the cache and apply using slicing.
    do iCell = 1, nCell
      coeffInd = 1
      do iAtom = 1, nAtom
        iSpecies = species(iAtom)
        print "(*(G0, 1X))", " -> Adding contribution of Atom no.", iAtom
        
        ! Load Alignment
        iChunk => chunkedIndices(:, iAtom, iCell)
        iMain => sliceIndicesMain(:, :, iAtom, iCell)
        iCache => sliceIndicesCache(:, :, iAtom, iCell)

        do iOrb = iStos(iSpecies), iStos(iSpecies + 1) - 1
          iL = angMoms(iOrb)
          do iM = -iL, iL
            cacheInd = cacheIndexMap(iM, iOrb)
            do iEig = 1, nPoints(4)
              if (tReal) then
                if (tAddDensities) then
                  ! Square the wavefunction
                  valueReal(iMain(1,1):iMain(1,2), iMain(2,1):iMain(2,2), iMain(3,1):iMain(3,2), iEig) = &
                      & valueReal(iMain(1,1):iMain(1,2), iMain(2,1):iMain(2,2), iMain(3,1):iMain(3,2), iEig) + &
                      & eigVecsReal(coeffInd, iEig) * wavefunctionCache(iCache(1,1):iCache(1,2), &
                      & cacheInd, iChunk(1), iCache(2,1):iCache(2,2), iChunk(2), &
                      & iCache(3,1):iCache(3,2), iChunk(3)) ** 2
                else
                  valueReal(iMain(1,1):iMain(1,2), iMain(2,1):iMain(2,2), iMain(3,1):iMain(3,2), iEig) = &
                      & valueReal(iMain(1,1):iMain(1,2), iMain(2,1):iMain(2,2), iMain(3,1):iMain(3,2), iEig) + &
                      & eigVecsReal(coeffInd, iEig) * wavefunctionCache(iCache(1,1):iCache(1,2), &
                      & cacheInd, iChunk(1), iCache(2,1):iCache(2,2), iChunk(2), &
                      & iCache(3,1):iCache(3,2), iChunk(3))
                end if
              else
                !TODO: Implement complex version
                stop "TODO: Complex not implemented yet"
              end if
            end do
            coeffInd = coeffInd + 1
          end do
        end do
        end do
      end do
  end subroutine localGetValue

end module waveplot_molorb
