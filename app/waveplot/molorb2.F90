!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

! Notes:
! Print cutoff size. We assume that maxval == minval.
! One would expect H to be smaller than e.g. Pb.
! https://github.com/dftbparams/matsci/releases has differing cutoffs including 5,6, and 8.
! That is not neglibible, storing the small orbitals in a large buffer would be a waste of memory.
! Idea: Attach the cache buffers to a derived type, one for each orbital.
! We would then have \sum_{iAtom} nOrbs_iAtom objects.
! This would allow different cutoffs since we no longer store everything in the same array.
! They can have an attribute isInitialised and an initialisation function that populates the array.
! We could hide the chunking by having a subroutine that returns a view (pointer) to the correct subgrid.
! Should we also move the indices-calculation into a typebound procedure?
! Or We could have a subroutine that directly adds the contribution of the wavefunction.
! But the orbitals need to multiplied by the eigenvectors, and occasionally squared.
! Can we even call those density function thingies wavefunctions? Probably not.
! How about OrbitalCache%cache?
! Todo: Figure out how waveplot is to be used as a library in the future.
! General Todo list: 
! Todo: Acquire a real-world example to base decisions on.
!        -> Figure out a sensible targetGridDistance
!        -> Base on available memory?
! Todo: Implement Complex Version
! Todo: Check if results align with unmodified version
!       -> Run both, compare
!       -> Investigate how subdivision affects accuracy
! Todo: Move targetGridDistance setting etc. to waveplot in hsd
! Use Fypp for complex/real distinction?
! What was the name of that fancy LLM page?
! Ask how compiler handels realTessY
! Does the compiler inline the call to getValueExplicit?
! Quit trying to apply premature optimisations
! Try to get the LSP working. Will require running fypp in the background.



#:include 'common.fypp'

!> Contains routines to calculate the value of one or more molecular orbitals composed from STOs on
!! an equidistant grid.
module waveplot_molorb2
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : imag
  use dftbp_dftb_boundarycond, only : TBoundaryConditions
  use dftbp_dftb_periodic, only : getCellTranslations
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_typegeometry, only : TGeometry
  use waveplot_slater, only : TSlaterOrbital, realTessY, getValue
  use dftbp_math_lapackroutines, only: gesv
  implicit none  
  
  public :: localGetValue


  ! Wavefunction cache overview:
  ! Wavefunctions are evaluated across a finer grid (subdivisionFactor as many points as the main one).
  ! When aligned, points could simply be read using a stride of subdivisionFactor.
  ! The X-coordinate is subdivided into <subdivisionFactor> chunked parts to ensure a contiguous memory layout.
  ! Y and Z have been subdivided as well, though only for convenience.
  ! Indices : [X, XChunked, Y, YChunked, Z, ZChunked,  cacheInd ]
  type :: TOrbitalCache
    private
    !> Cache for a single orbital.
    ! Layout: (X, XChunked, Y, YChunked, Z, ZChunked, iM)
    !TODO: Dont forget to deallocate me!
    real(dp), pointer :: cache(:,:,:,:,:,:,:) => null()

    !> Has the cache been initialised?
    ! (Can also be checked by looking at wavefunctionCache's allocation status)
    logical :: isInitialised = .false.

    !> The STO that is cached
    type(TSlaterOrbital) :: sto    

    !> Angular momentum (L) of the orbital
    integer :: angMom

    !> The grids basis vectors
    real(dp) :: gridVecs(3,3)

    !> How many shifted copies of the main grid are stored in the cache
    integer :: subdivisionFactor(3)

    !> Size of cache from center to cutoff
    integer, public :: nPointsHalved(3)
  contains
    procedure :: initialise => TOrbitalCache_initialise
    procedure :: access => TOrbitalCache_access
    procedure :: align => TOrbitalCache_align
  end type TOrbitalCache




contains

  !> Initialises the cache for a given STO and its angular momentum L.
  subroutine TOrbitalCache_initialise(this, sto, angMom, gridVecs, subdivisionFactor)
    class(TOrbitalCache), intent(inout) :: this
    type(TSlaterOrbital), intent(in) :: sto
    real(dp), intent(in) :: gridVecs(3,3)
    integer, intent(in) :: subdivisionFactor(3)
    integer, intent(in) :: angMom

    real(dp) :: curCoords(3,3), xx, val, xyz(3)
    integer :: i1, i2, i3, iM, iL, i1Chunked, i2Chunked, i3Chunked
    integer :: nPointsHalved(3)
    real(dp) :: expectedSizeMB
    !---------------------------------
    
    ! Skip initialisation if the cache is already populated. (Doesnt check if parameters changed)
    if (this%isInitialised) then
      print "(*(G0, 1X))", " -> Cache already initialised, skipping."
      return
    end if

    ! Todo: Sprinkle in some asserts
    this%sto = sto
    this%angMom = angMom
    this%gridVecs = gridVecs
    this%subdivisionFactor = subdivisionFactor
    iL = angMom

    ! Determine the size of the cache
    ! TODO: Currently assuming an orthogonal basis
    !       -> investigate using Gram matrix for correct cutoffs
    ! Alternative: Warn the user / stop if they provide a non-orthogonal basis

    nPointsHalved = ceiling(sto%cutoff / norm2(gridVecs, dim=1))
    this%nPointsHalved = nPointsHalved
    !print *, "Cutoff:", sto%cutoff
    !print *, "GridVec size:", norm2(gridVecs, dim=1)

    ! Allocate the cache
    expectedSizeMB = 0.0_dp
    expectedSizeMB = storage_size(1.0_dp) * product(nPointsHalved)
    expectedSizeMB = expectedSizeMB * 2**3 * product(subdivisionFactor) * (2 * iL + 1)
    expectedSizeMB = expectedSizeMB / 8.0 / 1024.0 / 1024.0
    print "(*(G0, 1X))", "Allocating Cache Grid of dimensions", nPointsHalved * 2
    print "(*(G0, 1X))", "Expected Cache Allocation Size", expectedSizeMB, "MB"
    if (expectedSizeMB > 8000) then
      stop "Expected cache array size exceeds 8GB"
    end if


    ! We define (0,0,0) to be the origin using Fortrans fancy arbitrary indices feature.
    allocate(this%cache(-nPointsHalved(1):nPointsHalved(1), 0:subdivisionFactor(1)-1, &
                      & -nPointsHalved(2):nPointsHalved(2), 0:subdivisionFactor(2)-1, &
                      & -nPointsHalved(3):nPointsHalved(3), 0:subdivisionFactor(3)-1, &
                      & -iL:iL), source=0.0_dp)
    !print "(*(G0, 1X))", " ->", size(wavefunctionCache), "elements,",  sizeof(wavefunctionCache) / 1000.0 / 1000.0, "MB"



    ! For every Point on the fine grid:
    lpI3Chunked : do i3Chunked = 0, subdivisionFactor(3)-1
      lpI3 : do i3 = -nPointsHalved(3), nPointsHalved(3)
        curCoords(:, 3) = real(i3 + i3Chunked / subdivisionFactor(3), dp) * gridVecs(:, 3)
        lpI2Chunked : do i2Chunked = 0, subdivisionFactor(2)-1
          lpI2 : do i2 = -nPointsHalved(2), nPointsHalved(2)
            curCoords(:, 2) = real(i2 + i2Chunked / subdivisionFactor(2), dp) * gridVecs(:, 2)
            lpI1Chunked : do i1Chunked = 0, subdivisionFactor(1)-1
              lpI1 : do i1 = -nPointsHalved(1), nPointsHalved(1)
                curCoords(:, 1) = real(i1 + i1Chunked / subdivisionFactor(1), dp) * gridVecs(:, 1)
                
                ! Calculate position
                xyz = sum(curCoords, dim=2)
                xx = norm2(xyz)

                if (xx <= sto%cutoff) then
                  ! Get radial dependence
                  call getValue(sto, xx, val)
                  lpIM : do iM = -iL, iL
                    ! Combine with angular dependence and add to cache
                    this%cache(i1, i1Chunked, i2, i2Chunked, i3, i3Chunked, iM) = val * realTessY(iL, iM, xyz, xx)
                  end do lpIM
                end if
              end do lpI1
            end do lpI1Chunked
          end do lpI2
        end do lpI2Chunked
      end do lpI3
    end do lpI3Chunked

    this%isInitialised = .true.

    ! Deallocate the radial dependence cache to free up some memory
    ! (Add a deallocation subroutine?)
    !deallocate(sto%gridValue)

  end subroutine TOrbitalCache_initialise

  !> Accesses the cache for a given angular momentum L and chunked indices.
  !! Returns a pointer to the requested subgrid.
  subroutine TOrbitalCache_access(this, cachePtr, iChunk, iM)
    class(TOrbitalCache), intent(in) :: this
    real(dp), pointer, intent(out)   :: cachePtr(:,:,:)
    integer, intent(in) :: iChunk(3)
    integer, intent(in) :: iM

    ! Access the cache using the chunked indices
    ! Fortran does not support simultaneous bound mapping and rank reduction :(
    cachePtr => this%cache(:, iChunk(1), :, iChunk(2), :, iChunk(3), iM)
    

  end subroutine TOrbitalCache_access

  !> Aligns the cache to the main grid.
  !! Returns the chunks to choose the correct subgrid,
  !! as well as the indices bounds for the main grid and an offset for the cache.
  subroutine TOrbitalCache_align(this, gridDims, shiftedPos, iOffset, iChunk, iMain)
      class(TOrbitalCache), intent(in) :: this
      integer, intent(in) :: gridDims(4)
      real(dp), intent(in) :: shiftedPos(:)
      integer, pointer, intent(out) :: iOffset(:)
      integer, pointer, intent(out) :: iChunk(:)
      integer, pointer, intent(out) :: iMain(:,:)

      ! Used for decomposing Atom position into cacheGrid Basis vecs using gesv
      real(dp) :: tmpBasis(3,3), pos(3,1)
      integer :: ii
      

      pos(:,1) = shiftedPos
      tmpBasis(:,:) = this%gridVecs
      ! Get the atom position in terms of the basis <gridVecs>, stored in pos
      call gesv(tmpBasis, pos)

      ! Cache Indices are derived from the main indices, and:
      ! -> offset by the atom position
      iOffset(:) = int(pos(:,1))
      ! -> shifted by (0...subdivisionFactor-1) to select the closest subgrid (chunk)
      ! TODO: By replacing int with nint we get the closest one, but would need to adjust the main indices.
      iChunk(:) = modulo(  int( pos(:,1) * this%subdivisionFactor(:) ), this%subdivisionFactor(:)   )

      ! Lower Main Indices need to include
      ! -> start of main grid (1)
      ! -> start of cache grid (atom offset - half cache size)
      iMain(:, 1) = max(1, iOffset(:) - this%nPointsHalved(:))
      ! Upper Main Indices need to include
      ! -> end of main grid (nPoints)
      ! -> end of cache grid (atom offset + half cache size)
      iMain(:, 2) = min(gridDims(:3), iOffset(:) + this%nPointsHalved(:))

      ! Quick Fix to counter lost bound mapping
      iOffset(:) =  this%nPointsHalved(:)- iOffset(:) 

      print "(*(G0, 1X))", " -> Atom Position in Basis vectors:", pos(:,1)
      print "(*(G0, 1X))", " -> Atom Position in Grid:", iOffset(:), "Chunk:", iChunk(:)
      print "(*(G0, 1X))", "Indices Mapping Main -> Cache:"
      do ii = 1, 3
        print "(*(G0, 1X))", " o", Char(ii +87), iMain(ii,1), ":" , &
                        & iMain(ii,2), "->",&
                        & iMain(ii,1) - iOffset(ii), ":", &
                        & iMain(ii,2) - iOffset(ii)
      end do
  end subroutine TOrbitalCache_align



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

    complex(dp) :: phases(nCell, size(kPoints, dim=2))
    integer :: nPoints(4)
    integer :: i1, i2, i3, iEig, iAtom, iOrb, iM, iSpecies, iL, iCell

    ! Cache variables

    ! Cache Indexing Arrays
    integer, allocatable, target, save:: iAtomOffsets(:, :, :)
    integer, allocatable, target, save:: iAtomChunks(:, :, :)
    integer, allocatable, target, save:: iMainIndices(:, :, :, :)
    integer, pointer :: iOffset(:), iChunk(:), iMain(:,:)

    real(dp), pointer :: cachePtr(:,:,:)

    ! Todo: This overrides targetGridDistance and subdivisionFactor.
    integer, parameter :: subdivisionFactorLead = 5

    ! Todo: Figure out what unit we are using here. A?
    real(dp) :: targetGridDistance(3)
    
    ! Subdivision factors calculated to meet the targetGridDistance requirements
    integer :: subdivisionFactor(3)
  
    real(dp) :: pos(3)
    
    integer :: nUniqueOrb, coeffInd

    real(dp) :: cacheValue
    real(dp) :: chunkSeparation(3)

    ! TOrbitalCache array, one for each unique orbital
    type(TOrbitalCache), allocatable :: orbitalCache(:)

    !---------------------------------

    ! TODO: Use Fypp for the complex/real distinction?
    if (.not. tReal) then
      stop "TODO: Complex not implemented yet"
    end if


    ! Determine subgrid count / subdivision Factor / Grid Resolution
    targetGridDistance = norm2(gridVecs, dim=1) / subdivisionFactorLead
    subdivisionFactor = ceiling(norm2(gridVecs, dim=1) / targetGridDistance)
    chunkSeparation = 1.0_dp / subdivisionFactor
    print "(*(G0, 1X))", "Subdivision Distance / chunk Separation:", chunkSeparation
    print "(*(G0, 1X))", "Subdivision Factors / chunk Count:", subdivisionFactor
   

    nUniqueOrb = SIZE(angMoms)
    ! One TOrbitalCache for each STO.
    ! Each TOrbitalCache has several dimensions for their respective angular Momenta.
    allocate(orbitalCache(nUniqueOrb))

    ! Go through all orbitals and initialise them (build the cache)
    do iOrb = 1, nUniqueOrb
      print "(*(G0, 1X))", " -> Caching orbital ", iOrb
      call orbitalCache(iOrb)%initialise(stos(iOrb), angMoms(iOrb), gridVecs, subdivisionFactor)
    end do

    ! Main grid size
    nPoints = shape(valueReal)
    print "(*(G0, 1X))", "Main Grid Dimensions", nPoints
    !print "(*(G0, 1X))", " ->", size(valueReal), "elements,",  sizeof(valueReal) / 1024.0 / 1024.0, "MB"


    ! allocate Index arrays
    allocate(iAtomOffsets(3, nAtom, nCell))
    allocate(iAtomChunks(3, nAtom, nCell))
    allocate(iMainIndices(3, 2, nAtom, nCell))
    
    ! Alignment
    ! Loop over all cells, atoms and their orbitals to select the appropriate cache subgrid, bounds, and offset
    do iCell = 1, nCell
      do iAtom = 1, nAtom
        iSpecies = species(iAtom)
        ! TODO: This currently overwrite because we loop over iOrb.
        ! The current examples have identical cutoffs, so this does not matter.
        !do iOrb = iStos(iSpecies), iStos(iSpecies + 1) - 1
          iOrb = iStos(iSpecies)
          WRITE (*, '(A,I0,A,I0,A,I0,A,I0)', ADVANCE='NO') CHAR(13) // "Aligning contribution from ",&
          &iAtom, " of ", nAtom, " in cell ", iCell, " of ", nCell

          iOffset => iAtomOffsets(:, iAtom, iCell)
          iChunk => iAtomChunks(:, iAtom, iCell)
          iMain => iMainIndices(:, :, iAtom, iCell)

          pos = coords(:, iAtom, iCell) - origin

          print "(*(G0, 1X))", " -> Atom Position: ", pos
          call orbitalCache(iOrb)%align(nPoints, pos, iOffset, iChunk, iMain)
        !end do
      end do
    end do
          

    print *, "Applying cached densities to grid"
    ! For each atom, determine the offsets, then align the cache and apply using slicing.
    do iCell = 1, nCell
      coeffInd = 1
      do iAtom = 1, nAtom
        iSpecies = species(iAtom)
        WRITE (*, '(A,I0,A,I0,A,I0,A,I0)', ADVANCE='NO') CHAR(13) // "Adding contribution from ",&
        &iAtom, " of ", nAtom, " in cell ", iCell, " of ", nCell

        ! Load Array Alignment boundarys
        iOffset => iAtomOffsets(:, iAtom, iCell)
        iChunk => iAtomChunks(:, iAtom, iCell)
        iMain => iMainIndices(:, :, iAtom, iCell)

        do iOrb = iStos(iSpecies), iStos(iSpecies + 1) - 1
          iL = angMoms(iOrb)
          do iM = -iL, iL
            ! Access the correct subgrid.
            ! Modifies the Offset to compensate for the lost bound mapping.
            call orbitalCache(iOrb)%access(cachePtr, iChunk, iM)
            
            ! Loop over the aligned gridpoints and add the contribution
            do i3 = iMain(3,1), iMain(3,2)
              do i2 = iMain(2,1), iMain(2,2)
                do iEig = 1, nPoints(4)
                  do i1 = iMain(1,1), iMain(1,2)


                    cacheValue = cachePtr(i1+iOffset(1), i2+iOffset(2), i3+iOffset(3))
                    if (tAddDensities) then
                      cacheValue = cacheValue * cacheValue
                    end if
                    valueReal(i1, i2, i3, iEig) = valueReal(i1, i2, i3, iEig) + eigVecsReal(coeffInd, iEig) * cacheValue
                  end do
                end do
              end do
            end do

            coeffInd = coeffInd + 1
          end do
        end do
      end do
    end do
    ! Print Newline
    print "(*(G0, 1X))", " -> Done!"
  end subroutine localGetValue

end module waveplot_molorb2
