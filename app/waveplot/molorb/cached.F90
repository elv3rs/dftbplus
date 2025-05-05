!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

! Notes:
! General Todo list: 
! Quit trying to apply premature optimisations!
! Always profile first.
!
! Near term todo:
! Investigate radial interpolation (-> Provide way to disable)
! Investigate periodic wrapping
! Compare results and execution time
! Readd pointwise molorb
! Add option to choose molorb
! Check complex version
! Look into repeated cache copy allocations


#:include 'common.fypp'

!> Contains routines to calculate the value of one or more molecular orbitals composed from STOs on
!! an equidistant grid.
module waveplot_molorb_cached
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
  
  public :: evaluateCached
  public :: gridBoxFromSphereRadius


  ! Wavefunctions are evaluated across a finer grid (subdivisionFactor as many points as the main one).
  ! The X-coordinate is subdivided into <subdivisionFactor> chunked parts to ensure a contiguous memory layout.
  ! Y and Z have been subdivided as well, though only for convenience.
  type :: TOrbitalCache
    !private / TODO: was set to public to write some tests, fix this
    !> Cache for a single orbital.
    !! Layout: (X, XChunked, Y, YChunked, Z, ZChunked, iM)
    real(dp), pointer :: cache(:,:,:,:,:,:,:) => null()

    !> How many shifted copies of the main grid are stored in the cache
    integer :: subdivisionFactor

    !> The grids basis vectors
    real(dp) :: gridVecs(3,3)

    !> Size of cache from center to cutoff
    integer :: nPointsHalved(3)
  contains
    procedure :: initialise => TOrbitalCache_initialise
    procedure :: access => TOrbitalCache_access
    procedure :: align => TOrbitalCache_align
    final :: TOrbitalCache_finalise
  end type TOrbitalCache




contains

  ! Finds the smallest integer array dimensions in terms of the basis gridVecs to fit
  ! the sphere defined by the cutoff radius.
  subroutine gridBoxFromSphereRadius(gridVecs, cutoff, nPointsHalved)
    !>  Basis vectors, where gridVecs(:,j) is the j-th vector.
    real(dp), intent(in) :: gridVecs(3,3)
    !>  Radius of sphere to fit within nPointsHalved.
    real(dp), intent(in) :: cutoff
    !>  Output array size 
    integer, intent(out) :: nPointsHalved(3)
    !--------------------
    !> Basis Vectors (copy of gridVecs)
    real(dp) :: B(3,3)
    !> Gram matrix G = B^T * B
    real(dp) :: G(3,3)
    !> Inverse Gram matrix
    real(dp) :: G_inv(3,3) 
    !> Loop index
    integer  :: i

    @:ASSERT(cutoff > 0.0_dp)

    B = gridVecs
    ! Gram Matrix G = B^T * B
    G = matmul(transpose(B), B)

    ! Copy and invert
    G_inv = G
    call invert33(G_inv)
    do i = 1, 3
       nPointsHalved(i) = ceiling(cutoff * sqrt(G_inv(i,i)))
    end do

  end subroutine gridBoxFromSphereRadius




    !> Initialises the cache at the provided subdivisionFactor for a given STO and the grid vectors.
    subroutine TOrbitalCache_initialise(this, sto, gridVecs, subdivisionFactor)
      class(TOrbitalCache), intent(inout) :: this
      !> The Slater orbital to be cached
      type(TSlaterOrbital), intent(in) :: sto
      !> Basis vectors of the main grid
      real(dp), intent(in) :: gridVecs(3,3)
      !> How many shifted subgrids are to be stored in the cache.
      !! Determines the accuracy of the approximation.
      integer, intent(in) :: subdivisionFactor
      !> Expected size of the cache in MB
      real(dp) :: expectedSizeMB

      integer :: nPointsHalved(3)
      !---------------------------------
      
      ! Skip initialisation if the cache is already populated.
      if (associated(this%cache)) then
        print "(*(G0, 1X))", " -> Cache already initialised, skipping."
        return
      end if

      ! Store parameters and determine cache dimensions
      this%gridVecs = gridVecs
      this%subdivisionFactor = subdivisionFactor
      call gridBoxFromSphereRadius(gridVecs, sto%cutoff, nPointsHalved)
      this%nPointsHalved = nPointsHalved

      ! Allocate and populate the cache
      expectedSizeMB = 0.0_dp
      expectedSizeMB = storage_size(1.0_dp) * product(this%nPointsHalved)
      expectedSizeMB = expectedSizeMB * 2**3 * this%subdivisionFactor**3 * (2 * sto%angMom + 1)
      expectedSizeMB = expectedSizeMB / 8.0 / 1024.0 / 1024.0
      print "(*(G0, 1X))", "Allocating Cache Grid of dimensions", this%nPointsHalved * 2
      print "(*(G0, 1X))", "Expected Cache Allocation Size", expectedSizeMB, "MB"
      if (expectedSizeMB > 8000) then
        stop "Expected cache array size exceeds 8GB"
      end if

      ! We define (0,0,0) to be the origin using Fortrans fancy arbitrary indices feature.
      ! When we access the cache via the pointer view, we need to correct the indices.
      allocate(this%cache(-this%nPointsHalved(1):this%nPointsHalved(1), 0:this%subdivisionFactor-1, &
                        & -this%nPointsHalved(2):this%nPointsHalved(2), 0:this%subdivisionFactor-1, &
                        & -this%nPointsHalved(3):this%nPointsHalved(3), 0:this%subdivisionFactor-1, &
                        & -sto%angMom:sto%angMom), source=0.0_dp)
      call populateCache(this, sto)

    end subroutine TOrbitalCache_initialise


    !> Populates the cache by evaluating the STO on the fine grid.
    subroutine populateCache(this, sto)
      class(TOrbitalCache), intent(in) :: this
      type(TSlaterOrbital), intent(in) :: sto
      !----------------
      real(dp) :: curCoords(3,3), xx, val, xyz(3)
      integer :: i1, i2, i3, iM,  i1Chunked, i2Chunked, i3Chunked
      real(dp) :: chunkFraction(3)

      ! For every point on the fine grid:
      lpI3Chunked : do i3Chunked = 0, this%subdivisionFactor-1
        chunkFraction(3) = real(i3Chunked, dp) / this%subdivisionFactor
        lpI3 : do i3 = -this%nPointsHalved(3), this%nPointsHalved(3)
          curCoords(:, 3) = (i3 + chunkFraction(3)) * this%gridVecs(:, 3)
          lpI2Chunked : do i2Chunked = 0, this%subdivisionFactor-1
            chunkFraction(2) = real(i2Chunked, dp) / this%subdivisionFactor
            lpI2 : do i2 = -this%nPointsHalved(2), this%nPointsHalved(2)
              curCoords(:, 2) = (i2 + chunkFraction(2)) * this%gridVecs(:, 2)
              lpI1Chunked : do i1Chunked = 0, this%subdivisionFactor-1
                chunkFraction(1) = real(i1Chunked, dp) / this%subdivisionFactor
                lpI1 : do i1 = -this%nPointsHalved(1), this%nPointsHalved(1)
                  curCoords(:, 1) = (i1 + chunkFraction(1)) * this%gridVecs(:, 1)
                  ! Calculate position
                  xyz = sum(curCoords, dim=2)
                  xx = norm2(xyz)
                  if (xx <= sto%cutoff) then
                    call sto%getRadial(xx, val)
                    lpIM : do iM = -sto%angMom, sto%angMom
                      ! Combine with angular dependence and add to cache
                      this%cache(i1, i1Chunked, i2, i2Chunked, i3, i3Chunked, iM) = val * realTessY(sto%angMom, iM, xyz, xx)

                       if(i1Chunked + i2Chunked + i3Chunked == 0) then
                         ! Only print the first chunk
                         !print  "(*(G0, 1X))", i1, i2, i3, "->", xx, "->", val, realTessY(sto%angMom, iM, xyz, xx)
                       end if

                    end do lpIM
                  end if
                end do lpI1
              end do lpI1Chunked
            end do lpI2
          end do lpI2Chunked
        end do lpI3
      end do lpI3Chunked
    end subroutine populateCache


  !> Provides a copy of the cache for a given angular momentum <iM> and
  !! the requested subgrid <iChunk>. 
  subroutine TOrbitalCache_access(this, cacheCopy, iChunk, iOffset, iM)
    class(TOrbitalCache), intent(in) :: this
    !> Contiguous copy of the requested cache memory
    real(dp), allocatable, intent(inout) :: cacheCopy(:,:,:)
    !> Which subgrid to access
    integer, intent(in) :: iChunk(3)
    !> Offset of the cache in the main grid, to be baked into the indices 
    integer, intent(in) :: iOffset(3)
    !> Which angular momentum to access
    integer, intent(in) :: iM
    !! Indices bounds
    integer :: nn(3,2)

    nn(:,1) = -this%nPointsHalved(:) - iOffset(:)
    nn(:,2) = this%nPointsHalved(:) - iOffset(:)

    ! Check if the cache is already allocated and large enough
    if(allocated(cacheCopy) &
      & .and. size(cacheCopy, 1) == nn(1,2) - nn(1,1) + 1 &
      & .and. size(cacheCopy, 2) == nn(2,2) - nn(2,1) + 1 &
      & .and. size(cacheCopy, 3) == nn(3,2) - nn(3,1) + 1) then
      ! Problem: Need to remap indices.
      deallocate(cacheCopy)
    else if (allocated(cacheCopy)) then
      deallocate(cacheCopy)
    end if

    if (.not. allocated(cacheCopy)) then
      allocate(cacheCopy(nn(1,1):nn(1,2), nn(2,1):nn(2,2), nn(3,1):nn(3,2)))
    end if

    cacheCopy(nn(1,1):nn(1,2), nn(2,1):nn(2,2), nn(3,1):nn(3,2)) &
          & = this%cache(:, iChunk(1), :, iChunk(2), :, iChunk(3), iM)

  end subroutine TOrbitalCache_access

  !> Aligns the cache to the main grid.
  !! Returns the chunks to choose the correct subgrid,
  !! as well as the indices bounds for the main grid and an offset for the cache.
  subroutine TOrbitalCache_align(this, gridDims, shiftedPos, iOffset, iChunk, iMain)
      class(TOrbitalCache), intent(in) :: this
      !> Extent of the main grid (x, y, z, eig)
      integer, intent(in) :: gridDims(4)
      !> Relative position of the atom in regular coordinates
      real(dp), intent(in) :: shiftedPos(:)
      !> Relevant Main grid coordinates
      integer, intent(out) :: iMain(3,2)
      !> Offset to align cache with main grid
      integer, intent(out) :: iOffset(3)
      !> Indices to access the closest cache subgrid
      integer, intent(out) :: iChunk(3)

      ! Used for decomposing Atom position into cacheGrid Basis vecs using gesv
      real(dp) :: tmpBasis(3,3), pos(3,1)
      ! Debug loop index
      integer :: ii
      !---------------------------------
      

      pos(:,1) = shiftedPos
      tmpBasis(:,:) = this%gridVecs
      ! Get the atom position in terms of the basis <gridVecs>, stored in pos
      call gesv(tmpBasis, pos)

      ! Cache Indices are derived from the main indices, and:
      ! -> offset by the atom position
      ! Added 1 since cache starts at 0 and main grid at 1
      iOffset(:) = int(pos(:,1)) + 1
      ! -> shifted by (0..subdivisionFactor-1)/subdivisionFactor to select the closest subgrid (chunk)
      ! TODO: By replacing int with nint we get the closest one, but would need to adjust the main indices.
      iChunk(:) = modulo(int(pos(:,1) * this%subdivisionFactor), this%subdivisionFactor)

      ! Lower Main Indices need to include
      ! -> start of main grid (1)
      ! -> start of cache grid (atom offset - half cache size)
      iMain(:, 1) = max(1, iOffset(:) - this%nPointsHalved(:))
      ! Upper Main Indices need to include
      ! -> end of main grid (nPoints)
      ! -> end of cache grid (atom offset + half cache size)
      iMain(:, 2) = min(gridDims(:3), iOffset(:) + this%nPointsHalved(:))


      !print "(*(G0, 1X))", " -> Atom Position in Basis vectors:", pos(:,1)
      !print "(*(G0, 1X))", " -> Atom Position in Grid:", iOffset(:), "Chunk:", iChunk(:)
      !print "(*(G0, 1X))", "Indices Mapping Main -> Cache:"
      
      ! Add the offset instead of subtracting it henceforth
      iOffset(:) = -iOffset(:)

      !do ii = 1, 3
      !  print "(*(G0, 1X))", " o", Char(ii + 87), iMain(ii,1), ":" , &
      !                  & iMain(ii,2), "->",&
      !                  & iMain(ii,1) + iOffset(ii), ":", &
      !                  & iMain(ii,2) + iOffset(ii)
      !end do
  end subroutine TOrbitalCache_align



  !> Deallocates the cache.
  subroutine TOrbitalCache_finalise(this)
    type(TOrbitalCache), intent(inout) :: this

    if (associated(this%cache)) then
       deallocate(this%cache)
    end if

  end subroutine TOrbitalCache_finalise


  !> Returns the values of several molecular orbitals on grids.
  !! Caveat: The flag tPeriodic decides if the complex or the real version is read/written for the
  !! various parameters.
  subroutine evaluateCached(origin, gridVecs, eigVecsReal, eigVecsCmpl, nAtom, nOrb, coords,&
      & species, iStos, stos, tPeriodic, tReal, latVecs, recVecs2p, kPoints,&
      & kIndexes, nCell, cellVec, tAddDensities, subdivisionFactor, valueReal, valueCmpl)

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

    !> Number of cached subgrids. Sets the accuracy of the approximation. 
    integer, intent(in) :: subdivisionFactor

    !> Contains the real grid on exit
    real(dp), intent(out) :: valueReal(:,:,:,:)

    !> Contains the complex grid on exit
    complex(dp), intent(out) :: valueCmpl(:,:,:,:)

    complex(dp) :: phases(nCell, size(kPoints, dim=2))
    integer :: nPoints(4)
    integer :: i1, i2, i3, iEig, iAtom, iOrb, iM, iSpecies, iL, iCell

    ! Cache variables

    integer :: iOffset(3), iChunk(3), iMain(3,2)

    real(dp), allocatable :: cacheCopy(:,:,:)

  
    real(dp) :: pos(3)
    
    integer :: coeffInd

    real(dp) :: val
    
    ! TOrbitalCache array, one for each unique orbital
    type(TOrbitalCache), allocatable, save :: orbitalCache(:)

    !> timing variables
    integer(dp) :: startTime, endTime, clockRate

    
    !---------------------------------

    call system_clock(count_rate=clockRate, count=startTime)

    print "(*(G0, 1X))", " -> Number of unique orbitals: ", size(stos)
    ! One TOrbitalCache for each STO.
    if(.not. allocated(orbitalCache)) then
      allocate(orbitalCache(size(stos)))
      ! Go through all orbitals and initialise them (build the cache)
      do iOrb = 1, size(stos)
        print "(*(G0, 1X))", " -> Caching orbital ", iOrb
        call orbitalCache(iOrb)%initialise(stos(iOrb), gridVecs, subdivisionFactor)
      end do
    end if

    call system_clock(count=endTime)
    print *, "InitTime:", real(endTime - startTime, dp) / clockRate



    ! Main grid size
    ! Zero the output array
    if (tReal) then
      nPoints = shape(valueReal)
      valueReal(:,:,:,:) = 0.0_dp
      print *, "Real"
    else
      nPoints = shape(valueCmpl)
      valueCmpl(:,:,:,:) = 0.0_dp
      print *, "Complex"
    end if
    print "(*(G0, 1X))", "Main Grid Dimensions", nPoints

    ! Phase factors for the periodic image cell. Note: This will be conjugated in the scalar product
    ! below. This is fine as, in contrast to what was published, DFTB+ implicitly uses exp(-ikr) as
    ! a phase factor, as the unpack routines assemble the lower triangular matrix with exp(ikr) as
    ! factor.
    phases(:,:) = exp(imag * matmul(transpose(cellVec), kPoints))

    ! For each atom, determine the offsets, then align the cache and apply using slicing.
    do iCell = 1, nCell
      coeffInd = 1
      do iAtom = 1, nAtom
        iSpecies = species(iAtom)
        print "(*(G0, 1X))",  "Applying contribution from ", iAtom, " of ", nAtom, " in cell ", iCell, " of ", nCell

        ! Where to shift the orbital origin to
        pos = coords(:, iAtom, iCell) - origin
        !print *, "."
        !print *, "Atom coords", coords(:, iAtom, iCell)
        !print *, "Origin", origin
        !print *, "Basis vector 1", gridVecs(:, 1)
        !print *, "Basis vector 2", gridVecs(:, 2)
        !print *, "Basis vector 3", gridVecs(:, 3)
        ! Todo: Periodic systems and correct unit cell mapping
        !if (tPeriodic) then
        !  frac(:) = matmul(pos, recVecs2p)
        !  pos(:) = matmul(latVecs, frac - real(floor(frac), dp))
        !end if

        do iOrb = iStos(iSpecies), iStos(iSpecies + 1) - 1
          ! Calculate alignment bounds and select the correct subgrid
          call orbitalCache(iOrb)%align(nPoints, pos, iOffset, iChunk, iMain)

          iL = stos(iOrb)%angMom
          do iM = -iL, iL
            call orbitalCache(iOrb)%access(cacheCopy, iChunk, iOffset, iM)

            ! Loop over the aligned gridpoints and add the contribution
            do i3 = iMain(3,1), iMain(3,2)
              do i2 = iMain(2,1), iMain(2,2)
                do iEig = 1, nPoints(4) 
                  do i1 = iMain(1,1), iMain(1,2)
                    val = cacheCopy(i1, i2, i3)
                    if (tReal) then
                      if (tAddDensities) then
                        val = val * val
                      end if
                      valueReal(i1, i2, i3, iEig) = valueReal(i1, i2, i3, iEig) + eigVecsReal(coeffInd, iEig) * val
                      !print "(*(G0, 1X))", i1, i2, i3, "->", valueReal(i1, i2, i3, iEig)
                    else
                      ! TODO: Verify this is correct
                      valueCmpl(i1, i2, i3, iEig) = valueCmpl(i1, i2, i3, iEig) + eigVecsCmpl(coeffInd, iEig) &
                          & * phases(iCell, kIndexes(iEig)) * val
                    end if
                  end do
                end do
              end do
            end do
            coeffInd = coeffInd + 1
          end do
        end do
      end do
    end do

    print "(*(G0, 1X))", " -> Done!"
  end subroutine evaluateCached

end module waveplot_molorb_cached
