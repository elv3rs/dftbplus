!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> A cache for calculating molecule orbitals on a grid.
!! This object is responsible for reading in the eigenvectors from a specified file and passing the
!! appropriate eigenvectors to the molecule orbital calculator.
module libwavegrid_gridcache
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : closeFile, openFile, TFileDescr
  use dftbp_common_globalenv, only : stdOut
#:if WITH_MPI
  use dftbp_common_schedule, only : getStartAndEndIndex
#:endif
  use libwavegrid_molorb, only : getValue, TMolecularOrbital
  use dftbp_io_message, only : error
  implicit none

  private
  public :: TGridCache

  !> Contains the data for a grid cache.
  type TGridCache

    !> Molecular orbital calculator
    type(TMolecularOrbital), pointer :: molorb

    !> Grid vectors
    real(dp) :: gridVec(3,3)

    !> Origin of the grid
    real(dp) :: origin(3)

    !> Real eigenvectors
    real(dp), allocatable :: eigenvecReal(:,:)

    !> Complex eigenvectors
    complex(dp), allocatable :: eigenvecCmpl(:,:)

    !> Parameters of the levels
    integer, allocatable :: levelIndex(:,:)

    !> Eigenvectors read so far
    integer :: nReadEigVec

    !> Grids processed so far
    integer :: iGrid

    !> Nr. of grids to process
    integer :: nGrid

    !> Position in the cache
    integer :: cachePos

    !> Are we ready?
    logical :: tFinished

    !> File descriptor for eigenvec
    type(TFileDescr) :: fdEigVec

    !> Size of the eigenvectors
    integer :: nOrb

    !> Levels in the eigenvec file
    integer :: nAllLevel

    !> The k-points in the eigenv. file
    integer :: nAllKPoint

    !> Spins in the eigenvec. file
    integer :: nAllSpin

    !> Verbose?
    logical :: tVerbose

    !> Nr. of cached grids
    integer :: nCached

    !> Subdivision factor
    integer :: subdivisionFactor

    !> Cache for real grids
    real(dp), allocatable :: gridCacheReal(:,:,:,:)

    !> Cache for complex grids
    complex(dp), allocatable :: gridCacheCmpl(:,:,:,:)

    !> KPoints
    real(dp), allocatable :: kPoints(:,:)

    !> Are eigenvectors real
    logical :: tReal

    !> Initialised?
    logical :: tInitialised = .false.
  contains
    procedure :: TGridCache_next_real
    procedure :: TGridCache_next_cmpl
    generic :: next => TGridCache_next_real, TGridCache_next_cmpl
    procedure :: init => TGridCache_init

  end type TGridCache


contains

  !> Initialises a GridCache instance.
  !! Caveat: Level index is not allowed to contain duplicate entries!
  subroutine TGridCache_init(sf, env, levelIndexAll, nOrb, nAllLevel, nAllKPoint, nAllSpin, nCached,&
      & nPoints, tVerbose, eigvecBin, gridVec, origin, kPointCoords, tReal, molorb, subdivisionFactor)

    !> Structure to initialise
    class(TgridCache), intent(out) :: sf

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Contains indexes (spin, kpoint, state) of the levels which must be calculated
    integer, intent(in) :: levelIndexAll(:,:)

    !> Nr. of orbitals (elements) in an eigenvectors
    integer, intent(in) :: nOrb

    !> Nr. of levels, for which eigenvector is provided in the eigenvector file
    integer, intent(in) :: nAllLevel

    !> Nr. of k-points for which eigenvectors are provided.
    integer, intent(in) :: nAllKPoint

    !> Nr. of all spins for which eigenvectors are provided.
    integer, intent(in) :: nAllSpin

    !> Nr. of cached grids
    integer, intent(in) :: nCached

    !> Nr. of grid points along the grid vectors
    integer, intent(in) :: nPoints(:)

    !> If verbosity should turned on
    logical, intent(in) :: tVerbose

    !> Name of the binary eigenvector file
    character(len=*), intent(in) :: eigvecBin

    !> Grid vectors
    real(dp), intent(in) :: gridVec(:,:)

    !> Origin of the grid
    real(dp), intent(in) :: origin(:)

    !> Coordinates of the k-points
    real(dp), intent(in) :: kPointCoords(:,:)

    !> If grids and eigenvectors are real
    logical, intent(in) :: tReal

    !> Molecular orbital calculator
    type(TMolecularOrbital), pointer, intent(in) :: molorb

    !> Subdivision factor for the cache grid
    integer, intent(in) :: subdivisionFactor

    !! Contains indexes (spin, kpoint, state) to be calculated by the current MPI process
    integer, allocatable :: levelIndex(:,:)

    !! Start and end level index for MPI parallelization, if applicable
    integer :: iLvlStart, iLvlEnd

    !! Auxiliary variables
    integer ::nAll
    integer :: iSpin, iKPoint, iLevel, ind, ii, iostat
    integer :: curVec(3)
    logical :: tFound
  
    integer :: startTime, endTime

    call system_clock(count=startTime)
  #:if WITH_MPI
    call getStartAndEndIndex(env, size(levelIndexAll, dim=2), iLvlStart, iLvlEnd)
  #:else
    iLvlStart = 1
    iLvlEnd = size(levelIndexAll, dim=2)
  #:endif

    levelIndex = levelIndexAll(:, iLvlStart:iLvlEnd)

    @:ASSERT(.not. sf%tInitialised)
    @:ASSERT(size(levelIndex, dim=1) == 3)
    @:ASSERT(size(levelIndex, dim=2) > 0)
    @:ASSERT(minval(levelIndex) > 0)
    @:ASSERT(maxval(levelIndex(1, :)) <= nAllLevel)
    @:ASSERT(maxval(levelIndex(2, :)) <= nAllKPoint)
    @:ASSERT(maxval(levelIndex(3, :)) <= nAllSpin)
    @:ASSERT(associated(molorb))
    @:ASSERT(all(shape(gridVec) == [3, 3]))
    @:ASSERT(size(origin) == 3)
    @:ASSERT(size(nPoints) == 3)
    @:ASSERT(all(shape(kPointCoords) == [3, nAllKPoint]))


    sf%molorb => molorb
    sf%gridVec(:,:) = gridVec
    sf%origin = origin
    sf%nOrb = nOrb
    sf%nAllLevel = nAllLevel
    sf%nAllKPoint = nAllKPoint
    sf%nAllSpin = nAllSpin
    sf%tVerbose = tVerbose
    allocate(sf%kPoints(3, nAllKPoint))
    sf%kPoints(:,:) = 2.0_dp * pi * kPointCoords
    sf%nCached = nCached
    sf%tReal = tReal
    sf%subdivisionFactor = subdivisionFactor

    if (sf%tReal) then
      allocate(sf%gridCacheReal(nPoints(1), nPoints(2), nPoints(3), nCached))
      allocate(sf%eigenvecReal(sf%nOrb, sf%nCached))
    else
      allocate(sf%gridCacheCmpl(nPoints(1), nPoints(2), nPoints(3), nCached))
      allocate(sf%eigenvecCmpl(sf%nOrb, sf%nCached))
    end if

    nAll = size(levelIndex, dim=2)
    allocate(sf%levelIndex(3, nAll))
    ! Make sure, entries are correctly sorted in the list
    ind = 1
    do iSpin = 1, nAllSpin
      do iKPoint = 1, nAllKPoint
        do iLevel = 1, nAllLevel
          curVec = [iLevel, iKPoint, iSpin]
          tFound = .false.
          lpLevelIndex: do ii = 1, size(levelIndex, dim=2)
            tFound = all(levelIndex(:, ii) == curVec)
            if (tFound) exit lpLevelIndex
          end do lpLevelIndex
          if (tFound) then
            sf%levelIndex(:, ind) = curVec
            ind = ind + 1
          end if
        end do
      end do
    end do
    @:ASSERT(ind == nAll + 1)
    sf%nGrid = nAll
    sf%iGrid = 1
    sf%cachePos = 1
    sf%nReadEigVec = 0
    sf%tFinished = .false.
    call openFile(sf%fdEigVec, eigvecBin, mode="rb", iostat=iostat)
    if (iostat /= 0) then
      call error("Can't open file '" // trim(eigvecBin) // "'.")
    end if
    read(sf%fdEigVec%unit) ii
    call system_clock(count=endTime)

    sf%tInitialised = .true.

    ! Timing
    print *,  "GridCache initialised in", endTime - startTime

  end subroutine TGridCache_init


  !> Returns the next entry from the cache (real version).
  subroutine TGridCache_next_real(sf, gridValReal, levelIndex, tFinished)

    !> Gridcache instance
    class(TgridCache), intent(inout) :: sf

    !> Contains the molecular orbital on the grid on exit
    real(dp), pointer :: gridValReal(:,:,:)

    !> Indices of the moleular orbital (spin, kpoint, level)
    integer, intent(out) :: levelIndex(:)

    !> If all orbitals had been processed.
    logical, intent(out) :: tFinished

    complex(dp), pointer, save :: gridValCmpl(:,:,:) => null()

    call localNext(sf, gridValReal, gridValCmpl, levelIndex, tFinished)

  end subroutine TGridCache_next_real


  !> Returns the next entry from the cache (complex version).
  subroutine TGridCache_next_cmpl(sf, gridValCmpl, levelIndex, tFinished)

    !> Gridcache instance
    class(TgridCache), intent(inout) :: sf

    !> Contains the molecular orbital on the grid on exit
    complex(dp), pointer :: gridValCmpl(:,:,:)

    !> Indices of the moleular orbital (spin, kpoint, level)
    integer, intent(out) :: levelIndex(:)

    !> If all orbitals had been processed.
    logical, intent(out) :: tFinished

    real(dp), pointer, save :: gridValReal(:,:,:) => null()

    call localNext(sf, gridValReal, gridValCmpl, levelIndex, tFinished)

  end subroutine TGridCache_next_cmpl


  !> Working subroutine for the TGridCache_next_* subroutines.
  subroutine localNext(sf, gridValReal, gridValCmpl, levelIndex, tFinished)

    !> Gridcache instance
    type(TgridCache), intent(inout), target :: sf

    !> Contains the real grid onexit
    real(dp), pointer :: gridValReal(:,:,:)

    !> Contains the complex grid on exit
    complex(dp), pointer :: gridValCmpl(:,:,:)

    !> Level indexes of the processed orbital on exit
    integer, intent(out) :: levelIndex(:)

    !> If all orbitals had been processed.
    logical, intent(out) :: tFinished

    integer :: iEnd, iStartAbs, iEndAbs, iLevel, iKPoint, iSpin
    integer :: ind, tmp
    real(dp), pointer :: eigReal(:,:)
    complex(dp), pointer :: eigCmpl(:,:)

    @:ASSERT(sf%tInitialised)
    @:ASSERT(.not. sf%tFinished)

    ! We passed back everything from the cache, fill it with new grids
    if (mod(sf%cachePos - 1, sf%nCached) == 0) then
      sf%cachePos = 1
      iStartAbs = sf%iGrid
      iEnd = min(sf%nGrid, sf%iGrid + sf%nCached - 1) - sf%iGrid + 1
      iEndAbs = sf%iGrid + iEnd - 1

      ind = 1
      do while (ind <= iEnd)
        if (sf%tReal) then
          read(sf%fdEigVec%unit) sf%eigenvecReal(:,ind)
        else
          read(sf%fdEigVec%unit) sf%eigenvecCmpl(:,ind)
        end if
        sf%nReadEigVec = sf%nReadEigVec + 1

        ! If eigenvec belongs to a level which must be plotted, keep it
        iSpin = (sf%nReadEigVec - 1) / (sf%nAllLevel * sf%nAllKPoint) + 1
        tmp = mod(sf%nReadEigVec - 1, sf%nAllLevel * sf%nAllKPoint)
        iKPoint = tmp / sf%nAllLevel + 1
        iLevel = mod(tmp, sf%nAllLevel) + 1
        if (all([iLevel, iKPoint, iSpin] == sf%levelIndex(:,iStartAbs+ind-1))) then
          ind = ind + 1
          if (sf%tVerbose) then
            write(stdout, "(I5,I7,I7,A8)") iSpin, iKPoint, iLevel, "read"
          end if
        end if
      end do

      ! Get molecular orbital for that eigenvector
      if (sf%tVerbose) then
        write(stdout, "(/,A,/)") "Calculating grid"
      end if
      if (sf%tReal) then
        eigReal => sf%eigenvecReal(:, :iEnd)
        call getValue(sf%molorb, sf%origin, sf%gridVec, eigReal, sf%subdivisionFactor, sf%gridCacheReal(:,:,:,:iEnd))
      else
        eigCmpl => sf%eigenvecCmpl(:, :iEnd)
        call getValue(sf%molorb, sf%origin, sf%gridVec, eigCmpl, sf%kPoints,&
            & sf%levelIndex(2, iStartAbs:iEndAbs), sf%subdivisionFactor, sf%gridCacheCmpl(:,:,:,:iEnd))
      end if
    end if

    ! Return the appropriate grid
    if (sf%tReal) then
      gridValReal => sf%gridCacheReal(:,:,:,sf%cachePos)
    else
      gridValCmpl => sf%gridCacheCmpl(:,:,:,sf%cachePos)
    end if
    levelIndex(:) = sf%levelIndex(:,sf%iGrid)

    ! Increase cache and grid counters
    sf%iGrid = sf%iGrid + 1
    sf%cachePos = sf%cachePos + 1
    if (sf%iGrid > sf%nGrid) then
      sf%tFinished = .true.
      call closeFile(sf%fdEigVec)
    end if
    tFinished = sf%tFinished

  end subroutine localNext

end module libwavegrid_gridcache
