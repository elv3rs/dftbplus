!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the routines for initialising Waveplot.
module waveplot_initwaveplot
  use waveplot_gridcache, only : TGridCache, TGridCache_init
  use waveplot_molorb, only : TMolecularOrbital, TMolecularOrbital_init, TSpeciesBasis
  use waveplot_slater, only : TSlaterOrbital_init
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : closeFile, openFile, setDefaultBinaryAccess, TFileDescr
  use dftbp_common_globalenv, only : stdOut, tIoProc
  use dftbp_common_release, only : releaseYear
  use dftbp_common_status, only : TStatus
  use dftbp_common_unitconversion, only : lengthUnits
  use dftbp_dftb_boundarycond, only : boundaryCondsEnum, TBoundaryConds,&
      & TBoundaryConds_init
  use dftbp_dftbplus_input_fileaccess, only : readBinaryAccessTypes
  use dftbp_io_charmanip, only : i2c, unquote
  use dftbp_io_formatout, only : printDftbHeader
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, getSelectedIndices
  use dftbp_io_unitconv, only : convertUnitHsd
  use hsd, only : hsd_warn_unprocessed, MAX_WARNING_LEN, hsd_error_t, hsd_rename_child,&
      & hsd_dump, hsd_clear_children, hsd_table_ptr, hsd_get_child_tables, hsd_get,&
      & hsd_get_or_set, hsd_set, hsd_get_table, hsd_get_attrib, HSD_STAT_OK, hsd_node,&
      & hsd_remove_child, hsd_get_matrix, &
      & hsd_get_name
  use hsd_data, only : data_load, DATA_FMT_AUTO, hsd_table, new_table
  use dftbp_io_message, only : error, warning
  use dftbp_math_simplealgebra, only : determinant33
  use dftbp_type_typegeometryhsd, only : readTGeometryGen, readTGeometryHSD, readTGeometryVasp,&
      & readTGeometryXyz, TGeometry, writeTGeometryHSD
  implicit none

  private
  public :: TProgramVariables, TProgramVariables_init


  !> Single variable-length integer array element (replaces TListIntR1)
  type :: TIntArray
    integer, allocatable :: data(:)
  end type


  !> Data type containing variables from detailed.xml.
  type TInput

    !> Geometry instance
    type(TGeometry) :: geo

    !> Identity of the run
    integer :: identity

    !> Nr. of orbitals per state
    integer :: nOrb

    !> True, if eigenvectors/hamiltonian is real-valued
    logical :: tRealHam

    !> Occupations
    real(dp), allocatable :: occupations(:,:,:)

  end type TInput


  !> Data type containing variables from the Option block.
  type TOption

    !> Nr. of grid points along 3 directions
    integer :: nPoints(3)

    !> Repeat box along 3 directions
    integer :: repeatBox(3)

    !> Levels to plot
    integer, allocatable :: plottedLevels(:)

    !> The k-points to plot
    integer, allocatable :: plottedKPoints(:)

    !> Spins to plot
    integer, allocatable :: plottedSpins(:)

    !> If box should filled with folded atoms
    logical :: tFillBox

    !> If coords should be folded to unit cell
    logical :: tFoldCoords

    !> If program should be verbose
    logical :: tVerbose

    !> If total charge should be plotted
    logical :: tPlotTotChrg

    !> If total charge should be calculated
    logical :: tCalcTotChrg

    !> If total spin pol. to be plotted
    logical :: tPlotTotSpin

    !> If total charge difference to be plotted
    logical :: tPlotTotDiff

    !> If atomic densities to be plotted
    logical :: tPlotAtomDens

    !> If atomic densities to be calculated
    logical :: tCalcAtomDens

    !> If charge for orbitals to be plotted
    logical :: tPlotChrg

    !> If charge difference for orbs. to be plotted
    logical :: tPlotChrgDiff

    !> If real part of the wfcs to plot.
    logical :: tPlotReal

    !> If imaginary part of the wfcs to plot
    logical :: tPlotImag

    !> Box vectors for the plotted region
    real(dp) :: boxVecs(3,3)

    !> Origin of the box
    real(dp) :: origin(3)

    !> Origin of the grid in the box
    real(dp) :: gridOrigin(3)

    !> List of levels to plot, whereby insignificant occupations were filtered out
    integer, allocatable :: levelIndex(:,:)

    !> File access types
    character(20) :: binaryAccessTypes(2)

  end type TOption


  !> Data type containing variables from the Basis block.
  type TBasis

    !> Definition of the wfcs
    type(TSpeciesBasis), allocatable :: basis(:)

    !> Resolution of the radial wfcs
    real(dp) :: basisResolution

  end type TBasis


  !> Data type containing variables from eigenvec.bin.
  type TEig

    !> Nr. of states
    integer :: nState

    !> Real eigenvectors
    real(dp), allocatable :: eigvecsReal(:,:)

    !> Complex eigenvectors
    complex(dp), allocatable :: eigvecsCplx(:,:)

  end type TEig


  !> Data type containing variables from the AtomicNumbers block.
  type TAtomicNumber

    !> Species-atomic nr. corresp.
    integer, allocatable :: atomicNumbers(:)

  end type TAtomicNumber


  !> Data type containing locally created variables.
  type TInternal

    !> Molecular orbital
    type(TMolecularOrbital), allocatable :: molOrb

    !> pointer to the orbital
    type(TMolecularOrbital), pointer :: pMolOrb

    !> Grid cache
    type(TGridCache) :: grid

    !> Grid vectors
    real(dp) :: gridVec(3,3)

    !> Volume of the grid
    real(dp) :: gridVol

    !> List of levels to plot
    integer, allocatable :: levelIndex(:,:)

  end type TInternal


  !> Data type containing program variables.
  type TProgramVariables

    !> Data of detailed.xml
    type(TInput) :: input

    !> Data of eigenvec.bin
    type(TEig) :: eig

    !> Data of Option block
    type(TOption) :: opt

    !> Boundary condition
    type(TBoundaryConds) :: boundaryCond

    !> Data of Basis block
    type(TBasis) :: basis

    !> Data of AtomicNumber block
    type(TAtomicNumber) :: aNr

    !> Locally created data
    type(TInternal) :: loc

  end type TProgramVariables


  !> Program version
  character(len=*), parameter :: version = "0.3"

  !> Root node name of the input tree
  character(len=*), parameter :: rootTag = "waveplot"

  !> Input file name
  character(len=*), parameter :: hsdInput = "waveplot_in.hsd"

  !> Parsed output name
  character(len=*), parameter :: hsdParsedInput = "waveplot_pin.hsd"

  !> Version of the input document
  integer, parameter :: parserVersion = 3

contains


  !> Initialises the program variables.
  subroutine TProgramVariables_init(this, env)

    !> Container of program variables
    type(TProgramVariables), intent(out), target :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !! Pointers to input nodes
    type(hsd_table), pointer :: root, tmp, detailed, hsdTree

    !! String buffer instance
    character(len=:), allocatable :: strBuffer

    !! Id of input/parser version
    integer :: inputVersion

    !! Nr. of cached grids
    integer :: nCached

    !! Nr. of K-points
    integer :: nKPoint

    !! Nr. of spins
    integer :: nSpin

    !! Wether to look for ground state occupations (True) or excited (False)
    logical :: tGroundState

    !! If grid should shifted by a half cell
    logical :: tShiftGrid

    !! K-points and weights
    real(dp), allocatable :: kPointsWeights(:,:)

    !! File with binary eigenvectors
    character(len=1024) :: eigVecBin

    !! Auxiliary variable
    integer :: ii

    !! Operation status, if an error needs to be returned
    type(TStatus) :: errStatus

    !! HSD parse error
    type(hsd_error_t), allocatable :: hsdError

    integer :: stat

    ! Write header
    call printDftbHeader('(WAVEPLOT '// version //')', releaseYear)

    ! Read in input file as HSD
    block
      type(hsd_table), pointer :: content
      allocate(content)
      call data_load(hsdInput, content, hsdError, fmt=DATA_FMT_AUTO)
      if (allocated(hsdError)) then
        call error("Error loading input file '" // trim(hsdInput) // "': " // trim(hsdError%message))
      end if
      content%name = rootTag
      allocate(hsdTree)
      call new_table(hsdTree, name="document")
      call hsdTree%add_child(content)
    end block
    call hsd_get_table(hsdTree, rootTag, root, stat, auto_wrap=.true.)
    if (.not. associated(root)) call error("Missing root block: '" // rootTag // "'")

    write(stdout, "(A)") "Interpreting input file '" // hsdInput // "'"

    ! Check if input version is the one, which we can handle
    call hsd_get_or_set(root, "InputVersion", inputVersion, parserVersion)
    if (inputVersion /= parserVersion) then
      call error("Version of input (" // i2c(inputVersion) // ") and parser (" &
          &// i2c(parserVersion) // ") do not match")
    end if

    call hsd_get_or_set(root, "GroundState", tGroundState, .true.)

    ! Read data from detailed output (supports .xml, .hsd, .json)
    call hsd_get(root, "DetailedXML", strBuffer, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(root, "Missing required value: 'DetailedXML'")
    allocate(tmp)
    call data_load(unquote(strBuffer), tmp, hsdError, fmt=DATA_FMT_AUTO)
    if (allocated(hsdError)) then
      call error("Error loading detailed output file '" // unquote(strBuffer) // "': "&
          & // trim(hsdError%message))
    end if
    call hsd_get_table(tmp, "detailedout", detailed, stat, auto_wrap=.true.)
    if (.not. associated(detailed)) call error("Missing 'detailedout' block in detailed output")
    call readDetailed(this, detailed, tGroundState, kPointsWeights)
    tmp => null()

    nKPoint = size(kPointsWeights, dim=2)
    nSpin = size(this%input%occupations, dim=3)
    this%eig%nState = size(this%input%occupations, dim=1)

    ! Read basis
    call hsd_get_table(root, "Basis", tmp, stat, auto_wrap=.true.)
    if (.not. associated(tmp)) call dftbp_error(root, "Missing required block: 'Basis'")
    call readBasis(this, tmp, this%input%geo%speciesNames)
    call hsd_get(root, "EigenvecBin", strBuffer, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(root, "Missing required value: 'EigenvecBin'")
    eigVecBin = unquote(strBuffer)

    ! Read options
    call hsd_get_table(root, "Options", tmp, stat, auto_wrap=.true.)
    if (.not. associated(tmp)) call dftbp_error(root, "Missing required block: 'Options'")
    call readOptions(this, tmp, this%eig%nState, nKPoint, nSpin, nCached, tShiftGrid)

    ! Issue warning about unprocessed nodes
    if (.not. .true.) then
      block
        character(len=MAX_WARNING_LEN), allocatable :: warnings(:)
        integer :: ii
        call hsd_warn_unprocessed(root, warnings)
        do ii = 1, size(warnings)
          call warning(trim(warnings(ii)))
        end do
      end block
    end if

    ! Finish parsing, dump parsed and processed input
    if (tIoProc) then
      ! TODO(Phase 4): call hsd_dump(hsdTree, hsdParsedInput)
    end if
    write(stdout, "(A)") "Processed input written as HSD to '" // hsdParsedInput &
        &//"'"
    write(stdout, "(A,/)") repeat("-", 80)
    hsdTree => null()

  #:if WITH_MPI
    call env%initMpi(1)
  #:endif

    if (this%input%geo%tPeriodic) then
      call TBoundaryConds_init(this%boundaryCond, boundaryCondsEnum%pbc3d, errStatus)
    else if (this%input%geo%tHelical) then
      call TBoundaryConds_init(this%boundaryCond, boundaryCondsEnum%helical, errStatus)
    else
      call TBoundaryConds_init(this%boundaryCond, boundaryCondsEnum%cluster, errStatus)
    end if
    if (errStatus%hasError()) call error(errStatus%message)

    ! Create grid vectors, shift them if necessary
    do ii = 1, 3
      this%loc%gridVec(:, ii) = this%opt%boxVecs(:, ii) / real(this%opt%nPoints(ii), dp)
    end do
    if (tShiftGrid) then
      this%opt%gridOrigin(:) = this%opt%origin(:) + 0.5_dp * sum(this%loc%gridVec, dim=2)
    else
      this%opt%gridOrigin(:) = this%opt%origin(:)
    end if
    this%loc%gridVol = abs(determinant33(this%loc%gridVec))

    write(stdout, "(A)") "Doing initialisation"

    ! Set the same access for readwrite as for write (we do not open any files in readwrite mode)
    call setDefaultBinaryAccess(this%opt%binaryAccessTypes(1), this%opt%binaryAccessTypes(2),&
        & this%opt%binaryAccessTypes(2))

    ! Check eigenvector id
    call checkEigenvecs(eigVecBin, this%input%identity)

    ! Initialize necessary (molecular orbital, grid) objects
    allocate(this%loc%molOrb)
    this%loc%pMolOrb => this%loc%molOrb
    call TMolecularOrbital_init(this%loc%molOrb, this%input%geo, this%boundaryCond,&
        & this%basis%basis)

    call TGridCache_init(this%loc%grid, env, this%loc%levelIndex, this%input%nOrb, this%eig%nState,&
        & nKPoint, nSpin, nCached, this%opt%nPoints, this%opt%tVerbose, eigVecBin,&
        & this%loc%gridVec, this%opt%gridOrigin, kPointsWeights(1:3, :), this%input%tRealHam,&
        & this%loc%pMolOrb)

  end subroutine TProgramVariables_init


  !> Interpret the information stored in detailed.xml.
  subroutine readDetailed(this, detailed, tGroundState, kPointsWeights)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Pointer to the node, containing the information
    type(hsd_table), pointer :: detailed

    !> Wether to look for ground state occupations (True) or excited (False)
    logical, intent(in) :: tGroundState

    !> The k-points and weights
    real(dp), intent(out), allocatable :: kPointsWeights(:,:)

    !! Pointers to input nodes
    type(hsd_table), pointer :: tmp, occ, spin

    !! Nr. of K-points
    integer :: nKPoint

    !! Nr. of spins
    integer :: nSpin

    !! Nr. of states
    integer :: nState

    !! Auxiliary variables
    integer :: iSpin, iKpoint
    integer :: stat

    call hsd_get(detailed, "Identity", this%input%identity, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(detailed, "Missing required value: 'Identity'")
    call hsd_get_table(detailed, "Geometry", tmp, stat, auto_wrap=.true.)
    if (.not. associated(tmp)) call dftbp_error(detailed, "Missing required block: 'Geometry'")
    call readGeometry(this%input%geo, tmp)

    call hsd_get(detailed, "Real", this%input%tRealHam, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(detailed, "Missing required value: 'Real'")
    call hsd_get(detailed, "NrOfKPoints", nKPoint, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(detailed, "Missing required value: 'NrOfKPoints'")
    call hsd_get(detailed, "NrOfSpins", nSpin, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(detailed, "Missing required value: 'NrOfSpins'")
    call hsd_get(detailed, "NrOfStates", nState, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(detailed, "Missing required value: 'NrOfStates'")
    call hsd_get(detailed, "NrOfOrbitals", this%input%nOrb, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(detailed, "Missing required value: 'NrOfOrbitals'")

    allocate(kPointsWeights(4, nKPoint))
    allocate(this%input%occupations(nState, nKPoint, nSpin))

    block
      integer :: nrows, ncols
      call hsd_get_matrix(detailed, "KPointsAndWeights", kPointsWeights, nrows, ncols, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(detailed, "Missing required value: 'KPointsAndWeights'")
    end block

    if (tGroundState) then
      call hsd_get_table(detailed, "Occupations", occ, stat, auto_wrap=.true.)
      if (.not. associated(occ)) call dftbp_error(detailed, "Missing required block: 'Occupations'")
      do iSpin = 1, nSpin
        call hsd_get_table(occ, "spin" // i2c(iSpin), spin, stat, auto_wrap=.true.)
        if (.not. associated(spin)) call dftbp_error(occ, "Missing spin block")
        do iKpoint = 1, nKPoint
          block
            real(dp), allocatable :: tmp_occ(:)
            call hsd_get(spin, "k" // i2c(iKpoint), tmp_occ, stat=stat)
            if (stat /= HSD_STAT_OK) call dftbp_error(spin, "Missing k-point occupations")
            this%input%occupations(:, iKpoint, iSpin) = tmp_occ
          end block
        end do
      end do
      do iKpoint = 1, nKPoint
        this%input%occupations(:, iKpoint, :) = this%input%occupations(:, iKpoint, :)&
            & * kPointsWeights(4, iKpoint)
      end do
    else
      call hsd_get_table(detailed, "ExcitedOccupations", occ, stat, auto_wrap=.true.)
      if (.not. associated(occ)) call dftbp_error(detailed, "Missing required block: 'ExcitedOccupations'")
      do iSpin = 1, nSpin
        call hsd_get_table(occ, "spin" // i2c(iSpin), spin, stat, auto_wrap=.true.)
        if (.not. associated(spin)) call dftbp_error(occ, "Missing spin block")
        do iKpoint = 1, nKPoint
          block
            real(dp), allocatable :: tmp_occ(:)
            call hsd_get(spin, "k" // i2c(iKpoint), tmp_occ, stat=stat)
            if (stat /= HSD_STAT_OK) call dftbp_error(spin, "Missing k-point occupations")
            this%input%occupations(:, iKpoint, iSpin) = tmp_occ
          end block
        end do
      end do
      do iKpoint = 1, nKPoint
        this%input%occupations(:, iKpoint, :) = this%input%occupations(:, iKpoint, :)&
            & * kPointsWeights(4, iKpoint)
      end do
    end if

  end subroutine readDetailed


  !> Read in the geometry stored as .xml in internal or .gen format.
  subroutine readGeometry(geo, geonode)

    !> Geometry instance
    type(TGeometry), intent(out) :: geo

    !> Node containing the geometry
    type(hsd_table), pointer :: geonode

    !! Pointers to input nodes
    type(hsd_table), pointer :: child

    !! String buffer instance
    character(len=:), allocatable :: buffer
    integer :: stat

    ! Get first table child for dispatch
    child => null()
    block
      class(hsd_node), pointer :: ch
      integer :: ii
      do ii = 1, geonode%num_children
        call geonode%get_child(ii, ch)
        select type (t => ch)
        type is (hsd_table)
          child => t
          exit
        end select
      end do
    end block
    call hsd_get_name(child, buffer, "#text")

    select case (buffer)
    case ("genformat")
      call readTGeometryGen(child, geo)
      call hsd_clear_children(geonode)
      call writeTGeometryHSD(geonode, geo)
    case ("xyzformat")
      call readTGeometryXyz(child, geo)
      call hsd_clear_children(geonode)
      call writeTGeometryHSD(geonode, geo)
    case ("vaspformat")
      call readTGeometryVasp(child, geo)
      call hsd_clear_children(geonode)
      call writeTGeometryHSD(geonode, geo)
    case default
      call readTGeometryHSD(geonode, geo)
    end select

  end subroutine readGeometry


  !> Interpret the options.
  subroutine readOptions(this, node, nLevel, nKPoint, nSpin, nCached, tShiftGrid)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Node containing the information
    type(hsd_table), pointer :: node

    !> Nr. of states in the calculation
    integer, intent(in) :: nLevel

    !> Nr. of K-points
    integer, intent(in) :: nKPoint

    !> Nr. of spins
    integer, intent(in) :: nSpin

    !> Nr. of cached grids
    integer, intent(out) :: nCached

    !> If grid should be shifted by half a cell
    logical, intent(out) :: tShiftGrid

    !! Pointer to the nodes, containing the information
    type(hsd_table), pointer :: subnode, field, value

    !! String buffer instances
    character(len=:), allocatable :: buffer, modifier

    !! Onedimensional integer-valued index list
    type(TIntArray), allocatable :: indexBuffer(:)

    !! Id of calculation at hand
    integer :: curId

    !! If current level is found be calculated explicitely
    logical :: tFound

    !! Warning issued, if the detailed.xml id does not match the eigenvector id
    character(len=63) :: warnId(3) = [&
        & "The external files you are providing differ from those provided",&
        & "when this input file was generated. The results you obtain with",&
        & "the current files could therefore be different.                "]

    !! Auxiliary variables
    integer :: ii, iLevel, iKPoint, iSpin, iAtom, iSpecies
    real(dp) :: tmpvec(3), minvals(3), maxvals(3)
    real(dp), allocatable :: mcutoffs(:)
    real(dp) :: minEdge
    integer :: stat

    ! Warning, if processed input is read in, but eigenvectors are different
    call hsd_get_or_set(node, "Identity", curId, this%input%identity)
    if (curId /= this%input%identity) then
      call warning(warnId)
    end if

    call hsd_get_or_set(node, "TotalChargeDensity", this%opt%tPlotTotChrg, .false.)

    if (nSpin == 2) then
      call hsd_rename_child(node, "TotalSpinPolarization", "TotalSpinPolarisation")
      call hsd_get_or_set(node, "TotalSpinPolarisation", this%opt%tPlotTotSpin, .false.)
    else
      this%opt%tPlotTotSpin = .false.
    end if

    call hsd_get_or_set(node, "TotalChargeDifference", this%opt%tPlotTotDiff, .false.)
    call hsd_get_or_set(node, "TotalAtomicDensity", this%opt%tPlotAtomDens, .false.)
    call hsd_get_or_set(node, "ChargeDensity", this%opt%tPlotChrg, .false.)
    call hsd_get_or_set(node, "ChargeDifference", this%opt%tPlotChrgDiff, .false.)

    this%opt%tCalcTotChrg = this%opt%tPlotTotChrg .or. this%opt%tPlotTotSpin&
        & .or. this%opt%tPlotTotDiff .or. this%opt%tPlotChrgDiff
    this%opt%tCalcAtomDens = this%opt%tPlotTotDiff .or. this%opt%tPlotChrgDiff&
        & .or. this%opt%tPlotAtomDens

    call hsd_get_or_set(node, "RealComponent", this%opt%tPlotReal, .false.)
    call hsd_get_or_set(node, "ImagComponent", this%opt%tPlotImag, .false.)
    call hsd_get_table(node, "ImagComponent", field)

    if (this%opt%tPlotImag .and. this%input%tRealHam) then
      call dftbp_warning(field, "Wave functions are real, no imaginary part will be plotted")
      this%opt%tPlotImag = .false.
    end if

    call hsd_get(node, "PlottedLevels", buffer, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "PlottedLevels must be present")
    call getSelectedIndices(node, buffer, [1, nLevel], this%opt%plottedLevels)

    if (this%input%geo%tPeriodic) then
      call hsd_get(node, "PlottedKPoints", buffer, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "PlottedKPoints must be present")
      call getSelectedIndices(node, buffer, [1, nKPoint], this%opt%plottedKPoints)
    else
      allocate(this%opt%plottedKPoints(1))
      this%opt%plottedKPoints(1) = 1
    end if

    call hsd_get(node, "PlottedSpins", buffer, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "PlottedSpins must be present")
    call getSelectedIndices(node, buffer, [1, nSpin], this%opt%plottedSpins)

    ! Create the list of the levels, which must be calculated explicitely
    allocate(indexBuffer(0))
    do iSpin = 1, nSpin
      do iKPoint = 1, nKPoint
        do iLevel = 1, nLevel
          tFound = any(this%opt%plottedLevels == iLevel)&
              & .and. any(this%opt%plottedKPoints == iKPoint)&
              & .and. any(this%opt%plottedSpins == iSpin)
          if ((.not. tFound) .and. this%opt%tCalcTotChrg) then
            tFound = this%input%occupations(iLevel, iKPoint, iSpin) > 1e-08_dp
          end if
          if (tFound) then
            indexBuffer = [indexBuffer, TIntArray([iLevel, iKPoint, iSpin])]
          end if
        end do
      end do
    end do

    if (size(indexBuffer) == 0) then
      call error("No levels specified for plotting")
    end if

    allocate(this%loc%levelIndex(3, size(indexBuffer)))
    do ii = 1, size(indexBuffer)
      this%loc%levelIndex(:, ii) = indexBuffer(ii)%data
    end do
    deallocate(indexBuffer)

    call hsd_get_or_set(node, "NrOfCachedGrids", nCached, 1)
    call hsd_get_table(node, "NrOfCachedGrids", field)

    if (nCached < 1 .and. nCached /= -1) then
      call dftbp_error(field, "Value must be -1 or greater than zero.")
    end if

    if (nCached == -1) then
      nCached = size(this%loc%levelIndex, dim=2)
    end if

    ! Plotted region: if last (and hopefully only) childnode is not an allowed method -> assume
    ! explicit setting, parse the node "PlottedRegion" for the appropriate children.
    call hsd_get_table(node, "PlottedRegion", subnode, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "PlottedRegion must be present")
    ! Get first table child for dispatch
    value => null()
    block
      class(hsd_node), pointer :: ch
      integer :: jj
      do jj = 1, subnode%num_children
        call subnode%get_child(jj, ch)
        select type (t => ch)
        type is (hsd_table)
          value => t
          exit
        end select
      end do
    end block
    call hsd_get_name(value, buffer, "#text")

    select case (buffer)
    case ("unitcell")
      !! Unit cell for the periodic case, smallest possible cuboid for cluster
      if (this%input%geo%tPeriodic) then
        this%opt%origin(:) = [0.0_dp, 0.0_dp, 0.0_dp]
        this%opt%boxVecs(:,:) = this%input%geo%latVecs(:,:)
      else
        call hsd_get_or_set(value, "MinEdgeLength", minEdge, 1.0_dp)
        call hsd_get_table(value, "MinEdgeLength", field)
        if (minEdge < 0.0_dp) then
          call dftbp_error(field, "Minimal edge length must be positive")
        end if
        this%opt%origin = minval(this%input%geo%coords, dim=2)
        tmpvec = maxval(this%input%geo%coords, dim=2) - this%opt%origin
        do ii = 1, 3
          if (tmpvec(ii) < minEdge) then
            this%opt%origin(ii) = this%opt%origin(ii) - 0.5_dp * (minEdge - tmpvec(ii))
            tmpvec(ii) = minEdge
          end if
        end do
        this%opt%boxVecs(:,:) = 0.0_dp
        do ii = 1, 3
          this%opt%boxVecs(ii, ii) = tmpvec(ii)
        end do
      end if

    case ("optimalcuboid")
      ! Determine optimal cuboid, so that no basis function leaks out
      call hsd_get_or_set(value, "MinEdgeLength", minEdge, 1.0_dp)
      call hsd_get_table(value, "MinEdgeLength", field)
      if (minEdge < 0.0_dp) then
        call dftbp_error(field, "Minimal edge length must be positive")
      end if
      allocate(mcutoffs(this%input%geo%nSpecies))
      do iSpecies = 1 , this%input%geo%nSpecies
        mcutoffs(iSpecies) = maxval(this%basis%basis(iSpecies)%cutoffs)
      end do
      minvals = this%input%geo%coords(:,1)
      maxvals = this%input%geo%coords(:,1)
      do iAtom = 1, this%input%geo%nAtom
        iSpecies = this%input%geo%species(iAtom)
        maxvals(:) = max(maxvals, this%input%geo%coords(:, iAtom) + mcutoffs(iSpecies))
        minvals(:) = min(minvals, this%input%geo%coords(:, iAtom) - mcutoffs(iSpecies))
      end do
      this%opt%origin(:) = minvals(:)
      tmpvec(:) = maxvals(:) - minvals(:)
      do ii = 1, 3
        if (tmpvec(ii) < minEdge) then
          this%opt%origin(ii) = this%opt%origin(ii) - 0.5_dp * (minEdge - tmpvec(ii))
          tmpvec(ii) = minEdge
        end if
      end do
      this%opt%boxVecs(:,:) = 0.0_dp
      do ii = 1, 3
        this%opt%boxVecs(ii, ii) = tmpvec(ii)
      end do

    case ("origin","box")
      ! Those nodes are part of an explicit specification -> explitic specif
      block
        real(dp), allocatable :: tmp_boxVecs(:,:)
        integer :: nrows, ncols
        call hsd_get_matrix(subnode, "Box", tmp_boxVecs, nrows, ncols, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(subnode, "Box must be present")
        this%opt%boxVecs(:,:) = tmp_boxVecs
      end block
      call hsd_get_table(subnode, "Box", field, auto_wrap=.true.)
      if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
      call convertUnitHsd(modifier, lengthUnits, field, this%opt%boxVecs)
      if (abs(determinant33(this%opt%boxVecs)) < 1e-08_dp) then
        call dftbp_error(field, "Vectors are linearly dependent")
      end if
      block
        real(dp), allocatable :: tmp_origin(:)
        call hsd_get(subnode, "Origin", tmp_origin, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(subnode, "Origin must be present")
        this%opt%origin(:) = tmp_origin
      end block
      call hsd_get_table(subnode, "Origin", field, auto_wrap=.true.)
      if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
      call convertUnitHsd(modifier, lengthUnits, field, this%opt%origin)

    case default
      ! Object with unknown name passed
      call dftbp_error(value, "Invalid element name")
    end select

    ! Replace existing PlottedRegion definition
    call hsd_remove_child(node, "PlottedRegion")
    block
      type(hsd_table) :: pr_table
      call new_table(pr_table, "PlottedRegion")
      call node%add_child(pr_table)
    end block
    call hsd_get_table(node, "PlottedRegion", field)
    call hsd_set(field, "Origin", this%opt%origin)
    call hsd_set(field, "Box", this%opt%boxVecs)

    block
      integer, allocatable :: tmp_nPoints(:)
      call hsd_get(node, "NrOfPoints", tmp_nPoints, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "NrOfPoints must be present")
      this%opt%nPoints(:) = tmp_nPoints
    end block
    call hsd_get_table(node, "NrOfPoints", field)

    if (any(this%opt%nPoints <= 0)) then
      call dftbp_error(field, "Specified numbers must be greater than zero")
    end if

    call hsd_get_or_set(node, "ShiftGrid", tShiftGrid, .true.)

    if (this%input%geo%tPeriodic) then
      call hsd_get_or_set(node, "FoldAtomsToUnitCell", this%opt%tFoldCoords, .false.)
      call hsd_get_or_set(node, "FillBoxWithAtoms", this%opt%tFillBox, .false.)
      this%opt%tFoldCoords = this%opt%tFoldCoords .or. this%opt%tFillBox
    else
      this%opt%tFillBox = .false.
      this%opt%tFoldCoords = .false.
    end if

    block
      integer, allocatable :: tmp_repeatBox(:)
      call hsd_get_or_set(node, "RepeatBox", tmp_repeatBox, [1, 1, 1])
      this%opt%repeatBox(:) = tmp_repeatBox
    end block
    call hsd_get_table(node, "RepeatBox", field)

    if (.not. all(this%opt%repeatBox > 0)) then
      call dftbp_error(field, "Indexes must be greater than zero")
    end if

    call hsd_get_or_set(node, "Verbose", this%opt%tVerbose, .false.)

    call readBinaryAccessTypes(node, this%opt%binaryAccessTypes)

  end subroutine readOptions


  !> Reads in the basis related informations.
  subroutine readBasis(this, node, speciesNames)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Node containing the basis definition
    type(hsd_table), pointer :: node

    !> Names of the species for which the basis should be read in
    character(len=*), intent(in) :: speciesNames(:)

    !! Name of current species
    character(len=len(speciesNames)) :: speciesName

    !! Input node instance, containing the information
    type(hsd_table), pointer :: speciesNode

    !! Total number of species in the system
    integer :: nSpecies

    !! Auxiliary variable
    integer :: ii
    integer :: stat

    nSpecies = size(speciesNames)

    @:ASSERT(nSpecies > 0)

    call hsd_get(node, "Resolution", this%basis%basisResolution, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Resolution must be present")

    allocate(this%basis%basis(nSpecies))
    allocate(this%aNr%atomicNumbers(nSpecies))

    do ii = 1, nSpecies
      speciesName = speciesNames(ii)
      call hsd_get_table(node, speciesName, speciesNode, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Species " // trim(speciesName) // " must be present")
      call readSpeciesBasis(speciesNode, this%basis%basisResolution, this%basis%basis(ii))
      this%aNr%atomicNumbers(ii) = this%basis%basis(ii)%atomicNumber
    end do

  end subroutine readBasis


  !> Read in basis function for a species.
  subroutine readSpeciesBasis(node, basisResolution, spBasis)

    !> Node containing the basis definition for a species
    type(hsd_table), pointer :: node

    !> Grid distance for discretising the basis functions
    real(dp), intent(in) :: basisResolution

    !> Contains the basis on return
    type(TSpeciesBasis), intent(out) :: spBasis

    !! Input node instances, containing the information
    type(hsd_table), pointer :: tmpNode, child

    !! Node list instance
    type(hsd_table_ptr), allocatable :: children(:)

    !! Basis coefficients and exponents
    real(dp), allocatable :: coeffs(:), exps(:)

    !! Auxiliary variable
    integer :: ii
    integer :: stat

    call hsd_get(node, "AtomicNumber", spBasis%atomicNumber, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "AtomicNumber must be present")
    call hsd_get_child_tables(node, "Orbital", children)
    spBasis%nOrb = size(children)

    if (spBasis%nOrb < 1) then
      call dftbp_error(node, "Missing orbital definitions")
    end if

    allocate(spBasis%angMoms(spBasis%nOrb))
    allocate(spBasis%occupations(spBasis%nOrb))
    allocate(spBasis%stos(spBasis%nOrb))
    allocate(spBasis%cutoffs(spBasis%nOrb))

    do ii = 1, spBasis%nOrb
      tmpNode => children(ii)%ptr
      call hsd_get(tmpNode, "AngularMomentum", spBasis%angMoms(ii), stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(tmpNode, "AngularMomentum must be present")
      call hsd_get(tmpNode, "Occupation", spBasis%occupations(ii), stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(tmpNode, "Occupation must be present")
      call hsd_get(tmpNode, "Cutoff", spBasis%cutoffs(ii), stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(tmpNode, "Cutoff must be present")
      call hsd_get(tmpNode, "Exponents", exps, stat=stat)
      if (stat /= HSD_STAT_OK .or. size(exps) == 0) then
        call dftbp_error(tmpNode, "Missing exponents")
      end if
      call hsd_get(tmpNode, "Coefficients", coeffs, stat=stat)
      if (stat /= HSD_STAT_OK .or. size(coeffs) == 0) then
        call dftbp_error(tmpNode, "Missing coefficients")
      end if
      if (mod(size(coeffs), size(exps)) /= 0) then
        call dftbp_error(tmpNode, "Number of coefficients incompatible with number of exponents")
      end if
      call TSlaterOrbital_init(spBasis%stos(ii), reshape(coeffs, [size(coeffs) / size(exps),&
          & size(exps)]), exps, ii - 1, basisResolution, spBasis%cutoffs(ii))
      deallocate(exps, coeffs)
    end do

  end subroutine readSpeciesBasis


  !> Checks, if the eigenvector file has the right identity number.
  subroutine checkEigenvecs(fileName, identity)

    !> File to check
    character(len=*), intent(in) :: fileName

    !> Identity number
    integer, intent(in) :: identity

    type(TFileDescr) :: fd
    integer :: id, iostat

    call openFile(fd, fileName, mode="rb", iostat=iostat)

    if (iostat /= 0) then
      call error("Can't open file '" // trim(fileName) // "'.")
    end if

    read(fd%unit) id

    if (id /= identity) then
      call error("Ids for eigenvectors ("// i2c(id) //") and xml-input ("// i2c(identity) // &
          & ") don't match.")
    end if

    call closeFile(fd)

  end subroutine checkEigenvecs

end module waveplot_initwaveplot
