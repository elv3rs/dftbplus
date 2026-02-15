!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module phonons_initphonons
  use dftbp_common_accuracy, only : dp, lc, mc
  use dftbp_common_atomicmass, only : getAtomicMass
  use dftbp_common_constants, only : amu__au
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : closeFile, openFile, TFileDescr
  use dftbp_common_globalenv, only : stdOut, tIoProc
  use dftbp_common_status, only : TStatus
  use dftbp_common_unitconversion, only : energyUnits, lengthUnits
  use dftbp_dftb_periodic, only : getCellTranslations, getNrOfNeighboursForAll, getSuperSampling,&
      & TNeighbourList, TNeighbourlist_init, updateNeighbourList
  use dftbp_io_charmanip, only : i2c, tolower, unquote
  use dftbp_io_hsdutils, only : dftbp_error, getSelectedAtomIndices, getSelectedIndices
  use dftbp_io_unitconv, only : convertUnitHsd
  use hsd, only : hsd_warn_unprocessed, MAX_WARNING_LEN, hsd_error_t, hsd_dump,&
      & hsd_table_ptr, hsd_get_child_tables, hsd_get, hsd_get_or_set, hsd_set,&
      & hsd_get_table, HSD_STAT_OK, hsd_node, hsd_get_matrix, &
      & hsd_get_name, hsd_get_inline_text
  use hsd_data, only : data_load, DATA_FMT_AUTO, hsd_table, new_table
  use dftbp_io_message, only : error, warning
  use dftbp_io_tokenreader, only : getNextToken, TOKEN_OK
  use dftbp_math_simplealgebra, only : determinant33
  use dftbp_transport_negfvars, only : ContactInfo, TNEGFtundos, TTransPar
  use dftbp_type_oldskdata, only : readFromFile, TOldSKData
  use dftbp_type_typegeometryhsd, only : readTGeometryGen, readTGeometryHSD, TGeometry
  use dftbp_type_wrappedintr, only : TWrappedInt1
  implicit none
  private

  character(len=*), parameter :: rootTag = "phonons"
  character(len=*), parameter :: autotestTag = "autotest.tag"
  character(len=*), parameter :: hsdInput = "phonons_in.hsd"
  character(len=*), parameter :: hsdParsedInput = "phonons_pin.hsd"
  character(len=*), parameter :: xmlInput = "phonons_in.xml"
  character(len=*), parameter :: xmlParsedInput = "phonons_pin.xml"

  public :: initProgramVariables, destructProgramVariables
  public :: TPdos, autotestTag

  type TPdos
    type(TWrappedInt1), allocatable :: iAtInRegion(:)
    character(lc), allocatable :: regionLabels(:)
  end type TPdos

  !> Wrapper type for an array of character(lc) strings
  type :: TCharLcArray
    character(lc), allocatable :: items(:)
  end type TCharLcArray

  !> Identity of the run
  integer, public :: identity

  !> Geometry container
  type(TGeometry), public :: geo

  !> Projected dos infos
  type(TPdos), public :: pdos

  !> Container of transport parameters
  type(TTransPar), public :: transpar

  !> Container of tunneling parameters
  type(TNEGFtundos), public :: tundos

  !> Verbose flag
  logical, public :: tVerbose

  !> Core Variables
  real(dp), allocatable, public :: atomicMasses(:)
  real(dp), allocatable, public :: dynMatrix(:,:)
  integer, allocatable, public :: iMovedAtoms(:)
  integer, public :: nMovedAtom, nAtomUnitCell

  !> Kpoints information
  real(dp), allocatable, public :: kPoint(:,:), kWeight(:)
  integer, public :: nKPoints

  !> Maps atom index in central cell
  integer, allocatable, public  :: Img2CentCell(:)

  !> Neighbor list
  type(TNeighbourList), public :: neighbourList

  !> Number of neighbors
  integer, allocatable, target, public :: nNeighbour(:)

  !> Cutoff for Hessian
  real(dp), public :: cutoff

  !> Temperature range
  real(dp), public :: TempMin, TempMax, TempStep

  !> Modes to analyze (e.g., longitudinal, transverse, in-plane, etc)
  integer, public :: selTypeModes

  !> Order=2 means harmonic, 3 is anharmonic 3rd order, etc.
  integer, public :: order

  !> Atomic temperature
  real(dp), public :: atTemperature

  !> Whether modes should be computed
  logical, public :: tCompModes

  !> Whether modes should be plotted
  logical, public :: tPlotModes

  !> Whether modes should be animated
  logical, public :: tAnimateModes

  !>
  logical, public :: tXmakeMol

  !>
  logical, public :: tTransport

  !> Whether phonon dispersions should be computed
  logical, public :: tPhonDispersion

  !> Number of repeated cells along lattice vectors
  integer, public :: nCells(3)

  !> Whether phonon dispersions should be computed
  character(4), public :: outputUnits

  !> Whether taggedoutput should be written
  logical, public :: tWriteTagged

  !> Which phonon modes to animate
  integer, allocatable, public :: modesToPlot(:)

  !> Number of phonon modes to animate
  integer, public :: nModesToPlot

  !> Number of cycles in mode animation
  integer, public :: nCycles

  !> Number of steps in mode animation
  integer, public, parameter :: nSteps = 10

  !> Version of the current parser
  integer, parameter :: parserVersion = 4

  !> Version of the oldest parser for which compatibility is maintained
  integer, parameter :: minVersion = 4


  !> Container type for parser related flags.
  type TParserFlags

    !> Stop after parsing?
    logical :: tStop

    !> Continue despite unprocessed nodes
    logical :: tIgnoreUnprocessed

    !> Write parsed HSD input
    logical :: tWriteHSD

    !> Write TaggedOutput

    logical :: tWriteTagged
  end type TParserFlags


  !> Constants parameters
  type TModeEnum
    !> Along x,y,z
    integer :: ALLMODES = 0
    !> Along x
    integer :: XX = 1
    !> Along y
    integer :: YY = 2
    !> Along z
    integer :: ZZ = 3
    !> Along z (assume transport along z)
    integer :: LONGITUDINAL = 4
    !> Along x,y
    integer :: TRANSVERSE = 5
    !> Along x,z (assume 2D structure in x-z)
    integer :: INPLANE = 6
    !> Along y
    integer :: OUTOFPLANE = 7
  end type TModeEnum

  type(TModeEnum), public, parameter :: modeEnum = TModeEnum()


contains

  !> Initialise program variables
  subroutine initProgramVariables(env)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    ! locals
    type(hsd_table), pointer :: hsdTree, root, node, tmp
    type(hsd_table), pointer :: child, value
    character(len=:), allocatable :: buffer, buffer2, modif
    integer :: inputVersion
    integer :: ii, iSp1, iAt
    logical :: tHSD, reqMass, tBadKPoints
    real(dp), allocatable :: speciesMass(:)
    integer :: nDerivs, nGroups
    type(TParserflags) :: parserFlags

    integer :: cubicType, quarticType
    integer :: stat

    write(stdOut, "(/, A)") "Starting initialization..."
    write(stdOut, "(A80)") repeat("-", 80)

    nGroups = 1
#:if WITH_MPI
    call env%initMpi(nGroups)
#:endif

    !! Read in input file as HSD or XML.
    write(stdOut, "(A)") "Interpreting input file '" // hsdInput // "'"
    write(stdOut, "(A)") repeat("-", 80)
    block
      type(hsd_error_t), allocatable :: hsdError
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
    call hsd_get_table(hsdTree, rootTag, root, stat=stat)
    if (stat /= HSD_STAT_OK) call error("Missing root tag '" // rootTag // "'")

    !! Check if input version is the one, which we can handle
    !! Handle parser options
    call hsd_get_table(root, "Options", child)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: opt_table
        call new_table(opt_table, "Options")
        call root%add_child(opt_table)
      end block
      call hsd_get_table(root, "Options", child)
    end if
    call readOptions(child, root, parserFlags)

    call hsd_get_table(root, "Geometry", tmp, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(root, "Geometry must be present")
    call readGeometry(tmp, geo)

    ! Read Transport block
    ! This defines system partitioning
    tTransport = .false.
    call hsd_get_table(root, "Transport", child)
    if (associated(child)) then
      tTransport = .true.
      call readTransportGeometry(child, geo, transpar)
    end if

    call hsd_get_or_set(root, "Atoms", buffer2, "1:-1")
    call hsd_get_table(root, "Atoms", child)
    call getSelectedAtomIndices(child, buffer2, geo%speciesNames, geo%species, iMovedAtoms)
    nMovedAtom = size(iMovedAtoms)

    tCompModes = .false.
    call hsd_get_table(root, "ComputeModes", node)
    if (associated(node)) then
      tCompModes = .true.
      tPlotModes = .false.
      nModesToPlot = 0
      tAnimateModes = .false.
      tXmakeMol = .false.
      call hsd_get_table(root, "DisplayModes", node)
      if (associated(node)) then
        tPlotModes = .true.
        call hsd_get_or_set(node, "PlotModes", buffer2, "1:-1")
        call hsd_get_table(node, "PlotModes", child)
        call getSelectedIndices(child, buffer2, [1, 3 * nMovedAtom], modesToPlot)
        nModesToPlot = size(modesToPlot)
        call hsd_get_or_set(node, "Animate", tAnimateModes, .true.)
        call hsd_get_or_set(node, "XMakeMol", tXmakeMol, .true.)
      end if
    end if

    if (tAnimateModes.and.tXmakeMol) then
      nCycles = 1
    else
      nCycles = 3
    end if

    ! Reading K-points for Phonon Dispersion calculation
    tPhonDispersion = .false.
    call hsd_get_table(root, "PhononDispersion", node)
    if  (associated(node))  then
      tPhonDispersion = .true.
      block
        integer, allocatable :: tmpCells(:)
        call hsd_get(node, "supercell", tmpCells, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing supercell data")
        nCells(:) = tmpCells
      end block
      nAtomUnitCell = geo%nAtom/(nCells(1)*nCells(2)*nCells(3))
      call hsd_get_or_set(node, "outputUnits", buffer, "H")
      select case(trim(buffer))
      case("H", "eV" , "meV", "THz", "cm")
        outputUnits=trim(buffer)
      case default
        call dftbp_error(node,"Unknown outputUnits "//trim(buffer))
      end select
      call readKPoints(node, geo, tBadKpoints)
    end if

    ! Read the atomic masses from SlaterKosterFiles or Masses
    allocate(speciesMass(geo%nSpecies))
    call hsd_get_table(root, "Masses", node, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(root, "Masses must be present")
    if ( associated(node) ) then
      call hsd_get_table(node, "SlaterKosterFiles", value)
      if ( associated(value) ) then
        call readSKfiles(value, geo, speciesMass)
      else
        call readMasses(node, geo, speciesMass)
      endif
    endif
    allocate(atomicMasses(nMovedAtom))
    do iAt = 1, nMovedAtom
      atomicMasses(iAt) = speciesMass(geo%species(iMovedAtoms(iAt)))
    end do
    deallocate(speciesMass)

    ! --------------------------------------------------------------------------------------
    ! Reading Hessian block parameters
    ! --------------------------------------------------------------------------------------
    call hsd_get_table(root, "Hessian", node, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(root, "Hessian must be present")
    ! cutoff used to cut out interactions
    call hsd_get_or_set(node, "Cutoff", cutoff, 9.45_dp)
    call hsd_get_table(node, "Cutoff", value, auto_wrap=.true.)
    if (allocated(value%attrib)) then; modif = value%attrib; else; modif = ""; end if
    call convertUnitHsd(modif, lengthUnits, value, cutoff)

    ! Reading the actual Hessian matrix
    call hsd_get_table(node, "Matrix", child, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Matrix must be present")
    ! Get first table child for dispatch
    value => null()
    block
      class(hsd_node), pointer :: ch
      integer :: jj
      do jj = 1, child%num_children
        call child%get_child(jj, ch)
        select type (t => ch)
        type is (hsd_table)
          value => t
          exit
        end select
      end do
    end block
    if (associated(value)) then
      call hsd_get_name(value, buffer, "#text")
    else
      buffer = "#text"
    end if
    select case(trim(buffer))
    case ("dftb")
      call readDftbHessian(value)
    case ("dynmatrix")
      call readDynMatrix(value)
    case ("cp2k")
      call readCp2kHessian(value)
    case default
      call dftbp_error(node,"Unknown Hessian type "//buffer)
    end select

    ! --------------------------------------------------------------------------------------

    ! --------------------------------------------------------------------------------------
    ! Reading cubic forces
    ! --------------------------------------------------------------------------------------
    order = 2
    call hsd_get_table(root, "Cubic", child)
    if (associated(child)) then
      order = 3
      call hsd_get(child, "Matrix", buffer, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child, "Matrix must be present")
      select case(trim(buffer))
      case ("gaussian")
        cubicType = 1
      case default
        call dftbp_error(root,"Unknown Cubic forces type "//buffer)
      end select
    end if

    call buildNeighbourList()

    call hsd_get_table(root, "Analysis", child)

    if (associated(child)) then
      if (tPhonDispersion) then
         call dftbp_error(root, "Analysis and PhononDispersion cannot coexist")
      end if
      call readAnalysis(child, geo, pdos, tundos, transpar, atTemperature)
    endif

    !! Issue warning about unprocessed nodes
    write(stdOut, "(/, A)") "check unprocessed nodes..."
    if (.not. parserFlags%tIgnoreUnprocessed) then
      block
        character(len=MAX_WARNING_LEN), allocatable :: warnings(:)
        integer :: ii
        call hsd_warn_unprocessed(root, warnings)
        do ii = 1, size(warnings)
          call warning(trim(warnings(ii)))
        end do
      end block
    end if

    !! Dump processed tree in HSD and XML format
    if (tIoProc .and. parserFlags%tWriteHSD) then
      block
        type(hsd_error_t), allocatable :: hsdErr
        call hsd_dump(hsdTree, hsdParsedInput, hsdErr)
        if (allocated(hsdErr)) call error("Error writing HSD file '" // hsdParsedInput // "'")
      end block
      write(stdOut, '(/,/,A)') "Processed input in HSD format written to '" &
          &// hsdParsedInput // "'"
    end if

    !! Stop, if only parsing is required
    if (parserFlags%tStop) then
      call error("Keyword 'StopAfterParsing' is set to Yes. Stopping.")
    end if

    write(stdOut, "(/, A)") "Initialization done..."

  end subroutine initProgramVariables

  !!* destruct the program variables created in initProgramVariables
  subroutine destructProgramVariables()
    deallocate(atomicMasses)
    deallocate(dynMatrix)
    if (allocated(iMovedAtoms)) then
      deallocate(iMovedAtoms)
    end if
    if (allocated(modesToPlot)) then
      deallocate(modesToPlot)
    end if
    write(stdOut, "(/,A)") repeat("=", 80)
  end subroutine destructProgramVariables


  !!* Read in parser options (options not passed to the main code)
  !!* @param node Node to get the information from
  !!* @param root Root of the entire tree (in the case it must be converted)
  !!* @param flags Contains parser flags on exit.
  subroutine readOptions(node, root, flags)
    type(hsd_table), pointer :: node
    type(hsd_table), pointer :: root
    type(TParserFlags), intent(out) :: flags

    integer :: inputVersion
    type(hsd_table), pointer :: child
    integer :: stat

    !! Check if input needs compatibility conversion.
    call hsd_get_or_set(node, "ParserVersion", inputVersion, parserVersion)
    call hsd_get_table(node, "ParserVersion", child)
    if (inputVersion < 1 .or. inputVersion > parserVersion) then
      call dftbp_error(child, "Invalid parser version (" // i2c(inputVersion) // ")")
    elseif (inputVersion < minVersion) then
      call dftbp_error(child, "Sorry, no compatibility mode for parser version "&
          & // i2c(inputVersion) // " (too old)")
    end if

    call hsd_get_or_set(node, "WriteAutotestTag", tWriteTagged, .false.)
    tWriteTagged = tWriteTagged .and. tIoProc
    call hsd_get_or_set(node, "WriteHSDInput", flags%tWriteHSD, .true.)
    call hsd_get_or_set(node, "StopAfterParsing", flags%tStop, .false.)

    call hsd_get_or_set(node, "IgnoreUnprocessedNodes", &
        &flags%tIgnoreUnprocessed, .false.)

  end subroutine readOptions


  !> Read in the geometry stored as xml in internal or gen format.
  subroutine readGeometry(geonode, geo)

    !> Node containing the geometry
    type(hsd_table), pointer :: geonode

    !> Contains the geometry information on exit
    type(TGeometry), intent(out) :: geo

    type(hsd_table), pointer :: child, value
    character(len=:), allocatable :: buffer

    ! Get first table child for dispatch
    child => geonode
    value => null()
    block
      class(hsd_node), pointer :: ch
      integer :: jj
      do jj = 1, geonode%num_children
        call geonode%get_child(jj, ch)
        select type (t => ch)
        type is (hsd_table)
          value => t
          exit
        end select
      end do
    end block
    if (associated(value)) then
      call hsd_get_name(value, buffer, "#text")
    else
      buffer = "#text"
    end if
    select case (buffer)
    case ("genformat")
      call readTGeometryGen(value, geo)
    case default
      if (associated(value)) value%processed = .false.
      call readTGeometryHSD(child, geo)
    end select

    if (geo%tPeriodic) then
      write(stdOut,*) 'supercell lattice vectors:'
      write(stdOut,*) 'a1:',geo%latVecs(1,:)
      write(stdOut,*) 'a2:',geo%latVecs(2,:)
      write(stdOut,*) 'a3:',geo%latVecs(3,:)
    end if

  end subroutine readGeometry


  !> Read geometry information for transport calculation
  subroutine readTransportGeometry(root, geom, tp)
    type(hsd_table), pointer :: root
    type(TGeometry), intent(inout) :: geom
    type(TTransPar), intent(inout) :: tp

    type(hsd_table), pointer :: pGeom, pDevice, pNode, pTask, pTaskType
    character(len=:), allocatable :: modif
    type(hsd_table), pointer :: pTmp, field
    type(hsd_table_ptr), allocatable :: pNodeList(:)
    integer :: ii, contact
    real(dp) :: acc, contactRange(2), sep
    integer :: stat

    tp%defined = .true.
    tp%tPeriodic1D = .not. geom%tPeriodic
    call hsd_get_table(root, "Device", pDevice, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(root, "Device must be present")
    block
      integer, allocatable :: tmp_idx(:)
      call hsd_get(pDevice, "AtomRange", tmp_idx, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(pDevice, "AtomRange must be present")
      tp%idxdevice(:) = tmp_idx
    end block
    call hsd_get_table(pDevice, "FirstLayerAtoms", pTmp)
    call readFirstLayerAtoms(pTmp, tp%PL, tp%nPLs, tp%idxdevice)
    if (.not.associated(pTmp)) then
      call hsd_set(pDevice, "FirstLayerAtoms", tp%PL)
    end if
    call hsd_get_child_tables(root, "Contact", pNodeList)
    tp%ncont = size(pNodeList)
    if (tp%ncont < 2) then
      call dftbp_error(pGeom, "At least two contacts must be defined")
    end if
    allocate(tp%contacts(tp%ncont))
    !! Parse contact geometry
    call readContacts(pNodeList, tp%contacts, geom, .true.)

  end subroutine readTransportGeometry


  !> Reads settings for the first layer atoms in principal layers
  subroutine readFirstLayerAtoms(pnode, pls, npl, idxdevice, check)

    type(hsd_table), pointer, intent(in) :: pnode

    !> Start atoms in the principal layers
    integer, allocatable, intent(out) :: pls(:)

    !> Number of principal layers
    integer, intent(out) :: npl

    !> Atoms range of the device
    integer, intent(in) :: idxdevice(2)

    !> Optional setting to turn on/off check (defaults to on if absent)
    logical, optional, intent(in) :: check


    integer :: stat
    logical :: checkidx

    checkidx = .true.
    if (present(check)) checkidx = check

    if (associated(pnode)) then
        call hsd_get(pnode, "#text", pls, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(pnode, "Missing first layer atom data")
        npl = size(pls)
        if (checkidx) then
          if (any(pls < idxdevice(1) .or. &
                  pls > idxdevice(2))) then
             call dftbp_error(pnode, "First layer atoms must be between " &
               &// i2c(idxdevice(1)) // " and " // i2c(idxdevice(2)) // ".")
          end if
        end if
      else
         npl = 1
         allocate(pls(npl))
         pls = (/ 1 /)
      end if

  end subroutine readFirstLayerAtoms


  !> Read bias information, used in Analysis and Green's function solver
  subroutine readContacts(pNodeList, contacts, geom, upload)
    type(hsd_table_ptr), intent(in) :: pNodeList(:)
    type(ContactInfo), allocatable, dimension(:), intent(inout) :: contacts
    type(TGeometry), intent(in) :: geom
    logical, intent(in) :: upload

    real(dp) :: contactLayerTol, vec(3)
    integer :: ii, jj
    type(hsd_table), pointer :: field, pNode, pTmp, pWide
    character(len=:), allocatable :: buffer, modif
    integer :: stat

    do ii = 1, size(contacts)

      contacts(ii)%wideBand = .false.
      contacts(ii)%wideBandDos = 0.0

      pNode => pNodeList(ii)%ptr
      call hsd_get(pNode, "Id", buffer, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "Id must be present")
      call hsd_get_table(pNode, "Id", pTmp)
      buffer = tolower(trim(unquote(buffer)))
      if (len(buffer) > mc) then
        call dftbp_error(pTmp, "Contact id may not be longer than " &
            &// i2c(mc) // " characters.")
      end if
      contacts(ii)%name = buffer
      if (any(contacts(1:ii-1)%name == contacts(ii)%name)) then
        call dftbp_error(pTmp, "Contact id '" // trim(contacts(ii)%name) &
            &//  "' already in use")
      end if

      call hsd_get_or_set(pNode, "PLShiftTolerance", contactLayerTol, 1e-5_dp)
      call hsd_get_table(pNode, "PLShiftTolerance", field, auto_wrap=.true.)
      if (allocated(field%attrib)) then; modif = field%attrib; else; modif = ""; end if
      call convertUnitHsd(modif, lengthUnits, field, contactLayerTol)

      block
        integer, allocatable :: tmp_range(:)
        call hsd_get(pNode, "AtomRange", tmp_range, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "AtomRange must be present")
        contacts(ii)%idxrange(:) = tmp_range
      end block
      call hsd_get_table(pNode, "AtomRange", pTmp)
      call getContactVectorII(contacts(ii)%idxrange, geom, ii, pTmp, contactLayerTol, &
                              & contacts(ii)%lattice, contacts(ii)%dir)
      contacts(ii)%length = sqrt(sum(contacts(ii)%lattice**2))

      ! Contact temperatures. Needed
      call hsd_get_or_set(pNode, "Temperature", contacts(ii)%kbT, 0.0_dp)
      call hsd_get_table(pNode, "Temperature", field, auto_wrap=.true.)
      if (allocated(field%attrib)) then; modif = field%attrib; else; modif = ""; end if
      call convertUnitHsd(modif, energyUnits, field, contacts(ii)%kbT)

      if (upload) then
        contacts(ii)%potential = 0.d0

        call hsd_get_or_set(pNode, "wideBand", contacts(ii)%wideBand, .false.)
        if (contacts(ii)%wideBand) then
          call hsd_get_or_set(pNode, 'LevelSpacing', contacts(ii)%wideBandDos, 0.735_dp)
          call hsd_get_table(pNode, 'LevelSpacing', field, auto_wrap=.true.)
          if (allocated(field%attrib)) then; modif = field%attrib; else; modif = ""; end if
          call convertUnitHsd(modif, energyUnits, field,&
                              &contacts(ii)%wideBandDos)
          !WideBandApproximation is defined as energy spacing between levels
          !In the code the inverse value (Density of states) is used
          !Convert the negf input value. Default is 20.e eV
          contacts(ii)%wideBandDos = 1.d0 / contacts(ii)%wideBandDos
        end if
        contacts(ii)%eFermi=0.d0
      end if

    enddo

  end subroutine readContacts


  !> Verification checking of atom ranges and returning contact vector and direction.
  subroutine getContactVectorII(atomrange, geom, id, pContact, plShiftTol, contactVec, contactDir)
    integer, intent(in) :: atomrange(2)
    type(TGeometry), intent(in) :: geom
    integer, intent(in) :: id
    type(hsd_table), pointer :: pContact
    real(dp), intent(in) :: plShiftTol
    real(dp), intent(out) :: contactVec(3)
    integer, intent(out) :: contactDir

    integer :: iStart, iStart2, iEnd
    logical :: mask(3)

    iStart = atomrange(1)
    iEnd = atomrange(2)
    !! Consistency check for the atom ranges
    if (iStart < 1 .or. iEnd < 1 .or. iStart > geom%nAtom &
        &.or. iEnd > geom%nAtom .or. iEnd < iStart) then
      call dftbp_error(pContact, "Invalid atom range '" // i2c(iStart) &
          &// " " // i2c(iEnd) // "', values should be between " // i2c(1) &
          &// " and " // i2c(geom%nAtom) // ".")
    end if
    if (mod(iEnd - iStart + 1, 2) /= 0) then
      call dftbp_error(pContact, "Nr. of atoms in the contact must be even")
    end if

    ! Determining contact vector
    iStart2 = iStart + (iEnd - iStart + 1) / 2
    contactVec = geom%coords(:,iStart) - geom%coords(:,iStart2)
    if (any(sqrt(sum((geom%coords(:,iStart:iStart2-1) - geom%coords(:,iStart2:iEnd) &
        &- spread(contactVec, dim=2, ncopies=iStart2-iStart))**2, dim=1)) > plShiftTol)) then
      write(stdOut,*) 'coords:', geom%coords(:,iStart)
      write(stdOut,*) 'coords:', geom%coords(:,iStart2)
      write(stdOut,*) 'Contact Vector:', contactVec(1:3)
      write(stdOut,*) iStart,iStart2,iEnd
      write(stdOut,*) 'X:'
      write(stdOut,*) ((geom%coords(1,iStart:iStart2-1) - geom%coords(1,iStart2:iEnd)&
          & - spread(contactVec(1), dim=1, ncopies=iStart2-iStart)))
      write(stdOut,*) 'Y:'
      write(stdOut,*) ((geom%coords(2,iStart:iStart2-1) - geom%coords(2,iStart2:iEnd) &
          & - spread(contactVec(2), dim=1, ncopies=iStart2-iStart)))
      write(stdOut,*) 'Z:'
      write(stdOut,*) ((geom%coords(3,iStart:iStart2-1) - geom%coords(3,iStart2:iEnd) &
          &- spread(contactVec(3), dim=1, ncopies=iStart2-iStart)))
      call error("Contact " // i2c(id) &
          &// " does not consist of two rigidly shifted layers."//new_line('a') &
          &// "Check structure or increase PLShiftTolerance.")
    end if

    contactDir = 0

  end subroutine getContactVectorII


  !> Used to read atomic masses from SK files
  subroutine readSKfiles(child, geo, speciesMass)
    type(hsd_table), pointer :: child
    type(TGeometry), intent(in) :: geo
    real(dp), dimension(:) :: speciesMass

    type(TOldSKData) :: skData
    type(TCharLcArray), allocatable :: skFiles(:)
    type(hsd_table), pointer :: value, child2
    character(len=:), allocatable :: buffer, buffer2
    character(lc) :: prefix, suffix, separator, elem1, elem2, strTmp, filename
    integer :: ii, iSp1, stat
    logical :: tLower, tExist

    write(stdOut, "(/, A)") "read atomic masses from sk files..."
    !! Slater-Koster files
    allocate(skFiles(geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
        allocate(skFiles(iSp1)%items(0))
    end do

    ! Get first table child for dispatch
    value => null()
    block
      class(hsd_node), pointer :: ch
      integer :: jj
      do jj = 1, child%num_children
        call child%get_child(jj, ch)
        select type (t => ch)
        type is (hsd_table)
          value => t
          exit
        end select
      end do
    end block
    if (associated(value)) then
      call hsd_get_name(value, buffer, "#text")
    else
      buffer = "#text"
    end if

    select case(buffer)
    case ("type2filenames")
      call hsd_get_or_set(value, "Prefix", buffer2, "")
      prefix = unquote(buffer2)
      call hsd_get_or_set(value, "Suffix", buffer2, "")
      suffix = unquote(buffer2)
      call hsd_get_or_set(value, "Separator", buffer2, "")
      separator = unquote(buffer2)
      call hsd_get_or_set(value, "LowerCaseTypeName", tLower, .false.)
      do iSp1 = 1, geo%nSpecies
        if (tLower) then
          elem1 = tolower(geo%speciesNames(iSp1))
        else
          elem1 = geo%speciesNames(iSp1)
        end if
        strTmp = trim(prefix) // trim(elem1) // trim(separator) &
            &// trim(elem1) // trim(suffix)
        call appendCharLc(skFiles(iSp1)%items, strTmp)
        inquire(file=strTmp, exist=tExist)
        if (.not. tExist) then
          call dftbp_error(value, "SK file with generated name '" &
              &// trim(strTmp) // "' does not exist.")
        end if
      end do
    case default
      if (associated(value)) value%processed = .false.
      do iSp1 = 1, geo%nSpecies
        strTmp = trim(geo%speciesNames(iSp1)) // "-" &
            &// trim(geo%speciesNames(iSp1))
        call hsd_get_table(child, trim(strTmp), child2)
        block
          character(:), allocatable :: skStrArr(:)
          call hsd_get(child, trim(strTmp), skStrArr, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing SK data for " // trim(strTmp))
          ! We can't handle selected shells here (also not needed I guess)
          if (size(skStrArr) /= 1) then
            call dftbp_error(child2, "Incorrect number of Slater-Koster &
                &files")
          end if
          do ii = 1, size(skStrArr)
            strTmp = skStrArr(ii)
            inquire(file=strTmp, exist=tExist)
            if (.not. tExist) then
              call dftbp_error(child2, "SK file '" // trim(strTmp) &
                  &// "' does not exist'")
            end if
            call appendCharLc(skFiles(iSp1)%items, strTmp)
          end do
        end block
      end do
    end select

    do iSp1 = 1, geo%nSpecies
      fileName = skFiles(iSp1)%items(1)
      call readFromFile(skData, fileName, .true.)
      deallocate(skData%skHam)
      deallocate(skData%skOver)
      speciesMass(iSp1) = skData%mass
    end do

    do iSp1 = 1, geo%nSpecies
      if (allocated(skFiles(iSp1)%items)) deallocate(skFiles(iSp1)%items)
    end do
    deallocate(skFiles)

  end subroutine readSKfiles


  subroutine readMasses(value, geo, speciesMass)
    type(hsd_table), pointer :: value
    type(TGeometry), intent(in) :: geo
    real(dp), dimension(:) :: speciesMass

    type(hsd_table), pointer :: child, child2
    character(len=:), allocatable :: modif
    integer :: iSp
    real(dp) :: mass, defmass

    write(stdOut, "(/, A)") "set atomic masses as IUPAC defaults ..."

    do iSp = 1, geo%nSpecies
      defmass = getAtomicMass(trim(geo%speciesNames(iSp)))
      call hsd_get_or_set(value, geo%speciesNames(iSp), mass, defmass)
      call hsd_get_table(value, geo%speciesNames(iSp), child2, auto_wrap=.true.)
      if (allocated(child2%attrib)) then; modif = child2%attrib; else; modif = ""; end if
      speciesMass(iSp) = mass
      write(stdOut,*) trim(geo%speciesNames(iSp)),": ", mass/amu__au, "amu", &
            &SpeciesMass(iSp),"a.u."
    end do

  end subroutine readMasses


  !> K-Points
  subroutine readKPoints(node, geo, tBadIntegratingKPoints)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Error check
    logical, intent(out) :: tBadIntegratingKPoints

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer, modifier
    integer :: ind, ii, jj, kk
    real(dp) :: coeffsAndShifts(3, 4)
    real(dp) :: rTmp3(3)
    integer, allocatable :: tmpI1(:)
    real(dp), allocatable :: kpts(:,:)
    character(lc) :: errorStr
    integer :: stat

    ! Assume SCC can has usual default number of steps if needed
    tBadIntegratingKPoints = .false.

    ! K-Points
    if (geo%tPeriodic) then
      call hsd_get_table(node, "KPointsAndWeights", child, stat=stat, auto_wrap=.true.)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "KPointsAndWeights must be present")
      if (allocated(child%attrib)) then; modifier = child%attrib; else; modifier = ""; end if
      ! Get first table child for dispatch
      value1 => null()
      block
        class(hsd_node), pointer :: ch
        integer :: jjj
        do jjj = 1, child%num_children
          call child%get_child(jjj, ch)
          select type (t => ch)
          type is (hsd_table)
            value1 => t
            exit
          end select
        end do
      end block
      if (associated(value1)) then
        call hsd_get_name(value1, buffer, "#text")
      else
        buffer = "#text"
      end if
      select case(buffer)

      case ("supercellfolding")
        tBadIntegratingKPoints = .false.
        if (len(modifier) > 0) then
          call dftbp_error(child, "No modifier is allowed, if the &
              &SupercellFolding scheme is used.")
        end if
        block
          real(dp), allocatable :: tmp_cs(:,:)
          integer :: nrows, ncols
          call hsd_get_matrix(value1, "", tmp_cs, nrows, ncols, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(value1, "SupercellFolding data must be present")
          coeffsAndShifts(:,:) = tmp_cs
        end block
        if (abs(determinant33(coeffsAndShifts(:,1:3))) - 1.0_dp < -1e-6_dp) then
          call dftbp_error(value1, "Determinant of the supercell matrix must &
              &be greater than 1")
        end if
        if (any(abs(modulo(coeffsAndShifts(:,1:3) + 0.5_dp, 1.0_dp) - 0.5_dp) &
            &> 1e-6_dp)) then
          call dftbp_error(value1, "The components of the supercell matrix &
              &must be integers.")
        end if
        call getSuperSampling(coeffsAndShifts(:,1:3), modulo(coeffsAndShifts(:,4), 1.0_dp),&
            & kPoint, kWeight, reduceByInversion=.true.)

        nKPoints = size(kPoint, dim=2)

      case ("klines")
        ! probably unable to integrate charge for SCC
        tBadIntegratingKPoints = .true.
        block
          character(len=:), allocatable :: text
          integer :: bufInt(1)
          real(dp) :: bufReal(3)
          integer :: iStart, iErr, nItem, nPairs

          call hsd_get(value1, "#text", text, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing k-line data")

          ! First pass: count entries
          nPairs = 0
          iStart = 1
          iErr = TOKEN_OK
          do while (iErr == TOKEN_OK)
            call getNextToken(text, bufInt, iStart, iErr, nItem)
            if (iErr /= TOKEN_OK) exit
            nPairs = nPairs + 1
            call getNextToken(text, bufReal, iStart, iErr, nItem)
          end do

          if (nPairs < 1) then
            call dftbp_error(value1, "At least one line must be specified.")
          end if

          ! Second pass: read values
          allocate(tmpI1(nPairs))
          allocate(kpts(3, 0:nPairs))
          iStart = 1
          do ii = 1, nPairs
            call getNextToken(text, bufInt, iStart, iErr, nItem)
            tmpI1(ii) = bufInt(1)
            call getNextToken(text, bufReal, iStart, iErr, nItem)
            kpts(:, ii) = bufReal
          end do
          kpts(:,0) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
        end block
        if (any(tmpI1 < 0)) then
          call dftbp_error(value1, "Interval steps must be greater equal to &
              &zero.")
        end if
        nKPoints = sum(tmpI1)
        if (nKPoints < 1) then
          call dftbp_error(value1, "Sum of the interval steps must be greater &
              &than zero.")
        end if
        ii = 1
        do while (tmpI1(ii) == 0)
          ii = ii + 1
        end do
        allocate(kPoint(3, nKPoints))
        allocate(kWeight(nKPoints))
        ind = 1
        do jj = ii, size(tmpI1)
          if (tmpI1(jj) == 0) then
            cycle
          end if
          rTmp3 = (kpts(:,jj) - kpts(:,jj-1)) / real(tmpI1(jj), dp)
          do kk = 1, tmpI1(jj)
            kPoint(:,ind) = kpts(:,jj-1) + real(kk, dp) * rTmp3
            ind = ind + 1
          end do
        end do
        kWeight(:) = 1.0_dp
        if (len(modifier) > 0) then
          select case (tolower(modifier))
          case ("relative")
          case ("absolute")
            kPoint(:,:) =  matmul(transpose(geo%latVecs), kPoint)
            kpts(:,:) = matmul(transpose(geo%latVecs), kpts)
          case default
            call dftbp_error(child, "Invalid modifier: '" // modifier &
                &// "'")
          end select
        end if
        deallocate(tmpI1)
        deallocate(kpts)

      case ("#text")

        ! no idea, but assume user knows what they are doing
        tBadIntegratingKPoints = .false.

        block
          real(dp), allocatable :: flatArr(:)
          call hsd_get(child, "#text", flatArr, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing k-point data")
          if (mod(size(flatArr), 4) /= 0) then
            call dftbp_error(child, "K-point data must have 4 values per point")
          end if
          nKPoints = size(flatArr) / 4
          if (nKPoints < 1) then
            call dftbp_error(child, "At least one k-point must be defined.")
          end if
          allocate(kpts(4, nKPoints))
          kpts = reshape(flatArr, [4, nKPoints])
        end block
        if (len(modifier) > 0) then
          select case (tolower(modifier))
          case ("relative")
            continue
          case ("absolute")
            kpts(1:3,:) =  matmul(transpose(geo%latVecs), kpts(1:3,:))
          case default
            call dftbp_error(child, "Invalid modifier: '" // modifier &
                &// "'")
          end select
        end if
        allocate(kPoint(3, nKPoints))
        allocate(kWeight(nKPoints))
        kPoint(:,:) = kpts(1:3, :)
        kWeight(:) = kpts(4, :)
        deallocate(kpts)
      case default
        if (associated(value1)) then
          call dftbp_error(value1, "Invalid K-point scheme")
        else
          call dftbp_error(child, "Invalid K-point scheme")
        end if
      end select
    end if

  end subroutine readKPoints


  subroutine  readKPointsFile(child)
    type(hsd_table),  pointer ::  child
    character(len=:), allocatable :: text

    call hsd_get_inline_text(child, text)
    call readKPointsFile_help(child, text)

  end subroutine  readKPointsFile


  subroutine readKPointsFile_help(child,text)
    type(hsd_table),  pointer ::  child
    character(len=*), intent(in) :: text
    integer :: iStart, iErr=0, ii, iOldStart
    real(dp), dimension(:), allocatable :: tmparray
    real(dp), dimension(:,:), allocatable :: kpts

    iStart = 1
    call getNextToken(text, nKPoints, iStart, iErr)

    allocate(tmparray(4))
    allocate(kpts(4, nKPoints))
    allocate(kPoint(3, nKPoints))
    allocate(kWeight(nKPoints))

    iErr = -2 !TOKEN_ERROR
    iOldStart = iStart
    iStart  = iOldStart

    do ii = 1, nKPoints
      call getNextToken(text, tmparray, iStart, iErr)
      kpts(:, ii) = tmparray(:)
    end do

    do ii = 1, nKPoints
        kPoint(1:3, ii)  = kpts(1:3, ii)
        kWeight(ii)  = kpts(4, ii)
    end do

  end subroutine readKPointsFile_help


  !>  Read DFTB hessian.
  !!
  !! The derivatives matrix must be stored in the following order:
  !!  For the x y z directions of atoms 1..n
  !!    d^2 E        d^2 E       d^2 E       d^2 E        d^2 E
  !!  ---------- + --------- + --------- + ---------- + ---------- +...
  !!  dx_1 dx_1    dy_1 dx_1   dz_1 dx_1   dx_2 dx_1    dy_2 dx_1
  !!
  subroutine readDftbHessian(child)
    type(hsd_table), pointer :: child

    integer :: iCount, jCount, ii, kk, jj, ll
    integer :: nDerivs

    type(TFileDescr) :: fd
    integer ::  n, j1, j2
    type(hsd_table), pointer :: child2
    character(len=:), allocatable :: filename
    logical :: texist
    character(lc) :: strTmp

    call hsd_get_or_set(child, "Filename", filename, "hessian.out")

    ! workaround for NAG7.1 Build 7148 in Debug build
    strTmp = filename
    inquire(file=strTmp, exist=texist )
    if (texist) then
      write(stdOut, "(/, A)") "read dftb hessian '"//trim(filename)//"'..."
    else
      call dftbp_error(child,"Hessian file "//trim(filename)//" does not exist")
    end if

    nDerivs = 3 * nMovedAtom
    allocate(dynMatrix(nDerivs,nDerivs))

    call openFile(fd, trim(filename))
    do ii = 1,  nDerivs
        read(fd%unit,'(4f16.10)') dynMatrix(1:nDerivs,ii)
    end do
    call closeFile(fd)

    ! Note: we read the transpose matrix to avoid temporary arrays (ifort warnings).
    ! It should be symmetric or could be symmetrized here
    dynMatrix = transpose(dynMatrix)
    ! mass weight the Hessian matrix to get the dynamical matrix
    iCount = 0
    do ii = 1, nMovedAtom
      do kk = 1, 3
        iCount = iCount + 1
        jCount = 0
        do jj = 1, nMovedAtom
          do ll = 1, 3
            jCount = jCount + 1
            dynMatrix(jCount,iCount) = dynMatrix(jCount,iCount) &
                & / (sqrt(atomicMasses(ii)) * sqrt(atomicMasses(jj)))
          end do
        end do
      end do
    end do

  end subroutine readDftbHessian


  subroutine readDynMatrix(child)
    type(hsd_table), pointer :: child

    integer :: nDerivs, stat

    nDerivs = 3 * nMovedAtom
    allocate(dynMatrix(nDerivs,nDerivs))

    !! The derivatives matrix must be stored as the following order:
    !!
    !! For the x y z directions of atoms 1..n
    !!   d^2 E        d^2 E       d^2 E       d^2 E        d^2 E
    !! ---------- + --------- + --------- + ---------- + ---------- +...
    !! dx_1 dx_1    dy_1 dx_1   dz_1 dx_1   dx_2 dx_1    dy_2 dx_1

    write(stdOut, "(/, A)") "read dynamical matrix..."

    block
      real(dp), allocatable :: flatArr(:)
      integer :: nRows
      call hsd_get(child, "#text", flatArr, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing derivatives data")
      if (mod(size(flatArr), nDerivs) /= 0) then
        call dftbp_error(child, "wrong number of derivatives supplied")
      end if
      nRows = size(flatArr) / nDerivs
      if (nRows /= nDerivs) then
        call dftbp_error(child, "wrong number of derivatives supplied:" &
            & // i2c(nRows) // " supplied, " &
            & // i2c(nDerivs) // " required.")
      end if
      dynMatrix = reshape(flatArr, [nDerivs, nDerivs])
    end block

  end subroutine readDynMatrix


  subroutine readCp2kHessian(child)
    type(hsd_table), pointer :: child

    integer :: iCount, jCount, ii, kk, jj, ll
    integer :: nDerivs, nBlocks

    type(TFileDescr) :: fd
    real, dimension(:,:), allocatable :: HessCp2k
    integer ::  n, j1, j2,  p,  q
    character(len=:), allocatable :: filename
    logical :: texist
    character(lc) :: strTmp

    nDerivs = 3 * nMovedAtom
    allocate(dynMatrix(nDerivs,nDerivs))

    call hsd_get_or_set(child, "Filename", filename, "hessian.cp2k")
    ! workaround for NAG7.1 Build 7148 in Debug build
    strTmp = filename
    inquire(file=strTmp, exist=texist )
    if (texist) then
      write(stdOut, "(/, A)") "read cp2k hessian '"//trim(filename)//"'..."
    else
      call dftbp_error(child,"Hessian file "//trim(filename)//" does not exist")
    end if

    !! The derivatives matrix must be stored as the following order:
    !!
    !! For the x y z directions of atoms 1..n
    !!   d^2 E        d^2 E       d^2 E       d^2 E        d^2 E
    !! ---------- + --------- + --------- + ---------- + ---------- +...
    !! dx_1 dx_1    dy_1 dx_1   dz_1 dx_1   dx_2 dx_1    dy_2 dx_1

    nBlocks = nDerivs/5.0
    allocate(HessCp2k(nDerivs*nBlocks,5))

    call openFile(fd, trim(filename))
    do  ii  = 1,  nDerivs*nBlocks
      read(fd%unit, *) HessCp2k(ii,1:5)
    end do
    call closeFile(fd)

    do ii = 1,  nBlocks
        do  jj  = 1, nDerivs
            p = 1+5*(ii-1)
            q = 5*ii
          dynMatrix(jj,p:q) =  HessCp2k(jj + nDerivs*(ii-1),1:5)
        end do
    end do

    ! mass weight the Hessian matrix to get the dynamical matrix
    iCount = 0
    do ii = 1, nMovedAtom
      do kk = 1, 3
        iCount = iCount + 1
        jCount = 0
        do jj = 1, nMovedAtom
          do ll = 1, 3
            jCount = jCount + 1
            dynMatrix(jCount,iCount) = dynMatrix(jCount,iCount) &
                & / (sqrt(atomicMasses(ii)) * sqrt(atomicMasses(jj)))
          end do
        end do
      end do
    end do


  end subroutine readCp2kHessian


  !> Subroutine removing entries in the Dynamical Matrix.
  !! Not used because identified as a wrong way
  subroutine selectModes()

    integer :: iCount, jCount, ii, jj, kk, ll

    select case ( selTypeModes )
    case(modeEnum%INPLANE)
      iCount = 0
      do ii = 1, nMovedAtom
        do kk = 1, 3
          iCount = iCount + 1
          jCount = 0
          do jj = 1, nMovedAtom
            do ll = 1, 3
              jCount = jCount + 1
              if (mod(iCount,3).eq.0 .or. mod(jCount,3).eq.0) then
                  dynMatrix(jCount,iCount) = 0.0
              end if
            end do
          end do
        end do
      end do
    case(modeEnum%OUTOFPLANE)
      iCount = 0
      do ii = 1, nMovedAtom
        do kk = 1, 3
          iCount = iCount + 1
          jCount = 0
          do jj = 1, nMovedAtom
            do ll = 1, 3
              jCount = jCount + 1
              if (mod(iCount,3).ne.0 .and. mod(jCount,3).ne.0) then
                  dynMatrix(jCount,iCount) = 0.0
              end if
            end do
          end do
        end do
      end do
    end select

  end subroutine selectModes


  !> Reads the Analysis block.
  subroutine readAnalysis(node, geo, pdos, tundos, transpar, atTemperature)
    type(hsd_table), pointer :: node, pnode
    type(TGeometry), intent(in) :: geo
    type(TPdos), intent(inout) :: pdos
    type(TNEGFTunDos), intent(inout) :: tundos
    type(TTransPar), intent(inout) :: transpar
    real(dp) :: atTemperature, TempRange(2)

    type(hsd_table), pointer :: val, child, field
    character(len=:), allocatable :: modif
    logical :: tBadKpoints
    integer :: stat

    call hsd_get_table(node, "TunnelingAndDOS", child)
    if (associated(child)) then
      if (.not.tTransport) then
        call dftbp_error(node, "Tunneling requires Transport block")
      end if
      call readTunAndDos(child, geo, tundos, transpar, maxval(transpar%contacts(:)%kbT) )
    endif

    !call readKPoints(node, geo, tBadKpoints)
    call hsd_get_table(node, "Conductance", child)
    if (associated(child)) then
      if (.not.tTransport) then
        call dftbp_error(node, "Conductance requires Transport block")
      end if
      block
        real(dp), allocatable :: tmp_range(:)
        call hsd_get(child, "TempRange", tmp_range, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child, "TempRange must be present")
        TempRange(:) = tmp_range
      end block
      call hsd_get_table(child, "TempRange", field, auto_wrap=.true.)
      if (allocated(field%attrib)) then; modif = field%attrib; else; modif = ""; end if
      call convertUnitHsd(modif, energyUnits, field, TempRange)

      call hsd_get(child, "TempStep", TempStep, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child, "TempStep must be present")
      call hsd_get_table(child, "TempStep", field, auto_wrap=.true.)
      if (allocated(field%attrib)) then; modif = field%attrib; else; modif = ""; end if
      call convertUnitHsd(modif, energyUnits, field, TempStep)

       TempMin = TempRange(1)
       TempMax = TempRange(2)
    else
       TempMin = 0.0_dp
       TempMax = 0.0_dp
       TempStep = 1.0_dp
    endif

  end subroutine readAnalysis


  subroutine readPDOSRegions(children, geo, iAtInregion, regionLabels)
    type(hsd_table_ptr), intent(in) :: children(:)
    type(TGeometry), intent(in) :: geo
    type(TWrappedInt1), allocatable, intent(out) :: iAtInRegion(:)
    character(lc), allocatable, intent(out) :: regionLabels(:)

    integer :: nReg, iReg
    integer, allocatable :: tmpI1(:)
    type(hsd_table), pointer :: child, child2
    character(len=:), allocatable :: buffer
    character(lc) :: strTmp
    integer :: stat

    nReg = size(children)
    allocate(regionLabels(nReg))
    allocate(iAtInRegion(nReg))
    do iReg = 1, nReg
      child => children(iReg)%ptr
      call hsd_get(child, "Atoms", buffer, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child, "Atoms must be present")
      call hsd_get_table(child, "Atoms", child2)
      call getSelectedAtomIndices(child2, buffer, geo%speciesNames, geo%species, tmpI1)
      iAtInRegion(iReg)%data = tmpI1
      write(strTmp, "('region',I0)") iReg
      call hsd_get_or_set(child, "Label", buffer, trim(strTmp))
      regionLabels(iReg) = unquote(buffer)
    end do

  end subroutine readPDOSRegions


  !> Read Tunneling and Dos options from analysis block
  subroutine readTunAndDos(root, geo, tundos, transpar, temperature)
    type(hsd_table), pointer :: root
    type(TGeometry), intent(in) :: geo

    !> The container to be filled
    type(TNEGFTunDos), intent(inout) :: tundos
    type(TTransPar), intent(inout) :: transpar
    real(dp), intent(in) :: temperature

    type(hsd_table), pointer :: pNode, pTmp, field
    type(hsd_table_ptr), allocatable :: pNodeList(:)
    integer :: ii, jj, ind, ncont, nKT
    real(dp) :: eRange(2), eRangeDefault(2)
    character(len=:), allocatable :: buffer, modif
    type(TWrappedInt1), allocatable :: iAtInRegion(:)
    logical, allocatable :: tDirectionResInRegion(:)
    character(lc), allocatable :: regionLabelPrefixes(:)
    integer :: stat

    tundos%defined = .true.
    ncont = transpar%ncont
    call hsd_get_or_set(root, "Verbosity", tundos%verbose, 51)
    call hsd_get_or_set(root, "WriteLDOS", tundos%writeLDOS, .true.)
    call hsd_get_or_set(root, "WriteTunn", tundos%writeTunn, .true.)

    ! Default meaningful: eRange= (0..10*kT]
    ! nKT is set to GreensFunction default, i.e. 10
    ! I avoid an explicit nKT option because I find it confusing here
    ! (it makes sense only out of equilibrium)
    ! What matters is w*[nB(w;T1)-nB(w;T2)] that is finite lim w->0
    nKT = 10
    eRangeDefault(1) = 0.0001_dp
    eRangeDefault(2) = nKT * temperature

    block
      real(dp), allocatable :: tmp_eRange(:)
      call hsd_get_or_set(root, "FreqRange", tmp_eRange, eRangeDefault)
      eRange(:) = tmp_eRange
    end block
    call hsd_get_table(root, "FreqRange", field, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modif = field%attrib; else; modif = ""; end if
    call convertUnitHsd(modif, energyUnits, field, eRange)
    tundos%emin = eRange(1)
    tundos%emax = eRange(2)

    if (eRange(1).le.0.d0) then
       call dftbp_error(root, "FreqRange must be > 0")
    end if
    if (eRange(2).lt.eRange(1)) then
       call dftbp_error(root, "Emax < Emin")
    end if

    call hsd_get_or_set(root, "FreqStep", tundos%estep, 1.0e-5_dp)
    call hsd_get_table(root, "FreqStep", field, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modif = field%attrib; else; modif = ""; end if

    call convertUnitHsd(modif, energyUnits, field, tundos%estep)

    ! Terminal currents
    call hsd_get_table(root, "TerminalCurrents", pTmp)

    if (associated(pTmp)) then
      call hsd_get_child_tables(pTmp, "EmitterCollector", pNodeList)
      allocate(tundos%ni(size(pNodeList)))
      allocate(tundos%nf(size(pNodeList)))
      do ii = 1, size(pNodeList)
        pNode => pNodeList(ii)%ptr
        call getEmitterCollectorByName(pNode, tundos%ni(ii),&
            & tundos%nf(ii), transpar%contacts(:)%name)
      end do
    else
      allocate(tundos%ni(ncont-1) )
      allocate(tundos%nf(ncont-1) )
      block
        type(hsd_table) :: tc_table
        call new_table(tc_table, "TerminalCurrents")
        call root%add_child(tc_table)
      end block
      call hsd_get_table(root, "TerminalCurrents", pTmp)
      ind = 1
      do ii = 1, 1
        do jj = ii + 1, ncont
          call hsd_set(pTmp, "EmitterCollector", &
              &[character(mc) :: transpar%contacts(ii)%name, transpar%contacts(jj)%name])
          tundos%ni(ind) = ii
          tundos%nf(ind) = jj
          ind = ind + 1
        end do
      end do
    end if

    call hsd_get_table(root, "DeltaModel", pNode, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(root, "DeltaModel must be present")
    call readDeltaModel(pNode, tundos)

    call hsd_get_or_set(root, "BroadeningDelta", tundos%broadeningDelta, 0.0_dp)
    call hsd_get_table(root, "BroadeningDelta", field, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modif = field%attrib; else; modif = ""; end if
    call convertUnitHsd(modif, energyUnits, field, &
        &tundos%broadeningDelta)

    call hsd_get_child_tables(root, "Region", pNodeList)
    call readPDOSRegions(pNodeList, geo, iAtInRegion, regionLabelPrefixes)

    call addAtomResolvedRegion(tundos%dosOrbitals, tundos%dosLabels)

    call setTypeOfModes(root, transpar)


    contains

      !! Adds one region with all the orbitals of the atoms in it.
      subroutine addAtomResolvedRegion(iOrbRegion, regionLabels)

        type(TWrappedInt1), allocatable, intent(out) :: iOrbRegion(:)
        character(lc), allocatable, intent(out) :: regionLabels(:)

        integer :: nRegion, nAtomInRegion, iReg, ind, ii, jj, iAt
        integer :: nIndices

        nRegion = size(iAtInRegion)
        allocate(iOrbRegion(nRegion))
        allocate(regionLabels(nRegion))

        do iReg = 1, nRegion
          nAtomInRegion = size(iAtInRegion(iReg)%data)
          nIndices = 3*nAtomInRegion
          allocate(iOrbRegion(iReg)%data(nIndices))
          ind = 1
          do ii = 1, nAtomInRegion
            iAt = iAtInRegion(iReg)%data(ii)
            do jj = 0, 2
              iOrbRegion(iReg)%data(ind) = 3*iAt - 2 + jj
              ind = ind + 1
            end do
          end do
          regionLabels(iReg) = regionLabelPrefixes(iReg)
        end do

      end subroutine addAtomResolvedRegion


  end subroutine readTunAndDos


  !> Get contacts for terminal currents by name
  subroutine getEmitterCollectorByName(pNode, emitter, collector, contactNames)
    type(hsd_table), pointer :: pNode
    integer, intent(out) :: emitter, collector
    character(len=*), intent(in) :: contactNames(:)

    character(:), allocatable :: strArr(:)
    character(len=mc) :: buffer
    integer :: ind, stat
    logical :: tFound

    call hsd_get(pNode, "#text", strArr, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "Missing contact data")
    if (size(strArr) /= 2) then
      call dftbp_error(pNode, "You must provide two contacts")
    end if
    buffer = strArr(1)
    emitter = getContactByName(contactNames, buffer, pNode)
    buffer = strArr(2)
    collector = getContactByName(contactNames, buffer, pNode)

  end subroutine getEmitterCollectorByName


  !> Getting the contact by name
  function getContactByName(contactNames, contName, pNode) result(contact)
    character(len=*), intent(in) :: contactNames(:)
    character(len=*), intent(in) :: contName
    type(hsd_table), pointer :: pNode
    integer :: contact

    logical :: tFound

    tFound = .false.
    do contact = 1, size(contactNames)
      tFound = (contactNames(contact) == contName)
      if (tFound) then
        exit
      end if
    end do
    if (.not. tFound) then
      call dftbp_error(pNode, "Invalid collector contact name '" &
          &// trim(contName) // "'")
    end if

  end function getContactByName


  !> Set model for w-dependent delta in G.F.
  subroutine readDeltaModel(root, tundos)
    type(hsd_table), pointer :: root
    type(TNEGFTunDos), intent(inout) :: tundos

    type(hsd_table), pointer :: pValue, pChild, field
    character(len=:), allocatable :: buffer, modif
    integer :: stat

    ! Get first table child for dispatch
    pValue => null()
    block
      class(hsd_node), pointer :: ch
      integer :: jj
      do jj = 1, root%num_children
        call root%get_child(jj, ch)
        select type (t => ch)
        type is (hsd_table)
          pValue => t
          exit
        end select
      end do
    end block
    pChild => root
    if (associated(pValue)) then
      call hsd_get_name(pValue, buffer, "#text")
    else
      buffer = "#text"
    end if
    ! Delta is repeated to allow different defaults if needed
    select case (trim(buffer))
    case("deltasquared")
      call hsd_get_or_set(pValue, "Delta", tundos%delta, 0.0001_dp)
      call hsd_get_table(pValue, "Delta", field, auto_wrap=.true.)
      if (allocated(field%attrib)) then; modif = field%attrib; else; modif = ""; end if
      call convertUnitHsd(modif, energyUnits, field, tundos%delta)
      tundos%deltaModel=0
    case("deltaomega")
      call hsd_get_or_set(pValue, "Delta", tundos%delta, 0.0001_dp)
      call hsd_get_table(pValue, "Delta", field, auto_wrap=.true.)
      if (allocated(field%attrib)) then; modif = field%attrib; else; modif = ""; end if
      call convertUnitHsd(modif, energyUnits, field, tundos%delta)
      tundos%deltaModel=1
    case("mingo")
      ! As in Numerical Heat transfer, Part B, 51:333, 2007, Taylor&Francis.
      ! Here Delta is just a dimensionless scaling factor
      call hsd_get_or_set(pValue, "Delta", tundos%delta, 0.0001_dp)
      ! We set a cutoff frequency of 2000 cm^-1.
      call hsd_get_or_set(pValue, "Wmax", tundos%wmax, 0.009_dp)
      call hsd_get_table(pValue, "Wmax", field, auto_wrap=.true.)
      if (allocated(field%attrib)) then; modif = field%attrib; else; modif = ""; end if
      call convertUnitHsd(modif, energyUnits, field, tundos%delta)
      tundos%deltaModel=2
      ! If Emax >> Wmax delta becomes negative
      if (tundos%Emax > tundos%Wmax) then
        call dftbp_error(pValue,"In Mingo model check Wmax <= Emax")
      end if
    case default
      call dftbp_error(pValue,"Unknown deltaModel "//trim(buffer))
    end select
  end subroutine ReadDeltaModel


  !> Build a simple neighbor list. Currently does not work for periodic systems.
  !! Have to fix this important point
  subroutine buildNeighbourList()

    integer ::  iAtom, jAtom, ii, jj, kk, PL1, PL2
    !* First guess for nr. of neighbors.
    integer, parameter :: nInitNeighbours = 100
    real :: disAtom, dd(3)
    integer :: nAllAtom
    real(dp) :: mCutoff
    real(dp), allocatable :: coords(:,:), cellVec(:,:), rCellVec(:,:)
    integer, allocatable :: iCellVec(:)
    type(TStatus) :: errStatus

    call TNeighbourlist_init(neighbourList, geo%nAtom, nInitNeighbours)

    mCutoff = 1.0_dp * cutoff

    if (geo%tPeriodic) then
      !! Make some guess for the nr. of all interacting atoms
      nAllAtom = int((real(geo%nAtom, dp)**(1.0_dp/3.0_dp) + 3.0_dp)**3)
      call getCellTranslations(cellVec, rCellVec, geo%latVecs, &
            &geo%recVecs2p, mCutoff)
    else
      nAllAtom = geo%nAtom
      allocate(rCellVec(3, 1))
      rCellVec(:, 1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
    end if

    allocate(coords(3, nAllAtom))
    allocate(img2CentCell(nAllAtom))
    allocate(iCellVec(nAllAtom))

    call updateNeighbourList(coords, img2CentCell, iCellVec, neighbourList, &
        &nAllAtom, geo%coords, mCutoff, rCellVec, errStatus, symmetric=.false.)
    if (errStatus%hasError()) then
      call error(errStatus%message)
    end if

    deallocate(coords)
    deallocate(iCellVec)
    deallocate(rCellVec)

    allocate(nNeighbour(geo%nAtom))
    nNeighbour(:) = 0

    call getNrOfNeighboursForAll(nNeighbour, neighbourList, mCutoff)

    ! Check PL size with neighbor list
    do iAtom = 1, transpar%idxdevice(2)
      PL1 = getPL(iAtom)
      do jj = 1, nNeighbour(iAtom)
        jAtom = img2CentCell(neighbourList%iNeighbour(jj,iAtom))
        if (jAtom > transpar%idxdevice(2)) cycle
        PL2 = getPL(jAtom)
        if (.not.(PL1.eq.PL2 .or. PL1.eq.PL2+1 .or. PL1.eq.PL2-1)) then
          write(stdOut,*) 'ERROR: PL size inconsistent with cutoff'
          stop
        end if
      end do
    end do

  end subroutine buildNeighbourList


  subroutine cutDynMatrix()

    integer :: iAtom, jAtom, jj
    real(dp), allocatable :: dynMat2(:,:)

    allocate(dynMat2(3*nMovedAtom, 3*nMovedAtom))
    dynMat2 = 0.0_dp

    do iAtom = 1, geo%nAtom
       do jj = 1, nNeighbour(iAtom)
          jAtom = img2CentCell(neighbourList%iNeighbour(jj, iAtom))
          if (neighbourList%neighDist2(jj,iAtom) .le. cutoff**2) then
            dynMat2(3*(iAtom-1)+1:3*(iAtom-1)+3, 3*(jAtom-1)+1:3*(jAtom-1)+3) = &
                dynMatrix(3*(iAtom-1)+1:3*(iAtom-1)+3, 3*(jAtom-1)+1:3*(jAtom-1)+3)
            dynMat2(3*(jAtom-1)+1:3*(jAtom-1)+3, 3*(iAtom-1)+1:3*(iAtom-1)+3) = &
                dynMatrix(3*(iAtom-1)+1:3*(iAtom-1)+3, 3*(jAtom-1)+1:3*(jAtom-1)+3)
          end if
       end do
    end do

    dynMatrix = dynMat2

    deallocate(dynMat2)

  end subroutine cutDynMatrix


  function getPL(iAt) result(PL)
    integer, intent(in) :: iAt
    integer :: PL

    integer :: ii

    do ii = 1, transpar%nPLs-1
      if (iAt>=transpar%PL(ii) .and. iAt<transpar%PL(ii+1)) then
        PL = ii
      end if
    end do
    if ( iAt>=transpar%PL(transpar%nPLs) .and. iAt<=transpar%idxdevice(2)) then
      PL = transpar%nPLs
    endif
    if (iAt > transpar%idxdevice(2)) then
      PL = transpar%nPLs + 1
    endif

  end function getPL


  !> Select family of modes to analyze and restrict transmission
  subroutine setTypeOfModes(root, transpar)
    type(hsd_table), pointer :: root
    type(TTransPar), intent(inout) :: transpar

    character(len=:), allocatable :: buffer

    !selecting the type of modes you want to analyze
    call hsd_get_or_set(root, "ModeType", buffer, "all")
    select case(trim(buffer))
    case("all")
      selTypeModes = modeEnum%ALLMODES
    case("along-x")
      selTypeModes = modeEnum%XX
    case("along-y")
      selTypeModes = modeEnum%YY
    case("along-z")
      selTypeModes = modeEnum%ZZ
    case("longitudinal")
      selTypeModes = modeEnum%LONGITUDINAL
      call checkTypeOfModes(root, transpar)
    case("transverse")
      selTypeModes = modeEnum%TRANSVERSE
      call checkTypeOfModes(root, transpar)
    case("in-plane")
      selTypeModes = modeEnum%INPLANE
      call checkTypeOfModes(root, transpar)
    case("out-of-plane")
      selTypeModes = modeEnum%OUTOFPLANE
      call checkTypeOfModes(root, transpar)
    case default
      call dftbp_error(root,"Unknown type of modes")
    end select

    transpar%typeModes = selTypeModes

  end subroutine setTypeOfModes


  !> Check that the geometry orientation is consistent with selTypeModes
  !! Currently only checks that transport direction is along z
  subroutine checkTypeOfModes(root, tp)
    type(hsd_table), pointer :: root
    type(TTransPar), intent(inout) :: tp

    select case (selTypeModes)
    case(modeEnum%LONGITUDINAL, modeEnum%TRANSVERSE)
      call checkAlongZ(root, tp)
    case(modeEnum%INPLANE, modeEnum%OUTOFPLANE)
      call checkAlongZ(root, tp)
    case default
    end select

  end subroutine checkTypeOfModes


  !> Check that transport direction is along z
  subroutine checkAlongZ(root, tp)
    type(hsd_table), pointer :: root
    type(TTransPar), intent(inout) :: tp

    real(dp) :: contactVec(3)
    logical :: mask(3)
    integer :: ii

    do ii = 1, size(tp%contacts)
      contactVec = tp%contacts(ii)%lattice
      ! Determine to which axis the contact vector is parallel.
      mask = (abs(abs(contactVec) - sqrt(sum(contactVec**2))) < 1.0e-8_dp)
      if (count(mask) /= 1 .or. .not.mask(3)) then
        call dftbp_error(root,"Transport direction is not along z")
      end if
    end do

  end subroutine checkAlongZ


  subroutine appendCharLc(arr, elem)
    character(lc), allocatable, intent(inout) :: arr(:)
    character(lc), intent(in) :: elem
    character(lc), allocatable :: tmp(:)
    integer :: nn

    nn = size(arr)
    allocate(tmp(nn + 1))
    if (nn > 0) tmp(:nn) = arr(:nn)
    tmp(nn + 1) = elem
    call move_alloc(tmp, arr)
  end subroutine appendCharLc


end module phonons_initphonons
