!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'


!> Fills the derived type with the input parameters from an HSD or an XML file.
module dftbp_dftbplus_parser
  use dftbp_common_globalenv, only : stdout
  use dftbp_common_status, only : TStatus
  use dftbp_dftbplus_inputdata, only : TInputData
  use dftbp_dftbplus_oldcompat, only : parserVersion
  use dftbp_io_hsdutils, only : dftbp_error
  use hsd, only : hsd_error_t, hsd_get_table, hsd_get_choice
  use hsd_data, only : hsd_table, data_load, DATA_FMT_AUTO, new_table
  use dftbp_io_message, only : error
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_dftbplus_parser_analysis, only : readAnalysis, readLaterAnalysis
  use dftbp_dftbplus_parser_driver, only : readDriver, readElecDynamics
  use dftbp_dftbplus_parser_excited, only : readExcited
  use dftbp_dftbplus_parser_general, only : readOptions
  use dftbp_dftbplus_parser_geometry, only : readGeometry
  use dftbp_dftbplus_parser_hamiltonian, only : readHamiltonian, readLaterHamiltonian
  use dftbp_dftbplus_parser_parallel, only : readParallel
  use dftbp_dftbplus_parser_parseroptions, only : TParserFlags, handleInputVersion,&
      & readParserOptions
  use dftbp_dftbplus_parser_reks, only : readReks
  use dftbp_dftbplus_parser_spin, only : readSpinConstants
#:if WITH_TRANSPORT
  use dftbp_dftbplus_parser_transport, only : readTransportGeometry, finalizeNegf
#:endif
  implicit none

  private
  public :: parserVersion, rootTag
  public :: TParserFlags
  public :: readHsdFile, parseHsdTree


  !> Tag at the head of the input document tree
  character(len=*), parameter :: rootTag = "dftbplusinput"


contains

  !> Reads input from a file (HSD, XML, JSON, or TOML format).
  !>
  !> The format is auto-detected from the file extension:
  !>   .hsd → HSD format
  !>   .xml → XML format
  !>   .json → JSON format
  !>   .toml → TOML format
  !>   other → defaults to HSD format
  !>
  !> The loaded content is wrapped in a rootTag child to match the tree
  !> structure expected by parseHsdTree (document → rootTag → children).
  subroutine readHsdFile(hsdFile, hsdTree)

    !> Name of the input file
    character(*), intent(in) :: hsdFile

    !> Data tree representation of the input (document → rootTag → content)
    type(hsd_table), pointer :: hsdTree

    type(hsd_table), pointer :: content, existingRoot
    type(hsd_error_t), allocatable :: hsdError
    integer :: stat

    ! Load raw content into a temporary table
    allocate(content)
    call data_load(hsdFile, content, hsdError, fmt=DATA_FMT_AUTO)
    if (allocated(hsdError)) then
      call error("Error loading input file '" // trim(hsdFile) // "': " // trim(hsdError%message))
    end if

    ! Check if content already has a rootTag child (e.g. from processed dftb_pin.hsd).
    ! If so, we must not wrap again — just use the content as the document node.
    existingRoot => null()
    call hsd_get_table(content, rootTag, existingRoot, stat)
    if (associated(existingRoot)) then
      ! Already wrapped (dftb_pin.hsd starts with "dftbplusinput { ... }")
      content%name = "document"
      hsdTree => content
    else
      ! Bare content — wrap in document structure: hsdTree → rootTag → {children}
      ! This matches the tree layout the legacy parseHSD() produced.
      content%name = rootTag
      allocate(hsdTree)
      call new_table(hsdTree, name="document")
      call hsdTree%add_child(content)
    end if

  end subroutine readHsdFile


  !> Parse input from an HSD/XML file
  subroutine parseHsdTree(hsdTree, input, parserFlags)

    !> Tree representation of the input
    type(hsd_table), pointer :: hsdTree

    !> Returns initialised input variables on exit
    type(TInputData), intent(out) :: input

    !> Special block containings parser related settings
    type(TParserFlags), intent(out) :: parserFlags

    type(TStatus) :: errStatus
    type(TOrbitals) :: orb
    type(hsd_table), pointer :: root, tmp, driverNode, hamNode, analysisNode, child, dummy
    logical :: tReadAnalysis
    integer, allocatable :: implicitParserVersion
    integer :: stat
    character(len=:), allocatable :: choiceName

    write(stdout, '(A,1X,I0,/)') 'Parser version:', parserVersion
    write(stdout, "(A)") repeat("-", 80)

    call hsd_get_table(hsdTree, rootTag, root, stat, auto_wrap=.true.)
    if (.not. associated(root)) call error("Missing root node: '" // rootTag // "'")

    call handleInputVersion(root, implicitParserVersion)
    call hsd_get_table(root, "ParserOptions", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="parseroptions")
        call root%add_child(tmpTbl)
      end block
      call hsd_get_table(root, "ParserOptions", child, stat)
    end if
    call readParserOptions(child, root, parserFlags, implicitParserVersion)

    ! Read the geometry unless the list of atoms has been provided through the API
    if (.not. allocated(input%geom%coords)) then
      call hsd_get_table(root, "Geometry", tmp, stat, auto_wrap=.true.)
      if (.not. associated(tmp)) call dftbp_error(root, "Missing required block: 'Geometry'")
      call readGeometry(tmp, input)
    end if
    input%geom%areContactsPresent = .false.

    ! Hamiltonian settings that need to know settings from the REKS block
    call hsd_get_table(root, "Reks", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: defContainer, defChild
        call new_table(defContainer, name="reks")
        call new_table(defChild, name="None")
        call defContainer%add_child(defChild)
        call root%add_child(defContainer)
      end block
      call hsd_get_table(root, "Reks", child, stat)
    end if
    call hsd_get_choice(child, "", choiceName, dummy, stat)
    call readReks(child, dummy, input%ctrl, input%geom)

    call hsd_get_table(root, "Transport", child, stat, auto_wrap=.true.)

  #:if WITH_TRANSPORT

    ! Read in transport and modify geometry if it is only a contact calculation
    if (associated(child)) then
      call readTransportGeometry(child, input%geom, input%transpar)
    else
      input%transpar%ncont=0
      allocate(input%transpar%contacts(0))
      ! set range of atoms in the 'device', as there are no contacts
      input%transpar%idxdevice(1) = 1
      input%transpar%idxdevice(2) = input%geom%nAtom
    end if

    call hsd_get_table(root, "Dephasing", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      call dftbp_error(child, "Be patient... Dephasing feature will be available soon!")
      !call readDephasing(child, input%slako%orb, input%geom, input%transpar, input%ginfo%tundos)
    end if

    ! electronic Hamiltonian
    call hsd_get_choice(root, "Hamiltonian", choiceName, hamNode, stat)
    if (.not. associated(hamNode)) call dftbp_error(root, "Missing required block: 'Hamiltonian'")
    call readHamiltonian(hamNode, input%ctrl, input%geom, input%slako, input%transpar,&
        & input%ginfo%greendens, input%poisson, errStatus)

  #:else

    if (associated(child)) then
      call dftbp_error(child, "Program has been compiled without transport enabled")
    end if

    ! electronic Hamiltonian
    call hsd_get_choice(root, "Hamiltonian", choiceName, hamNode, stat)
    if (.not. associated(hamNode)) call dftbp_error(root, "Missing required block: 'Hamiltonian'")
    call readHamiltonian(hamNode, input%ctrl, input%geom, input%slako, input%poisson, errStatus)

  #:endif

    if (errStatus%hasError()) then
      call error(errStatus%message)
    end if

    call hsd_get_table(root, "Driver", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="driver")
        call root%add_child(tmpTbl)
      end block
      call hsd_get_table(root, "Driver", child, stat)
    end if
    call hsd_get_choice(child, "", choiceName, driverNode, stat)
  #:if WITH_TRANSPORT
    call readDriver(driverNode, child, input%geom, input%ctrl, input%transpar)
  #:else
    call readDriver(driverNode, child, input%geom, input%ctrl)
  #:endif

    tReadAnalysis = .true.
    call hsd_get_table(root, "ElectronDynamics", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      allocate(input%ctrl%elecDynInp)
      call readElecDynamics(child, input%ctrl%elecDynInp, input%geom, input%ctrl%masses)
      if (input%ctrl%elecDynInp%tReadRestart .and. .not.input%ctrl%elecDynInp%tPopulations) then
        tReadAnalysis = .false.
      end if

    end if

    if (tReadAnalysis) then
      ! Analysis of properties
      call hsd_get_table(root, "Analysis", analysisNode, stat, auto_wrap=.true.)
      if (.not. associated(analysisNode)) then
        block
          type(hsd_table) :: tmpTbl
          call new_table(tmpTbl, name="analysis")
          call root%add_child(tmpTbl)
        end block
        call hsd_get_table(root, "Analysis", analysisNode, stat)
      end if

    #:if WITH_TRANSPORT
      call readAnalysis(analysisNode, input%ctrl, input%geom, input%slako%orb, input%transpar, &
          & input%ginfo%tundos)

      call finalizeNegf(input)
    #:else
      call readAnalysis(analysisNode, input%ctrl, input%geom)
    #:endif

    end if

    call hsd_get_table(root, "ExcitedState", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="excitedstate")
        call root%add_child(tmpTbl)
      end block
      call hsd_get_table(root, "ExcitedState", child, stat)
    end if
    call readExcited(child, input%geom, input%ctrl)

    ! Hamiltonian settings that need to know about settings from the blocks above
    call readLaterHamiltonian(hamNode, input%ctrl, driverNode, input%geom)

    call hsd_get_table(root, "Options", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="options")
        call root%add_child(tmpTbl)
      end block
      call hsd_get_table(root, "Options", child, stat)
    end if
    call readOptions(child, input%ctrl, input%geom)

    ! W values if needed by Hamiltonian or excited state calculation
    if (allocated(input%ctrl%tbliteInp)) then
      call input%ctrl%tbliteInp%setupOrbitals(input%geom%species, orb)
      call readSpinConstants(hamNode, input%geom, orb, input%ctrl)
    else
      call readSpinConstants(hamNode, input%geom, input%slako%orb, input%ctrl)
    end if

    ! analysis settings that need to know settings from the options block
    if (tReadAnalysis) then
      call readLaterAnalysis(analysisNode, input%ctrl)
    end if

    ! read parallel calculation settings
    call readParallel(root, input)

    ! input data strucutre has been initialised
    input%tInitialized = .true.

  end subroutine parseHsdTree


end module dftbp_dftbplus_parser
