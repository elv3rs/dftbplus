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
  use dftbp_common_accuracy, only : distFudge, distFudgeOld, dp, lc, mc, minTemp, sc
  use dftbp_common_constants, only : Bohr__AA, boltzmann, maxL, pi, shellNames, symbolToNumber
  use dftbp_common_file, only : closeFile, openFile, TFileDescr
  use dftbp_common_filesystem, only : findFile, getParamSearchPaths, joinPathsPrettyErr
  use dftbp_common_globalenv, only : abortProgram, stdout, withMpi, withScalapack
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_status, only : TStatus
  use dftbp_common_unitconversion, only : EFieldUnits, energyUnits, freqUnits, lengthUnits
  use dftbp_dftb_coordnumber, only : cnType, getD3Radius, getElectronegativity, TCNInput
  use dftbp_dftb_dftbplusu, only : plusUFunctionals
  use dftbp_dftb_dftd4param, only : getEeqChi, getEeqGam, getEeqKcn, getEeqRad
  use dftbp_dftb_dispersions, only : getUffValues, TDispDftD4Inp, TDispersionInp, TDispSlaKirkInp,&
      & TDispUffInp, TSimpleDftD3Input
  use dftbp_dftb_elecconstraints, only : readElecConstraintInput
  use dftbp_dftb_encharges, only : TEeqInput
  use dftbp_dftb_etemp, only : fillingTypes
  use dftbp_dftb_halogenx, only : halogenXSpecies1, halogenXSpecies2
  use dftbp_dftb_hybridxc, only : checkSupercellFoldingMatrix, hybridXcAlgo, hybridXcFunc,&
      & hybridXcGammaTypes, THybridXcSKTag
  use dftbp_dftb_nonscc, only : diffTypes
  use dftbp_dftb_periodic, only : getCellTranslations, getSuperSampling, TNeighbourList,&
      & TNeighbourlist_init, updateNeighbourList
  use dftbp_dftb_repulsive_chimesrep, only : TChimesRepInp
  use dftbp_dftb_repulsive_polyrep, only : TPolyRep, TPolyRepInp
  use dftbp_dftb_repulsive_splinerep, only : TSplineRep, TSplineRepInp
  use dftbp_dftb_slakocont, only : addTable, init
  use dftbp_dftb_slakoeqgrid, only : init, skEqGridNew, skEqGridOld, TSlakoEqGrid
  use dftbp_dftb_mdftb, only : TMdftbAtomicIntegrals
  use dftbp_dftbplus_forcetypes, only : forceTypes
  use dftbp_dftbplus_input_fileaccess, only : readBinaryAccessTypes
  use dftbp_dftbplus_inputdata, only : TBlacsOpts, TControl, THybridXcInp, TInputData, TSlater
  use dftbp_dftbplus_oldcompat, only : convertOldHSD, minVersion, parserVersion, versionMaps
  use dftbp_dftbplus_specieslist, only : readSpeciesList
  use dftbp_elecsolvers_elecsolvers, only : electronicSolverTypes, providesEigenvalues
  use dftbp_extlibs_arpack, only : withArpack
  use dftbp_extlibs_elsiiface, only : withELSI, withPEXSI
  use dftbp_extlibs_poisson, only : TPoissonInfo, withPoisson
  use dftbp_extlibs_sdftd3, only : dampingFunction, TSDFTD3Input
  use dftbp_extlibs_tblite, only : tbliteMethod
  use dftbp_io_charmanip, only : i2c, newline, tolower, unquote
  use dftbp_io_hsdcompat, only : hsd_table, hsd_child_list, hsd_error_t, textNodeName, &
      & detailedError, detailedWarning, getChild, getChildren, getChildValue, &
      & getSelectedAtomIndices, setChild, setChildValue, getNodeName, getNodeHSDName, &
      & getNodeName2, getLength, getItem1, destroyNodeList, removeChild, &
      & convertUnitHsd, hsd_rename_child, setUnprocessed, new_table, hasInlineData, &
      & hsd_get_table
  use dftbp_extlibs_hsddata, only : data_load, data_detect_format, DATA_FMT_AUTO
  use dftbp_io_message, only : error, warning
  use dftbp_math_simplealgebra, only : cross3, determinant33, diagonal
  use dftbp_md_thermostats, only : thermostatTypes
  use dftbp_reks_reks, only : reksTypes
  use dftbp_solvation_solvparser, only : readCM5, readSolvation
  use dftbp_timedep_linresptypes, only : linRespSolverTypes
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_linkedlist, only : append, asArray, asVector, destruct, get, init, intoArray, len,&
      & TListCharLc, TListInt, TListIntR1, TListReal, TListRealR1, TListRealR2, TListString
  use dftbp_type_oldskdata, only : parseHybridXcTag, readFromFile, TOldSKData
  use dftbp_type_orbitals, only : getShellnames
  use dftbp_type_typegeometry, only : reduce, setLattice, TGeometry
  use dftbp_type_typegeometryhsd, only : readTGeometryGen, readTGeometryHsd, readTGeometryLammps,&
      & readTGeometryVasp, readTGeometryXyz
  use dftbp_type_wrappedintr, only : TWrappedInt1
#:if WITH_MBD
  use dftbp_dftb_dispmbd, only : TDispMbdInp
#:endif

#:if WITH_TRANSPORT
  use dftbp_transport_negfvars, only : ContactInfo, TElPh, TNEGFGreenDensInfo, TNEGFTunDos
#:endif
  use dftbp_transport_negfvars, only : TTransPar
  use dftbp_dftbplus_parser_analysis, only : readAnalysis, readLaterAnalysis
  use dftbp_dftbplus_parser_dispersion, only : readDispersion
  use dftbp_dftbplus_parser_driver, only : readDriver, readElecDynamics
  use dftbp_dftbplus_parser_excited, only : readExcited
  use dftbp_dftbplus_parser_hybrid, only : parseHybridBlock, parseChimes
  use dftbp_dftbplus_parser_kpoints, only : readKPoints, maxSelfConsIterations
  use dftbp_dftbplus_parser_parallel, only : readParallel
  use dftbp_dftbplus_parser_reks, only : readReks
  use dftbp_dftbplus_parser_skfiles, only : readSKFiles
#:if WITH_TRANSPORT
  use dftbp_dftbplus_parser_transport, only : readTransportGeometry, readGreensFunction, &
      & readTunAndDos, finalizeNegf
#:endif
#:if WITH_POISSON
  use dftbp_dftbplus_parser_transport, only : readPoisson
#:endif
  implicit none

  private
  public :: parserVersion, rootTag
  public :: TParserFlags
  public :: readHsdFile, parseHsdTree


  !> Tag at the head of the input document tree
  character(len=*), parameter :: rootTag = "dftbplusinput"


  !> Container type for parser related flags.
  type TParserFlags

    !> stop after parsing?
    logical :: tStop

    !> Continue despite unprocessed nodes
    logical :: tIgnoreUnprocessed

    !> HSD output?
    logical :: tWriteHSD
  end type TParserFlags


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

    write(stdout, '(A,1X,I0,/)') 'Parser version:', parserVersion
    write(stdout, "(A)") repeat("-", 80)

    call getChild(hsdTree, rootTag, root)

    call handleInputVersion(root, implicitParserVersion)
    call getChildValue(root, "ParserOptions", dummy, "", child=child, list=.true.,&
        & allowEmptyValue=.true., dummyValue=.true.)
    call readParserOptions(child, root, parserFlags, implicitParserVersion)

    ! Read the geometry unless the list of atoms has been provided through the API
    if (.not. allocated(input%geom%coords)) then
      call getChild(root, "Geometry", tmp)
      call readGeometry(tmp, input)
    end if
    input%geom%areContactsPresent = .false.

    ! Hamiltonian settings that need to know settings from the REKS block
    call getChildValue(root, "Reks", dummy, "None", child=child)
    call readReks(child, dummy, input%ctrl, input%geom)

    call getChild(root, "Transport", child, requested=.false.)

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

    call getChild(root, "Dephasing", child, requested=.false.)
    if (associated(child)) then
      call detailedError(child, "Be patient... Dephasing feature will be available soon!")
      !call readDephasing(child, input%slako%orb, input%geom, input%transpar, input%ginfo%tundos)
    end if

    ! electronic Hamiltonian
    call getChildValue(root, "Hamiltonian", hamNode)
    call readHamiltonian(hamNode, input%ctrl, input%geom, input%slako, input%transpar,&
        & input%ginfo%greendens, input%poisson, errStatus)

  #:else

    if (associated(child)) then
      call detailedError(child, "Program has been compiled without transport enabled")
    end if

    ! electronic Hamiltonian
    call getChildValue(root, "Hamiltonian", hamNode)
    call readHamiltonian(hamNode, input%ctrl, input%geom, input%slako, input%poisson, errStatus)

  #:endif

    if (errStatus%hasError()) then
      call error(errStatus%message)
    end if

    call getChildValue(root, "Driver", driverNode, "", child=child, allowEmptyValue=.true.)
  #:if WITH_TRANSPORT
    call readDriver(driverNode, child, input%geom, input%ctrl, input%transpar)
  #:else
    call readDriver(driverNode, child, input%geom, input%ctrl)
  #:endif

    tReadAnalysis = .true.
    call getChild(root, "ElectronDynamics", child=child, requested=.false.)
    if (associated(child)) then
      allocate(input%ctrl%elecDynInp)
      call readElecDynamics(child, input%ctrl%elecDynInp, input%geom, input%ctrl%masses)
      if (input%ctrl%elecDynInp%tReadRestart .and. .not.input%ctrl%elecDynInp%tPopulations) then
        tReadAnalysis = .false.
      end if

    end if

    if (tReadAnalysis) then
      ! Analysis of properties
      call getChildValue(root, "Analysis", dummy, "", child=analysisNode, list=.true.,&
          & allowEmptyValue=.true., dummyValue=.true.)

    #:if WITH_TRANSPORT
      call readAnalysis(analysisNode, input%ctrl, input%geom, input%slako%orb, input%transpar, &
          & input%ginfo%tundos)

      call finalizeNegf(input)
    #:else
      call readAnalysis(analysisNode, input%ctrl, input%geom)
    #:endif

    end if

    call getChildValue(root, "ExcitedState", dummy, "", child=child, list=.true.,&
        & allowEmptyValue=.true., dummyValue=.true.)
    call readExcited(child, input%geom, input%ctrl)

    ! Hamiltonian settings that need to know about settings from the blocks above
    call readLaterHamiltonian(hamNode, input%ctrl, driverNode, input%geom)

    call getChildValue(root, "Options", dummy, "", child=child, list=.true.,&
        & allowEmptyValue=.true., dummyValue=.true.)
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


  !> Converts input version to parser version and removes InputVersion node if present.
  subroutine handleInputVersion(root, implicitParserVersion)

    !> Root eventually containing InputVersion
    type(hsd_table), pointer, intent(in) :: root

    !> Parser version corresponding to input version, or unallocated if none has been found
    integer, allocatable, intent(out) :: implicitParserVersion

    type(hsd_table), pointer :: child, dummy
    character(len=:), allocatable :: versionString

    call getChild(root, "InputVersion", child, requested=.false.)
    if (associated(child)) then
      call getChildValue(child, "", versionString)
      implicitParserVersion = parserVersionFromInputVersion(trim(unquote(versionString)))
      if (implicitParserVersion == 0) then
        call detailedError(child, "Input version '" // trim(unquote(versionString))&
            & // "' not recognized")
      end if
    end if

  end subroutine handleInputVersion


  !> Read in parser options (options not passed to the main code)
  subroutine readParserOptions(node, root, flags, implicitVersion)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Root of the entire tree (in case it needs to be converted, for example because of
    !> compatibility options)
    type(hsd_table), pointer :: root

    !> Contains parser flags on exit.
    type(TParserFlags), intent(out) :: flags

    !> Parser version implied by version number
    integer, intent(in), optional :: implicitVersion

    integer :: inputVersion
    type(hsd_table), pointer :: child

    call getChild(node, "ParserVersion", child, requested=.false.)
    if (.not. associated(child) .and. .not. present(implicitVersion)) then
      call detailedWarning(root, "Input containing neither InputVersion nor ParserVersion is&
          & DEPRECATED!(!!) Specify the InputVersion keyword in your input to ensure that future&
          & versions of DFTB+ can also parse it.")
      inputVersion = parserVersion
      call setChildValue(node, "ParserVersion", inputVersion)
    else if (.not. associated(child) .and. present(implicitVersion)) then
      inputVersion = implicitVersion
    else if (associated(child) .and. .not. present(implicitVersion)) then
      call getChildValue(child, "", inputVersion)
    else
      call getChildValue(child, "", inputVersion)
      if (inputVersion /= implicitVersion) then
        call detailedError(child, "Parser version deduced from InputVersion ("&
            & // i2c(implicitVersion) // ") differs from version explicitely set in&
            & ParserVersion (" // i2c(inputVersion) // ")")
      end if
    end if

    if (inputVersion < 1 .or. inputVersion > parserVersion) then
      call detailedError(child, "Invalid parser version (" // i2c(inputVersion) // ")")
    else if (inputVersion < minVersion) then
      call detailedError(child, &
          & "Sorry, no compatibility mode for parser version " // i2c(inputVersion)&
          & // " (too old)")
    else if (inputVersion /= parserVersion) then
      write(stdout, "(A,I2,A,I2,A)") "***  Converting input from parser version ",&
          & inputVersion, " to parser version ", parserVersion, " ..."
      call convertOldHSD(root, inputVersion, parserVersion)
      write(stdout, "(A,/)") "***  Done."
    end if

    call getChildValue(node, "WriteHSDInput", flags%tWriteHSD, .true.)
    if (.not. flags%tWriteHSD) then
      call detailedWarning(node, "WriteHSDInput turned off. You are not guaranteed" // newline // &
          &" to able to obtain the same results with a later version of the code!" // newline // &
          & "(the dftb_pin.hsd file DOES guarantee this)")
    end if
    call getChildValue(node, "StopAfterParsing", flags%tStop, .false.)

    call getChildValue(node, "IgnoreUnprocessedNodes", &
        &flags%tIgnoreUnprocessed, .false.)

  end subroutine readParserOptions


  !> Read in Geometry
  subroutine readGeometry(node, input)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Input structure to be filled
    type(TInputData), intent(inout) :: input

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer

    call getChildValue(node, "", value1, child=child)
    call getNodeName(value1, buffer)
    input%geom%tPeriodic = .false.
    input%geom%tHelical = .false.
    select case (buffer)
    case ("genformat")
      call readTGeometryGen(value1, input%geom)
    case ("xyzformat")
      call readTGeometryXyz(value1, input%geom)
    case ("vaspformat")
      call readTGeometryVasp(value1, input%geom)
    case ("lammpsformat")
      call readTGeometryLammps(value1, input%geom)
    case default
      call setUnprocessed(value1)
      call readTGeometryHSD(child, input%geom)
    end select

  end subroutine readGeometry


  !> Reads Hamiltonian
#:if WITH_TRANSPORT
  subroutine readHamiltonian(node, ctrl, geo, slako, tp, greendens, poisson, errStatus)
#:else
  subroutine readHamiltonian(node, ctrl, geo, slako, poisson, errStatus)
#:endif

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    !> Slater-Koster structure to be filled
    type(TSlater), intent(inout) :: slako

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens
  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    character(len=:), allocatable :: buffer
    type(hsd_table), pointer :: child

    call getNodeName(node, buffer)
    select case (buffer)
    case ("dftb")
  #:if WITH_TRANSPORT
      call readDFTBHam(node, ctrl, geo, slako, tp, greendens, poisson, errStatus)
  #:else
      call readDFTBHam(node, ctrl, geo, slako, poisson, errStatus)
  #:endif
      @:PROPAGATE_ERROR(errStatus)
    case ("xtb")
  #:if WITH_TRANSPORT
      call readXTBHam(node, ctrl, geo, tp, greendens, poisson, errStatus)
  #:else
      call readXTBHam(node, ctrl, geo, poisson, errStatus)
  #:endif
      @:PROPAGATE_ERROR(errStatus)
    case default
      call detailedError(node, "Invalid Hamiltonian")
    end select

  #:if WITH_API
    call getChild(node, "ASI", child, requested=.false.)
    if (associated(child)) then
      ctrl%isASICallbackEnabled = .true.
    else
      ctrl%isASICallbackEnabled = .false.
    end if
  #:endif

  end subroutine readHamiltonian


  !> Reads DFTB-Hamiltonian
#:if WITH_TRANSPORT
  subroutine readDFTBHam(node, ctrl, geo, slako, tp, greendens, poisson, errStatus)
#:else
  subroutine readDFTBHam(node, ctrl, geo, slako, poisson, errStatus)
#:endif

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Slater-Koster structure to be filled
    type(TSlater), intent(inout) :: slako

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(hsd_table), pointer :: value1, child, child2, child3
    type(hsd_child_list), pointer :: children
    character(len=:), allocatable :: buffer, buffer2, modifier
    type(TListInt) :: li
    type(TListInt), allocatable :: liN(:)
    type(TListIntR1), allocatable :: li1N(:)
    type(TListReal), allocatable :: lrN(:)
    type(TListCharLc), allocatable :: skFiles(:,:)
    type(TListString) :: lStr
    type(TListIntR1), allocatable :: angShells(:)
    logical, allocatable :: repPoly(:,:)
    integer :: iSp1, iSp2, ii
    character(lc) :: prefix, suffix, separator, elem1, elem2, strTmp, str2Tmp
    character(lc) :: errorStr
    logical :: tLower
    integer, allocatable :: pTmpI1(:)
    real(dp) :: rTmp
    integer, allocatable :: iTmpN(:)
    integer :: nShell, skInterMeth
    real(dp) :: rSKCutOff
    character(len=:), allocatable :: searchPath(:)
    character(len=:), allocatable :: strOut, strJoin
    logical :: isHalogenXCorr

    !> For hybrid functional calculations
    type(THybridXcSKTag) :: hybridXcSK

    ctrl%hamiltonian = hamiltonianTypes%dftb

    call readMaxAngularMomentum(node, geo, angShells)

    ! Orbitals and angular momenta for the given shells (once the SK files contain the full
    ! information about the basis, this will be moved to the SK reading routine).
    allocate(slako%orb)
    call setupOrbitals(slako%orb, geo, angShells)

    ! Slater-Koster files
    call getParamSearchPaths(searchPath)
    strJoin = joinPathsPrettyErr(searchPath)
    allocate(skFiles(geo%nSpecies, geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
      do iSp2 = 1, geo%nSpecies
        call init(skFiles(iSp2, iSp1))
      end do
    end do
    call getChildValue(node, "SlaterKosterFiles", value1, child=child)
    call getNodeName(value1, buffer)
    select case(buffer)
    case ("type2filenames")
      call getChildValue(value1, "Prefix", buffer2, "")
      prefix = unquote(buffer2)
      call getChildValue(value1, "Suffix", buffer2, "")
      suffix = unquote(buffer2)
      call getChildValue(value1, "Separator", buffer2, "")
      separator = unquote(buffer2)
      call getChildValue(value1, "LowerCaseTypeName", tLower, .false.)
      do iSp1 = 1, geo%nSpecies
        if (tLower) then
          elem1 = tolower(geo%speciesNames(iSp1))
        else
          elem1 = geo%speciesNames(iSp1)
        end if
        do iSp2 = 1, geo%nSpecies
          if (tLower) then
            elem2 = tolower(geo%speciesNames(iSp2))
          else
            elem2 = geo%speciesNames(iSp2)
          end if
          strTmp = trim(prefix) // trim(elem1) // trim(separator) // trim(elem2) // trim(suffix)
          call findFile(searchPath, strTmp, strOut)
          if (.not. allocated(strOut)) then
            call detailedError(value1, "SK file with generated name '" // trim(strTmp)&
                & // "' not found." // newline // "   (search path(s): " // strJoin // ").")
          end if
          strTmp = strOut
          call append(skFiles(iSp2, iSp1), strTmp)
        end do
      end do
    case default
      call setUnprocessed(value1)
      call getChildValue(child, "Prefix", buffer2, "")
      prefix = unquote(buffer2)
      call getChild(child, "Suffix", child2, requested=.false.)
      if (associated(child2)) then
        call detailedError(child2, "Keyword requires SlaterKosterFiles = Type2Filenames {")
      end if
      call getChild(child, "Separator", child2, requested=.false.)
      if (associated(child2)) then
        call detailedError(child2, "Keyword requires SlaterKosterFiles = Type2Filenames {")
      end if
      call getChild(child, "LowerCaseTypeName", child2, requested=.false.)
      if (associated(child2)) then
        call detailedError(child2, "Keyword requires SlaterKosterFiles = Type2Filenames {")
      end if
      do iSp1 = 1, geo%nSpecies
        do iSp2 = 1, geo%nSpecies
          strTmp = trim(geo%speciesNames(iSp1)) // "-" // trim(geo%speciesNames(iSp2))
          call init(lStr)
          call getChildValue(child, trim(strTmp), lStr, child=child2)
          if (len(lStr) /= len(angShells(iSp1)) * len(angShells(iSp2))) then
            write(errorStr, "(A,I0,A,I0)") "Incorrect number of Slater-Koster files for " //&
                & trim(strTmp) // ", expected ", len(angShells(iSp1)) * len(angShells(iSp2)),&
                & " but received ", len(lStr)
            call detailedError(child2, errorStr)
          end if
          do ii = 1, len(lStr)
            call get(lStr, str2Tmp, ii)
            strTmp = trim(prefix) // str2Tmp
            call findFile(searchPath, strTmp, strOut)
            if (.not. allocated(strOut)) then
              call detailedError(child2, "SK file '" // trim(strTmp) // "' not found." // newline&
                  & // "   (search path(s): " // strJoin // ").")
            end if
            strTmp = strOut
            call append(skFiles(iSp2, iSp1), strTmp)
          end do
          call destruct(lStr)
        end do
      end do
    end select

    ! Which repulsive is defined by polynomial? (Default: None)
    allocate(repPoly(geo%nSpecies, geo%nSpecies))
    call getChildValue(node, "PolynomialRepulsive", value1, "", child=child, list=.true.,&
        & allowEmptyValue=.true., dummyValue=.true.)
    call getNodeName2(value1, buffer)
    if (buffer == "" .and. hasInlineData(child)) buffer = textNodeName
    select case (buffer)
    case ("")
      repPoly(:,:) = .false.
    case("setforall")
      call getChildValue(value1, "", repPoly(1,1))
      repPoly(:,:) = repPoly(1,1)
    case default
      do iSp1 = 1, geo%nSpecies
        do iSp2 = 1, geo%nSpecies
          strTmp = trim(geo%speciesNames(iSp1)) // "-" &
              &// trim(geo%speciesNames(iSp2))
          call getChildValue(child, trim(strTmp), repPoly(iSp2, iSp1), .false.)
        end do
      end do
      if (.not. all(repPoly .eqv. transpose(repPoly))) then
        call detailedError(value1, "Asymmetric definition (both A-B and B-A must&
            & be defined for using polynomial repulsive)")
      end if
    end select

    call parseChimes(node, ctrl%chimesRepInput)

    ! SCC
    call getChildValue(node, "SCC", ctrl%tSCC, .false.)

    call parseHybridBlock(node, ctrl%hybridXcInp, ctrl, geo, skFiles)

    if (allocated(ctrl%hybridXcInp)) then
      if (.not.ctrl%tSCC) then
        call detailedError(node, "Hybrid calculations require SCC = Yes")
      end if
    end if

    if (ctrl%tSCC) then
      call getChildValue(node, "ShellResolvedSCC", ctrl%tShellResolved, .false.)
    else
      ctrl%tShellResolved = .false.
    end if

    call getChildValue(node, "OldSKInterpolation", ctrl%oldSKInter, .false.)
    if (ctrl%oldSKInter) then
      skInterMeth = skEqGridOld
    else
      skInterMeth = skEqGridNew
    end if

    if (.not. allocated(ctrl%hybridXcInp)) then
      call getChild(node, "TruncateSKRange", child, requested=.false.)
      if (associated(child)) then
        call warning("Artificially truncating the SK table, this is normally a bad idea!")
        call SKTruncations(child, rSKCutOff, skInterMeth)
        call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tShellResolved,&
            & skInterMeth, repPoly, rSKCutOff)
      else
        rSKCutOff = 0.0_dp
        call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tShellResolved,&
            & skInterMeth, repPoly)
      end if
    else
      call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tShellResolved,&
          & skInterMeth, repPoly, hybridXcSK=hybridXcSK)
      ctrl%hybridXcInp%omega = hybridXcSK%omega
      ctrl%hybridXcInp%camAlpha = hybridXcSK%camAlpha
      ctrl%hybridXcInp%camBeta = hybridXcSK%camBeta
    end if

    do iSp1 = 1, geo%nSpecies
      call destruct(angShells(iSp1))
      do iSp2 = 1, geo%nSpecies
        call destruct(skFiles(iSp2, iSp1))
      end do
    end do
    deallocate(angShells)
    deallocate(skFiles)
    deallocate(repPoly)

    ! SCC parameters
    ifSCC: if (ctrl%tSCC) then

      ! get charge mixing options
      call readSccOptions(node, ctrl, geo)

      ! DFTB hydrogen bond corrections
      call readHCorrection(node, geo, ctrl)

      !> TI-DFTB varibles for Delta DFTB
      call getChild(node, "NonAufbau", child, requested=.false.)
      if (associated(child)) then
        ctrl%isNonAufbau = .true.
        call getChildValue(child, "SpinPurify", ctrl%isSpinPurify, .true.)
        call getChildValue(child, "GroundGuess", ctrl%isGroundGuess, .false.)
        ctrl%tSpin = .true.
        ctrl%t2Component = .false.
        ctrl%nrSpinPol = 0.0_dp
        ctrl%tSpinSharedEf = .false.
      else
        ctrl%isNonAufbau = .false.
      end if

    end if ifSCC

    ! Customize the reference atomic charges for virtual doping
    call readCustomReferenceOcc(node, slako%orb, slako%skOcc, geo, &
        & ctrl%customOccAtoms, ctrl%customOccFillings)

    ! Spin calculation
    if (ctrl%reksInp%reksAlg == reksTypes%noReks  .and. .not.ctrl%isNonAufbau) then
    #:if WITH_TRANSPORT
      call readSpinPolarisation(node, ctrl, geo, tp)
    #:else
      call readSpinPolarisation(node, ctrl, geo)
    #:endif
    end if

    ! temporararily removed until debugged
    !if (.not. ctrl%tscc) then
    !  !! In a non-SCC calculation it is possible to upload charge shifts
    !  !! This is useful if the calculation can jump directly to the Analysis block
    !  call getChildValue(node, "ReadShifts", ctrl%tReadShifts, .false.)
    !end if
    ctrl%tReadShifts = .false.

    ! External fields and potentials
    call readExternal(node, ctrl, geo)

    ! Non-self-consistent spin-orbit coupling
    call readSpinOrbit(node, ctrl, geo, slako%orb)

    ! Electronic solver
  #:if WITH_TRANSPORT
    call readSolver(node, ctrl, geo, tp, greendens, poisson)

    if (tp%taskUpload) then
      ! Initialise variable, but unused
      ctrl%nrChrg =  0.0_dp
    else
      ! Charge
      call getChildValue(node, "Charge", ctrl%nrChrg, 0.0_dp)
    end if
  #:else
    call readSolver(node, ctrl, geo, poisson)

    ! Charge
    call getChildValue(node, "Charge", ctrl%nrChrg, 0.0_dp)
  #:endif

    ! K-Points
    call readKPoints(node, ctrl, geo, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    if (ctrl%tscc) then

      call getChild(node, "OrbitalPotential", child, requested=.false.)
      if (associated(child)) then
        allocate(ctrl%dftbUInp)
        call getChildValue(child, "Functional", buffer, "fll")
        select case(tolower(buffer))
        case ("fll")
          ctrl%dftbUInp%iFunctional = plusUFunctionals%fll
        case ("psic")
          ctrl%dftbUInp%iFunctional = plusUFunctionals%pSic
        case default
          call detailedError(child,"Unknown orbital functional :"// buffer)
        end select

        allocate(ctrl%dftbUInp%nUJ(geo%nSpecies))
        ctrl%dftbUInp%nUJ(:) = 0

        ! to hold list of U-J values for each atom
        allocate(lrN(geo%nSpecies))
        ! to hold count of U-J values for each atom
        allocate(liN(geo%nSpecies))
        ! to hold list of shells for each U-J block of values
        allocate(li1N(geo%nSpecies))

        do iSp1 = 1, geo%nSpecies
          call init(lrN(iSp1))
          call init(liN(iSp1))
          call init(li1N(iSp1))
          call getChildren(child, trim(geo%speciesNames(iSp1)), children)
          ctrl%dftbUInp%nUJ(iSp1) = getLength(children)
          do ii = 1, ctrl%dftbUInp%nUJ(iSp1)
            call getItem1(children, ii, child2)

            call init(li)
            call getChildValue(child2,"Shells",li)
            allocate(pTmpI1(len(li)))
            call asArray(li,pTmpI1)
            call append(li1N(iSp1),pTmpI1)
            call append(liN(iSp1),size(pTmpI1))
            deallocate(pTmpI1)
            call destruct(li)
            call getChildValue(child2, "uj", rTmp, 0.0_dp, modifier=modifier, &
                & child=child3)
            call convertUnitHsd(modifier, energyUnits, child3, rTmp)
            if (rTmp < 0.0_dp) then
              write(errorStr,"(F12.8)")rTmp
              call detailedError(child2,"Negative value of U-J:"//errorStr)
            end if
            if (rTmp <= 1.0E-10_dp) then
              write(errorStr,"(F12.8)")rTmp
              call detailedError(child2,"Invalid value of U-J, too small: " &
                  & //errorStr)
            end if
            call append(lrN(iSp1),rTmp)
          end do
          call destroyNodeList(children)
        end do

        do iSp1 = 1, geo%nSpecies
          ctrl%dftbUInp%nUJ(iSp1) = len(lrN(iSp1))
        end do
        allocate(ctrl%dftbUInp%UJ(maxval(ctrl%dftbUInp%nUJ),geo%nSpecies))
        ctrl%dftbUInp%UJ(:,:) = 0.0_dp
        allocate(ctrl%dftbUInp%niUJ(maxval(ctrl%dftbUInp%nUJ),geo%nSpecies))
        ctrl%dftbUInp%niUJ(:,:) = 0
        do iSp1 = 1, geo%nSpecies
          call asArray(lrN(iSp1),ctrl%dftbUInp%UJ(1:len(lrN(iSp1)),iSp1))
          allocate(iTmpN(len(liN(iSp1))))
          call asArray(liN(iSp1),iTmpN)
          ctrl%dftbUInp%niUJ(1:len(liN(iSp1)),iSp1) = iTmpN(:)
          deallocate(iTmpN)
          call destruct(lrN(iSp1))
          call destruct(liN(iSp1))
        end do
        allocate(ctrl%dftbUInp%iUJ(maxval(ctrl%dftbUInp%niUJ),&
            & maxval(ctrl%dftbUInp%nUJ),geo%nSpecies))
        ctrl%dftbUInp%iUJ(:,:,:) = 0
        do iSp1 = 1, geo%nSpecies
          do ii = 1, ctrl%dftbUInp%nUJ(iSp1)
            allocate(iTmpN(ctrl%dftbUInp%niUJ(ii,iSp1)))
            call get(li1N(iSp1),iTmpN,ii)
            ctrl%dftbUInp%iUJ(1:ctrl%dftbUInp%niUJ(ii,iSp1),ii,iSp1) = iTmpN(:)
            deallocate(iTmpN)
          end do
          call destruct(li1N(iSp1))
        end do

        deallocate(li1N)
        deallocate(lrN)
        deallocate(liN)

        ! check input values
        allocate(iTmpN(slako%orb%mShell))
        do iSp1 = 1, geo%nSpecies
          iTmpN = 0
          ! loop over number of blocks for that species
          do ii = 1, ctrl%dftbUInp%nUJ(iSp1)
            iTmpN(ctrl%dftbUInp%iUJ(1:ctrl%dftbUInp%niUJ(ii,iSp1),ii,iSp1)) = &
                & iTmpN(ctrl%dftbUInp%iUJ(1:ctrl%dftbUInp%niUJ(ii,iSp1),ii,iSp1)) + 1
          end do
          if (any(iTmpN(:)>1)) then
            write(stdout, *)'Multiple copies of shells present in OrbitalPotential!'
            write(stdout, "(A,A3,A,I2)") &
                & 'The count for the occurrence of shells of species ', &
                & trim(geo%speciesNames(iSp1)),' are:'
            write(stdout, *)iTmpN(1:slako%orb%nShell(iSp1))
            call abortProgram()
          end if
        end do
        deallocate(iTmpN)

      end if

      ! On-site
      call getChildValue(node, "OnSiteCorrection", value1, "", child=child, allowEmptyValue=.true.,&
          & dummyValue=.true.)
      if (associated(value1)) then
        allocate(ctrl%onSiteElements(slako%orb%mShell, slako%orb%mShell, 2, geo%nSpecies))
        do iSp1 = 1, geo%nSpecies
          call getChildValue(child, trim(geo%speciesNames(iSp1))//"uu",&
              & ctrl%onSiteElements(:slako%orb%nShell(iSp1), :slako%orb%nShell(iSp1), 1, iSp1))
          call getChildValue(child, trim(geo%speciesNames(iSp1))//"ud",&
              & ctrl%onSiteElements(:slako%orb%nShell(iSp1), :slako%orb%nShell(iSp1), 2, iSp1))
        end do
      end if

    end if

    ! Dispersion
    call getChildValue(node, "Dispersion", value1, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%dispInp)
      call readDispersion(child, geo, ctrl%dispInp, ctrl%nrChrg, ctrl%tSCC)
    end if

    ! Solvation
    call getChildValue(node, "Solvation", value1, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%solvInp)
      call readSolvation(child, geo, ctrl%solvInp)
      call getChildValue(value1, "RescaleSolvatedFields", ctrl%isSolvatedFieldRescaled, .true.)
    end if

    ! Electronic constraints
    call getChildValue(node, "ElectronicConstraints", value1, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true., list=.true.)
    if (associated(value1)) then
      allocate(ctrl%elecConstraintInp)
      call readElecConstraintInput(child, geo, ctrl%tSpin, ctrl%t2Component, ctrl%elecConstraintInp)
      if (.not. allocated(ctrl%elecConstraintInp%mullikenConstrs)) then
        call detailedWarning(child, "No electronic constraint specified")
        deallocate(ctrl%elecConstraintInp)
      end if
    end if

    if (ctrl%tLatOpt .and. .not. geo%tPeriodic) then
      call error("Lattice optimisation only applies for periodic structures.")
    end if

    if (ctrl%tSCC) then
    #:if WITH_TRANSPORT
      call readElectrostatics(node, ctrl, geo, tp, poisson)
    #:else
      call readElectrostatics(node, ctrl, geo, poisson)
    #:endif
    end if

    ! Multipole expansion
    ctrl%isMdftb = .false.
    call readMdftb(node, ctrl, geo)

    ! Third order stuff
    ctrl%t3rd = .false.
    ctrl%t3rdFull = .false.
    if (ctrl%tSCC) then
      call getChildValue(node, "ThirdOrder", ctrl%t3rd, .false.)
      call getChildValue(node, "ThirdOrderFull", ctrl%t3rdFull, .false.)
      if (ctrl%t3rd .and. ctrl%t3rdFull) then
        call detailedError(node, "You must choose either ThirdOrder or&
            & ThirdOrderFull")
      end if
      if (ctrl%t3rd .and. ctrl%tShellResolved) then
        call error("Only full third-order DFTB is compatible with orbital&
            & resolved SCC")
      end if
      if (ctrl%t3rd .or. ctrl%t3rdFull) then
        call getChild(node, 'HubbardDerivs', child, requested=.true.)
        allocate(ctrl%HubDerivs(slako%orb%mShell, geo%nSpecies))
        ctrl%hubDerivs(:,:) = 0.0_dp
        do iSp1 = 1, geo%nSpecies
          nShell = slako%orb%nShell(iSp1)
          if (ctrl%tShellResolved) then
            call getChildValue(child, geo%speciesNames(iSp1),&
                & ctrl%hubDerivs(1:nShell, iSp1))
          else
            call getChildValue(child, geo%speciesNames(iSp1),&
                & ctrl%hubDerivs(1, iSp1))
            ctrl%hubDerivs(2:nShell, iSp1) = ctrl%hubDerivs(1, iSp1)
          end if
        end do
        if (ctrl%t3rd) then
          allocate(ctrl%thirdOrderOn(geo%nAtom, 2))
          ctrl%thirdOrderOn(:,1) = 0.0_dp
          ctrl%thirdOrderOn(:,2) = ctrl%hubDerivs(1, geo%species)
        end if

        ! Halogen correction to the DFTB3 model
        isHalogenXCorr =&
            & any([(any(halogenXSpecies1(ii) == geo%speciesNames), ii=1, size(halogenXSpecies1))])&
            & .and. &
            & any([(any(halogenXSpecies2(ii) == geo%speciesNames), ii=1, size(halogenXSpecies2))])

        if (isHalogenXCorr) then
          call getChildValue(node, "HalogenXCorr", ctrl%tHalogenX, .false.)
        end if

      end if
    end if

    call readDifferentiation(node, ctrl)

    if (ctrl%tSCC) then
      ! Force type
      call readForceOptions(node, ctrl)
    else
      ctrl%forceType = forceTypes%orig
    end if

    call readCustomisedHubbards(node, geo, slako%orb, ctrl%tShellResolved, ctrl%hubbU)

  end subroutine readDFTBHam


  !> Reads xTB-Hamiltonian
#:if WITH_TRANSPORT
  subroutine readXTBHam(node, ctrl, geo, tp, greendens, poisson, errStatus)
#:else
  subroutine readXTBHam(node, ctrl, geo, poisson, errStatus)
#:endif

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer
    character(len=:), allocatable :: searchPath(:)
    integer :: method
    character(len=:), allocatable :: paramFile, paramTmp
    type(TOrbitals) :: orb

    ctrl%hamiltonian = hamiltonianTypes%xtb

    allocate(ctrl%tbliteInp)
    call ctrl%tbliteInp%setupGeometry(geo%nAtom, geo%species, geo%coords, geo%speciesNames,&
        & geo%latVecs)

    call getChild(node, "Method", child, requested=.false.)
    if (associated(child)) then
      call getChildValue(child, "", buffer)
      select case(unquote(buffer))
      case default
        call detailedError(child, "Unknown method "//buffer//" for xTB Hamiltonian")
      case("GFN1-xTB")
        method = tbliteMethod%gfn1xtb
      case("GFN2-xTB")
        method = tbliteMethod%gfn2xtb
      case("IPEA1-xTB")
        method = tbliteMethod%ipea1xtb
      end select
      call ctrl%tbliteInp%setupCalculator(method)
      ctrl%tbliteInp%info%name = trim(unquote(buffer))
    else
      call getChildValue(node, "ParameterFile", value1, "", child=child, allowEmptyValue=.true.,&
          & dummyValue=.true.)
      if (associated(value1)) then
        call getChildValue(child, "", buffer)
        paramFile = trim(unquote(buffer))
        call getParamSearchPaths(searchPath)
        call findFile(searchPath, paramFile, paramTmp)
        if (allocated(paramTmp)) call move_alloc(paramTmp, paramFile)
        write(stdOut, '(a)') "Using parameter file '"//paramFile//"' for xTB Hamiltonian"
        call ctrl%tbliteInp%setupCalculator(paramFile)
      else
        call detailedError(node, "Either a Method or ParameterFile must be specified for xTB")
      end if
    end if

    call getChildValue(node, "ShellResolvedSCC", ctrl%tShellResolved, .true.)

    ! SCC parameters
    call getChildValue(node, "SCC", ctrl%tSCC, .true.)
    ifSCC: if (ctrl%tSCC) then

      ! get charge mixing options etc.
      call readSccOptions(node, ctrl, geo)

      !> TI-DFTB varibles for Delta DFTB
      call getChild(node, "NonAufbau", child, requested=.false.)
      if (associated(child)) then
        ctrl%isNonAufbau = .true.
        call getChildValue(child, "SpinPurify", ctrl%isSpinPurify, .true.)
        call getChildValue(child, "GroundGuess", ctrl%isGroundGuess, .false.)
        ctrl%tSpin = .true.
        ctrl%t2Component = .false.
        ctrl%nrSpinPol = 0.0_dp
        ctrl%tSpinSharedEf = .false.
      else
        ctrl%isNonAufbau = .false.
      end if

    end if ifSCC

    ! Spin calculation
    if (ctrl%reksInp%reksAlg == reksTypes%noReks .and. .not.ctrl%isNonAufbau .and. ctrl%tSCC) then
    #:if WITH_TRANSPORT
      call readSpinPolarisation(node, ctrl, geo, tp)
    #:else
      call readSpinPolarisation(node, ctrl, geo)
    #:endif
    end if

    ! temporararily removed until debugged
    !if (.not. ctrl%tscc) then
    !  !! In a non-SCC calculation it is possible to upload charge shifts
    !  !! This is useful if the calculation can jump directly to the Analysis block
    !  call getChildValue(node, "ReadShifts", ctrl%tReadShifts, .false.)
    !end if
    ctrl%tReadShifts = .false.

    ! External fields and potentials
    call readExternal(node, ctrl, geo)

    ! Non-self-consistent spin-orbit coupling
    call ctrl%tbliteInp%setupOrbitals(geo%species, orb)
    call readSpinOrbit(node, ctrl, geo, orb)

    ! Electronic solver
  #:if WITH_TRANSPORT
    call readSolver(node, ctrl, geo, tp, greendens, poisson)

    if (tp%taskUpload) then
      ! Initialise, but unused
      ctrl%nrChrg =  0.0_dp
    else
      ! Charge
      call getChildValue(node, "Charge", ctrl%nrChrg, 0.0_dp)
    end if
  #:else
    call readSolver(node, ctrl, geo, poisson)

    ! Charge
    call getChildValue(node, "Charge", ctrl%nrChrg, 0.0_dp)
  #:endif

    ! K-Points
    call readKPoints(node, ctrl, geo, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Dispersion
    call getChildValue(node, "Dispersion", value1, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%dispInp)
      call readDispersion(child, geo, ctrl%dispInp, ctrl%nrChrg, ctrl%tSCC)
    end if

    ! Solvation
    call getChildValue(node, "Solvation", value1, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%solvInp)
      call readSolvation(child, geo, ctrl%solvInp)
      call getChildValue(value1, "RescaleSolvatedFields", ctrl%isSolvatedFieldRescaled, .true.)
    end if

    if (ctrl%tLatOpt .and. .not. geo%tPeriodic) then
      call error("Lattice optimisation only applies for periodic structures.")
    end if

  #:if WITH_TRANSPORT
    call readElectrostatics(node, ctrl, geo, tp, poisson)
  #:else
    call readElectrostatics(node, ctrl, geo, poisson)
  #:endif

    ! Third order stuff
    ctrl%t3rd = .true.
    ctrl%t3rdFull = .false.

    call readDifferentiation(node, ctrl)

    if (ctrl%tSCC) then
      ! Force type
      call readForceOptions(node, ctrl)
    else
      ctrl%forceType = forceTypes%orig
    end if

    ! Electronic constraints
    call getChildValue(node, "ElectronicConstraints", value1, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true., list=.true.)
    if (associated(value1)) then
      allocate(ctrl%elecConstraintInp)
      call readElecConstraintInput(child, geo, ctrl%tSpin, ctrl%t2Component, ctrl%elecConstraintInp)
    end if

  end subroutine readXTBHam


  !> Reads in settings for spin orbit enabled calculations
  subroutine readSpinOrbit(node, ctrl, geo, orb)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Information about the orbitals of the species/atoms in the system
    class(TOrbitals), intent(in) :: orb

    type(hsd_table), pointer :: child, child2
    character(len=:), allocatable :: modifier
    integer :: iSp

    call getChild(node, "SpinOrbit", child, requested=.false.)
    if (.not. associated(child)) then
      ctrl%tSpinOrbit = .false.
      allocate(ctrl%xi(0,0))
    else
      if (ctrl%tSpin .and. .not. ctrl%t2Component) then
        call error("Spin-orbit coupling incompatible with collinear spin.")
      end if

      ctrl%tSpinOrbit = .true.
      ctrl%t2Component = .true.

      call getChildValue(child, "Dual", ctrl%tDualSpinOrbit, .true.)

      allocate(ctrl%xi(orb%mShell,geo%nSpecies), source = 0.0_dp)
      do iSp = 1, geo%nSpecies
        call getChildValue(child, geo%speciesNames(iSp), &
            & ctrl%xi(:orb%nShell(iSp),iSp), modifier=modifier, child=child2 )
        call convertUnitHsd(modifier, energyUnits, child2,&
            & ctrl%xi(:orb%nShell(iSp),iSp))
      end do
    end if

  end subroutine readSpinOrbit


  !> Read in maximal angular momenta or selected shells
  subroutine readMaxAngularMomentum(node, geo, angShells)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> List containing the angular momenta of the shells
    type(TListIntR1), allocatable, intent(out) :: angShells(:)

    type(hsd_table), pointer :: value1, child, child2
    character(len=:), allocatable :: buffer
    integer :: iSp1, ii, jj, kk
    character(lc) :: strTmp
    integer :: nShell
    character(1) :: tmpCh
    type(TListString) :: lStr
    logical :: tShellIncl(4), tFound
    integer :: angShell(maxL+1), angShellOrdered(maxL+1)

    do ii = 1, maxL+1
      angShellOrdered(ii) = ii - 1
    end do
    call getChild(node, "MaxAngularMomentum", child)
    allocate(angShells(geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
      call init(angShells(iSp1))
      call getChildValue(child, geo%speciesNames(iSp1), value1, child=child2)
      call getNodeName(value1, buffer)
      select case(buffer)
      case("selectedshells")
        call init(lStr)
        call getChildValue(value1, "", lStr)
        do ii = 1, len(lStr)
          call get(lStr, strTmp, ii)
          strTmp = tolower(unquote(trim(strTmp)))
          if (len_trim(strTmp) > 4 .or. len_trim(strTmp) < 1) then
            call detailedError(value1, "Invalid shell selection '" &
                &// trim(strTmp) &
                &// "'. Nr. of selected shells must be between 1 and 4.")
          end if
          tShellIncl(:) = .false.
          nShell = len_trim(strTmp)
          do jj = 1, nShell
            tmpCh = strTmp(jj:jj)
            tFound = .false.
            do kk = 1, size(shellNames)
              if (tmpCh == trim(shellNames(kk))) then
                if (tShellIncl(kk)) then
                  call detailedError(value1, "Double selection of the same shell&
                      & '" // tmpCh // "' in shell selection block '" &
                      &// trim(strTmp) // "'")
                end if
                tShellIncl(kk) = .true.
                angShell(jj) = kk - 1
                tFound = .true.
                exit
              end if
            end do
            if (.not. tFound) then
              call detailedError(value1, "Invalid shell name '" // tmpCh // "'")
            end if
          end do
          call append(angShells(iSp1), angShell(1:nShell))
        end do
        call destruct(lStr)

      case(textNodeName)
        call getChildValue(child2, "", buffer)
        strTmp = unquote(buffer)
        do jj = 1, size(shellNames)
          if (trim(strTmp) == trim(shellNames(jj))) then
            call append(angShells(iSp1), angShellOrdered(:jj))
          end if
        end do
        if (len(angShells(iSp1)) < 1) then
          call detailedError(child2, "Invalid orbital name '" // &
              &trim(strTmp) // "'")
        end if

      case default
        call getNodeHSDName(value1, buffer)
        call detailedError(child2, "Invalid shell specification method '" //&
            & buffer // "'")
      end select
    end do

  end subroutine readMaxAngularMomentum


  !> Setup information about the orbitals of the species/atoms from angShell lists
  subroutine setupOrbitals(orb, geo, angShells)

    !> Information about the orbitals of the species/atoms in the system
    class(TOrbitals), intent(out) :: orb

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> List containing the angular momenta of the shells,
    !> must be inout, since intoArray requires inout arguments
    type(TListIntR1), intent(inout) :: angShells(:)

    integer :: nShell, iSp1, iSh1, ii, jj, ind
    integer :: angShell(maxL+1)

    allocate(orb%nShell(geo%nSpecies))
    allocate(orb%nOrbSpecies(geo%nSpecies))
    allocate(orb%nOrbAtom(geo%nAtom))
    orb%mOrb = 0
    orb%mShell = 0
    do iSp1 = 1, geo%nSpecies
      orb%nShell(iSp1) = 0
      orb%nOrbSpecies(iSp1) = 0
      do ii = 1, len(angShells(iSp1))
        call intoArray(angShells(iSp1), angShell, nShell, ii)
        orb%nShell(iSp1) = orb%nShell(iSp1) + nShell
        do jj = 1, nShell
          orb%nOrbSpecies(iSp1) = orb%nOrbSpecies(iSp1) &
              &+ 2 * angShell(jj) + 1
        end do
      end do
    end do
    orb%mShell = maxval(orb%nShell)
    orb%mOrb = maxval(orb%nOrbSpecies)
    orb%nOrbAtom(:) = orb%nOrbSpecies(geo%species(:))
    orb%nOrb = sum(orb%nOrbAtom)

    allocate(orb%angShell(orb%mShell, geo%nSpecies))
    allocate(orb%iShellOrb(orb%mOrb, geo%nSpecies))
    allocate(orb%posShell(orb%mShell+1, geo%nSpecies))
    orb%angShell(:,:) = 0
    do iSp1 = 1, geo%nSpecies
      ind = 1
      iSh1 = 1
      do ii = 1, len(angShells(iSp1))
        call intoArray(angShells(iSp1), angShell, nShell, ii)
        do jj = 1, nShell
          orb%posShell(iSh1, iSp1) = ind
          orb%angShell(iSh1, iSp1) = angShell(jj)
          orb%iShellOrb(ind:ind+2*angShell(jj), iSp1) = iSh1
          ind = ind + 2 * angShell(jj) + 1
          iSh1 = iSh1 + 1
        end do
        orb%posShell(iSh1, iSp1) = ind
      end do
    end do

  end subroutine setupOrbitals


#:if WITH_TRANSPORT
  subroutine readElectrostatics(node, ctrl, geo, tp, poisson)
#:else
  subroutine readElectrostatics(node, ctrl, geo, poisson)
#:endif

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp
  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer

    ctrl%tPoisson = .false.

    ! Read in which kind of electrostatics method to use.
    call getChildValue(node, "Electrostatics", value1, "GammaFunctional", child=child)
    call getNodeName(value1, buffer)

    select case (buffer)

    case ("gammafunctional")
    #:if WITH_TRANSPORT
      if (tp%taskUpload .and. ctrl%tSCC) then
        call detailedError(value1, "GammaFunctional not available, if you upload contacts in an SCC&
            & calculation.")
      end if
    #:endif

    case ("poisson")
      if (.not. withPoisson) then
        call detailedError(value1, "Poisson not available as binary was built without the Poisson&
            &-solver")
      end if
      #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
        ctrl%tPoisson = .true.
        #:if WITH_TRANSPORT
          call readPoisson(value1, poisson, geo%tPeriodic, tp, geo%latVecs, ctrl%updateSccAfterDiag)
        #:else
          call readPoisson(value1, poisson, geo%tPeriodic, geo%latVecs, ctrl%updateSccAfterDiag)
        #:endif
      #:endblock

    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Unknown electrostatics '" // buffer // "'")
    end select

  end subroutine readElectrostatics


  !> Read in the mdftb parameters
  subroutine readMdftb(node, ctrl, geo)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    type(hsd_table), pointer :: value1, child, child2
    character(len=:), allocatable :: buffer
    integer :: iSp1

    ctrl%isMdftb = .false.
    if (ctrl%tSCC) then
      call getChildValue(node, "Mdftb", value1, "None", child=child, allowEmptyValue=.true.,&
          & dummyValue=.false.)
      if (associated(value1)) then
        call getNodeName(value1, buffer)
        select case(buffer)
        case("onecenterapproximation")
          ctrl%isMdftb = .true.
          allocate(ctrl%mdftbAtomicIntegrals)
          allocate(ctrl%mdftbAtomicIntegrals%DScaling(geo%nSpecies), source=1.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%QScaling(geo%nSpecies), source=1.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%SXPx(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%PxXDxxyy(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%PxXDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%PyYDxxyy(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%PzZDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%SXXS(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%PxXXPx(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%PyXXPy(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%SXXDxxyy(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%SXXDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%SYYDxxyy(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%SZZDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DxyXXDxy(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DyzXXDyz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DxxyyXXDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DzzXXDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DxxyyYYDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DzzZZDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DxzXZDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DyzYZDxxyy(geo%nSpecies), source=0.0_dp)

          call getChild(value1, 'AtomDIntegralScalings', child2, requested=.false.)
          if (associated(child2)) then
            do iSp1 = 1, geo%nSpecies
              call getChildValue(child2, trim(geo%speciesNames(iSp1)),&
                  & ctrl%mdftbAtomicIntegrals%DScaling(iSp1), 1.0_dp)
            end do
          end if

          call getChild(value1, 'AtomQIntegralScalings', child2, requested=.false.)
          if (associated(child2)) then
            do iSp1 = 1, geo%nSpecies
              call getChildValue(child2, trim(geo%speciesNames(iSp1)),&
                  & ctrl%mdftbAtomicIntegrals%QScaling(iSp1), 1.0_dp)
            end do
          end if

          call getChild(value1, 'OneCenterAtomIntegrals', child2, requested=.true.)
          do iSp1 = 1, geo%nSpecies
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":S|X|Px",&
                & ctrl%mdftbAtomicIntegrals%SXPx(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Px|X|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%PxXDxxyy (iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Px|X|Dzz",&
                & ctrl%mdftbAtomicIntegrals%PxXDzz(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Py|Y|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%PyYDxxyy(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Pz|Z|Dzz",&
                & ctrl%mdftbAtomicIntegrals%PzZDzz(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":S|XX|S",&
                & ctrl%mdftbAtomicIntegrals%SXXS(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Px|XX|Px",&
                & ctrl%mdftbAtomicIntegrals%PxXXPx(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Py|XX|Py",&
                & ctrl%mdftbAtomicIntegrals%PyXXPy(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":S|XX|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%SXXDxxyy(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":S|XX|Dzz",&
                & ctrl%mdftbAtomicIntegrals%SXXDzz(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":S|YY|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%SYYDxxyy(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":S|ZZ|Dzz",&
                & ctrl%mdftbAtomicIntegrals%SZZDzz(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Dxy|XX|Dxy",&
                & ctrl%mdftbAtomicIntegrals%DxyXXDxy(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Dyz|XX|Dyz",&
                & ctrl%mdftbAtomicIntegrals%DyzXXDyz(iSp1), 0.0_dp)
            !call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Dxx-yy|XX|Dzz",&
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Dzz|XX|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%DxxyyXXDzz(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Dzz|XX|Dzz",&
                & ctrl%mdftbAtomicIntegrals%DzzXXDzz(iSp1), 0.0_dp)
            !call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Dxx-yy|YY|Dzz",&
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Dzz|YY|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%DxxyyYYDzz(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Dzz|ZZ|Dzz",&
                & ctrl%mdftbAtomicIntegrals%DzzZZDzz(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Dxz|XZ|Dzz",&
                & ctrl%mdftbAtomicIntegrals%DxzXZDzz(iSp1), 0.0_dp)
            call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Dyz|YZ|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%DyzYZDxxyy(iSp1), 0.0_dp)
          end do
        case("none")
          ctrl%isMdftb = .false.
        case default
          call detailedError(child,"Unknown functions :"// buffer)
        end select
      end if
    end if

  end subroutine readMdftb


  !> Spin calculation
#:if WITH_TRANSPORT
  subroutine readSpinPolarisation(node, ctrl, geo, tp)
#:else
  subroutine readSpinPolarisation(node, ctrl, geo)
#:endif

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp
  #:endif

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer

    call hsd_rename_child(node, "SpinPolarization", "SpinPolarisation")
    call getChildValue(node, "SpinPolarisation", value1, "", child=child, allowEmptyValue=.true.)
    call getNodeName2(value1, buffer)
    select case(buffer)
    case ("")
      ctrl%tSpin = .false.
      ctrl%t2Component = .false.
      ctrl%nrSpinPol = 0.0_dp

    case ("colinear", "collinear")
      ctrl%tSpin = .true.
      ctrl%t2Component = .false.
      call getChildValue(value1, 'UnpairedElectrons', ctrl%nrSpinPol, 0.0_dp)
      call getChildValue(value1, 'RelaxTotalSpin', ctrl%tSpinSharedEf, .false.)
      if (.not. ctrl%tReadChrg) then
        call getInitialSpins(value1, geo, 1, ctrl%initialSpins)
      end if

    case ("noncolinear", "noncollinear")
      ctrl%tSpin = .true.
      ctrl%t2Component = .true.
      if (.not. ctrl%tReadChrg) then
        call getInitialSpins(value1, geo, 3, ctrl%initialSpins)
      end if

    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Invalid spin polarisation type '" //&
          & buffer // "'")
    end select

  end subroutine readSpinPolarisation


  ! External field(s) and potential(s)
  subroutine readExternal(node, ctrl, geo)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    type(hsd_table), pointer :: value1, child, child2, child3
    type(hsd_child_list), pointer :: children
    character(len=:), allocatable :: modifier, buffer, buffer2
    real(dp) :: rTmp
    type(TFileDescr) :: file
    integer :: ind, ii, iErr, nElem
    real(dp), allocatable :: tmpR1(:), tmpR2(:,:)
    type(TListRealR2) :: lCharges
    type(TListRealR1) :: lBlurs, lr1
    type(TListReal) :: lr
    type(TListInt) :: li

    call getChildValue(node, "ElectricField", value1, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! external applied field
    call getChild(child, "External", child2, requested=.false.)
    if (associated(child2)) then
      allocate(ctrl%electricField)
      ctrl%tMulliken = .true.
      call getChildValue(child2, "Strength", ctrl%electricField%EFieldStrength, modifier=modifier,&
          & child=child3)
      call convertUnitHsd(modifier, EFieldUnits, child3, ctrl%electricField%EFieldStrength)
      call getChildValue(child2, "Direction", ctrl%electricField%EfieldVector)
      if (sum(ctrl%electricField%EfieldVector**2) < 1e-8_dp) then
        call detailedError(child2,"Vector too small")
      else
        ctrl%electricField%EfieldVector = ctrl%electricField%EfieldVector&
            & / sqrt(sum(ctrl%electricField%EfieldVector**2))
      end if
      call getChildValue(child2, "Frequency", ctrl%electricField%EFieldOmega, 0.0_dp, &
          & modifier=modifier, child=child3)
      call convertUnitHsd(modifier, freqUnits, child3, ctrl%electricField%EFieldOmega)
      if (ctrl%electricField%EFieldOmega > 0.0) then
        ! angular frequency
        ctrl%electricField%EFieldOmega = 2.0_dp * pi * ctrl%electricField%EFieldOmega
        ctrl%electricField%isTDEfield = .true.
      else
        ctrl%electricField%isTDEfield = .false.
        ctrl%electricField%EFieldOmega = 0.0_dp
      end if
      ctrl%electricField%EfieldPhase = 0
      if (ctrl%electricField%isTDEfield) then
        call getChildValue(child2, "Phase", ctrl%electricField%EfieldPhase, 0)
      end if
    end if

    ctrl%nExtChrg = 0
    if (ctrl%hamiltonian == hamiltonianTypes%dftb) then

      call getChildren(child, "PointCharges", children)
      if (getLength(children) > 0) then
        ! Point charges present
        if (.not.ctrl%tSCC) then
          call error("External charges can only be used in an SCC calculation")
        end if
        call init(lCharges)
        call init(lBlurs)
        ctrl%nExtChrg = 0
        do ii = 1, getLength(children)
          call getItem1(children, ii, child2)
          call getChildValue(child2, "CoordsAndCharges", value1, modifier=modifier, child=child3)
          call getNodeName(value1, buffer)
          select case(buffer)
          case (textNodeName)
            call init(lr1)
            call getChildValue(child3, "", 4, lr1, modifier=modifier)
            allocate(tmpR2(4, len(lr1)))
            call asArray(lr1, tmpR2)
            ctrl%nExtChrg = ctrl%nExtChrg + len(lr1)
            call destruct(lr1)
          case ("directread")
            call getChildValue(value1, "Records", ind)
            call getChildValue(value1, "File", buffer2)
            allocate(tmpR2(4, ind))
            call openFile(file, unquote(buffer2), mode="r", iostat=iErr)
            if (iErr /= 0) then
              call detailedError(value1, "Could not open file '"&
                  & // trim(unquote(buffer2)) // "' for direct reading" )
            end if
            read(file%unit, *, iostat=iErr) tmpR2
            if (iErr /= 0) then
              call detailedError(value1, "Error during direct reading '"&
                  & // trim(unquote(buffer2)) // "'")
            end if
            call closeFile(file)
            ctrl%nExtChrg = ctrl%nExtChrg + ind
          case default
            call detailedError(value1, "Invalid block name")
          end select
          call convertUnitHsd(modifier, lengthUnits, child3, tmpR2(1:3,:))
          call append(lCharges, tmpR2)
          call getChildValue(child2, "GaussianBlurWidth", rTmp, 0.0_dp, modifier=modifier,&
              & child=child3)
          if (rTmp < 0.0_dp) then
            call detailedError(child3, "Gaussian blur width may not be negative")
          end if
          call convertUnitHsd(modifier, lengthUnits, child3, rTmp)
          allocate(tmpR1(size(tmpR2, dim=2)))
          tmpR1(:) = rTmp
          call append(lBlurs, tmpR1)
          deallocate(tmpR1)
          deallocate(tmpR2)
        end do

        allocate(ctrl%extChrg(4, ctrl%nExtChrg))
        ind = 1
        do ii = 1, len(lCharges)
          call intoArray(lCharges, ctrl%extChrg(:, ind:), nElem, ii)
          ind = ind + nElem
        end do
        call destruct(lCharges)

        allocate(ctrl%extChrgBlurWidth(ctrl%nExtChrg))
        ind = 1
        do ii = 1, len(lBlurs)
          call intoArray(lBlurs, ctrl%extChrgBlurWidth(ind:), nElem, ii)
          ind = ind + nElem
        end do
        call destruct(lBlurs)
        call destroyNodeList(children)
      end if

    else

      call getChildren(child, "PointCharges", children)
      if (getLength(children) > 0) then
        call detailedError(child, "External charges are not currently supported for this model")
      end if

    end if

    call getChild(node, "AtomSitePotential", child, requested=.false.)
    if (associated(child)) then
      allocate(ctrl%atomicExtPotential)

      call getChild(child, "Net", child2, requested=.false.)
      if (associated(child2)) then
        ! onsites
        ctrl%tNetAtomCharges = .true.
        call init(li)
        call init(lr)
        call getChildValue(child2, "Atoms", li)
        call getChildValue(child2, "Vext", lr, modifier=modifier, child=child3)
        if (len(li) /= len(lr)) then
          call detailedError(child2, "Mismatch in number of sites and potentials")
        end if
        allocate(ctrl%atomicExtPotential%iAtOnSite(len(li)))
        call asArray(li, ctrl%atomicExtPotential%iAtOnSite)
        allocate(ctrl%atomicExtPotential%VextOnSite(len(lr)))
        call asArray(lr, ctrl%atomicExtPotential%VextOnSite)
        call convertUnitHsd(modifier, energyUnits, child3, ctrl%atomicExtPotential%VextOnSite)
        call destruct(li)
        call destruct(lr)
      end if

      call getChild(child, "Gross", child2, requested=.false.)
      if (associated(child2)) then
        ! atomic
        call init(li)
        call init(lr)
        call getChildValue(child2, "Atoms", li)
        call getChildValue(child2, "Vext", lr, modifier=modifier, child=child3)
        if (len(li) /= len(lr)) then
          call detailedError(child2, "Mismatch in number of sites and potentials")
        end if
        allocate(ctrl%atomicExtPotential%iAt(len(li)))
        call asArray(li, ctrl%atomicExtPotential%iAt)
        allocate(ctrl%atomicExtPotential%Vext(len(lr)))
        call asArray(lr, ctrl%atomicExtPotential%Vext)
        call convertUnitHsd(modifier, energyUnits, child3, ctrl%atomicExtPotential%Vext)
        call destruct(li)
        call destruct(lr)
      end if

      if (.not.allocated(ctrl%atomicExtPotential%iAt)&
          & .and. .not.allocated(ctrl%atomicExtPotential%iAtOnSite)) then
        call detailedError(child, "No atomic potentials specified")
      end if

    end if

  end subroutine readExternal


  !> Filling of electronic levels
  subroutine readFilling(node, ctrl, geo, temperatureDefault)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to test for periodicity
    type(TGeometry), intent(in) :: geo

    !> Default temperature for filling
    real(dp), intent(in) :: temperatureDefault

    type(hsd_table), pointer :: value1, child, child2, child3, field
    character(len=:), allocatable :: buffer, modifier
    character(lc) :: errorStr

    call getChildValue(node, "Filling", value1, "Fermi", child=child)
    call getNodeName(value1, buffer)

    select case (buffer)
    case ("fermi")
      ctrl%iDistribFn = fillingTypes%Fermi ! Fermi function
    case ("gaussian")
      ctrl%iDistribFn = fillingTypes%Methfessel ! Gauss function broadening of levels (0th order MP)
    case ("methfesselpaxton")
      ! Set the order of the Methfessel-Paxton step function approximation, defaulting to 1st order
      call getChildValue(value1, "Order", ctrl%iDistribFn, 1)
      if (ctrl%iDistribFn < 1) then
        call getNodeHSDName(value1, buffer)
        select case(ctrl%iDistribFn)
        case (0)
          write(errorStr, "(A)")"Methfessel-Paxton filling order 0 is equivalent to gaussian&
              & smearing"
          call detailedWarning(child, errorStr)
        case default
          write(errorStr, "(A,A,A,I4)")"Filling order must be above zero '", buffer,"' :",&
              &ctrl%iDistribFn
          call detailedError(child, errorStr)
        end select
      end if
      ctrl%iDistribFn = ctrl%iDistribFn + fillingTypes%Methfessel
    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Invalid filling method '" //buffer// "'")
    end select

    if (.not. ctrl%tSetFillingTemp) then
      call getChildValue(value1, "Temperature", ctrl%tempElec, temperatureDefault, &
          &modifier=modifier, child=field)
      call convertUnitHsd(modifier, energyUnits, field, ctrl%tempElec)
      if (ctrl%tempElec < minTemp) then
        ctrl%tempElec = minTemp
      end if
    end if

    call getChild(value1, "FixedFermiLevel", child=child2, modifier=modifier, requested=.false.)
    ctrl%tFixEf = associated(child2)
    if (ctrl%tFixEf) then
      if (ctrl%tSpin .and. .not.ctrl%t2Component) then
        allocate(ctrl%Ef(2))
      else
        allocate(ctrl%Ef(1))
      end if
      call getChildValue(child2, "", ctrl%Ef, modifier=modifier, child=child3)
      call convertUnitHsd(modifier, energyUnits, child3, ctrl%Ef)
    end if

    if (geo%tPeriodic .and. .not.ctrl%tFixEf) then
      call getChildValue(value1, "IndependentKFilling", ctrl%tFillKSep, .false.)
    end if

  end subroutine readFilling


  !> Electronic Solver
#:if WITH_TRANSPORT
  subroutine readSolver(node, ctrl, geo, tp, greendens, poisson)
#:else
  subroutine readSolver(node, ctrl, geo, poisson)
#:endif

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer, modifier

    integer :: iTmp

    ! Electronic solver
    call getChildValue(node, "Solver", value1, "RelativelyRobust")
    call getNodeName(value1, buffer)

    select case(buffer)

    case ("qr")
      ctrl%solver%isolver = electronicSolverTypes%qr

    case ("divideandconquer")
      ctrl%solver%isolver = electronicSolverTypes%divideandconquer

    case ("relativelyrobust")
      ctrl%solver%isolver = electronicSolverTypes%relativelyrobust

    case ("magma")
  #:if WITH_MAGMA
      ctrl%solver%isolver = electronicSolverTypes%magmaGvd
      call getChildValue(value1, "DensityMatrixGPU", ctrl%isDmOnGpu, .true.)
  #:else
      call detailedError(node, "DFTB+ must be compiled with MAGMA support in order to enable&
          & this solver")
  #:endif

    case ("elpa")
      allocate(ctrl%solver%elsi)
      call getChildValue(value1, "Sparse", ctrl%solver%elsi%elsiCsr, .false.)
      if (ctrl%solver%elsi%elsiCsr) then
        ctrl%solver%isolver = electronicSolverTypes%elpadm
      else
        ctrl%solver%isolver = electronicSolverTypes%elpa
      end if
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      call getChildValue(value1, "Mode", ctrl%solver%elsi%elpaSolver, 2)
      call getChildValue(value1, "Autotune", ctrl%solver%elsi%elpaAutotune, .false.)
      call getChildValue(value1, "Gpu", ctrl%solver%elsi%elpaGpu, .false., child=child)
      #:if not WITH_GPU
        if (ctrl%solver%elsi%elpaGpu) then
          call detailedError(child, "DFTB+ must be compiled with GPU support in order to enable&
              & the GPU acceleration for the ELPA solver")
        end if
      #:endif

    case ("omm")
      ctrl%solver%isolver = electronicSolverTypes%omm
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      call getChildValue(value1, "nIterationsELPA", ctrl%solver%elsi%ommIterationsElpa, 5)
      call getChildValue(value1, "Tolerance", ctrl%solver%elsi%ommTolerance, 1.0E-10_dp)
      call getChildValue(value1, "Choleskii", ctrl%solver%elsi%ommCholesky, .true.)

    case ("pexsi")
      ctrl%solver%isolver = electronicSolverTypes%pexsi
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
    #:if ELSI_VERSION > 2.5
      call getChildValue(value1, "Method", ctrl%solver%elsi%pexsiMethod, 3)
    #:else
      call getChildValue(value1, "Method", ctrl%solver%elsi%pexsiMethod, 2)
    #:endif
      select case(ctrl%solver%elsi%pexsiMethod)
      case(1)
        iTmp = 60
      case(2)
        iTmp = 20
      case(3)
        iTmp = 30
      end select
      call getChildValue(value1, "Poles", ctrl%solver%elsi%pexsiNPole, iTmp)
      if (ctrl%solver%elsi%pexsiNPole < 10) then
        call detailedError(value1, "Too few PEXSI poles")
      end if
      select case(ctrl%solver%elsi%pexsiMethod)
      case(1)
        if (mod(ctrl%solver%elsi%pexsiNPole,10) /= 0 .or. ctrl%solver%elsi%pexsiNPole > 120) then
          call detailedError(value1, "Unsupported number of PEXSI poles for method 1")
        end if
      case(2,3)
        if (mod(ctrl%solver%elsi%pexsiNPole,5) /= 0 .or. ctrl%solver%elsi%pexsiNPole > 40) then
          call detailedError(value1, "Unsupported number of PEXSI poles for this method")
        end if
      end select
      call getChildValue(value1, "ProcsPerPole", ctrl%solver%elsi%pexsiNpPerPole, 1)
      call getChildValue(value1, "muPoints", ctrl%solver%elsi%pexsiNMu, 2)
      call getChildValue(value1, "SymbolicFactorProcs", ctrl%solver%elsi%pexsiNpSymbo, 1)
      call getChildValue(value1, "SpectralRadius", ctrl%solver%elsi%pexsiDeltaE, 10.0_dp,&
          & modifier=modifier, child=child)
      call convertUnitHsd(modifier, energyUnits, child, ctrl%solver%elsi%pexsiDeltaE)

    case ("ntpoly")
      ctrl%solver%isolver = electronicSolverTypes%ntpoly
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      if (ctrl%tSpin) then
        call detailedError(value1, "Solver does not currently support spin polarisation")
      end if
      call getChildValue(value1, "PurificationMethod", ctrl%solver%elsi%ntpolyMethod, 2)
      call getChildValue(value1, "Tolerance", ctrl%solver%elsi%ntpolyTolerance, 1.0E-5_dp)
      call getChildValue(value1, "Truncation", ctrl%solver%elsi%ntpolyTruncation, 1.0E-10_dp)

  #:if WITH_TRANSPORT
    case ("greensfunction")
      ctrl%solver%isolver = electronicSolverTypes%GF
      ! need electronic temperature to be read for this solver:
      call readElectronicFilling(node, ctrl, geo)
      if (tp%defined .and. .not.tp%taskUpload) then
        call detailederror(node, "greensfunction solver cannot be used "// &
            &  "when task = contactHamiltonian")
      end if
      call readGreensFunction(value1, greendens, tp, ctrl%tempElec)
      ! fixEf also avoids checks of total charge later on in the run
      ctrl%tFixEf = .true.
    case ("transportonly")
      if (tp%defined .and. .not.tp%taskUpload) then
        call detailederror(node, "transportonly cannot be used when task = contactHamiltonian")
      end if
      call readGreensFunction(value1, greendens, tp, ctrl%tempElec)
      ctrl%solver%isolver = electronicSolverTypes%OnlyTransport
      ctrl%tFixEf = .true.
  #:endif

    case default
      call detailedError(value1, "Unknown electronic solver")

    end select

    if ((ctrl%solver%isolver == electronicSolverTypes%omm .or.&
        & ctrl%solver%isolver == electronicSolverTypes%pexsi ) .and. .not.ctrl%tSpinSharedEf&
        & .and. ctrl%tSpin .and. .not. ctrl%t2Component) then
      call detailedError(value1, "This solver currently requires spin values to be relaxed")
    end if
    if (ctrl%solver%isolver == electronicSolverTypes%pexsi .and. .not.withPEXSI) then
      call error("Not compiled with PEXSI support via ELSI")
    end if
    if (any(ctrl%solver%isolver == [electronicSolverTypes%elpa, electronicSolverTypes%omm,&
        & electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly])) then
      if (.not.withELSI) then
        call error("Not compiled with ELSI supported solvers")
      end if
    end if

    if (any(ctrl%solver%isolver == [electronicSolverTypes%omm, electronicSolverTypes%pexsi,&
        & electronicSolverTypes%ntpoly])) then
      call getChildValue(value1, "Sparse", ctrl%solver%elsi%elsiCsr, .true.)
      if (.not.ctrl%solver%elsi%elsiCsr) then
        if (any(ctrl%solver%isolver == [electronicSolverTypes%pexsi,electronicSolverTypes%ntpoly]))&
            & then
          call getChildValue(value1, "Threshold", ctrl%solver%elsi%elsi_zero_def, 1.0E-15_dp)
        end if
      end if
    end if

  #:if WITH_TRANSPORT
    if (all(ctrl%solver%isolver /= [electronicSolverTypes%GF,electronicSolverTypes%OnlyTransport])&
        & .and. tp%taskUpload) then
      call detailedError(value1, "Eigensolver incompatible with transport calculation&
          & (GreensFunction or TransportOnly required)")
    end if
  #:endif

  end subroutine readSolver


  !> SCC options that are need for different hamiltonian choices
  subroutine readSccOptions(node, ctrl, geo)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    ctrl%tMulliken = .true.

    call getChildValue(node, "ReadInitialCharges", ctrl%tReadChrg, .false.)
    if (.not. ctrl%tReadChrg) then
      call getInitialCharges(node, geo, ctrl%initialCharges)
    end if

    call getChildValue(node, "SCCTolerance", ctrl%sccTol, 1.0e-5_dp)

    ! temporarily removed until debugged
    ! call getChildValue(node, "WriteShifts", ctrl%tWriteShifts, .false.)
    ctrl%tWriteShifts = .false.

    if (geo%tPeriodic) then
      call getChildValue(node, "EwaldParameter", ctrl%ewaldAlpha, 0.0_dp)
      call getChildValue(node, "EwaldTolerance", ctrl%tolEwald, 1.0e-9_dp)
    end if

    if (geo%tHelical) then
      ! Tolerance for k-points being commensurate with C_n rotation
      call getChildValue(node, "HelicalSymmetryTol", ctrl%helicalSymTol, 1.0E-6_dp)
    end if

    ! self consistency required or not to proceed
    call getChildValue(node, "ConvergentSCCOnly", ctrl%isSccConvRequired, .true.)

  end subroutine readSccOptions


  !> Force evaluation options that are need for different hamiltonian choices
  subroutine readForceOptions(node, ctrl)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    type(hsd_table), pointer :: child
    character(len=:), allocatable :: buffer

    call getChildValue(node, "ForceEvaluation", buffer, "Traditional", child=child)
    select case (tolower(unquote(buffer)))
    case("traditional")
      ctrl%forceType = forceTypes%orig
    case("dynamicst0")
      ctrl%forceType = forceTypes%dynamicT0
    case("dynamics")
      ctrl%forceType = forceTypes%dynamicTFinite
    case default
      call detailedError(child, "Invalid force evaluation method.")
    end select

  end subroutine readForceOptions


  !> Options for truncation of the SK data sets at a fixed distance
  subroutine SKTruncations(node, truncationCutOff, skInterMeth)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> This is the resulting cutoff distance
    real(dp), intent(out) :: truncationCutOff

    !> Method of the sk interpolation
    integer, intent(in) :: skInterMeth

    logical :: tHardCutOff
    type(hsd_table), pointer :: field
    character(len=:), allocatable :: modifier

    ! Artificially truncate the SK table
    call getChildValue(node, "SKMaxDistance", truncationCutOff, modifier=modifier, child=field)
    call convertUnitHsd(modifier, lengthUnits, field, truncationCutOff)

    call getChildValue(node, "HardCutOff", tHardCutOff, .true.)
    if (tHardCutOff) then
      ! Adjust by the length of the tail appended to the cutoff
      select case(skInterMeth)
      case(skEqGridOld)
        truncationCutOff = truncationCutOff - distFudgeOld
      case(skEqGridNew)
        truncationCutOff = truncationCutOff - distFudge
      end select
    end if
    if (truncationCutOff < epsilon(0.0_dp)) then
      call detailedError(field, "Truncation is shorter than the minimum distance over which SK data&
          & goes to 0")
    end if

  end subroutine SKTruncations


  !> Reads initial charges
  subroutine getInitialCharges(node, geo, initCharges)

    !> relevant node in input tree
    type(hsd_table), pointer :: node

    !> geometry, including atomic type information
    type(TGeometry), intent(in) :: geo

    !> initial atomic charges
    real(dp), allocatable :: initCharges(:)

    type(hsd_table), pointer :: child, child2, child3, val
    type(hsd_child_list), pointer :: children
    integer, allocatable :: pTmpI1(:)
    character(len=:), allocatable :: buffer
    real(dp) :: rTmp
    integer :: ii, jj, iAt

    call getChildValue(node, "InitialCharges", val, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! Read either all atom charges, or individual atom specifications
    call getChild(child, "AllAtomCharges", child2, requested=.false.)
    if (associated(child2)) then
      allocate(initCharges(geo%nAtom))
      call getChildValue(child2, "", initCharges)
    else
      call getChildren(child, "AtomCharge", children)
      if (getLength(children) > 0) then
        allocate(initCharges(geo%nAtom))
        initCharges = 0.0_dp
      end if
      do ii = 1, getLength(children)
        call getItem1(children, ii, child2)
        call getChildValue(child2, "Atoms", buffer, child=child3, multiple=.true.)
        call getSelectedAtomIndices(child3, buffer, geo%speciesNames, geo%species, pTmpI1)
        call getChildValue(child2, "ChargePerAtom", rTmp)
        do jj = 1, size(pTmpI1)
          iAt = pTmpI1(jj)
          if (initCharges(iAt) /= 0.0_dp) then
            call detailedWarning(child3, "Previous setting for the charge &
                &of atom" // i2c(iAt) // " overwritten")
          end if
          initCharges(iAt) = rTmp
        end do
        deallocate(pTmpI1)
      end do
      call destroyNodeList(children)
    end if

  end subroutine getInitialCharges


  !> Reads initial spins
  subroutine getInitialSpins(node, geo, nSpin, initSpins)

    !> relevant node in input data
    type(hsd_table), pointer :: node

    !> geometry, including atomic information
    type(TGeometry), intent(in) :: geo

    !> number of spin channels
    integer, intent(in) :: nSpin

    !> initial spins on return
    real(dp), allocatable :: initSpins(:,:)

    type(hsd_table), pointer :: child, child2, child3, val
    type(hsd_child_list), pointer :: children
    integer, allocatable :: pTmpI1(:)
    character(len=:), allocatable :: buffer
    real(dp), allocatable :: rTmp(:)
    integer :: ii, jj, iAt

    @:ASSERT(nSpin == 1 .or. nSpin == 3)

    call getChildValue(node, "InitialSpins", val, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! Read either all atom spins, or individual spin specifications
    call getChild(child, "AllAtomSpins", child2, requested=.false.)
    if (associated(child2)) then
      allocate(initSpins(nSpin, geo%nAtom))
      call getChildValue(child2, "", initSpins)
    else
      call getChildren(child, "AtomSpin", children)
      if (getLength(children) > 0) then
        allocate(initSpins(nSpin, geo%nAtom))
        initSpins = 0.0_dp
      end if
      allocate(rTmp(nSpin))
      do ii = 1, getLength(children)
        call getItem1(children, ii, child2)
        call getChildValue(child2, "Atoms", buffer, child=child3, multiple=.true.)
        call getSelectedAtomIndices(child3, buffer, geo%speciesNames, geo%species, pTmpI1)
        call getChildValue(child2, "SpinPerAtom", rTmp)
        do jj = 1, size(pTmpI1)
          iAt = pTmpI1(jj)
          if (any(initSpins(:,iAt) /= 0.0_dp)) then
            call detailedWarning(child3, "Previous setting for the spin of atom" // i2c(iAt) //&
                & " overwritten")
          end if
          initSpins(:,iAt) = rTmp
        end do
        deallocate(pTmpI1)
      end do
      deallocate(rTmp)
      call destroyNodeList(children)
    end if

  end subroutine getInitialSpins


  !> Reads numerical differentiation method to be used
  subroutine readDifferentiation(node, ctrl)

    !> relevant node in input tree
    type(hsd_table), pointer, intent(in) :: node

    !> control structure to fill
    type(TControl), intent(inout) :: ctrl


    !> default of a reasonable choice for round off when using a second order finite difference
    !> formula
    real(dp), parameter :: defDelta = epsilon(1.0_dp)**0.25_dp

    character(len=:), allocatable :: buffer, modifier
    type(hsd_table), pointer :: val, child

    call getChildValue(node, "Differentiation", val, "FiniteDiff",&
        & child=child)
    call getNodeName(val, buffer)
    select case (buffer)
    case ("finitediff")
      ctrl%iDerivMethod = diffTypes%finiteDiff
      call getChildValue(val, "Delta", ctrl%deriv1stDelta, defDelta,&
          & modifier=modifier, child=child)
      call convertUnitHsd(modifier, lengthUnits, child,&
          & ctrl%deriv1stDelta)
    case ("richardson")
      ctrl%iDerivMethod = diffTypes%richardson
    case default
      call getNodeHSDName(val, buffer)
      call detailedError(child, "Invalid derivative calculation '" &
          & // buffer // "'")
    end select

  end subroutine readDifferentiation


  !> Reads the H corrections (H5, Damp)
  subroutine readHCorrection(node, geo, ctrl)

    !> Node containing the h-bond correction sub-block.
    type(hsd_table), pointer, intent(in) :: node

    !> Geometry.
    type(TGeometry), intent(in) :: geo

    !> Control structure
    type(TControl), intent(inout) :: ctrl

    type(hsd_table), pointer :: value1, child, child2
    character(len=:), allocatable :: buffer
    real(dp) :: h5ScalingDef
    integer :: iSp

    ! X-H interaction corrections including H5 and damping
    ctrl%tDampH = .false.
    call getChildValue(node, "HCorrection", value1, "None", child=child)
    call getNodeName(value1, buffer)

    select case (buffer)

    case ("none")
      ! nothing to do

    case ("damping")
      ! Switch the correction on
      ctrl%tDampH = .true.
      call getChildValue(value1, "Exponent", ctrl%dampExp)

    case ("h5")
      allocate(ctrl%h5Input)
      associate (h5Input => ctrl%h5Input)
        call getChildValue(value1, "RScaling", h5Input%rScale, 0.714_dp)
        call getChildValue(value1, "WScaling", h5Input%wScale, 0.25_dp)
        allocate(h5Input%elementParams(geo%nSpecies))
        call getChild(value1, "H5Scaling", child2, requested=.false., emptyIfMissing=.true.)
        do iSp = 1, geo%nSpecies
          select case (geo%speciesNames(iSp))
          case ("O")
            h5ScalingDef = 0.06_dp
          case ("N")
            h5ScalingDef = 0.18_dp
          case ("S")
            h5ScalingDef = 0.21_dp
          case default
            ! Default value is -1, this indicates that the element should be ignored
            h5ScalingDef = -1.0_dp
          end select
          call getChildValue(child2, geo%speciesNames(iSp), h5Input%elementParams(iSp),&
              & h5ScalingDef)
        end do
        h5Input%speciesNames = geo%speciesNames
      end associate

    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Invalid HCorrection '" // buffer // "'")
    end select

  end subroutine readHCorrection



  !> Reads the option block
  subroutine readOptions(node, ctrl, geom)

    !> Node to parse
    type(hsd_table), pointer :: node

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure
    type(TGeometry), intent(in) :: geom

    type(hsd_table), pointer :: child, value1
    logical :: tWriteDetailedOutDef

  #:if WITH_SOCKETS
    tWriteDetailedOutDef = .not. allocated(ctrl%socketInput)
  #:else
    tWriteDetailedOutDef = .true.
  #:endif
    call getChildValue(node, "WriteDetailedOut", ctrl%tWriteDetailedOut, tWriteDetailedOutDef)

    call getChildValue(node, "WriteAutotestTag", ctrl%tWriteTagged, .false.)
    call getChildValue(node, "WriteDetailedXML", ctrl%tWriteDetailedXML, .false.)
    block
      character(len=:), allocatable :: tmpFmt
      call getChildValue(node, "DetailedOutputFormat", tmpFmt, "xml")
      ctrl%detailedOutputFormat = tmpFmt
    end block
    call getChildValue(node, "WriteResultsTag", ctrl%tWriteResultsTag, .false.)

    if (.not.(ctrl%tMD.or.ctrl%isGeoOpt.or.allocated(ctrl%geoOpt))) then
      if (ctrl%tSCC) then
        call getChildValue(node, "RestartFrequency", ctrl%restartFreq, 20)
      else
        ctrl%restartFreq = 0
      end if
    end if
    if (ctrl%tMD) then
      allocate(ctrl%mdOutput)
      call getChildValue(node, "MDOutput", value1, "", child=child, allowEmptyValue=.true.,&
          & dummyValue=.true.)
      if (associated(value1)) then
        if (providesEigenvalues(ctrl%solver%isolver)) then
          call getChildValue(child, "AppendBandOut", ctrl%mdOutput%bandStructure, .false.)
        end if
        if (ctrl%tPrintForces) then
          call getChildValue(child, "Derivatives", ctrl%mdOutput%printForces, .true.)
        end if
        if (ctrl%tPrintMulliken) then
          call getChildValue(child, "Charges", ctrl%mdOutput%printCharges, .true.)
        end if
        if (ctrl%tAtomicEnergy) then
          call getChildValue(child, "AtomEnergies", ctrl%mdOutput%printAtomEnergies, .true.)
        end if
      end if
    end if
    call getChildValue(node, "RandomSeed", ctrl%iSeed, 0, child=child)
    if (ctrl%iSeed < 0) then
      call detailedError(child, "Random seed must be greater or equal zero")
    end if
    call getChildValue(node, "WriteHS", ctrl%tWriteHS, .false.)
    call getChildValue(node, "WriteRealHS", ctrl%tWriteRealHS, .false.)
    call hsd_rename_child(node, "MinimizeMemoryUsage", "MinimiseMemoryUsage")
    call getChildValue(node, "MinimiseMemoryUsage", ctrl%tMinMemory, .false., child=child)
    if (ctrl%tMinMemory) then
      call detailedWarning(child, "Memory minimisation is not working currently, normal calculation&
          & will be used instead")
    end if
    call getChildValue(node, "ShowFoldedCoords", ctrl%tShowFoldedCoord, .false.)
  #:if DEBUG > 0
    call getChildValue(node, "TimingVerbosity", ctrl%timingLevel, -1)
  #:else
    call getChildValue(node, "TimingVerbosity", ctrl%timingLevel, 1)
  #:endif

    if (ctrl%tReadChrg) then
      call getChildValue(node, "ReadChargesAsText", ctrl%tReadChrgAscii, .false.)
    end if

    call getChildValue(node, "WriteCharges", ctrl%tWriteCharges, .true.)
    if (ctrl%tWriteCharges) then
      call getChildValue(node, "WriteChargesAsText", ctrl%tWriteChrgAscii, .false.)
    end if

    ctrl%tSkipChrgChecksum = .false.
    if (.not. ctrl%tFixEf .and. ctrl%tReadChrg) then
      call getChildValue(node, "SkipChargeTest", ctrl%tSkipChrgChecksum, .false.)
    end if

    if (geom%tPeriodic .or. geom%tHelical .or. geom%areContactsPresent) then
      call getChildValue(node, "WriteAllAtomGeometry", ctrl%areAllAtomsPrinted, .false.)
    end if

    call readBinaryAccessTypes(node, ctrl%binaryAccessTypes)

  end subroutine readOptions


  !> Read in hamiltonian settings that are influenced by those read from REKS{}, electronDynamics{}
  subroutine readLaterHamiltonian(hamNode, ctrl, driverNode, geo)

    !> Hamiltonian node to parse
    type(hsd_table), pointer :: hamNode

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> Geometry driver node to parse
    type(hsd_table), pointer :: driverNode

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    type(hsd_table), pointer :: value1, value2, child, child2
    character(len=:), allocatable :: buffer, buffer2
    type(TListRealR1) :: lr1

    if (ctrl%reksInp%reksAlg == reksTypes%noReks) then

      if (ctrl%tSCC) then

        call getChildValue(hamNode, "Mixer", value1, "Broyden", child=child)
        call getNodeName(value1, buffer)
        select case(buffer)

        case ("broyden")

          allocate(ctrl%mixerInp%broydenMixerInp)
          associate (inp => ctrl%mixerInp%broydenMixerInp)
            call getChildValue(value1, "MixingParameter", inp%mixParam, 0.2_dp)
            call getChildValue(value1, "InverseJacobiWeight", inp%omega0, 0.01_dp)
            call getChildValue(value1, "MinimalWeight", inp%minWeight, 1.0_dp)
            call getChildValue(value1, "MaximalWeight", inp%maxWeight, 1.0e5_dp)
            call getChildValue(value1, "WeightFactor", inp%weightFac, 1.0e-2_dp)
          end associate

        case ("anderson")

          allocate(ctrl%mixerInp%andersonMixerInp)
          associate (inp => ctrl%mixerInp%andersonMixerInp)
            call getChildValue(value1, "MixingParameter", inp%mixParam, 0.05_dp)
            call getChildValue(value1, "Generations", inp%iGenerations, 4)
            call getChildValue(value1, "InitMixingParameter", inp%initMixParam, 0.01_dp)
            call getChildValue(value1, "DynMixingParameters", value2, "", child=child,&
                & allowEmptyValue=.true.)
            call getNodeName2(value2, buffer2)
            if (buffer2 == "" .and. .not. hasInlineData(child)) then
              inp%nConvMixParam = 0
            else
              call init(lr1)
              call getChildValue(child, "", 2, lr1, child=child2)
              if (len(lr1) < 1) then
                call detailedError(child2, "At least one dynamic mixing parameter must be defined.")
              end if
              inp%nConvMixParam = len(lr1)
              allocate(inp%convMixParam(2, inp%nConvMixParam))
              call asArray(lr1, inp%convMixParam)
              call destruct(lr1)
            end if
            call getChildValue(value1, "DiagonalRescaling", inp%omega0, 1.0e-2_dp)
          end associate

        case ("simple")

          allocate(ctrl%mixerInp%simpleMixerInp)
          associate (inp => ctrl%mixerInp%simpleMixerInp)
            call getChildValue(value1, "MixingParameter", inp%mixParam, 0.05_dp)
          end associate

        case ("diis")

          allocate(ctrl%mixerInp%diisMixerInp)
          associate (inp => ctrl%mixerInp%diisMixerInp)
            call getChildValue(value1, "InitMixingParameter", inp%initMixParam, 0.2_dp)
            call getChildValue(value1, "Generations", inp%iGenerations, 6)
            call getChildValue(value1, "UseFromStart", inp%tFromStart, .true.)
          end associate

        case default

          call getNodeHSDName(value1, buffer)
          call detailedError(child, "Invalid mixer '" // buffer // "'")

        end select

      end if

      if (ctrl%tMD) then
        if (ctrl%thermostatInp%thermostatType /= thermostatTypes%dummy) then
          call getChildValue(driverNode, "Thermostat", child, child=child2)
          if (ctrl%reksInp%reksAlg == reksTypes%noReks) then
            call getChildValue(child, "AdaptFillingTemp", ctrl%tSetFillingTemp, .false.)
          end if
        end if
      end if

    end if

    hamNeedsT: if (ctrl%reksInp%reksAlg == reksTypes%noReks) then

      if (allocated(ctrl%elecDynInp)) then
        if (ctrl%elecDynInp%tReadRestart .and. .not.ctrl%elecDynInp%tPopulations) then
          exit hamNeedsT
        end if
      end if

      if (ctrl%solver%isolver /= electronicSolverTypes%GF) then
        call readElectronicFilling(hamNode, ctrl, geo)
      end if

    end if hamNeedsT

  end subroutine readLaterHamiltonian


  !> Parses for electronic filling temperature (should only read if not either REKS or electron
  !> dynamics from a supplied density matrix)
  subroutine readElectronicFilling(hamNode, ctrl, geo)

    !> Relevant node in input tree
    type(hsd_table), pointer :: hamNode

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to test for periodicity
    type(TGeometry), intent(in) :: geo

    select case(ctrl%hamiltonian)
    case(hamiltonianTypes%xtb)
      call readFilling(hamNode, ctrl, geo, 300.0_dp*Boltzmann)
    case(hamiltonianTypes%dftb)
      call readFilling(hamNode, ctrl, geo, 0.0_dp)
    end select

  end subroutine readElectronicFilling


  !> Reads W values if required by settings in the Hamiltonian or the excited state
  subroutine readSpinConstants(hamNode, geo, orb, ctrl)

    !> node for Hamiltonian data
    type(hsd_table), pointer :: hamNode

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> control structure
    type(TControl), intent(inout) :: ctrl

    type(hsd_table), pointer :: child
    logical :: tLRNeedsSpinConstants, tShellResolvedW
    integer :: iSp1, nConstants
    type(TListReal) :: realBuffer
    character(lc) :: strTmp
    real(dp) :: rWork(maxval(orb%nShell)**2)

    tLRNeedsSpinConstants = .false.

    if (allocated(ctrl%lrespini)) then
      select case (ctrl%lrespini%sym)
      case ("T", "B", " ")
        tLRNeedsSpinConstants = .true.
      case ("S")
        tLRNeedsSpinConstants = .false.
      case default
      end select
    end if

    if (tLRNeedsSpinConstants .or. ctrl%tSpin .or. &
        & ctrl%reksInp%reksAlg /= reksTypes%noReks) then
      allocate(ctrl%spinW(orb%mShell, orb%mShell, geo%nSpecies))
      ctrl%spinW(:,:,:) = 0.0_dp

      call getChild(hamNode, "SpinConstants", child)
      if (ctrl%hamiltonian == hamiltonianTypes%xtb) then
        call getChildValue(child, "ShellResolvedSpin", tShellResolvedW, .true.)
      else
        if (.not.ctrl%tShellResolved) then
          call getChildValue(child, "ShellResolvedSpin", tShellResolvedW, .false.)
        else
          tShellResolvedW = .true.
        end if
      end if

      do iSp1 = 1, geo%nSpecies
        call init(realBuffer)
        call getChildValue(child, geo%speciesNames(iSp1), realBuffer)
        nConstants = len(realBuffer)
        if (tShellResolvedW) then
          if (nConstants == orb%nShell(iSp1)**2) then
            call asArray(realBuffer, rWork(:orb%nShell(iSp1)**2))
            ctrl%spinW(:orb%nShell(iSp1), :orb%nShell(iSp1), iSp1) =&
                & reshape(rWork(:orb%nShell(iSp1)**2), [orb%nShell(iSp1), orb%nShell(iSp1)])
          else
            write(strTmp, "(A,I0,A,I0,A,A,A)")'Expecting a ', orb%nShell(iSp1), ' x ',&
                & orb%nShell(iSp1), ' spin constant matrix for "', trim(geo%speciesNames(iSp1)),&
                & '", as ShellResolvedSpin enabled.'
            call detailedError(child, trim(strTmp))
          end if
        else
          if (nConstants == 1) then
            call asArray(realBuffer, rWork(:1))
            ! only one value for all atom spin constants
            ctrl%spinW(:orb%nShell(iSp1), :orb%nShell(iSp1), iSp1) = rWork(1)
          else
            write(strTmp, "(A,A,A)")'Expecting a single spin constant for "',&
                & trim(geo%speciesNames(iSp1)),'", as ShellResolvedSpin not enabled.'
            call detailedError(child, trim(strTmp))
          end if
        end if
        call destruct(realBuffer)
      end do
    end if

  end subroutine readSpinConstants


  !> Reads customised Hubbard U values that over-ride the SK file values
  subroutine readCustomisedHubbards(node, geo, orb, tShellResolvedScc, hubbU)

    !> input data to parse
    type(hsd_table), pointer, intent(in) :: node

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> is this a shell resolved calculation, or only one U value per atom
    logical, intent(in) :: tShellResolvedScc

    !> hubbard U values on exit
    real(dp), allocatable, intent(out) :: hubbU(:,:)

    type(hsd_table), pointer :: child, child2
    integer :: iSp1

    call hsd_rename_child(node, "CustomizedHubbards", "CustomisedHubbards")
    call getChild(node, "CustomisedHubbards", child, requested=.false.)
    if (associated(child)) then
      allocate(hubbU(orb%mShell, geo%nSpecies))
      hubbU(:,:) = 0.0_dp
      do iSp1 = 1, geo%nSpecies
        call getChild(child, geo%speciesNames(iSp1), child2, requested=.false.)
        if (.not. associated(child2)) then
          cycle
        end if
        if (tShellResolvedScc) then
          call getChildValue(child2, "", hubbU(:orb%nShell(iSp1), iSp1))
        else
          call getChildValue(child2, "", hubbU(1, iSp1))
          hubbU(:orb%nShell(iSp1), iSp1) = hubbU(1, iSp1)
        end if
      end do
    end if

  end subroutine readCustomisedHubbards


  !> This subroutine overrides the neutral (reference) atom electronic occupation
  subroutine readCustomReferenceOcc(root, orb, referenceOcc, geo, iAtInRegion, customOcc)

    !> Node to be parsed
    type(hsd_table), pointer, intent(in) :: root

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Default reference occupations
    real(dp), intent(in) :: referenceOcc(:,:)

    !> Geometry information
    type(TGeometry), intent(in) :: geo

    !> Atom indices corresponding to user defined reference atomic charges
    type(TWrappedInt1), allocatable, intent(out) :: iAtInRegion(:)

    !> User-defined reference atomic charges
    real(dp), allocatable, intent(out) :: customOcc(:,:)

    type(hsd_table), pointer :: node, container, child
    type(hsd_child_list), pointer :: nodes
    character(len=:), allocatable :: buffer
    integer :: nCustomOcc, iCustomOcc, iShell, iSpecies, nAtom
    character(sc), allocatable :: shellNamesTmp(:)
    logical, allocatable :: atomOverriden(:)

    call hsd_rename_child(root, "CustomizedOccupations", "CustomisedOccupations")
    call getChild(root, "CustomisedOccupations", container, requested=.false.)
    if (.not. associated(container)) then
      return
    end if

    call getChildren(container, "ReferenceOccupation", nodes)
    nCustomOcc = getLength(nodes)
    nAtom = size(geo%species)
    allocate(iAtInRegion(nCustomOcc))
    allocate(customOcc(orb%mShell, nCustomOcc))
    allocate(atomOverriden(nAtom))
    atomOverriden(:) = .false.
    customOcc(:,:) = 0.0_dp

    do iCustomOcc = 1, nCustomOcc
      call getItem1(nodes, iCustomOcc, node)
      call getChildValue(node, "Atoms", buffer, child=child, multiple=.true.)
      call getSelectedAtomIndices(child, buffer, geo%speciesNames, geo%species,&
          & iAtInRegion(iCustomOcc)%data)
      if (any(atomOverriden(iAtInRegion(iCustomOcc)%data))) then
        call detailedError(child, "Atom region contains atom(s) which have already been overridden")
      end if
      atomOverriden(iAtInRegion(iCustomOcc)%data) = .true.
      iSpecies = geo%species(iAtInRegion(iCustomOcc)%data(1))
      if (any(geo%species(iAtInRegion(iCustomOcc)%data) /= iSpecies)) then
        call detailedError(child, "All atoms in a ReferenceOccupation declaration must have the&
            & same type.")
      end if
      call getShellNames(iSpecies, orb, shellNamesTmp)
      do iShell = 1, orb%nShell(iSpecies)
          call getChildValue(node, shellNamesTmp(iShell), customOcc(iShell, iCustomOcc), &
            & default=referenceOcc(iShell, iSpecies))
      end do
      deallocate(shellNamesTmp)
    end do
    call destroyNodeList(nodes)

  end subroutine readCustomReferenceOcc
  function is_numeric(string) result(is)
    character(len=*), intent(in) :: string
    logical :: is

    real :: x
    integer :: err

    read(string,*,iostat=err) x
    is = (err == 0)
  end function is_numeric


  !> Returns parser version for a given input version or throws an error if not possible.
  function parserVersionFromInputVersion(versionString) result(parserVersion)

    !> Input version string
    character(len=*), intent(in) :: versionString

    !> Corresponding parser version, or 0 if no corresponding parser version was found
    integer :: parserVersion

    integer :: ii

    parserVersion = 0
    do ii = 1, size(versionMaps)
      if (versionMaps(ii)%inputVersion == versionString) then
        parserVersion = versionMaps(ii)%parserVersion
        return
      end if
    end do

  end function parserVersionFromInputVersion

end module dftbp_dftbplus_parser
