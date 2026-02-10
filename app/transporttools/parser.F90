!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Fills the derived type with the input parameters from an HSD or an XML file.
module transporttools_parser
  use transporttools_helpsetupgeom, only : setupGeometry
  use transporttools_inputdata, only : TInputData
  use dftbp_common_accuracy, only : distFudge, distFudgeOld, dp, lc, mc
  use dftbp_common_constants, only : Bohr__AA
  use dftbp_common_globalenv, only : stdOut, tIoProc
  use dftbp_common_unitconversion, only : lengthUnits
  use dftbp_dftb_slakoeqgrid, only : skEqGridNew, skEqGridOld
  use dftbp_dftbplus_oldcompat, only : convertOldHsd
  use dftbp_io_charmanip, only : i2c, newline, tolower, unquote
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, getSelectedAtomIndices
  use dftbp_io_unitconv, only : convertUnitHsd
  use hsd, only : hsd_warn_unprocessed, MAX_WARNING_LEN, hsd_error_t, hsd_dump,&
      & hsd_table_ptr, hsd_get_child_tables, hsd_get, hsd_get_or_set, hsd_get_table,&
      & hsd_get_choice, hsd_set, hsd_get_attrib, HSD_STAT_OK, &
      & hsd_get_name
  use hsd_data, only : hsd_table, data_load, DATA_FMT_AUTO, new_table
  use dftbp_io_message, only : error, warning
  use dftbp_transport_negfvars, only : ContactInfo, TTransPar
  use dftbp_type_oldskdata, only : readFromFile, TOldSKData
  use dftbp_type_typegeometryhsd, only : readTGeometryGen, readTGeometryHsd, TGeometry
  use dftbp_type_wrappedintr, only : TWrappedInt1
  implicit none

  private
  public :: parseHsdInput, parserVersion


  ! Default file names

  !> Main HSD input file
  character(len=*), parameter :: hsdInputName = "setup_in.hsd"

  !> XML input file
  character(len=*), parameter :: xmlInputName = "setup_in.xml"

  !> Processed HSD input
  character(len=*), parameter :: hsdProcInputName = "setup_pin.hsd"

  !> Processed  XML input
  character(len=*), parameter :: xmlProcInputName = "setup_pin.xml"

  !> Tag at the head of the input document tree
  character(len=*), parameter :: rootTag = "setup_in"

  !> Wrapper type for an array of character(lc) strings
  type :: TCharLcArray
    character(lc), allocatable :: items(:)
  end type TCharLcArray


  !> Version of the current parser
  integer, parameter :: parserVersion = 6


  !> Version of the oldest parser for which compatibility is still maintained
  integer, parameter :: minVersion = 1


  !> Container type for parser related flags.
  type TParserFlags

    !> stop after parsing?
    logical :: tStop

    !> Continue despite unprocessed nodes
    logical :: tIgnoreUnprocessed

    !> XML output?
    logical :: tWriteXML

    !> HSD output?
    logical :: tWriteHSD
  end type TParserFlags


contains


  !> Parse input from an HSD/XML file
  subroutine parseHsdInput(input)

    !> Returns initialised input variables on exit
    type(TInputData), intent(out) :: input

    type(hsd_table), pointer :: hsdTree
    type(hsd_table), pointer :: root, tmp, child
    type(TParserflags) :: parserFlags
    integer :: stat

    write(stdOut, "(/, A, /)") "***  Parsing and initializing"

    ! Read in the input
    block
      type(hsd_error_t), allocatable :: hsdError
      type(hsd_table), pointer :: content
      allocate(content)
      call data_load(hsdInputName, content, hsdError, fmt=DATA_FMT_AUTO)
      if (allocated(hsdError)) then
        call error("Error loading input file '" // trim(hsdInputName) // "': " &
            & // trim(hsdError%message))
      end if
      content%name = rootTag
      allocate(hsdTree)
      call new_table(hsdTree, name="document")
      call hsdTree%add_child(content)
    end block
    call hsd_get_table(hsdTree, rootTag, root, stat, auto_wrap=.true.)
    if (.not. associated(root)) call dftbp_error(hsdTree, "Missing required block: '" // rootTag // "'")

    write(stdout, '(A,1X,I0,/)') 'Parser version:', parserVersion
    write(stdout, "(A)") "Interpreting input file '" // hsdInputName // "'"
    write(stdout, "(A)") repeat("-", 80)

    ! Handle parser options
    call hsd_get_table(root, "ParserOptions", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      call hsd_set(root, "ParserOptions", "")
      call hsd_get_table(root, "ParserOptions", child, stat, auto_wrap=.true.)
    end if
    call readParserOptions(child, root, parserFlags)

    ! Read in the different blocks

    ! Atomic geometry and boundary conditions
    call hsd_get_table(root, "Geometry", tmp, stat, auto_wrap=.true.)
    if (.not. associated(tmp)) call dftbp_error(root, "Missing required block: 'Geometry'")
    call readGeometry(tmp, input)

    call hsd_get_table(root, "Transport", child, stat, auto_wrap=.true.)

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
    ! input data strucutre has been initialised
    input%tInitialized = .true.

    ! Issue warning about unprocessed nodes
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

    ! Dump processed tree in HSD and XML format
    if (tIoProc .and. parserFlags%tWriteHSD) then
      block
        type(hsd_error_t), allocatable :: hsdErr
        call hsd_dump(hsdTree, hsdProcInputName, hsdErr)
        if (allocated(hsdErr)) call error("Error writing HSD file '" // hsdProcInputName // "'")
      end block
      write(stdout, '(/,/,A)') "Processed input in HSD format written to '" &
          &// hsdProcInputName // "'"
    end if

    ! Stop, if only parsing is required
    if (parserFlags%tStop) then
      call error("Keyword 'StopAfterParsing' is set to Yes. Stopping.")
    end if

    hsdTree => null()

    write(stdout,*) 'Geometry processed. Job finished'

  end subroutine parseHsdInput


  !> Read in parser options (options not passed to the main code)
  subroutine readParserOptions(node, root, flags)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Root of the entire tree (in case it needs to be converted, for example because dftbp_of
    !> compatibility options)
    type(hsd_table), pointer :: root

    !> Contains parser flags on exit.
    type(TParserFlags), intent(out) :: flags

    integer :: inputVersion
    integer :: stat
    type(hsd_table), pointer :: child

    ! Check if input needs compatibility conversion.
    call hsd_get_or_set(node, "ParserVersion", inputVersion, parserVersion, &
        &child=child)
    if (inputVersion < 1 .or. inputVersion > parserVersion) then
      call dftbp_error(child, "Invalid parser version (" // i2c(inputVersion)&
          &// ")")
    elseif (inputVersion < minVersion) then
      call dftbp_error(child, &
          &"Sorry, no compatibility mode for parser version " &
          &// i2c(inputVersion) // " (too old)")
    elseif (inputVersion /= parserVersion) then
      write(stdout, "(A,I2,A,I2,A)") "***  Converting input from version ", &
          &inputVersion, " to version ", parserVersion, " ..."
      call convertOldHSD(root, inputVersion, parserVersion)
      write(stdout, "(A,/)") "***  Done."
    end if

    call hsd_get_or_set(node, "WriteHSDInput", flags%tWriteHSD, .true.)
    call hsd_get_or_set(node, "WriteXMLInput", flags%tWriteXML, .false.)
    if (.not. (flags%tWriteHSD .or. flags%tWriteXML)) then
      call dftbp_warning(node, &
          &"WriteHSDInput and WriteXMLInput both turned off. You are not&
          & guaranteed" &
          &// newline // &
          &" to able to obtain the same results with a later version of the&
          & code!")
    end if
    call hsd_get_or_set(node, "StopAfterParsing", flags%tStop, .false.)

    call hsd_get_or_set(node, "IgnoreUnprocessedNodes", &
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
    integer :: stat

    call hsd_get_choice(node, "", buffer, value1, stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value in Geometry block")
    child => node
    select case (buffer)
    case ("genformat")
      call readTGeometryGen(value1, input%geom)
    case default
      value1%processed = .false.
      call readTGeometryHSD(child, input%geom)
    end select

  end subroutine readGeometry


  !> Read geometry information for transport calculation
  subroutine readTransportGeometry(root, geom, transpar)

    !> Root node containing the current block
    type(hsd_table), pointer :: root

    !> geometry of the system, which may be modified for some types of calculation
    type(TGeometry), intent(inout) :: geom

    !> Parameters of the transport calculation
    type(TTransPar), intent(inout) :: transpar

    type(hsd_table), pointer :: pDevice, pTask, pTaskType
    character(len=:), allocatable :: buffer
    type(hsd_table_ptr), allocatable :: pNodeList(:)
    real(dp) :: skCutoff
    type(TWrappedInt1), allocatable :: iAtInRegion(:)
    integer, allocatable :: nPLs(:)
    logical :: printDebug
    integer :: stat

    transpar%defined = .true.
    transpar%tPeriodic1D = .not. geom%tPeriodic

    !! Note: we parse first the task because dftbp_we need to know it to define the
    !! mandatory contact entries. On the other hand we need to wait that
    !! contacts are parsed to resolve the name of the contact for task =
    !! contacthamiltonian
    call hsd_get_table(root, "Task", pTaskType, stat, auto_wrap=.true.)
    if (associated(pTaskType)) then
      call hsd_get_choice(root, "Task", buffer, pTask, stat)
      if (.not. associated(pTask)) buffer = 'uploadcontacts'
    else
      buffer = 'uploadcontacts'
      pTask => null()
      pTaskType => root
    end if

    if (buffer.ne."setupgeometry") then
      call hsd_get_table(root, "Device", pDevice, stat, auto_wrap=.true.)
      if (.not. associated(pDevice)) call dftbp_error(root, "Missing required block: 'Device'")
      block
        integer, allocatable :: tmpRange(:)
        call hsd_get(pDevice, "AtomRange", tmpRange, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(pDevice, "Missing required value: 'AtomRange'")
        transpar%idxdevice = tmpRange
      end block
    end if

    call hsd_get_child_tables(root, "Contact", pNodeList)
    transpar%ncont = size(pNodeList)
    if (transpar%ncont < 2) then
      call dftbp_error(root, "At least two contacts must be defined")
    end if
    allocate(transpar%contacts(transpar%ncont))

    select case (buffer)
    case ("setupgeometry")

      call readContacts(pNodeList, transpar%contacts, geom, buffer, iAtInRegion, nPLs)
      call getSKcutoff(pTask, geom, skCutoff)
      write(stdOut,*) 'Maximum SK cutoff:', SKcutoff*Bohr__AA,'(A)'
      call hsd_get_or_set(pTask, "printInfo", printDebug, .false.)
      call setupGeometry(geom, iAtInRegion, transpar%contacts, skCutoff, nPLs, printDebug)

    case default

      call hsd_get_name(pTask, buffer)
      call dftbp_error(pTaskType, "Invalid task '" // buffer // "'")

   end select


  end subroutine readTransportGeometry


  !> Read bias information, used in Analysis and Green's function eigensolver
  subroutine readContacts(pNodeList, contacts, geom, task, iAtInRegion, nPLs)
    type(hsd_table_ptr), intent(in) :: pNodeList(:)
    type(ContactInfo), allocatable, dimension(:), intent(inout) :: contacts
    type(TGeometry), intent(in) :: geom
    character(*), intent(in) :: task
    type(TWrappedInt1), allocatable, intent(out) :: iAtInRegion(:)
    integer, intent(out), allocatable :: nPLs(:)

    real(dp) :: contactLayerTol, vec(3)
    integer :: selectionRange(2)
    integer :: ii, ishift, stat
    type(hsd_table), pointer :: field, pNode, pTmp
    character(len=:), allocatable :: buffer, modif

    allocate(iAtInRegion(size(contacts)+1))
    allocate(nPLs(size(contacts)))

    do ii = 1, size(contacts)

      contacts(ii)%wideBand = .false.
      contacts(ii)%wideBandDos = 0.0_dp

      pNode => pNodeList(ii)%ptr
      call hsd_get(pNode, "Id", buffer, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "Missing required value: 'Id'")
      call hsd_get_table(pNode, "Id", pTmp, stat, auto_wrap=.true.)
      buffer = tolower(trim(unquote(buffer)))
      if (len(buffer) > mc) then
        call dftbp_error(pTmp, "Contact id may not be longer than " // i2c(mc) // " characters.")
      end if
      contacts(ii)%name = buffer
      if (any(contacts(1:ii-1)%name == contacts(ii)%name)) then
        call dftbp_error(pTmp, "Contact id '" // trim(contacts(ii)%name) &
            &//  "' already in use")
      end if

      call hsd_get_or_set(pNode, "PLShiftTolerance", contactLayerTol, 1e-5_dp, child=field)
      call hsd_get_attrib(pNode, "PLShiftTolerance", modif, stat)
      if (stat /= HSD_STAT_OK) modif = ""
      call convertUnitHsd(modif, lengthUnits, field, contactLayerTol)

      if (task .eq. "setupgeometry") then
        call hsd_get(pNode, "PLsDefined", nPLs(ii), stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "Missing required value: 'PLsDefined'")
        call hsd_get(pNode, "Atoms", buffer, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "Missing required value: 'Atoms'")
        call hsd_get_table(pNode, "Atoms", pTmp, stat, auto_wrap=.true.)
        call hsd_get_attrib(pNode, "Atoms", modif, stat)
        if (stat /= HSD_STAT_OK) modif = ""
        if (isZeroBased(modif)) then
          selectionRange(:) = [0, size(geom%species) - 1]
          ishift = 1
        else
          selectionRange(:) = [1, size(geom%species)]
          ishift = 0
        end if
        call getSelectedAtomIndices(pTmp, buffer, geom%speciesNames, geom%species, &
            & iAtInRegion(ii)%data, selectionRange=selectionRange)
        iAtInRegion(ii)%data = iAtInRegion(ii)%data + ishift
        block
          real(dp), allocatable :: vecArr(:)
          call hsd_get(pNode, "ContactVector", vecArr, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(pNode, &
              & "Missing required value: 'ContactVector'")
          call hsd_get_attrib(pNode, "ContactVector", modif, stat)
          if (stat /= HSD_STAT_OK) modif = ""
          if (size(vecArr) == 3) then
             vec = vecArr
             call convertUnitHsd(modif, lengthUnits, pNode, vec)
             ! check vector is along x y or z
             if (count(vec == 0.0_dp) < 2 ) then
               call error("ContactVector must be along either x, y or z")
             end if
             contacts(ii)%lattice = vec
             contacts(ii)%shiftAccuracy = contactLayerTol
          else
             call error("ContactVector must define three entries")
          end if
        end block
      else
        call error("Invalid task for setpugeometry tool")
      end if

    end do

    contains

      function isZeroBased(chr)
        character(*), intent(in) :: chr
        logical :: isZeroBased

        if (trim(chr) == "" .or. tolower(trim(chr)) == "onebased") then
          isZeroBased = .false.
        else if (tolower(trim(chr)) .eq. "zerobased") then
          isZeroBased = .true.
        else
          call error("Modifier in Atoms " // trim(chr) // " not recongnized")
        end if

      end function isZeroBased

  end subroutine readContacts


  subroutine getSKcutoff(node, geo, mSKCutoff)
    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Maximum cutoff distance from sk files
    real(dp), intent(out) :: mSKCutoff

    ! Locals
    type(hsd_table), pointer :: child
    integer :: skInterMeth, stat
    logical :: oldSKInter

    call hsd_get_or_set(node, "OldSKInterpolation", oldSKInter, .false.)
    if (oldSKInter) then
      skInterMeth = skEqGridOld
    else
      skInterMeth = skEqGridNew
    end if

    call hsd_get_table(node, "TruncateSKRange", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      call warning("Artificially truncating the SK table, this is normally a bad idea!")
      call SKTruncations(child, mSKCutOff, skInterMeth)
    else
      call readSKFiles(node, geo%nSpecies, geo%speciesNames, mSKCutOff)
    end if
    ! The fudge distance is added to get complete cutoff
    select case(skInterMeth)
    case(skEqGridOld)
      mSKCutOff = mSKCutOff + distFudgeOld
    case(skEqGridNew)
      mSKCutOff = mSKCutOff + distFudge
    end select

  end subroutine getSKcutoff

  !> Reads Slater-Koster files
  !> Should be replaced with a more sophisticated routine, once the new SK-format has been
  !> established
  subroutine readSKFiles(node, nSpecies, speciesNames, maxSKcutoff)
    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Nr. of species in the system
    integer, intent(in) :: nSpecies

    !> Array with specie names
    character(mc), intent(in) :: speciesNames(:)

    !> Maximum SK cutoff distance obtained from SK files
    real(dp), intent(out) :: maxSKcutoff

    type(hsd_table), pointer :: value1, child, child2
    character(len=:), allocatable :: buffer, buffer2
    type(TCharLcArray), allocatable :: skFiles(:,:)
    type(TOldSKData) :: skData
    integer :: iSp1, iSp2, ii, stat
    character(lc) :: prefix, suffix, separator, elem1, elem2, strTmp
    character(lc) :: fileName
    logical :: tLower, tExist

    ! Slater-Koster files
    allocate(skFiles(nSpecies, nSpecies))
    do iSp1 = 1, nSpecies
      do iSp2 = 1, nSpecies
        allocate(skFiles(iSp2, iSp1)%items(0))
      end do
    end do
    call hsd_get_choice(node, "SlaterKosterFiles", buffer, value1, stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'SlaterKosterFiles'")
    call hsd_get_table(node, "SlaterKosterFiles", child, stat, auto_wrap=.true.)
    select case(buffer)
    case ("type2filenames")
      call hsd_get_or_set(value1, "Prefix", buffer2, "")
      prefix = unquote(buffer2)
      call hsd_get_or_set(value1, "Suffix", buffer2, "")
      suffix = unquote(buffer2)
      call hsd_get_or_set(value1, "Separator", buffer2, "")
      separator = unquote(buffer2)
      call hsd_get_or_set(value1, "LowerCaseTypeName", tLower, .false.)
      do iSp1 = 1, nSpecies
        if (tLower) then
          elem1 = tolower(speciesNames(iSp1))
        else
          elem1 = speciesNames(iSp1)
        end if
        do iSp2 = 1, nSpecies
          if (tLower) then
            elem2 = tolower(speciesNames(iSp2))
          else
            elem2 = speciesNames(iSp2)
          end if
          strTmp = trim(prefix) // trim(elem1) // trim(separator) &
              &// trim(elem2) // trim(suffix)
          skFiles(iSp2, iSp1)%items = [skFiles(iSp2, iSp1)%items, strTmp]
          inquire(file=strTmp, exist=tExist)
          if (.not. tExist) then
            call dftbp_error(value1, "SK file with generated name '" &
                &// trim(strTmp) // "' does not exist.")
          end if
        end do
      end do
    case default
      value1%processed = .false.
      do iSp1 = 1, nSpecies
        do iSp2 = 1, nSpecies
          strTmp = trim(speciesNames(iSp1)) // "-" &
              &// trim(speciesNames(iSp2))
          block
            character(len=:), allocatable :: strArr(:)
            call hsd_get(child, trim(strTmp), strArr, stat=stat)
            if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing required value: '" &
                & // trim(strTmp) // "'")
            call hsd_get_table(child, trim(strTmp), child2, stat, auto_wrap=.true.)
            do ii = 1, size(strArr)
              strTmp = strArr(ii)
              inquire(file=strTmp, exist=tExist)
              if (.not. tExist) then
                call dftbp_error(child2, "SK file '" // trim(strTmp) &
                    &// "' does not exist'")
              end if
              skFiles(iSp2, iSp1)%items = [skFiles(iSp2, iSp1)%items, strTmp]
            end do
          end block
        end do
      end do
    end select

    write(stdout, "(A)") "Reading SK-files:"
    do iSp1 = 1, nSpecies
      do iSp2 = iSp1, nSpecies
        fileName = skFiles(iSp2, iSp1)%items(1)
        write(stdout,*) trim(fileName)
        call readFromFile(skData, fileName, (iSp1 == iSp2))
        maxSKcutoff = max(maxSKcutoff, skData%dist * size(skData%skHam,1))
      end do
    end do
    write(stdout, "(A)") "Done."
    write(stdout, *)

    do iSp1 = 1, nSpecies
      do iSp2 = 1, nSpecies
        if (allocated(skFiles(iSp2, iSp1)%items)) deallocate(skFiles(iSp2, iSp1)%items)
      end do
    end do
    deallocate(skFiles)

  end subroutine readSKFiles


  !> Options for truncation of the SK data sets at a fixed distance
  subroutine SKTruncations(node, truncationCutOff, skInterMeth)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> This is the resulting cutoff distance
    real(dp), intent(out) :: truncationCutOff

    !> Method of the sk interpolation
    integer, intent(in) :: skInterMeth

    logical :: tHardCutOff
    integer :: stat
    type(hsd_table), pointer :: field
    character(len=:), allocatable :: modifier

    ! Artificially truncate the SK table
    call hsd_get(node, "SKMaxDistance", truncationCutOff, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'SKMaxDistance'")
    call hsd_get_attrib(node, "SKMaxDistance", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    call hsd_get_table(node, "SKMaxDistance", field, stat, auto_wrap=.true.)
    call convertUnitHsd(modifier, lengthUnits, field, truncationCutOff)

    call hsd_get_or_set(node, "HardCutOff", tHardCutOff, .true.)
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
      call dftbp_error(field, "Truncation is shorter than the minimum distance over which SK data&
          & goes to 0")
    end if

  end subroutine SKTruncations

end module transporttools_parser
