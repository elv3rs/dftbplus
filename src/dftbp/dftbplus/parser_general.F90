#:include 'common.fypp'
#:include 'error.fypp'

!> Parser routines for the Options block.
module dftbp_dftbplus_parser_general
  use dftbp_dftbplus_input_fileaccess, only : readBinaryAccessTypes
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_elecsolvers_elecsolvers, only : providesEigenvalues
  use hsd, only : hsd_rename_child, hsd_get_or_set, hsd_get_table
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning
  use hsd_data, only : hsd_table, new_table
  use dftbp_type_typegeometry, only : TGeometry
  implicit none

  private
  public :: readOptions

contains


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
    integer :: stat

  #:if WITH_SOCKETS
    tWriteDetailedOutDef = .not. allocated(ctrl%socketInput)
  #:else
    tWriteDetailedOutDef = .true.
  #:endif
    call hsd_get_or_set(node, "WriteDetailedOut", ctrl%tWriteDetailedOut, tWriteDetailedOutDef)

    call hsd_get_or_set(node, "WriteAutotestTag", ctrl%tWriteTagged, .false.)
    call hsd_get_or_set(node, "WriteDetailedXML", ctrl%tWriteDetailedXML, .false.)
    block
      character(len=:), allocatable :: tmpFmt
      call hsd_get_or_set(node, "DetailedOutputFormat", tmpFmt, "xml")
      ctrl%detailedOutputFormat = tmpFmt
    end block
    call hsd_get_or_set(node, "WriteResultsTag", ctrl%tWriteResultsTag, .false.)
    block
      character(len=:), allocatable :: tmpFmt
      call hsd_get_or_set(node, "ResultsOutputFormat", tmpFmt, "tag")
      ctrl%resultsOutputFormat = tmpFmt
    end block

    if (.not.(ctrl%tMD.or.ctrl%isGeoOpt.or.allocated(ctrl%geoOpt))) then
      if (ctrl%tSCC) then
        call hsd_get_or_set(node, "RestartFrequency", ctrl%restartFreq, 20)
      else
        ctrl%restartFreq = 0
      end if
    end if
    if (ctrl%tMD) then
      allocate(ctrl%mdOutput)
      call hsd_get_table(node, "MDOutput", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) then
        block
          type(hsd_table) :: tmpTbl
          call new_table(tmpTbl, name="mdoutput")
          call node%add_child(tmpTbl)
        end block
        call hsd_get_table(node, "MDOutput", child, stat)
      end if
      if (associated(child)) then
        if (providesEigenvalues(ctrl%solver%isolver)) then
          call hsd_get_or_set(child, "AppendBandOut", ctrl%mdOutput%bandStructure, .false.)
        end if
        if (ctrl%tPrintForces) then
          call hsd_get_or_set(child, "Derivatives", ctrl%mdOutput%printForces, .true.)
        end if
        if (ctrl%tPrintMulliken) then
          call hsd_get_or_set(child, "Charges", ctrl%mdOutput%printCharges, .true.)
        end if
        if (ctrl%tAtomicEnergy) then
          call hsd_get_or_set(child, "AtomEnergies", ctrl%mdOutput%printAtomEnergies, .true.)
        end if
      end if
    end if
    call hsd_get_or_set(node, "RandomSeed", ctrl%iSeed, 0, child=child)
    if (ctrl%iSeed < 0) then
      call dftbp_error(child, "Random seed must be greater or equal zero")
    end if
    call hsd_get_or_set(node, "WriteHS", ctrl%tWriteHS, .false.)
    call hsd_get_or_set(node, "WriteRealHS", ctrl%tWriteRealHS, .false.)
    call hsd_rename_child(node, "MinimizeMemoryUsage", "MinimiseMemoryUsage")
    call hsd_get_or_set(node, "MinimiseMemoryUsage", ctrl%tMinMemory, .false., child=child)
    if (ctrl%tMinMemory) then
      call dftbp_warning(child, "Memory minimisation is not working currently, normal calculation&
          & will be used instead")
    end if
    call hsd_get_or_set(node, "ShowFoldedCoords", ctrl%tShowFoldedCoord, .false.)
  #:if DEBUG > 0
    call hsd_get_or_set(node, "TimingVerbosity", ctrl%timingLevel, -1)
  #:else
    call hsd_get_or_set(node, "TimingVerbosity", ctrl%timingLevel, 1)
  #:endif

    if (ctrl%tReadChrg) then
      call hsd_get_or_set(node, "ReadChargesAsText", ctrl%tReadChrgAscii, .false.)
    end if

    call hsd_get_or_set(node, "WriteCharges", ctrl%tWriteCharges, .true.)
    if (ctrl%tWriteCharges) then
      call hsd_get_or_set(node, "WriteChargesAsText", ctrl%tWriteChrgAscii, .false.)
    end if

    ctrl%tSkipChrgChecksum = .false.
    if (.not. ctrl%tFixEf .and. ctrl%tReadChrg) then
      call hsd_get_or_set(node, "SkipChargeTest", ctrl%tSkipChrgChecksum, .false.)
    end if

    if (geom%tPeriodic .or. geom%tHelical .or. geom%areContactsPresent) then
      call hsd_get_or_set(node, "WriteAllAtomGeometry", ctrl%areAllAtomsPrinted, .false.)
    end if

    call readBinaryAccessTypes(node, ctrl%binaryAccessTypes)

  end subroutine readOptions

end module dftbp_dftbplus_parser_general
