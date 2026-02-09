#:include 'common.fypp'
#:include 'error.fypp'

!> Parser routines for the Options block.
module dftbp_dftbplus_parser_general
  use dftbp_dftbplus_input_fileaccess, only : readBinaryAccessTypes
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_elecsolvers_elecsolvers, only : providesEigenvalues
  use hsd, only : hsd_rename_child
  use dftbp_io_hsdutils, only : getChildValue
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning
  use hsd_data, only : hsd_table
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
    block
      character(len=:), allocatable :: tmpFmt
      call getChildValue(node, "ResultsOutputFormat", tmpFmt, "tag")
      ctrl%resultsOutputFormat = tmpFmt
    end block

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
      call dftbp_error(child, "Random seed must be greater or equal zero")
    end if
    call getChildValue(node, "WriteHS", ctrl%tWriteHS, .false.)
    call getChildValue(node, "WriteRealHS", ctrl%tWriteRealHS, .false.)
    call hsd_rename_child(node, "MinimizeMemoryUsage", "MinimiseMemoryUsage")
    call getChildValue(node, "MinimiseMemoryUsage", ctrl%tMinMemory, .false., child=child)
    if (ctrl%tMinMemory) then
      call dftbp_warning(child, "Memory minimisation is not working currently, normal calculation&
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

end module dftbp_dftbplus_parser_general
