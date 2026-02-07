!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> HSD-parsing related helper routines.
!>
!> Supports multiple input formats: HSD, XML, JSON, TOML. The format is auto-detected
!> from the file extension. The default input file is "dftb_in.hsd" but the code will
!> also search for dftb_in.xml, dftb_in.json, and dftb_in.toml if the HSD file is not found.
module dftbp_dftbplus_hsdhelpers
  use dftbp_common_globalenv, only : stdOut, tIoProc
  use dftbp_dftbplus_inputdata, only : TInputData
  use dftbp_dftbplus_parser, only : parseHsdTree, readHsdFile, rootTag, TParserFlags
  use dftbp_io_hsdcompat, only : hsd_table, destroyNode, dumpHsd, getChild, warnUnprocessedNodes
  use dftbp_io_message, only : error
  implicit none

  private
  public :: parseHsdInput, doPostParseJobs


  !> Supported input file base name (without extension)
  character(*), parameter :: inputBaseName = "dftb_in"

  !> Name of the DFTB+ input file (default HSD format)
  character(*), parameter :: hsdFileName = inputBaseName // ".hsd"

  !> Supported input formats with their extensions (searched in order)
  integer, parameter :: nFormats = 4
  character(len=5), parameter :: inputExtensions(nFormats) = [character(len=5) :: &
      & ".hsd", ".xml", ".json", ".toml"]

  !> Name of the DFTB+ processed input file
  character(*), parameter :: hsdProcFileName = "dftb_pin.hsd"

contains

  !> Parses input file and returns initialised input structure.
  !>
  !> Searches for input files in the following order: dftb_in.hsd, dftb_in.xml,
  !> dftb_in.json, dftb_in.toml. Uses the first one found.
  subroutine parseHsdInput(input)

    !> Input data parsed from the input file
    type(TInputData), intent(out) :: input

    type(hsd_table), pointer :: hsdTree
    type(TParserFlags) :: parserFlags
    character(len=:), allocatable :: inputFile

    call findInputFile_(inputFile)
    write(stdout, "(A)") "Reading input file '" // inputFile // "'"
    call readHsdFile(inputFile, hsdTree)
    call parseHsdTree(hsdTree, input, parserFlags)
    call doPostParseJobs(hsdTree, parserFlags)
    call destroyNode(hsdTree)

  end subroutine parseHsdInput


  !> Searches for the input file across supported formats and returns the first found.
  !>
  !> Searches for: dftb_in.hsd, dftb_in.xml, dftb_in.json, dftb_in.toml
  !> Stops with an error if none is found.
  subroutine findInputFile_(inputFile)

    !> Path of the found input file
    character(len=:), allocatable, intent(out) :: inputFile

    logical :: tExist
    integer :: ii

    do ii = 1, nFormats
      inputFile = inputBaseName // trim(inputExtensions(ii))
      inquire(file=inputFile, exist=tExist)
      if (tExist) return
    end do

    call error("No input file found. Searched for: " // &
        & "dftb_in.hsd, dftb_in.xml, dftb_in.json, dftb_in.toml")

  end subroutine findInputFile_


  !> Execute parser related tasks (warning, processed input dumping) needed after parsing
  subroutine doPostParseJobs(hsdTree, parserFlags)

    !> Tree representation of the HSD input
    type(hsd_table), pointer, intent(in) :: hsdTree

    !> Parser specific settings in the output
    type(TParserFlags), intent(in) :: parserFlags

    type(hsd_table), pointer :: root

    call getChild(hsdTree, rootTag, root)

    ! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, parserFlags%tIgnoreUnprocessed)

    ! Dump processed tree in HSD format
    if (tIoProc .and. parserFlags%tWriteHSD) then
      call dumpHsd(hsdTree, hsdProcFileName)
      write(stdout, '(/,/,A)') "Processed input in HSD format written to '" // hsdProcFileName&
          & // "'"
    end if

    ! Stop, if only parsing is required
    if (parserFlags%tStop) then
      call error("Keyword 'StopAfterParsing' is set to Yes. Stopping.")
    end if

  end subroutine doPostParseJobs


end module dftbp_dftbplus_hsdhelpers
