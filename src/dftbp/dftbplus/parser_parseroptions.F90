#:include 'common.fypp'
#:include 'error.fypp'

!> Parser routines for the ParserOptions block and input version handling.
module dftbp_dftbplus_parser_parseroptions
  use dftbp_common_globalenv, only : stdout
  use dftbp_dftbplus_oldcompat, only : convertOldHSD, minVersion, parserVersion, versionMaps
  use dftbp_io_charmanip, only : i2c, newline, unquote
  use hsd, only : hsd_get_or_set, hsd_get, hsd_get_table, hsd_set, HSD_STAT_OK
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning
  use hsd_data, only : hsd_table
  implicit none

  private
  public :: TParserFlags
  public :: handleInputVersion, readParserOptions, parserVersionFromInputVersion


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


  !> Converts input version to parser version and removes InputVersion node if present.
  subroutine handleInputVersion(root, implicitParserVersion)

    !> Root eventually containing InputVersion
    type(hsd_table), pointer, intent(in) :: root

    !> Parser version corresponding to input version, or unallocated if none has been found
    integer, allocatable, intent(out) :: implicitParserVersion

    type(hsd_table), pointer :: child, dummy
    character(len=:), allocatable :: versionString
    integer :: stat

    call hsd_get_table(root, "InputVersion", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      call hsd_get(child, "#text", versionString, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing required value in 'InputVersion'")
      implicitParserVersion = parserVersionFromInputVersion(trim(unquote(versionString)))
      if (implicitParserVersion == 0) then
        call dftbp_error(child, "Input version '" // trim(unquote(versionString))&
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

    integer :: inputVersion, stat
    type(hsd_table), pointer :: child

    call hsd_get_table(node, "ParserVersion", child, stat, auto_wrap=.true.)
    if (.not. associated(child) .and. .not. present(implicitVersion)) then
      call dftbp_warning(root, "Input containing neither InputVersion nor ParserVersion is&
          & DEPRECATED!(!!) Specify the InputVersion keyword in your input to ensure that future&
          & versions of DFTB+ can also parse it.")
      inputVersion = parserVersion
      call hsd_set(node, "ParserVersion", inputVersion)
    else if (.not. associated(child) .and. present(implicitVersion)) then
      inputVersion = implicitVersion
    else if (associated(child) .and. .not. present(implicitVersion)) then
      call hsd_get(child, "#text", inputVersion, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing required value in 'ParserVersion'")
    else
      call hsd_get(child, "#text", inputVersion, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing required value in 'ParserVersion'")
      if (inputVersion /= implicitVersion) then
        call dftbp_error(child, "Parser version deduced from InputVersion ("&
            & // i2c(implicitVersion) // ") differs from version explicitely set in&
            & ParserVersion (" // i2c(inputVersion) // ")")
      end if
    end if

    if (inputVersion < 1 .or. inputVersion > parserVersion) then
      call dftbp_error(child, "Invalid parser version (" // i2c(inputVersion) // ")")
    else if (inputVersion < minVersion) then
      call dftbp_error(child, &
          & "Sorry, no compatibility mode for parser version " // i2c(inputVersion)&
          & // " (too old)")
    else if (inputVersion /= parserVersion) then
      write(stdout, "(A,I2,A,I2,A)") "***  Converting input from parser version ",&
          & inputVersion, " to parser version ", parserVersion, " ..."
      call convertOldHSD(root, inputVersion, parserVersion)
      write(stdout, "(A,/)") "***  Done."
    end if

    call hsd_get_or_set(node, "WriteHSDInput", flags%tWriteHSD, .true.)
    if (.not. flags%tWriteHSD) then
      call dftbp_warning(node, "WriteHSDInput turned off. You are not guaranteed" // newline // &
          &" to able to obtain the same results with a later version of the code!" // newline // &
          & "(the dftb_pin.hsd file DOES guarantee this)")
    end if
    call hsd_get_or_set(node, "StopAfterParsing", flags%tStop, .false.)

    call hsd_get_or_set(node, "IgnoreUnprocessedNodes", &
        &flags%tIgnoreUnprocessed, .false.)

  end subroutine readParserOptions


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

end module dftbp_dftbplus_parser_parseroptions
