!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_io_hsdcompat
  use fortuno_serial, only : suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp
  use dftbp_common_unitconversion, only : TUnit
  use dftbp_extlibs_hsddata, only : hsd_table, hsd_error_t, hsd_get, hsd_set, new_table, &
      & data_load_string, DATA_FMT_HSD, HSD_STAT_OK, hsd_has_child, hsd_get_attrib, &
      & hsd_set_attrib, hsd_get_table
  use dftbp_io_hsdcompat, only : getChildValue, getChild, setChildValue, &
      & detailedWarning, convertUnitHsd, getNodeName2, warnUnprocessedNodes
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

  !> Simple unit for testing convertUnitHsd: 1 Angstrom = 1.8897259886 Bohr
  real(dp), parameter :: angstrom_to_bohr = 1.8897259886_dp

contains


  $:TEST("getChildValue_int_required", label="hsdcompat")
    !! Get a required integer value
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    integer :: val

    call data_load_string("NAtoms = 5", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "NAtoms", val)
    @:ASSERT(val == 5)
  $:END_TEST()


  $:TEST("getChildValue_int_default", label="hsdcompat")
    !! Get integer with default (absent key → returns default AND writes it back)
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    integer :: val, check

    call data_load_string("Existing = 10", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    ! Key doesn't exist → should return default 42
    call getChildValue(root, "Missing", val, 42)
    @:ASSERT(val == 42)
    ! The default should have been written back to the tree
    call hsd_get(root, "Missing", check, stat=val)
    @:ASSERT(val == HSD_STAT_OK)
    @:ASSERT(check == 42)
  $:END_TEST()


  $:TEST("getChildValue_int_existing", label="hsdcompat")
    !! Get integer with default when key exists → returns actual value
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    integer :: val

    call data_load_string("NAtoms = 5", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "NAtoms", val, 99)
    @:ASSERT(val == 5)
  $:END_TEST()


  $:TEST("getChildValue_real_required", label="hsdcompat")
    !! Get a required real value
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    real(dp) :: val

    call data_load_string("Temp = 300.0", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "Temp", val)
    @:ASSERT(abs(val - 300.0_dp) < 1.0e-10_dp)
  $:END_TEST()


  $:TEST("getChildValue_real_default", label="hsdcompat")
    !! Get real with default
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    real(dp) :: val

    call data_load_string("X = 1.0", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "MissingReal", val, 3.14_dp)
    @:ASSERT(abs(val - 3.14_dp) < 1.0e-10_dp)
  $:END_TEST()


  $:TEST("getChildValue_logical_required", label="hsdcompat")
    !! Get a required logical value
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    logical :: val

    call data_load_string("Periodic = Yes", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "Periodic", val)
    @:ASSERT(val)
  $:END_TEST()


  $:TEST("getChildValue_logical_default", label="hsdcompat")
    !! Get logical with default
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    logical :: val

    call data_load_string("X = 1", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "MissingLogical", val, .false.)
    @:ASSERT(.not. val)
  $:END_TEST()


  $:TEST("getChildValue_string", label="hsdcompat")
    !! Get a required string value
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    character(len=:), allocatable :: val

    call data_load_string('Name = "Hello"', root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "Name", val)
    @:ASSERT(val == "Hello")
  $:END_TEST()


  $:TEST("getChildValue_string_default", label="hsdcompat")
    !! Get string with default
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    character(len=:), allocatable :: val

    call data_load_string("X = 1", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "Name", val, "DefaultName")
    @:ASSERT(val == "DefaultName")
  $:END_TEST()


  $:TEST("getChild_found", label="hsdcompat")
    !! Get a child table that exists
    type(hsd_table) :: root
    type(hsd_table), pointer :: child
    type(hsd_error_t), allocatable :: err
    character(len=*), parameter :: input = &
        & "Options {" // new_line('a') // &
        & "  Debug = Yes" // new_line('a') // &
        & "}"

    call data_load_string(input, root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call getChild(root, "Options", child)
    @:ASSERT(associated(child))
  $:END_TEST()


  $:TEST("getChild_optional_missing", label="hsdcompat")
    !! Get an optional child table that doesn't exist → null pointer
    type(hsd_table) :: root
    type(hsd_table), pointer :: child
    type(hsd_error_t), allocatable :: err

    call data_load_string("X = 1", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call getChild(root, "Missing", child, requested=.false.)
    @:ASSERT(.not. associated(child))
  $:END_TEST()


  $:TEST("setChildValue_int", label="hsdcompat")
    !! Set an integer value and read it back
    type(hsd_table) :: root
    integer :: val, stat

    call new_table(root, "root")
    call setChildValue(root, "Count", 42)
    call hsd_get(root, "Count", val, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(val == 42)
  $:END_TEST()


  $:TEST("setChildValue_real", label="hsdcompat")
    !! Set a real value and read it back
    type(hsd_table) :: root
    real(dp) :: val
    integer :: stat

    call new_table(root, "root")
    call setChildValue(root, "Energy", -1.5_dp)
    call hsd_get(root, "Energy", val, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(abs(val + 1.5_dp) < 1.0e-10_dp)
  $:END_TEST()


  $:TEST("setChildValue_char", label="hsdcompat")
    !! Set a character value and read it back
    type(hsd_table) :: root
    character(len=:), allocatable :: val
    integer :: stat

    call new_table(root, "root")
    call setChildValue(root, "Method", "DFTB")
    call hsd_get(root, "Method", val, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(val == "DFTB")
  $:END_TEST()


  $:TEST("getNodeName2_test", label="hsdcompat")
    !! Get the name of a table node
    type(hsd_table) :: root
    character(len=:), allocatable :: name

    call new_table(root, "MyTable")
    call getNodeName2(root, name)
    @:ASSERT(name == "MyTable")
  $:END_TEST()


  $:TEST("getChild_modifier", label="hsdcompat")
    !! Get a child table and its modifier (attribute)
    type(hsd_table) :: root
    type(hsd_table), pointer :: child
    type(hsd_error_t), allocatable :: err
    character(len=:), allocatable :: modifier
    character(len=*), parameter :: input = &
        & 'Solver [ELPA] {' // new_line('a') // &
        & '  X = 1' // new_line('a') // &
        & '}'

    call data_load_string(input, root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call getChild(root, "Solver", child, modifier=modifier)
    @:ASSERT(associated(child))
    @:ASSERT(modifier == "ELPA")
  $:END_TEST()


  $:TEST("getChildValue_dispatch", label="hsdcompat")
    !! Get a child table for dispatch (legacy getChildValue(node, "", child) pattern)
    type(hsd_table) :: root
    type(hsd_table), pointer :: child
    type(hsd_error_t), allocatable :: err
    character(len=:), allocatable :: name
    character(len=*), parameter :: input = &
        & 'Driver {' // new_line('a') // &
        & '  ConjugateGradient {' // new_line('a') // &
        & '    MaxSteps = 100' // new_line('a') // &
        & '  }' // new_line('a') // &
        & '}'

    call data_load_string(input, root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    ! Get the Driver table first
    call getChild(root, "Driver", child)
    @:ASSERT(associated(child))
    ! Now get the dispatch child (first table child, name = "" pattern)
    block
      type(hsd_table), pointer :: dispatch
      call getChildValue(child, "", dispatch)
      @:ASSERT(associated(dispatch))
      call getNodeName2(dispatch, name)
      @:ASSERT(name == "ConjugateGradient")
    end block
  $:END_TEST()


  $:TEST("warnUnprocessedNodes_noop", label="hsdcompat")
    !! warnUnprocessedNodes is a no-op stub — just verify it doesn't crash
    type(hsd_table) :: root
    call new_table(root, "root")
    call warnUnprocessedNodes(root)
    call warnUnprocessedNodes(root, tIgnoreUnprocessed=.true.)
    @:ASSERT(.true.)
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("hsdcompat", test_list([&
            $:TEST_ITEMS(label="hsdcompat")
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_io_hsdcompat
