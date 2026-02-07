!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_io_hsddata
  use fortuno_serial, only : suite => serial_suite_item, test_list
  use dftbp_extlibs_hsddata, only : hsd_table, hsd_error_t, hsd_get, hsd_set, new_table, &
      & data_load_string, data_dump_to_string, DATA_FMT_HSD, DATA_FMT_JSON, hsd_has_child, &
      & HSD_STAT_OK, dp
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains


  $:TEST("load_string_get_integer", label="hsddata")
    !! Load a simple HSD string and retrieve an integer value
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: error
    integer :: val
    integer :: stat

    call data_load_string("Answer = 42", root, DATA_FMT_HSD, error)
    @:ASSERT(.not. allocated(error))
    call hsd_get(root, "Answer", val, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(val == 42)
  $:END_TEST()


  $:TEST("load_string_get_real", label="hsddata")
    !! Load a simple HSD string and retrieve a real value
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: error
    real(dp) :: val
    integer :: stat

    call data_load_string("Pi = 3.14159", root, DATA_FMT_HSD, error)
    @:ASSERT(.not. allocated(error))
    call hsd_get(root, "Pi", val, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(abs(val - 3.14159_dp) < 1.0e-10_dp)
  $:END_TEST()


  $:TEST("load_string_get_string", label="hsddata")
    !! Load a simple HSD string and retrieve a string value
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: error
    character(len=:), allocatable :: val
    integer :: stat

    call data_load_string('Name = "DFTB+"', root, DATA_FMT_HSD, error)
    @:ASSERT(.not. allocated(error))
    call hsd_get(root, "Name", val, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(val == "DFTB+")
  $:END_TEST()


  $:TEST("load_string_nested_table", label="hsddata")
    !! Load an HSD string with a nested block and navigate the tree
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: error
    integer :: val
    integer :: stat
    character(len=*), parameter :: input = &
        & "Geometry {" // new_line('a') // &
        & "  NAtoms = 2" // new_line('a') // &
        & "}"

    call data_load_string(input, root, DATA_FMT_HSD, error)
    @:ASSERT(.not. allocated(error))
    @:ASSERT(hsd_has_child(root, "Geometry"))
    call hsd_get(root, "Geometry/NAtoms", val, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(val == 2)
  $:END_TEST()


  $:TEST("build_tree_and_dump", label="hsddata")
    !! Build an HSD tree programmatically and dump to string
    type(hsd_table) :: root
    character(len=:), allocatable :: output

    call new_table(root, "root")
    call hsd_set(root, "Answer", 42)
    call data_dump_to_string(root, output, DATA_FMT_HSD)
    @:ASSERT(allocated(output))
    @:ASSERT(len_trim(output) > 0)
  $:END_TEST()


  $:TEST("roundtrip_hsd", label="hsddata")
    !! Round-trip: build tree -> dump to HSD string -> load back -> verify
    type(hsd_table) :: root, loaded
    type(hsd_error_t), allocatable :: error
    character(len=:), allocatable :: hsd_str
    integer :: val
    integer :: stat

    call new_table(root, "root")
    call hsd_set(root, "Value", 99)
    call data_dump_to_string(root, hsd_str, DATA_FMT_HSD)
    @:ASSERT(allocated(hsd_str))

    call data_load_string(hsd_str, loaded, DATA_FMT_HSD, error)
    @:ASSERT(.not. allocated(error))
    call hsd_get(loaded, "Value", val, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(val == 99)
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("hsddata", test_list([&
            $:TEST_ITEMS(label="hsddata")
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_io_hsddata
