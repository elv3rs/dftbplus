!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                              !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                        !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

!> Smoke tests verifying hsd-fortran integration into DFTB+.
!>
!> These tests confirm that hsd-fortran is correctly linked and usable
!> from dftbplus code.  They exercise both the low-level hsd API (via
!> dftbp_extlibs_hsd) and the compatibility bridge (dftbp_io_hsdcompat).
module test_io_hsd_smoke
  use fortuno_serial, only : suite => serial_suite_item, test_list
  use dftbp_extlibs_hsd, only : hsd_table, hsd_error_t, hsd_get, hsd_set, &
      & hsd_load_string, hsd_dump_to_string, hsd_has_child, &
      & new_table, HSD_STAT_OK, dp
  use dftbp_io_hsdcompat, only : getChildValue, setChildValue, detailedError, &
      & detailedWarning, getChild, getNodeName2, warnUnprocessedNodes
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains

  ! ---------- Low-level hsd API tests ----------

  $:TEST("parse_get_integer", label="hsd_smoke")
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: error
    integer :: val
    integer :: stat

    call hsd_load_string("Answer = 42", root, error)
    @:ASSERT(.not. allocated(error))
    call hsd_get(root, "Answer", val, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(val == 42)
    call root%destroy()
  $:END_TEST()

  $:TEST("parse_get_real", label="hsd_smoke")
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: error
    real(dp) :: val
    integer :: stat

    call hsd_load_string("Pi = 3.14159", root, error)
    @:ASSERT(.not. allocated(error))
    call hsd_get(root, "Pi", val, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(abs(val - 3.14159_dp) < 1.0e-10_dp)
    call root%destroy()
  $:END_TEST()

  $:TEST("parse_get_string", label="hsd_smoke")
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: error
    character(len=:), allocatable :: val
    integer :: stat

    call hsd_load_string('Name = "DFTB+"', root, error)
    @:ASSERT(.not. allocated(error))
    call hsd_get(root, "Name", val, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(val == "DFTB+")
    call root%destroy()
  $:END_TEST()

  $:TEST("nested_table", label="hsd_smoke")
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: error
    integer :: val
    integer :: stat

    call hsd_load_string("Hamiltonian { MaxSccIterations = 100 }", root, error)
    @:ASSERT(.not. allocated(error))
    call hsd_get(root, "Hamiltonian/MaxSccIterations", val, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(val == 100)
    call root%destroy()
  $:END_TEST()

  $:TEST("build_and_dump", label="hsd_smoke")
    type(hsd_table) :: root
    character(len=:), allocatable :: output

    call new_table(root, "root")
    call hsd_set(root, "Value", 7)
    call hsd_dump_to_string(root, output)
    @:ASSERT(len_trim(output) > 0)
    call root%destroy()
  $:END_TEST()

  $:TEST("round_trip", label="hsd_smoke")
    type(hsd_table) :: root1, root2
    type(hsd_error_t), allocatable :: error
    character(len=:), allocatable :: dump1
    integer :: val1, val2
    integer :: stat

    call hsd_load_string("X = 99", root1, error)
    @:ASSERT(.not. allocated(error))
    call hsd_dump_to_string(root1, dump1)
    call hsd_load_string(dump1, root2, error)
    @:ASSERT(.not. allocated(error))
    call hsd_get(root1, "X", val1, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    call hsd_get(root2, "X", val2, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(val1 == val2)
    call root1%destroy()
    call root2%destroy()
  $:END_TEST()

  ! ---------- Compatibility module tests ----------

  $:TEST("compat_getChildValue_int", label="hsd_smoke")
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: error
    integer :: val

    call hsd_load_string("Count = 5", root, error)
    @:ASSERT(.not. allocated(error))
    call getChildValue(root, "Count", val)
    @:ASSERT(val == 5)
    call root%destroy()
  $:END_TEST()

  $:TEST("compat_setChildValue_real", label="hsd_smoke")
    type(hsd_table) :: root
    real(dp) :: val

    call new_table(root, "root")
    call setChildValue(root, "Energy", 1.5_dp)
    call getChildValue(root, "Energy", val)
    @:ASSERT(abs(val - 1.5_dp) < 1.0e-12_dp)
    call root%destroy()
  $:END_TEST()

  ! ---------- Test collection ----------

  function tests()
    type(test_list) :: tests
    tests = test_list([&
        suite("hsd_smoke", test_list([&
            $:TEST_ITEMS(label="hsd_smoke")
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()
  end function tests

end module test_io_hsd_smoke
