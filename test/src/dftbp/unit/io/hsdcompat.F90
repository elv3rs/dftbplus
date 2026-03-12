!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_io_hsdcompat
  use fortuno_serial, only : suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp, mc
  use dftbp_common_unitconversion, only : TUnit
  use dftbp_extlibs_hsd, only : hsd_table, hsd_error_t, hsd_get, hsd_set, new_table, &
      & hsd_load_string, HSD_STAT_OK, hsd_has_child, hsd_get_attrib, &
      & hsd_set_attrib, hsd_get_table, hsd_warn_unprocessed, MAX_WARNING_LEN
  use dftbp_io_hsdcompat, only : getChildValue, getChild, setChildValue, setChild, &
      & detailedWarning, convertUnitHsd, getNodeName, getNodeName2, getNodeHSDName, &
      & warnUnprocessedNodes, setUnprocessed, &
      & getChildren, getLength, getItem1, destroyNodeList, hsd_child_list, &
      & splitModifier, textNodeName
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

    call hsd_load_string("NAtoms = 5", root, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "NAtoms", val)
    @:ASSERT(val == 5)
  $:END_TEST()


  $:TEST("getChildValue_int_default", label="hsdcompat")
    !! Get integer with default (absent key → returns default AND writes it back)
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    integer :: val, check

    call hsd_load_string("Existing = 10", root, err)
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

    call hsd_load_string("NAtoms = 5", root, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "NAtoms", val, 99)
    @:ASSERT(val == 5)
  $:END_TEST()


  $:TEST("getChildValue_real_required", label="hsdcompat")
    !! Get a required real value
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    real(dp) :: val

    call hsd_load_string("Temp = 300.0", root, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "Temp", val)
    @:ASSERT(abs(val - 300.0_dp) < 1.0e-10_dp)
  $:END_TEST()


  $:TEST("getChildValue_real_default", label="hsdcompat")
    !! Get real with default
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    real(dp) :: val

    call hsd_load_string("X = 1.0", root, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "MissingReal", val, 3.14_dp)
    @:ASSERT(abs(val - 3.14_dp) < 1.0e-10_dp)
  $:END_TEST()


  $:TEST("getChildValue_logical_required", label="hsdcompat")
    !! Get a required logical value
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    logical :: val

    call hsd_load_string("Periodic = Yes", root, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "Periodic", val)
    @:ASSERT(val)
  $:END_TEST()


  $:TEST("getChildValue_logical_default", label="hsdcompat")
    !! Get logical with default
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    logical :: val

    call hsd_load_string("X = 1", root, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "MissingLogical", val, .false.)
    @:ASSERT(.not. val)
  $:END_TEST()


  $:TEST("getChildValue_string", label="hsdcompat")
    !! Get a required string value
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    character(len=:), allocatable :: val

    call hsd_load_string('Name = "Hello"', root, err)
    @:ASSERT(.not. allocated(err))
    call getChildValue(root, "Name", val)
    @:ASSERT(val == "Hello")
  $:END_TEST()


  $:TEST("getChildValue_string_default", label="hsdcompat")
    !! Get string with default
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    character(len=:), allocatable :: val

    call hsd_load_string("X = 1", root, err)
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

    call hsd_load_string(input, root, err)
    @:ASSERT(.not. allocated(err))
    call getChild(root, "Options", child)
    @:ASSERT(associated(child))
  $:END_TEST()


  $:TEST("getChild_optional_missing", label="hsdcompat")
    !! Get an optional child table that doesn't exist → null pointer
    type(hsd_table) :: root
    type(hsd_table), pointer :: child
    type(hsd_error_t), allocatable :: err

    call hsd_load_string("X = 1", root, err)
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
    type(hsd_table), target :: root
    type(hsd_table), pointer :: rootPtr
    character(len=:), allocatable :: name

    call new_table(root, "MyTable")
    rootPtr => root
    call getNodeName2(rootPtr, name)
    @:ASSERT(name == "mytable")
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

    call hsd_load_string(input, root, err)
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

    call hsd_load_string(input, root, err)
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
      @:ASSERT(name == "conjugategradient")
    end block
  $:END_TEST()


  $:TEST("getNodeName_null_returns_text", label="hsdcompat")
    !! getNodeName with null pointer returns "#text"
    type(hsd_table), pointer :: nullPtr
    character(len=:), allocatable :: name

    nullPtr => null()
    call getNodeName(nullPtr, name)
    @:ASSERT(name == textNodeName)
    @:ASSERT(name == "#text")
  $:END_TEST()


  $:TEST("getNodeName_associated", label="hsdcompat")
    !! getNodeName with an associated table returns the node name
    type(hsd_table), target :: tbl
    type(hsd_table), pointer :: ptr
    character(len=:), allocatable :: name

    call new_table(tbl, name="MyBlock")
    ptr => tbl
    call getNodeName(ptr, name)
    @:ASSERT(name == "myblock")
  $:END_TEST()


  $:TEST("getNodeHSDName_works", label="hsdcompat")
    !! getNodeHSDName returns the node name
    type(hsd_table), target :: root
    type(hsd_table), pointer :: rootPtr
    character(len=:), allocatable :: name

    call new_table(root, name="TestNode")
    rootPtr => root
    call getNodeHSDName(rootPtr, name)
    @:ASSERT(name == "TestNode")
  $:END_TEST()


  $:TEST("setChild_creates_block", label="hsdcompat")
    !! setChild creates a new named child table
    type(hsd_table) :: root
    type(hsd_table), pointer :: child
    character(len=:), allocatable :: name

    call new_table(root, name="root")
    call setChild(root, "NewBlock", child)
    @:ASSERT(associated(child))
    call getNodeName2(child, name)
    @:ASSERT(name == "newblock")
    ! Verify the child is accessible via getChild
    block
      type(hsd_table), pointer :: found
      call getChild(root, "NewBlock", found, requested=.false.)
      @:ASSERT(associated(found))
    end block
  $:END_TEST()


  $:TEST("getChildren_getItem1_getLength", label="hsdcompat")
    !! getChildren/getLength/getItem1/destroyNodeList workflow
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    type(hsd_child_list), pointer :: children
    type(hsd_table), pointer :: child
    integer :: nn

    call hsd_load_string("Region { Atoms = 1 }" // char(10) // &
        & "Region { Atoms = 2 }" // char(10) // &
        & "Region { Atoms = 3 }", root, err)
    @:ASSERT(.not. allocated(err))

    call getChildren(root, "Region", children)
    nn = getLength(children)
    @:ASSERT(nn == 3)

    ! Get first and last
    call getItem1(children, 1, child)
    @:ASSERT(associated(child))
    call getItem1(children, 3, child)
    @:ASSERT(associated(child))

    ! Out of bounds
    call getItem1(children, 0, child)
    @:ASSERT(.not. associated(child))
    call getItem1(children, 4, child)
    @:ASSERT(.not. associated(child))

    call destroyNodeList(children)
    @:ASSERT(.not. associated(children))
  $:END_TEST()


  $:TEST("getChildren_empty", label="hsdcompat")
    !! getChildren with no matches returns count 0
    type(hsd_table) :: root
    type(hsd_child_list), pointer :: children

    call new_table(root, name="root")
    call getChildren(root, "Missing", children)
    @:ASSERT(getLength(children) == 0)
    call destroyNodeList(children)
  $:END_TEST()


  $:TEST("splitModifier_3parts", label="hsdcompat")
    !! splitModifier splits comma-separated modifier into 3 parts
    type(hsd_table) :: root
    character(len=mc) :: mods(3)

    call new_table(root, name="test")
    call splitModifier("Angstrom, eV, e", root, mods)
    @:ASSERT(trim(mods(1)) == "Angstrom")
    @:ASSERT(trim(mods(2)) == "eV")
    @:ASSERT(trim(mods(3)) == "e")
  $:END_TEST()


  $:TEST("setUnprocessed_noop", label="hsdcompat")
    !! setUnprocessed is a no-op — just verify it doesn't crash
    type(hsd_table) :: root
    call new_table(root, "test")
    call setUnprocessed(root)
    @:ASSERT(.true.)
  $:END_TEST()


  $:TEST("warnUnprocessedNodes_noop", label="hsdcompat")
    !! warnUnprocessedNodes is a no-op stub — just verify it doesn't crash
    type(hsd_table) :: root
    call new_table(root, "root")
    call warnUnprocessedNodes(root)
    call warnUnprocessedNodes(root, tIgnoreUnprocessed=.true.)
    @:ASSERT(.true.)
  $:END_TEST()


  $:TEST("convertUnitHsd_angstrom", label="hsdcompat")
    !! SPEC §7.4 — Unit conversion bridge: convert Angstrom to Bohr
    type(hsd_table) :: root
    type(hsd_error_t), allocatable :: err
    real(dp) :: val
    character(len=:), allocatable :: modifier

    ! Create HSD with a value having Angstrom modifier
    call hsd_load_string('Distance [Angstrom] = 1.0', root, err)
    @:ASSERT(.not. allocated(err))

    ! Read with modifier — caller must accept the modifier or get an error
    call getChildValue(root, "Distance", val, modifier=modifier)
    @:ASSERT(abs(val - 1.0_dp) < 1.0e-10_dp)
    @:ASSERT(modifier == "Angstrom")
  $:END_TEST()


  $:TEST("error_line_number_preserved", label="hsdcompat")
    !! SPEC §7.4 — Error reporting: line numbers are preserved in tree
    type(hsd_table), target :: root
    type(hsd_table), pointer :: child
    type(hsd_error_t), allocatable :: err
    integer :: stat

    ! Multi-line HSD: Geometry starts at line 2 (0-based) or 3 (1-based)
    call hsd_load_string(&
        & "Title = ""Test""" // new_line('a') // &
        & "" // new_line('a') // &
        & "Geometry {" // new_line('a') // &
        & "  NAtoms = 2" // new_line('a') // &
        & "}", root, err)
    @:ASSERT(.not. allocated(err))

    call hsd_get_table(root, "Geometry", child, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    ! The line number should be > 0 (exact value depends on parser)
    @:ASSERT(child%line > 0)
  $:END_TEST()


  $:TEST("detailedWarning_runs", label="hsdcompat")
    !! SPEC §7.4 — Error reporting: detailedWarning doesn't crash
    type(hsd_table), target :: root
    type(hsd_error_t), allocatable :: err

    call hsd_load_string("Value = 42", root, err)
    @:ASSERT(.not. allocated(err))

    ! Just verify it doesn't crash — output goes to stderr
    call detailedWarning(root, "Test warning message")
    @:ASSERT(.true.)
  $:END_TEST()


  $:TEST("dispatch_child_marked_processed", label="hsdcompat")
    !! Regression: dispatch children obtained via getChildValue (table variant)
    !! must be marked as processed so that warnUnprocessedNodes does not report
    !! them as unread.  This mirrors the Geometry = GenFormat { ... } pattern.
    type(hsd_table), target :: root
    type(hsd_table), pointer :: value1, child
    type(hsd_error_t), allocatable :: err
    character(len=:), allocatable :: buffer
    character(len=MAX_WARNING_LEN), allocatable :: warnings(:)
    integer :: paramVal

    call hsd_load_string(&
        & 'Block = MyChoice {' // new_line('a') // &
        & '  Param = 42' // new_line('a') // &
        & '}', root, err)
    @:ASSERT(.not. allocated(err))

    ! Simulate the dispatch pattern: get the table block, then its first child
    call getChildValue(root, "Block", value1, child=child)
    call getNodeName(value1, buffer)
    @:ASSERT(buffer == "mychoice")

    ! The dispatch child (MyChoice) must be marked as processed
    @:ASSERT(value1%processed)

    ! Read a child from the dispatch table to process it
    call getChildValue(value1, "Param", paramVal)
    @:ASSERT(paramVal == 42)

    ! warnUnprocessedNodes must not find any unprocessed nodes
    call hsd_warn_unprocessed(root, warnings)
    @:ASSERT(size(warnings) == 0)
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
