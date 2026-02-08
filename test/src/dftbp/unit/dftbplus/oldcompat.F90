!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

!> Unit tests for oldcompat version conversion functions.
module test_dftbplus_oldcompat
  use fortuno_serial, only : suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp
  use hsd, only : hsd_table, hsd_error_t, hsd_get, hsd_set, HSD_STAT_OK, hsd_has_child, &
      & hsd_get_table
  use hsd_data, only : new_table, data_load_string, DATA_FMT_HSD
  use dftbp_dftbplus_oldcompat, only : convertOldHSD
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains


  $:TEST("convert_1_2_rename_speciesnames", label="oldcompat")
    !! Version 1→2: SpeciesNames → TypeNames
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot, geo
    type(hsd_error_t), allocatable :: err
    integer :: stat

    call data_load_string(&
        & "Geometry {" // new_line("a") // &
        & "  SpeciesNames = { ""Si"" ""C"" }" // new_line("a") // &
        & "}" // new_line("a") // &
        & "ParserOptions {}", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    pRoot => root
    call convertOldHSD(pRoot, 1, 2)

    ! SpeciesNames should be gone, TypeNames should exist
    call hsd_get_table(root, "Geometry", geo, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(hsd_has_child(geo, "TypeNames", .true.))
    @:ASSERT(.not. hsd_has_child(geo, "SpeciesNames", .true.))
  $:END_TEST()


  $:TEST("convert_6_7_rename_scc_solver", label="oldcompat")
    !! Version 6→7: OrbitalResolvedSCC → ShellResolvedSCC, Eigensolver → Solver
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot, ham, dftb
    type(hsd_error_t), allocatable :: err
    integer :: stat

    call data_load_string(&
        & "Hamiltonian {" // new_line("a") // &
        & "  DFTB {" // new_line("a") // &
        & "    OrbitalResolvedSCC = Yes" // new_line("a") // &
        & "    Eigensolver {" // new_line("a") // &
        & "      DivideAndConquer {}" // new_line("a") // &
        & "    }" // new_line("a") // &
        & "  }" // new_line("a") // &
        & "}" // new_line("a") // &
        & "ParserOptions {}", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    pRoot => root
    call convertOldHSD(pRoot, 6, 7)

    call hsd_get_table(root, "Hamiltonian", ham, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    call hsd_get_table(ham, "DFTB", dftb, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(hsd_has_child(dftb, "ShellResolvedSCC", .true.))
    @:ASSERT(.not. hsd_has_child(dftb, "OrbitalResolvedSCC", .true.))
    @:ASSERT(hsd_has_child(dftb, "Solver", .true.))
    @:ASSERT(.not. hsd_has_child(dftb, "Eigensolver", .true.))
  $:END_TEST()


  $:TEST("convert_7_8_rename_txt_remove_xml", label="oldcompat")
    !! Version 7→8: EigenvectorsAsTxt → EigenvectorsAsText, WriteXMLInput removed
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot, analysis, parserOpts
    type(hsd_error_t), allocatable :: err
    integer :: stat

    call data_load_string(&
        & "Analysis {" // new_line("a") // &
        & "  EigenvectorsAsTxt = Yes" // new_line("a") // &
        & "}" // new_line("a") // &
        & "ParserOptions {" // new_line("a") // &
        & "  WriteXMLInput = No" // new_line("a") // &
        & "}", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    pRoot => root
    call convertOldHSD(pRoot, 7, 8)

    call hsd_get_table(root, "Analysis", analysis, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(hsd_has_child(analysis, "EigenvectorsAsText", .true.))
    @:ASSERT(.not. hsd_has_child(analysis, "EigenvectorsAsTxt", .true.))

    call hsd_get_table(root, "ParserOptions", parserOpts, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(.not. hsd_has_child(parserOpts, "WriteXMLInput", .true.))
  $:END_TEST()


  $:TEST("convert_8_9_lbfgs_linesearch", label="oldcompat")
    !! Version 8→9: lBFGS gets LineSearch=Yes and oldLineSearch=Yes
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot, driver, lbfgs
    type(hsd_error_t), allocatable :: err
    integer :: stat
    logical :: lineSearch, oldLineSearch

    call data_load_string(&
        & "Driver {" // new_line("a") // &
        & "  lBFGS {" // new_line("a") // &
        & "    MaxSteps = 100" // new_line("a") // &
        & "  }" // new_line("a") // &
        & "}" // new_line("a") // &
        & "ParserOptions {}", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    pRoot => root
    call convertOldHSD(pRoot, 8, 9)

    call hsd_get_table(root, "Driver", driver, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    call hsd_get_table(driver, "lBFGS", lbfgs, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(hsd_has_child(lbfgs, "LineSearch", .true.))
    @:ASSERT(hsd_has_child(lbfgs, "oldLineSearch", .true.))
    call hsd_get(lbfgs, "LineSearch", lineSearch, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK .and. lineSearch)
    call hsd_get(lbfgs, "oldLineSearch", oldLineSearch, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK .and. oldLineSearch)
  $:END_TEST()


  $:TEST("convert_10_11_solvation_rescale", label="oldcompat")
    !! Version 10→11: GeneralizedBorn gets RescaleSolvatedFields = No
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot, ham, dftb, solv, gb
    type(hsd_error_t), allocatable :: err
    integer :: stat
    logical :: rescale

    call data_load_string(&
        & "Hamiltonian {" // new_line("a") // &
        & "  DFTB {" // new_line("a") // &
        & "    Solvation {" // new_line("a") // &
        & "      GeneralizedBorn {}" // new_line("a") // &
        & "    }" // new_line("a") // &
        & "  }" // new_line("a") // &
        & "}" // new_line("a") // &
        & "ParserOptions {}", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    pRoot => root
    call convertOldHSD(pRoot, 10, 11)

    call hsd_get_table(root, "Hamiltonian", ham, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    call hsd_get_table(ham, "DFTB", dftb, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    call hsd_get_table(dftb, "Solvation", solv, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    call hsd_get_table(solv, "GeneralizedBorn", gb, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(hsd_has_child(gb, "RescaleSolvatedFields", .true.))
    call hsd_get(gb, "RescaleSolvatedFields", rescale, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK .and. (.not. rescale))
  $:END_TEST()


  $:TEST("convert_13_14_rename_forces_hybrid", label="oldcompat")
    !! Version 13→14: CalculateForces → PrintForces, Rangeseparated → Hybrid
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot, analysis
    type(hsd_error_t), allocatable :: err
    integer :: stat

    call data_load_string(&
        & "Analysis {" // new_line("a") // &
        & "  CalculateForces = Yes" // new_line("a") // &
        & "}" // new_line("a") // &
        & "ParserOptions {}", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    pRoot => root
    call convertOldHSD(pRoot, 13, 14)

    call hsd_get_table(root, "Analysis", analysis, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(hsd_has_child(analysis, "PrintForces", .true.))
    @:ASSERT(.not. hsd_has_child(analysis, "CalculateForces", .true.))
  $:END_TEST()


  $:TEST("convert_multi_version", label="oldcompat")
    !! Multi-version upgrade: 1 → 14 (full chain)
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot, geo, parserOpts
    type(hsd_error_t), allocatable :: err
    integer :: stat, ver

    ! Start with version 1 format: Geometry/SpeciesNames
    call data_load_string(&
        & "Geometry {" // new_line("a") // &
        & "  SpeciesNames = { ""Si"" }" // new_line("a") // &
        & "}" // new_line("a") // &
        & "ParserOptions {}", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    pRoot => root

    ! Convert from version 1 all the way to 14
    call convertOldHSD(pRoot, 1, 14)

    ! Verify version 1→2 transformation applied
    call hsd_get_table(root, "Geometry", geo, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    @:ASSERT(hsd_has_child(geo, "TypeNames", .true.))
    @:ASSERT(.not. hsd_has_child(geo, "SpeciesNames", .true.))

    ! Verify parser version was updated
    call hsd_get_table(root, "ParserOptions", parserOpts, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK)
    call hsd_get(parserOpts, "ParserVersion", ver, stat=stat)
    @:ASSERT(stat == HSD_STAT_OK .and. ver == 14)
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("oldcompat", test_list([&
            $:TEST_ITEMS(label="oldcompat")
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()
  end function tests


end module test_dftbplus_oldcompat
