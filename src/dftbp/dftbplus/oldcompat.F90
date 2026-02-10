!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to convert HSD input for old parser to the current format.
!> Note: parserVersion is set in parser.F90
module dftbp_dftbplus_oldcompat
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_release, only : TVersionMap
  use dftbp_io_charmanip, only : i2c, newline, tolower
  use hsd_data, only : hsd_table, new_table
  use hsd, only : hsd_get, hsd_has_child, hsd_remove_child, hsd_get_table, HSD_STAT_OK,&
      & hsd_clear_children, hsd_table_ptr, hsd_get_child_tables, hsd_set, hsd_get_or_set,&
      & hsd_get_choice, hsd_rename_child, hsd_get_name
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning
  use dftbp_io_message, only : error
  implicit none

  private
  public :: convertOldHSD, minVersion, parserVersion, versionMaps


  !> Actual input version <-> parser version maps (must be updated at every public release)
  type(TVersionMap), parameter :: versionMaps(*) = [&
      & TVersionMap("25.1", 14),&
      & TVersionMap("24.1", 14), TVersionMap("23.1", 13), TVersionMap("22.2", 12),&
      & TVersionMap("22.1", 11), TVersionMap("21.2", 10), TVersionMap("21.1", 9),&
      & TVersionMap("20.2", 9), TVersionMap("20.1", 8), TVersionMap("19.1", 7),&
      & TVersionMap("18.2", 6), TVersionMap("18.1", 5), TVersionMap("17.1", 5)]

  !> Version of the oldest parser for which compatibility is still maintained
  integer, parameter :: minVersion = 1

  !> Version of the current parser (as latest version)
  integer, parameter :: parserVersion = maxval(versionMaps(:)%parserVersion)


contains


  !> Converts an HSD input for an older parser to the current format
  subroutine convertOldHSD(root, oldVersion, curVersion)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    !> Version number of the old parser
    integer, intent(in) :: oldVersion

    !> Version number of the current parser
    integer, intent(in) :: curVersion

    integer :: version, stat
    type(hsd_table), pointer :: ch1, par

    version = oldVersion
    do while (version < curVersion)
      select case(version)
      case(1)
        call convert_1_2(root)
        version = 2
      case(2)
        call convert_2_3(root)
        version = 3
      case (3)
        call convert_3_4(root)
        version = 4
      case (4)
        call convert_4_5(root)
        version = 5
      case (5)
        call convert_5_6(root)
        version = 6
      case (6)
        call convert_6_7(root)
        version = 7
      case (7)
        call convert_7_8(root)
        version = 8
      case (8)
        call convert_8_9(root)
        version = 9
      case (9)
        call convert_9_10(root)
        version = 10
      case (10)
        call convert_10_11(root)
        version = 11
      case (11)
        call convert_11_12(root)
        version = 12
      case (12)
        call convert_12_13(root)
        version = 13
      case (13)
        call convert_13_14(root)
        version = 14
      end select
    end do

    ! increase the parser version number in the tree - since the resulting dftb_pin would not work
    ! with the old parser as the options have changed to the new parser by now
    call hsd_get_table(root, "ParserOptions", par, stat, auto_wrap=.true.)
    if (.not. associated(par)) then
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="parseroptions")
        call root%add_child(tmpTbl)
      end block
      call hsd_get_table(root, "ParserOptions", par, stat)
    end if
    call hsd_set(par, "ParserVersion", version)

  end subroutine convertOldHSD


  !> Converts input from version 1 to 2. (Version 2 introduced in August 2006)
  subroutine convert_1_2(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table), pointer :: child1
    integer :: stat

    call hsd_get_table(root, "Geometry", child1, stat, auto_wrap=.true.)
    if (associated(child1)) then
      if (associated(child1)) child1%processed = .false.
      call hsd_rename_child(child1, "speciesnames", "typenames", stat, case_insensitive=.true.)
    end if

  end subroutine convert_1_2


  !> Converts input from version 2 to 3. (Version 3 introduced in Nov. 2006)
  subroutine convert_2_3(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table), pointer :: ch1, ch2, par
    logical :: tValue
    integer :: stat

    call hsd_rename_child(root, &
        &"driver/velocityverlet/thermostat/andersen/rescalingprobability", &
        &"reselectprobability", stat, case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, &
        &"Keyword renamed to 'ReselectProbability'.")

    call hsd_rename_child(root, &
        &"driver/velocityverlet/thermostat/andersen/rescaleindividually", &
        &"reselectindividually", stat, case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, &
        &"Keyword renamed to 'ReselectIndividually'.")

    call hsd_get_table(root, "hamiltonian/dftb/variational", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tValue, stat=stat)
      if (associated(ch1)) ch1%processed = .false.
      if (.not. tValue) then
        call dftbp_error(ch1, "Sorry, non-variational energy calculation &
            &is not supported any more!")
      else
        call dftbp_warning(ch1, "Energy calculation is made only variational, option removed.")
        ch1 => null()
      end if
    end if

    call hsd_get_table(root, "hamiltonian/dftb", par, stat, auto_wrap=.true.)
    call hsd_get_table(root, "hamiltonian/dftb/scc", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tValue, stat=stat)
      if (associated(ch1)) ch1%processed = .false.
      if (tValue) then
        call hsd_set(par, "OrbitalResolvedSCC", .true.)
        call hsd_get_table(par, "OrbitalResolvedSCC", ch2, stat, auto_wrap=.true.)
        if (associated(ch2)) ch2%processed = .false.
        call dftbp_warning(ch2, "Calculations are not orbital resolved &
            &per default any more. Keyword 'OrbitalResolvedSCC' added.")
      end if
    end if

    call hsd_rename_child(root, "options/printeigenvectors", "writeeigenvectors", stat, &
        &case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, &
        &"Keyword converted to 'WriteEigenvectors'")

    call hsd_rename_child(root, "options/writetaggedout", "writeautotesttag", stat, &
        &case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, &
        &"Keyword converted to 'WriteAutotestTag'. &
          &Output file name changed to 'autotest.out'")

    call hsd_rename_child(root, "options/writebanddat", "writebandout", stat, &
        &case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, &
        &"Keyword converted to 'WriteBandOut'. &
          &Output file name changed to 'band.out'")

  end subroutine convert_2_3


  !> Converts input from version 3 to 4. (Version 4 introduced in Mar. 2010)
  subroutine convert_3_4(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table),pointer :: node, node2, node3, par
    type(hsd_table_ptr), allocatable :: children(:)
    integer :: ii, stat

    ! Replace range operator with short start:end syntax
    call hsd_get_table(root, "Driver/SteepestDescent/MovedAtoms", node)
    call replaceRange(node)
    call hsd_get_table(root, "Driver/ConjugateGradient/MovedAtoms", node)
    call replaceRange(node)
    call hsd_get_table(root, "Driver/SecondDerivatives/Atoms", node)
    call replaceRange(node)
    call hsd_get_table(root, "Driver/VelocityVerlet/MovedAtoms", node)
    call replaceRange(node)
    call hsd_get_table(root, "Hamiltonian/DFTB/SpinPolarisation/Colinear&
        &/InitialSpin", node)
    if (associated(node)) then
      call hsd_get_child_tables(node, "AtomSpin", children)
      do ii = 1, size(children)
        node2 => children(ii)%ptr
        call hsd_get_table(node2, "Atoms", node3, stat, auto_wrap=.true.)
        if (.not. associated(node3)) call dftbp_error(node2, "Missing required block: 'Atoms'")
        call replaceRange(node3)
      end do
    end if

    call hsd_rename_child(root, "hamiltonian/dftb/spinpolarisation/colinear&
        &/initialspin", "initialspins", stat, case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, "Keyword renamed to 'InitialSpins'.")

  end subroutine convert_3_4

  !> Helper function for Range keyword in convert_3_4
  subroutine replaceRange(node)

    !> node to process
    type(hsd_table), pointer :: node

    type(hsd_table), pointer :: node2
    integer :: bounds(2)
    integer :: stat

    if (associated(node)) then
      call hsd_get_table(node, "Range", node2, stat, auto_wrap=.true.)
      if (associated(node2)) then
        block
          integer, allocatable :: tmpBounds(:)
          call hsd_get(node2, "#text", tmpBounds, stat=stat)
          bounds = tmpBounds
        end block
        call hsd_clear_children(node)
        call hsd_set(node, "#text", &
            &i2c(bounds(1)) // ":" // i2c(bounds(2)))
        call dftbp_warning(node, "Specification 'Range { start end }' &
            &not supported any more, using 'start:end' instead")
      end if
    end if

  end subroutine replaceRange


  !> Converts input from version 4 to 5. (Version 5 introduced in Dec. 2014)
  subroutine convert_4_5(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table), pointer :: ch1, ch2, ch3
    logical :: tVal
    integer :: stat

    call hsd_rename_child(root, "hamiltonian/dftb/eigensolver/standard", &
        &"qr", stat, case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, "Keyword renamed to 'QR'.")

    call hsd_get_table(root, "options/mullikenanalysis", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal, stat=stat)
      call dftbp_warning(ch1, "Keyword moved to Analysis block.")
      call hsd_remove_child(root, "options/mullikenanalysis")
      ch1 => null()
      call hsd_get_table(root, "Analysis", ch1, stat, auto_wrap=.true.)
      if (.not.associated(ch1)) then
        block
          type(hsd_table) :: tmpTbl
          call new_table(tmpTbl, name="analysis")
          call root%add_child(tmpTbl)
        end block
        call hsd_get_table(root, "Analysis", ch1, stat)
      end if
      call hsd_set(ch1, "MullikenAnalysis", tVal)
      if (associated(ch1)) ch1%processed = .false.
    end if

    call hsd_get_table(root, "options/atomresolvedenergies", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal, stat=stat)
      call dftbp_warning(ch1, "Keyword moved to Analysis block.")
      call hsd_remove_child(root, "options/atomresolvedenergies")
      ch1 => null()
      call hsd_get_table(root, "Analysis", ch1, stat, auto_wrap=.true.)
      if (.not.associated(ch1)) then
        block
          type(hsd_table) :: tmpTbl
          call new_table(tmpTbl, name="analysis")
          call root%add_child(tmpTbl)
        end block
        call hsd_get_table(root, "Analysis", ch1, stat)
      end if
      call hsd_set(ch1, "AtomResolvedEnergies", tVal)
      if (associated(ch1)) ch1%processed = .false.
    end if

    call hsd_get_table(root, "options/writeeigenvectors", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal, stat=stat)
      call dftbp_warning(ch1, "Keyword moved to Analysis block.")
      call hsd_remove_child(root, "options/writeeigenvectors")
      ch1 => null()
      call hsd_get_table(root, "Analysis", ch1, stat, auto_wrap=.true.)
      if (.not.associated(ch1)) then
        block
          type(hsd_table) :: tmpTbl
          call new_table(tmpTbl, name="analysis")
          call root%add_child(tmpTbl)
        end block
        call hsd_get_table(root, "Analysis", ch1, stat)
      end if
      call hsd_set(ch1, "WriteEigenvectors", tVal)
      if (associated(ch1)) ch1%processed = .false.
    end if

    call hsd_get_table(root, "options/writebandout", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal, stat=stat)
      call dftbp_warning(ch1, "Keyword moved to Analysis block.")
      call hsd_remove_child(root, "options/writebandout")
      ch1 => null()
      call hsd_get_table(root, "Analysis", ch1, stat, auto_wrap=.true.)
      if (.not.associated(ch1)) then
        block
          type(hsd_table) :: tmpTbl
          call new_table(tmpTbl, name="analysis")
          call root%add_child(tmpTbl)
        end block
        call hsd_get_table(root, "Analysis", ch1, stat)
      end if
      call hsd_set(ch1, "WriteBandOut", tVal)
      if (associated(ch1)) ch1%processed = .false.
    end if

    call hsd_get_table(root, "options/calculateforces", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal, stat=stat)
      call dftbp_warning(ch1, "Keyword moved to Analysis block.")
      call hsd_remove_child(root, "options/calculateforces")
      ch1 => null()
      call hsd_get_table(root, "Analysis", ch1, stat, auto_wrap=.true.)
      if (.not.associated(ch1)) then
        block
          type(hsd_table) :: tmpTbl
          call new_table(tmpTbl, name="analysis")
          call root%add_child(tmpTbl)
        end block
        call hsd_get_table(root, "Analysis", ch1, stat)
      end if
      call hsd_set(ch1, "CalculateForces", tVal)
      if (associated(ch1)) ch1%processed = .false.
    end if

    call hsd_get_table(root, "hamiltonian/dftb", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="differentiation")
        call ch1%add_child(tmpTbl)
      end block
      call hsd_get_table(ch1, "Differentiation", ch2, stat)
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="finitediff")
        call ch2%add_child(tmpTbl)
      end block
      call hsd_get_table(ch2, "FiniteDiff", ch3, stat)
      call hsd_set(ch3, "Delta", 1.0e-2_dp)
      call dftbp_warning(ch2, "Adding legacy step size for finite difference&
          & differentiation")
    end if

    call hsd_get_table(root, "hamiltonian/dftb/spinconstants", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_set(ch1, "ShellResolvedSpin", .true.)
    end if

  end subroutine convert_4_5

  !> Converts input from version 5 to 6. (Version 6 introduced in May. 2018)
  subroutine convert_5_6(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table), pointer :: ch1, ch2, ch3, ch4
    logical :: tVal
    real(dp) :: rTmp
    integer :: stat

    call hsd_rename_child(root, "analysis/localise/pipekmezey/tollerance", &
        &"tolerance", stat, case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, "Keyword converted to 'Tolerance'.")

    call hsd_rename_child(root, "analysis/localise/pipekmezey/sparsetollerances", &
        &"sparsetolerances", stat, case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, "Keyword converted to 'SparseTolerances'.")

    call hsd_get_table(root, "hamiltonian/dftb/dampxh", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal, stat=stat)
      call hsd_get_table(root, "hamiltonian/dftb/dampxhexponent", ch2, stat, auto_wrap=.true.)
      if (tVal .neqv. associated(ch2)) then
        call error("Incompatible combinaton of DampXH and DampXHExponent")
      end if
      if (associated(ch2)) then
        call hsd_get(ch2, "#text", rTmp, stat=stat)
      end if
      call dftbp_warning(ch1, "Keyword DampXH moved to HCorrection block")
      call hsd_remove_child(root, "hamiltonian/dftb/dampxh")
      ch1 => null()
      call hsd_remove_child(root, "hamiltonian/dftb/dampxhexponent")
      ch2 => null()

      ! clean out any HCorrection entry
      call hsd_get_table(root, "hamiltonian/dftb/hcorrection", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) then
        call dftbp_error(ch2, "HCorrection already present.")
      end if

      call hsd_get_table(root, "hamiltonian/dftb", ch2, stat, auto_wrap=.true.)
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="hcorrection")
        call ch2%add_child(tmpTbl)
      end block
      call hsd_get_table(ch2, "HCorrection", ch3, stat)
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="damping")
        call ch3%add_child(tmpTbl)
      end block
      call hsd_get_table(ch3, "Damping", ch4, stat)
      call hsd_set(ch4, "Exponent", rTmp)
      call dftbp_warning(ch3, "Adding Damping to HCorrection")
    end if

  end subroutine convert_5_6


  !> Converts input from version 6 to 7. (Version 7 introduced in April 2019)
  subroutine convert_6_7(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    integer :: stat

    call hsd_rename_child(root, "hamiltonian/dftb/orbitalresolvedscc", &
        &"shellresolvedscc", stat, case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, "Keyword converted to 'ShellResolvedSCC'.")
    call handleD3Defaults(root)

    call hsd_rename_child(root, "hamiltonian/dftb/eigensolver", &
        &"solver", stat, case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, "Keyword renamed to 'Solver'.")

  end subroutine convert_6_7


  !> Converts input from version 7 to 8. (Version 8 introduced in October 2019)
  subroutine convert_7_8(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table), pointer :: ch1, ch2
    logical :: tVal
    integer :: stat
    type(hsd_table), pointer :: pTaskType
    character(len=:), allocatable :: buffer

    call hsd_rename_child(root, "analysis/eigenvectorsastxt", &
        &"eigenvectorsastext", stat, case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, &
        &"Keyword converted to 'EigenvectorsAsText'.")

    call hsd_get_table(root, "transport", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get_table(ch1, "Task", ch2)
      if (.not. associated(ch2)) then
        call hsd_set(ch1, "readBinaryContact", .false.)
        call hsd_get_table(ch1, "readBinaryContact", ch2, stat, auto_wrap=.true.)
      else
        call hsd_get_choice(ch2, "", buffer, pTaskType, stat)
        call hsd_get_name(pTaskType, buffer, "#text")
        select case (buffer)
        case ("contacthamiltonian")
          call hsd_set(ch1, "writeBinaryContact", .false.)
          call hsd_get_table(ch1, "writeBinaryContact", ch2, stat, auto_wrap=.true.)
        case ("uploadcontacts")
          call hsd_set(ch1, "readBinaryContact", .false.)
          call hsd_get_table(ch1, "readBinaryContact", ch2, stat, auto_wrap=.true.)
        end select
      end if
    end if

    call hsd_get_table(root, "ParserOptions", ch1)
    if (associated(ch1) .and. hsd_has_child(ch1, "WriteXMLInput", .true.)) then
      call hsd_get(ch1, "WriteXMLInput", tVal, stat=stat)
      if (stat == HSD_STAT_OK .and. tVal) then
        call dftbp_warning(ch1, "Sorry, XML export of the dftb_in.hsd is not supported any more&
            & so is removed")
      else
        call dftbp_warning(ch1, "XML export option is removed.")
      end if
      call hsd_remove_child(ch1, "WriteXMLInput", stat, case_insensitive=.true.)
    end if

  end subroutine convert_7_8


  !> Converts input from version 8 to 9. (Version 9 introduced in August 2020)
  subroutine convert_8_9(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table), pointer :: ch1, ch2
    logical :: tVal1, tVal2
    integer :: stat

    ! If this is an electron dynamics restart, then remove keywords for the (un-needed) ground state
    ! calculation (unless the eigenvectors are required)
    call hsd_get_table(root, "electrondynamics/restart", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal1, stat=stat)
      if (associated(ch1)) ch1%processed = .false.
      tVal2 = .false.
      ! Population projection requires eigenvectors, which are not currently stored in the restart
      ! file.
      call hsd_get_table(root, "electrondynamics/populations", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) then
        call hsd_get(ch2, "#text", tVal2, stat=stat)
        if (associated(ch2)) ch2%processed = .false.
      end if
      if (tVal1 .and. .not.tVal2) then
        call hsd_get_table(root, "Hamiltonian/DFTB/Filling", ch1)
        if (associated(ch1)) then
          call dftbp_warning(ch1, "Restarted electronDynamics does not require Filling{}&
              & settings unless projected onto ground state")
          ch1 => null()
        end if
        call hsd_get_table(root, "Analysis", ch1)
        if (associated(ch1)) then
          call dftbp_warning(ch1, "Restarted electronDynamics does not use the Analysis{} block")
          ch1 => null()
        end if
      end if
    end if

    call hsd_get_table(root, "Driver/lBFGS", ch1)
    if (associated(ch1)) then
      call hsd_set(ch1, "LineSearch", .true.)
      call hsd_get_table(ch1, "LineSearch", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) ch2%processed = .false.
      call dftbp_warning(ch2, "Set 'LineSearch = Yes'")
      call hsd_set(ch1, "oldLineSearch", .true.)
      call hsd_get_table(ch1, "oldLineSearch", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) ch2%processed = .false.
      call dftbp_warning(ch2, "Set 'oldLineSearch = Yes'")
    end if

  end subroutine convert_8_9


  !> Converts input from version 9 to 10. (Version 10 introduced in November 2021)
  subroutine convert_9_10(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table), pointer :: ch1, ch2, ch3, ch4
    logical :: tVal1, tVal2
    integer :: stat

    call hsd_get_table(root, "ExcitedState/Casida", ch1)
    if (associated(ch1)) then
      call hsd_get_or_set(ch1, "WriteStatusArnoldi", tVal1, .false.)
      call hsd_get_table(ch1, "WriteStatusArnoldi", ch2, stat, auto_wrap=.true.)
      call hsd_remove_child(ch1, ch2%name)
      ch2 => null()
      call hsd_get_or_set(ch1, "TestArnoldi", tVal2, .false.)
      call hsd_get_table(ch1, "TestArnoldi", ch2, stat, auto_wrap=.true.)
      call hsd_remove_child(ch1, ch2%name)
      ch2 => null()
      call dftbp_warning(ch1, "Keyword moved to Diagonaliser block.")
      if (associated(ch1)) ch1%processed = .false.
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="diagonaliser")
        call ch1%add_child(tmpTbl)
      end block
      call hsd_get_table(ch1, "Diagonaliser", ch2, stat)
      if (associated(ch2)) ch2%processed = .false.
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="arpack")
        call ch2%add_child(tmpTbl)
      end block
      call hsd_get_table(ch2, "Arpack", ch3, stat)
      if (associated(ch3)) ch3%processed = .false.
      call hsd_set(ch3, "WriteStatusArnoldi", tVal1)
      call hsd_get_table(ch3, "WriteStatusArnoldi", ch4, stat, auto_wrap=.true.)
      if (associated(ch4)) ch4%processed = .false.
      call hsd_set(ch3, "TestArnoldi", tVal2)
      call hsd_get_table(ch3, "TestArnoldi", ch4, stat, auto_wrap=.true.)
      if (associated(ch4)) ch4%processed = .false.
    end if

    ! move ConvergentSccOnly and ConvergentForces into a common keyword

    call hsd_get_table(root, "hamiltonian/dispersion/ts/convergentscconly", ch1, stat, &
        &auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal1, stat=stat)
      call dftbp_warning(ch1, "Keyword Moved to Hamiltonian {}.")
      call hsd_remove_child(root, "hamiltonian/dispersion/ts/convergentscconly")
      ch1 => null()
      call hsd_get_table(root, "Hamiltonian", ch1)
      call hsd_set(ch1, "ConvergentSCCOnly", tVal1)
      call hsd_get_table(ch1, "ConvergentSCCOnly", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) ch2%processed = .false.
    end if

    call hsd_get_table(root, "hamiltonian/dispersion/mbd/convergentscconly", ch1, stat, &
        &auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal1, stat=stat)
      call dftbp_warning(ch1, "Keyword Moved to Hamiltonian {}.")
      call hsd_remove_child(root, "hamiltonian/dispersion/mbd/convergentscconly")
      ch1 => null()
      call hsd_get_table(root, "hamiltonian/convergentscconly", ch3, stat, auto_wrap=.true.)
      if (associated(ch3)) then
        call dftbp_error(ch3, "ConvergentSCCOnly already present.")
      end if
      call hsd_get_table(root, "Hamiltonian", ch1)
      call hsd_set(ch1, "ConvergentSCCOnly", tVal1)
      call hsd_get_table(ch1, "ConvergentSCCOnly", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) ch2%processed = .false.
    end if

    call hsd_get_table(root, "driver/conjugategradient/convergentforcesonly", ch1, stat, &
        &auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal1, stat=stat)
      call dftbp_warning(ch1, "Keyword Moved to Hamiltonian {}.")
      call hsd_remove_child(root, "driver/conjugategradient/convergentforcesonly")
      ch1 => null()
      call hsd_get_table(root, "hamiltonian/convergentscconly", ch3, stat, auto_wrap=.true.)
      if (associated(ch3)) then
        call dftbp_error(ch3, "ConvergentSCCOnly already present.")
      end if
      call hsd_get_table(root, "Hamiltonian", ch1)
      call hsd_set(ch1, "ConvergentSCCOnly", tVal1)
      call hsd_get_table(ch1, "ConvergentSCCOnly", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) ch2%processed = .false.
    end if

    call hsd_get_table(root, "driver/velocityverlet/convergentforcesonly", ch1, stat, &
        &auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal1, stat=stat)
      call dftbp_warning(ch1, "Keyword Moved to Hamiltonian {}.")
      call hsd_remove_child(root, "driver/velocityverlet/convergentforcesonly")
      ch1 => null()
      call hsd_get_table(root, "hamiltonian/convergentscconly", ch3, stat, auto_wrap=.true.)
      if (associated(ch3)) then
        call dftbp_error(ch3, "ConvergentSCCOnly already present.")
      end if
      call hsd_get_table(root, "Hamiltonian", ch1)
      call hsd_set(ch1, "ConvergentSCCOnly", tVal1)
      call hsd_get_table(ch1, "ConvergentSCCOnly", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) ch2%processed = .false.
    end if

    call hsd_get_table(root, "driver/steepestdescent/convergentforcesonly", ch1, stat, &
        &auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal1, stat=stat)
      call dftbp_warning(ch1, "Keyword Moved to Hamiltonian {}.")
      call hsd_remove_child(root, "driver/steepestdescent/convergentforcesonly")
      ch1 => null()
      call hsd_get_table(root, "hamiltonian/convergentscconly", ch3, stat, auto_wrap=.true.)
      if (associated(ch3)) then
        call dftbp_error(ch3, "ConvergentSCCOnly already present.")
      end if
      call hsd_get_table(root, "Hamiltonian", ch1)
      call hsd_set(ch1, "ConvergentSCCOnly", tVal1)
      call hsd_get_table(ch1, "ConvergentSCCOnly", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) ch2%processed = .false.
    end if

    call hsd_get_table(root, "driver/gdiis/convergentforcesonly", ch1, stat, &
        &auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal1, stat=stat)
      call dftbp_warning(ch1, "Keyword Moved to Hamiltonian {}.")
      call hsd_remove_child(root, "driver/gdiis/convergentforcesonly")
      ch1 => null()
      call hsd_get_table(root, "hamiltonian/convergentscconly", ch3, stat, auto_wrap=.true.)
      if (associated(ch3)) then
        call dftbp_error(ch3, "ConvergentSCCOnly already present.")
      end if
      call hsd_get_table(root, "Hamiltonian", ch1)
      call hsd_set(ch1, "ConvergentSCCOnly", tVal1)
      call hsd_get_table(ch1, "ConvergentSCCOnly", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) ch2%processed = .false.
    end if

    call hsd_get_table(root, "driver/lbfgs/convergentforcesonly", ch1, stat, &
        &auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal1, stat=stat)
      call dftbp_warning(ch1, "Keyword Moved to Hamiltonian {}.")
      call hsd_remove_child(root, "driver/lbfgs/convergentforcesonly")
      ch1 => null()
      call hsd_get_table(root, "hamiltonian/convergentscconly", ch3, stat, auto_wrap=.true.)
      if (associated(ch3)) then
        call dftbp_error(ch3, "ConvergentSCCOnly already present.")
      end if
      call hsd_get_table(root, "Hamiltonian", ch1)
      call hsd_set(ch1, "ConvergentSCCOnly", tVal1)
      call hsd_get_table(ch1, "ConvergentSCCOnly", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) ch2%processed = .false.
    end if

    call hsd_get_table(root, "driver/fire/convergentforcesonly", ch1, stat, &
        &auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal1, stat=stat)
      call dftbp_warning(ch1, "Keyword Moved to Hamiltonian {}.")
      call hsd_remove_child(root, "driver/fire/convergentforcesonly")
      ch1 => null()
      call hsd_get_table(root, "hamiltonian/convergentscconly", ch3, stat, auto_wrap=.true.)
      if (associated(ch3)) then
        call dftbp_error(ch3, "ConvergentSCCOnly already present.")
      end if
      call hsd_get_table(root, "Hamiltonian", ch1)
      call hsd_set(ch1, "ConvergentSCCOnly", tVal1)
      call hsd_get_table(ch1, "ConvergentSCCOnly", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) ch2%processed = .false.
    end if

    call hsd_get_table(root, "driver/secondderivatives/convergentforcesonly", ch1, stat, &
        &auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal1, stat=stat)
      call dftbp_warning(ch1, "Keyword Moved to Hamiltonian {}.")
      call hsd_remove_child(root, "driver/secondderivatives/convergentforcesonly")
      ch1 => null()
      call hsd_get_table(root, "hamiltonian/convergentscconly", ch3, stat, auto_wrap=.true.)
      if (associated(ch3)) then
        call dftbp_error(ch3, "ConvergentSCCOnly already present.")
      end if
      call hsd_get_table(root, "Hamiltonian", ch1)
      call hsd_set(ch1, "ConvergentSCCOnly", tVal1)
      call hsd_get_table(ch1, "ConvergentSCCOnly", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) ch2%processed = .false.
    end if

    call hsd_get_table(root, "driver/socket/convergentforcesonly", ch1, stat, &
        &auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", tVal1, stat=stat)
      call dftbp_warning(ch1, "Keyword Moved to Hamiltonian {}.")
      call hsd_remove_child(root, "driver/socket/convergentforcesonly")
      ch1 => null()
      call hsd_get_table(root, "hamiltonian/convergentscconly", ch3, stat, auto_wrap=.true.)
      if (associated(ch3)) then
        call dftbp_error(ch3, "ConvergentSCCOnly already present.")
      end if
      call hsd_get_table(root, "Hamiltonian", ch1)
      call hsd_set(ch1, "ConvergentSCCOnly", tVal1)
      call hsd_get_table(ch1, "ConvergentSCCOnly", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) ch2%processed = .false.
    end if

  end subroutine convert_9_10


  !> Converts input from version 10 to 11. (Version 11 introduced in April 2022)
  subroutine convert_10_11(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table), pointer :: ch1, ch2
    integer :: stat

    call hsd_get_table(root, "Hamiltonian/DFTB/Solvation/GeneralizedBorn", ch1)
    if (associated(ch1)) then
      call dftbp_warning(ch1, "Set solvated field scaling (RescaleSolvatedFields) to No.")
      call hsd_set(ch1, "RescaleSolvatedFields", .false.)
      call hsd_get_table(ch1, "RescaleSolvatedFields", ch2, stat, auto_wrap=.true.)
    end if

    call hsd_get_table(root, "Hamiltonian/DFTB/Solvation/Cosmo", ch1)
    if (associated(ch1)) then
      call dftbp_warning(ch1, "Set solvated field scaling (RescaleSolvatedFields) to No.")
      call hsd_set(ch1, "RescaleSolvatedFields", .false.)
      call hsd_get_table(ch1, "RescaleSolvatedFields", ch2, stat, auto_wrap=.true.)
    end if

    call hsd_get_table(root, "Hamiltonian/DFTB/Solvation/Sasa", ch1)
    if (associated(ch1)) then
      call dftbp_warning(ch1, "Set solvated field scaling (RescaleSolvatedFields) to No.")
      call hsd_set(ch1, "RescaleSolvatedFields", .false.)
      call hsd_get_table(ch1, "RescaleSolvatedFields", ch2, stat, auto_wrap=.true.)
    end if

  end subroutine convert_10_11


  !> Converts input from version 11 to 12. (Version 12 introduced in June 2022)
  subroutine convert_11_12(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table), pointer :: ch1, ch2
    character(len=:), allocatable :: buffer
    integer :: stat

    call hsd_get_table(root, "Transport", ch1)
    if (associated(ch1)) then
      call hsd_get_table(root, "Transport/Task", ch2, stat, auto_wrap=.true.)
      if (.not. associated(ch2)) then
        block
          type(hsd_table) :: tmpTbl, defChild
          call new_table(tmpTbl, name="task")
          call new_table(defChild, name="uploadcontacts")
          call tmpTbl%add_child(defChild)
          call ch1%add_child(tmpTbl)
        end block
        call hsd_get_table(root, "Transport/Task", ch2, stat)
      end if
      call hsd_get_choice(ch2, "", buffer, ch1, stat)
      call hsd_get_name(ch1, buffer, "#text")
      if (buffer /= "contacthamiltonian") then
      #:for LABEL in [("xtb"), ("dftb")]
        call hsd_get_table(root, "hamiltonian/${LABEL}$/charge", ch1, stat, auto_wrap=.true.)
        if (associated(ch1)) then
          if (associated(ch1)) ch1%processed = .false.
          call dftbp_warning(ch1, "Device region charge cannot be set if contacts are present.")
          ch1 => null()
        end if
      #:endfor
      end if
    end if

  end subroutine convert_11_12


  !> Converts input from version 12 to 13. (Version 13 introduced in February 2023)
  subroutine convert_12_13(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table), pointer :: ch1, ch2
    integer :: maxIter
    logical :: isPerturb, isConvRequired
    real(dp) :: sccTol
    integer :: stat

    call hsd_get_table(root, "Analysis/Polarisability", ch1)
    isPerturb = associated(ch1)
    if (.not.isPerturb) then
      call hsd_get_table(root, "Analysis/ResponseKernel", ch1)
      isPerturb = associated(ch1)
    end if

    if (isPerturb) then

      call hsd_rename_child(root, "analysis/eta", "perturbeta", stat, case_insensitive=.true.)
      if (stat == HSD_STAT_OK) call dftbp_warning(root, &
          &"Keyword renamed to 'PerturbEta'.")

      call hsd_rename_child(root, "analysis/degeneracytolerance", "pertubdegentol", stat, &
          &case_insensitive=.true.)
      if (stat == HSD_STAT_OK) call dftbp_warning(root, &
          &"Keyword renamed to 'PertubDegenTol'.")

      call hsd_get_table(root, "hamiltonian/dftb/maxscciterations", ch1, stat, auto_wrap=.true.)
      if (associated(ch1)) then
        call hsd_get(ch1, "#text", maxIter, stat=stat)
        call hsd_get_table(root, "Analysis", ch1)
        call hsd_set(ch1, "MaxPerturbIter", maxIter)
        call hsd_get_table(ch1, "MaxPerturbIter", ch2, stat, auto_wrap=.true.)
        if (associated(ch2)) ch2%processed = .false.
      end if

      call hsd_get_table(root, "hamiltonian/dftb/convergentscconly", ch1, stat, &
          &auto_wrap=.true.)
      if (associated(ch1)) then
        call hsd_get(ch1, "#text", isConvRequired, stat=stat)
        call hsd_get_table(root, "Analysis", ch1)
        call hsd_set(ch1, "ConvergedPerturb", isConvRequired)
        call hsd_get_table(ch1, "ConvergedPerturb", ch2, stat, auto_wrap=.true.)
        if (associated(ch2)) ch2%processed = .false.
      end if

      call hsd_get_table(root, "hamiltonian/dftb/scctolerance", ch1, stat, auto_wrap=.true.)
      if (associated(ch1)) then
        call hsd_get(ch1, "#text", sccTol, stat=stat)
        call hsd_get_table(root, "Analysis", ch1)
        call hsd_set(ch1, "PerturbSccTol", sccTol)
        call hsd_get_table(ch1, "PerturbSccTol", ch2, stat, auto_wrap=.true.)
        if (associated(ch2)) ch2%processed = .false.
      end if

    end if

  end subroutine convert_12_13


  !> Converts input from version 13 to 14. (Version 14 introduced in August 2023)
  subroutine convert_13_14(root)

    !> Root tag of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table), pointer :: ch1, ch2, ch3, par, hybridAlgorithm
    character(len=:), allocatable :: buffer
    logical :: isScc, isNoneAlgorithm
    integer :: iOrder
    character(lc) :: strTmp
    real(dp) :: rTol
    integer :: stat

    call hsd_rename_child(root, "analysis/calculateforces", "printforces", stat, &
        &case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, &
        &"Keyword renamed to 'PrintForces'.")

    call hsd_rename_child(root, "hamiltonian/dftb/rangeseparated", "hybrid", stat, &
        &case_insensitive=.true.)
    if (stat == HSD_STAT_OK) call dftbp_warning(root, &
        &"'Hamiltonian/DFTB/Rangeseparated' block renamed to 'Hamiltonian/DFTB/Hybrid'.")

    call hsd_get_table(root, "Hamiltonian/DFTB/Filling/MethfesselPaxton", ch1)
    if (.not.associated(ch1)) then
      call hsd_get_table(root, "Hamiltonian/xTB/Filling/MethfesselPaxton", ch1)
    end if
    if (associated(ch1)) then
      call hsd_get_table(ch1, "order", ch2, stat, auto_wrap=.true.)
      if (associated(ch2)) then
        call hsd_get(ch2, "#text", iOrder, stat=stat)
        ch3 => ch2
        if (iOrder > 1) then
          write(strtmp,"(A,I0,A,I0,A)")"Older Methfessel-Paxton requested order of ", iOrder,&
              & " is now equivalent to ", (iOrder -1), " from parser version 14."
          call dftbp_warning(ch2, strTmp)
        elseif (iOrder == 1) then
          write(strtmp,"(A)")"Older Methfessel-Paxton requested order of 1 is now&
              & equivalent to Gaussian smearing (order 0) from parser version 14." // newline //&
              & "   Please test your calculation carefully, due to a (corrected) error for array&
              & bounds in this case. "
          call dftbp_warning(ch2, strTmp)
        else
          write(strtmp,"(A,I0,A,I0,A)")"Older Methfessel-Paxton requested order of ", iOrder,&
              & " is now equivalent to a negative order of ", (iOrder -1), " from parser version 14&
              & and is incorrect."
          call dftbp_error(ch2, strTmp)
        end if
        call hsd_remove_child(ch1, ch2%name)
        ch2 => null()
        iOrder = iOrder - 1
        call hsd_set(ch1, "Order", iOrder)
        call hsd_get_table(ch1, "Order", ch2, stat, auto_wrap=.true.)
      else
        call dftbp_warning(ch1, "The default (i.e. unspecified) Methfessel-Paxton order in old&
            & code versions was equivalent to Order=1." // newline //&
            & "   This order is now the current default order from parser version 14.")
      end if

    end if

    call hsd_get_table(root, "Hamiltonian/DFTB/Hybrid/LC", ch1)
    if (associated(ch1)) then
      call hsd_get_table(ch1, "Screening", ch2)
      if (.not. associated(ch2)) then
        block
          type(hsd_table) :: tmpTbl
          call new_table(tmpTbl, name="screening")
          call ch1%add_child(tmpTbl)
        end block
        call hsd_get_table(ch1, "Screening", ch2, stat)
        block
          type(hsd_table) :: tmpTbl
          call new_table(tmpTbl, name="thresholded")
          call ch2%add_child(tmpTbl)
        end block
        call hsd_get_table(ch2, "Thresholded", ch3, stat)
        if (associated(ch1)) ch1%processed = .false.
        if (associated(ch2)) ch2%processed = .false.
      end if
    end if

    call hsd_get_table(root, "analysis", par, stat, auto_wrap=.true.)
    call hsd_get_table(root, "analysis/pertubdegentol", ch1, stat, auto_wrap=.true.)
    if (associated(ch1)) then
      call hsd_get(ch1, "#text", rTol, stat=stat)
      if (rTol < 1.0_dp) then
        call dftbp_error(ch1, "Perturbation degeneracy tolerance must be above 1x")
      end if
      call hsd_remove_child(root, "analysis/pertubdegentol")
      ch1 => null()
      call hsd_set(par, "PerturbDegenTol", rTol * epsilon(0.0_dp))
      call hsd_get_table(par, "PerturbDegenTol", ch1, stat, auto_wrap=.true.)
      call dftbp_warning(par, "Keyword renamed to 'PerturbDegenTol'.")
    end if

  end subroutine convert_13_14


  !> Update values in the DftD3 block to match behaviour of v6 parser
  subroutine handleD3Defaults(root)

    !> Root node of the HSD-tree
    type(hsd_table), pointer :: root

    type(hsd_table), pointer :: pD3, pDampMethod, pChild
    character(len=:), allocatable :: buffer
    integer :: stat

    call hsd_get_table(root, "Hamiltonian/DFTB/Dispersion/DftD3", pD3)
    if (.not. associated(pD3)) then
      return
    end if

    call useDftb3Default(pD3, "s6", 1.0_dp)
    call useDftb3Default(pD3, "s8", 0.5883_dp)

    call hsd_get_table(pD3, "Damping", pChild, stat, auto_wrap=.true.)
    if (.not. associated(pChild)) then
      block
        type(hsd_table) :: tmpTbl, defChild
        call new_table(tmpTbl, name="damping")
        call new_table(defChild, name="beckejohnson")
        call tmpTbl%add_child(defChild)
        call pD3%add_child(tmpTbl)
      end block
      call hsd_get_table(pD3, "Damping", pChild, stat)
    end if
    if (associated(pChild)) pChild%processed = .false.
    call hsd_get_choice(pChild, "", buffer, pDampMethod, stat)
    if (associated(pDampMethod)) pDampMethod%processed = .false.
    call hsd_get_name(pDampMethod, buffer, "#text")

    select case (buffer)
    case ("beckejohnson")
      call useDftb3Default(pDampMethod, "a1", 0.5719_dp)
      call useDftb3Default(pDampMethod, "a2", 3.6017_dp)
    end select

  end subroutine handleD3Defaults


  !> Helper routine to update values in the DftD3 block to match behaviour of v6 parser
  subroutine useDftb3Default(root, option, default)

    !> Root node of the HSD-tree
    type(hsd_table), pointer, intent(in) :: root

    !> Name of option inside the DftD3 block
    character(*), intent(in) :: option

    !> Default value to set
    real(dp), intent(in) :: default

    type(hsd_table), pointer :: pChild
    integer :: stat

    call hsd_get_table(root, option, pChild, stat, auto_wrap=.true.)
    if (.not. associated(pChild)) then
      call hsd_set(root, option, default)
      call hsd_get_table(root, option, pChild, stat, auto_wrap=.true.)
      call dftbp_warning(pChild, "Using DFTB3 optimised default value for parameter " // option)
    end if
    if (associated(pChild)) pChild%processed = .false.

  end subroutine useDftb3Default


end module dftbp_dftbplus_oldcompat
