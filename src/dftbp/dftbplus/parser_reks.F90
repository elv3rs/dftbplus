!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Reads in REKS (Real-time Electron Kohn-Sham) related settings from the HSD input.
module dftbp_dftbplus_parser_reks
  use dftbp_common_accuracy, only : dp
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_io_charmanip, only : i2c, tolower
  use hsd_data, only : hsd_table
  use dftbp_io_hsdutils, only : dftbp_error, getNodeName, getNodeHSDName, getNodeName2,&
      & hasInlineData
  use hsd, only : hsd_get, hsd_get_or_set, hsd_get_table, hsd_get_choice, &
      & HSD_STAT_OK, new_table
  use dftbp_reks_reks, only : reksTypes
  use dftbp_type_typegeometry, only : TGeometry
  implicit none

  private
  public :: readReks

contains


  !> Reads the REKS block
  subroutine readReks(node, dummy, ctrl, geo)

    !> Node to parse
    type(hsd_table), pointer, intent(in) :: node

    !> Node to parse
    type(hsd_table), pointer, intent(in) :: dummy

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    character(len=:), allocatable :: buffer

    ! SSR(2,2) or SSR(4,4) stuff
    call getNodeName(dummy, buffer)

    select case (buffer)
    case ("none")
      ctrl%reksInp%reksAlg = reksTypes%noReks
    case ("ssr22")
      ctrl%reksInp%reksAlg = reksTypes%ssr22
      call readSSR22(dummy, ctrl, geo)
    case ("ssr44")
      ctrl%reksInp%reksAlg = reksTypes%ssr44
      call dftbp_error(node, "SSR(4,4) is not implemented yet.")
    case default
      call getNodeHSDName(dummy, buffer)
      call dftbp_error(node, "Invalid Algorithm '" // buffer // "'")
    end select

  end subroutine readReks


  !> Reads the SSR(2,2) block
  subroutine readSSR22(node, ctrl, geo)

    !> Node to parse
    type(hsd_table), pointer, intent(in) :: node

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    type(hsd_table), pointer :: child1, value2, child2
    character(len=:), allocatable :: buffer2
    character(len=:), allocatable :: tmpFunc(:)
    integer :: ii, nFunc, stat
    logical :: tFunc = .true.


    !> Read 'Energy' block
    call hsd_get_table(node, "Energy", child1, stat, auto_wrap=.true.)
    if (.not. associated(child1)) call dftbp_error(node, "Missing required block: 'Energy'")

    !> Read 'Functional' block in 'Energy' block
    call hsd_get(child1, "Functional", tmpFunc, stat=stat)
    if (stat /= 0) call dftbp_error(child1, "Error reading Functional")

    !> Decide the energy functionals to be included in SA-REKS(2,2)
    nFunc = size(tmpFunc, dim=1)
    if (nFunc == 1) then
      if (trim(tmpFunc(1)) == "PPS") then
        !> Minimized energy functional : PPS
        ctrl%reksInp%Efunction = 1
      else
        tFunc = .false.
      end if
    else if (nFunc == 2) then
      if (trim(tmpFunc(1)) == "PPS" .and. trim(tmpFunc(2)) == "OSS") then
        !> Minimized energy functional : (PPS+OSS)/2
        ctrl%reksInp%Efunction = 2
      else
        tFunc = .false.
      end if
    else
      tFunc = .false.
    end if

    if (.not. tFunc) then
      write(stdOut,'(A)',advance="no") "Current Functional : "
      do ii = 1, nFunc
        if (ii == nFunc) then
          write(stdOut,'(A)') "'" // trim(tmpFunc(ii)) // "'"
        else
          write(stdOut,'(A)',advance="no") "'" // trim(tmpFunc(ii)) // "' "
        end if
      end do
      call dftbp_error(child1, "Invalid Functional")
    end if

    !> Decide the energy states in SA-REKS
    !> If true, it includes all possible states in current active space
    !> If false, it includes the states used in minimized energy functional
    call hsd_get_or_set(child1, "IncludeAllStates", ctrl%reksInp%tAllStates, .false.)
    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    call hsd_get_or_set(child1, "StateInteractions", ctrl%reksInp%tSSR, .false.)


    !> Target SSR state
    call hsd_get_or_set(node, "TargetState", ctrl%reksInp%rstate, 1)
    !> Target microstate
    call hsd_get_or_set(node, "TargetMicrostate", ctrl%reksInp%Lstate, 0)

    !> Read initial guess for eigenvectors in REKS
    !> If true, initial eigenvectors are obtained from 'eigenvec.bin'
    !> If false, initial eigenvectors are obtained from diagonalisation of H0
    call hsd_get_or_set(node, "ReadEigenvectors", ctrl%reksInp%tReadMO, .false.)
    !> Maximum iteration used in FON optimisation
    call hsd_get_or_set(node, "FonMaxIter", ctrl%reksInp%FonMaxIter, 20)
    !> Shift value in SCC cycle
    call hsd_get_or_set(node, "Shift", ctrl%reksInp%shift, 0.3_dp)

    !> Read "SpinTuning" block with 'nType' elements
    call readSpinTuning(node, ctrl, geo%nSpecies)

    !> Calculate transition dipole moments
    call hsd_get_or_set(node, "TransitionDipole", ctrl%reksInp%tTDP, .false.)


    !> Read 'Gradient' block
    !> Algorithms to calculate analytical gradients
    call hsd_get_table(node, "Gradient", child2, stat, auto_wrap=.true.)
    if (.not. associated(child2)) then
      block
        type(hsd_table) :: defContainer, defChild
        call new_table(defContainer, name="gradient")
        call new_table(defChild, name="ConjugateGradient")
        call defContainer%add_child(defChild)
        call node%add_child(defContainer)
      end block
      call hsd_get_table(node, "Gradient", child2, stat, auto_wrap=.true.)
    end if
    call hsd_get_choice(child2, "", buffer2, value2, stat)
    if (.not. associated(value2) .and. len_trim(buffer2) == 0) buffer2 = "conjugategradient"

    select case (tolower(buffer2))
    case ("conjugategradient")
      !> Maximum iteration used in calculation of gradient with PCG and CG
      call hsd_get_or_set(value2, "CGmaxIter", ctrl%reksInp%CGmaxIter, 20)
      !> Tolerance used in calculation of gradient with PCG and CG
      call hsd_get_or_set(value2, "Tolerance", ctrl%reksInp%Glimit, 1.0E-8_dp)
      !> Use preconditioner for conjugate gradient algorithm
      call hsd_get_or_set(value2, "Preconditioner", ctrl%reksInp%tPrecond, .false.)
      !> Save 'A' and 'Hxc' to memory in gradient calculation
      call hsd_get_or_set(value2, "SaveMemory", ctrl%reksInp%tSaveMem, .false.)
      if (ctrl%reksInp%tPrecond) then
        !> 1: preconditioned conjugate gradient (PCG)
        ctrl%reksInp%Glevel = 1
      else
        !> 2: conjugate gradient (CG)
        ctrl%reksInp%Glevel = 2
      end if
    case ("direct")
      !> 3: direct inverse-matrix multiplication
      ctrl%reksInp%Glevel = 3
    case default
      call dftbp_error(child2, "Invalid Algorithm '" // buffer2 // "'")
    end select

    !> Calculate relaxed density of SSR or SA-REKS state
    call hsd_get_or_set(node, "RelaxedDensity", ctrl%reksInp%tRD, .false.)
    !> Calculate nonadiabatic coupling vectors
    call hsd_get_or_set(node, "NonAdiabaticCoupling", ctrl%reksInp%tNAC, .false.)

    !> Print level in standard output file
    call hsd_get_or_set(node, "VerbosityLevel", ctrl%reksInp%Plevel, 1)

  end subroutine readSSR22


  !> Reads SpinTuning block in REKS input
  subroutine readSpinTuning(node, ctrl, nType)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Number of types for atoms
    integer, intent(in) :: nType

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer, modifier
    integer :: nAtom, iType, nrows, ncols, stat
    real(dp), allocatable :: tmpTuning(:)

    call hsd_get_table(node, "SpinTuning", child, stat, auto_wrap=.true.)
    if (.not. associated(child) .or. .not. hasInlineData(child)) then
      ! no 'SpinTuning' block in REKS input
      allocate(ctrl%reksInp%Tuning(nType))
      do iType = 1, nType
        ctrl%reksInp%Tuning(iType) = 1.0_dp
      end do
    else
      ! 'SpinTuning' block in REKS input
      call hsd_get(node, "SpinTuning", tmpTuning, stat=stat)
      if (stat /= 0) call dftbp_error(node, "Error reading SpinTuning")
      if (size(tmpTuning) /= nType) then
        call dftbp_error(node, "Incorrect number of 'SpinTuning' block: " &
            & // i2c(size(tmpTuning)) // " supplied, " &
            & // i2c(nType) // " required.")
      end if
      allocate(ctrl%reksInp%Tuning(nType))
      ctrl%reksInp%Tuning(:) = tmpTuning(:)
    end if

  end subroutine readSpinTuning


end module dftbp_dftbplus_parser_reks
