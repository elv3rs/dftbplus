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
  use dftbp_common_accuracy, only : dp, sc
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_io_charmanip, only : i2c
  use dftbp_io_hsdcompat, only : hsd_table, detailedError, getChild, getChildValue, &
      & getNodeName, getNodeHSDName, getNodeName2, hasInlineData
  use dftbp_reks_reks, only : reksTypes
  use dftbp_type_linkedlist, only : asArray, destruct, init, len, TListRealR1, TListString
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
      call detailedError(node, "SSR(4,4) is not implemented yet.")
    case default
      call getNodeHSDName(dummy, buffer)
      call detailedError(node, "Invalid Algorithm '" // buffer // "'")
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
    type(TListString) :: strBuffer
    character(len=:), allocatable :: buffer2
    character(sc), allocatable :: tmpFunc(:)
    integer :: ii, nFunc
    logical :: tFunc = .true.


    !> Read 'Energy' block
    call getChild(node, "Energy", child=child1)

    !> Read 'Functional' block in 'Energy' block
    call init(strBuffer)
    call getChildValue(child1, "Functional", strBuffer)
    allocate(tmpFunc(len(strBuffer)))
    call asArray(strBuffer, tmpFunc)
    call destruct(strBuffer)

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
      call detailedError(child1, "Invalid Functional")
    end if

    !> Decide the energy states in SA-REKS
    !> If true, it includes all possible states in current active space
    !> If false, it includes the states used in minimized energy functional
    call getChildValue(child1, "IncludeAllStates", ctrl%reksInp%tAllStates, default=.false.)
    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    call getChildValue(child1, "StateInteractions", ctrl%reksInp%tSSR, default=.false.)


    !> Target SSR state
    call getChildValue(node, "TargetState", ctrl%reksInp%rstate, default=1)
    !> Target microstate
    call getChildValue(node, "TargetMicrostate", ctrl%reksInp%Lstate, default=0)

    !> Read initial guess for eigenvectors in REKS
    !> If true, initial eigenvectors are obtained from 'eigenvec.bin'
    !> If false, initial eigenvectors are obtained from diagonalisation of H0
    call getChildValue(node, "ReadEigenvectors", ctrl%reksInp%tReadMO, default=.false.)
    !> Maximum iteration used in FON optimisation
    call getChildValue(node, "FonMaxIter", ctrl%reksInp%FonMaxIter, default=20)
    !> Shift value in SCC cycle
    call getChildValue(node, "Shift", ctrl%reksInp%shift, default=0.3_dp)

    !> Read "SpinTuning" block with 'nType' elements
    call readSpinTuning(node, ctrl, geo%nSpecies)

    !> Calculate transition dipole moments
    call getChildValue(node, "TransitionDipole", ctrl%reksInp%tTDP, default=.false.)


    !> Read 'Gradient' block
    !> Algorithms to calculate analytical gradients
    call getChildValue(node, "Gradient", value2, "ConjugateGradient", child=child2)
    call getNodeName(value2, buffer2)

    select case (buffer2)
    case ("conjugategradient")
      !> Maximum iteration used in calculation of gradient with PCG and CG
      call getChildValue(value2, "CGmaxIter", ctrl%reksInp%CGmaxIter, default=20)
      !> Tolerance used in calculation of gradient with PCG and CG
      call getChildValue(value2, "Tolerance", ctrl%reksInp%Glimit, default=1.0E-8_dp)
      !> Use preconditioner for conjugate gradient algorithm
      call getChildValue(value2, "Preconditioner", ctrl%reksInp%tPrecond, default=.false.)
      !> Save 'A' and 'Hxc' to memory in gradient calculation
      call getChildValue(value2, "SaveMemory", ctrl%reksInp%tSaveMem, default=.false.)
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
      call getNodeHSDName(value2, buffer2)
      call detailedError(child2, "Invalid Algorithm '" // buffer2 // "'")
    end select

    !> Calculate relaxed density of SSR or SA-REKS state
    call getChildValue(node, "RelaxedDensity", ctrl%reksInp%tRD, default=.false.)
    !> Calculate nonadiabatic coupling vectors
    call getChildValue(node, "NonAdiabaticCoupling", ctrl%reksInp%tNAC, default=.false.)

    !> Print level in standard output file
    call getChildValue(node, "VerbosityLevel", ctrl%reksInp%Plevel, default=1)

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
    type(TListRealR1) :: realBuffer
    integer :: nAtom, iType
    real(dp), allocatable :: tmpTuning(:,:)

    call getChildValue(node, "SpinTuning", value1, "", child=child, modifier=modifier,&
        & allowEmptyValue=.true.)
    call getNodeName2(value1, buffer)
    if (buffer == "" .and. .not. hasInlineData(child)) then
      ! no 'SpinTuning' block in REKS input
      allocate(ctrl%reksInp%Tuning(nType))
      do iType = 1, nType
        ctrl%reksInp%Tuning(iType) = 1.0_dp
      end do
    else
      ! 'SpinTuning' block in REKS input
      call init(realBuffer)
      call getChildValue(child, "", 1, realBuffer, modifier=modifier)
      nAtom = len(realBuffer)
      if (nAtom /= nType) then
        call detailedError(node, "Incorrect number of 'SpinTuning' block: " &
            & // i2c(nAtom) // " supplied, " &
            & // i2c(nType) // " required.")
      end if
      allocate(tmpTuning(1,nAtom))
      call asArray(realBuffer, tmpTuning)
      call destruct(realBuffer)
      allocate(ctrl%reksInp%Tuning(nType))
      ctrl%reksInp%Tuning(:) = tmpTuning(1,:)
    end if

  end subroutine readSpinTuning


end module dftbp_dftbplus_parser_reks
