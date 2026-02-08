!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Reads the excited state data block from the HSD input.
module dftbp_dftbplus_parser_excited
  use dftbp_common_accuracy, only : dp
  use dftbp_common_unitconversion, only : dipoleUnits, energyUnits
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_dftbplus_specieslist, only : readSpeciesList
  use dftbp_extlibs_arpack, only : withArpack
  use dftbp_io_charmanip, only : tolower, unquote
  use hsd, only : hsd_rename_child
  use hsd_data, only : hsd_table
  use dftbp_io_hsdutils, only : getChild, getChildValue, setChildValue
  use dftbp_io_hsdutils, only : dftbp_error, getNodeName
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_timedep_linresptypes, only : linRespSolverTypes
  use dftbp_type_typegeometry, only : TGeometry
  implicit none

  private
  public :: readExcited

contains


  !> Reads the excited state data block
  subroutine readExcited(node, geo, ctrl)

    !> Node to parse
    type(hsd_table), pointer :: node

    !> geometry object, which contains atomic species information
    type(TGeometry), intent(in) :: geo

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    type(hsd_table), pointer :: child
    type(hsd_table), pointer :: child2, child3
    type(hsd_table), pointer :: value
    character(len=:), allocatable :: buffer, modifier

    ! Linear response stuff
    call getChild(node, "Casida", child, requested=.false.)

    if (associated(child)) then

      allocate(ctrl%lrespini)
      ctrl%lrespini%tPrintEigVecs = .false.

      if (ctrl%tSpin) then
        ctrl%lrespini%sym = ' '
      else
        call getChildValue(child, "Symmetry", buffer, child=child2)
        select case (unquote(buffer))
        case ("Singlet" , "singlet")
          ctrl%lrespini%sym = 'S'
        case ("Triplet" , "triplet")
          ctrl%lrespini%sym = 'T'
        case ("Both" , "both")
          ctrl%lrespini%sym = 'B'
        case default
          call dftbp_error(child2, "Invalid symmetry value '"  // buffer // &
              & "' (must be 'Singlet', 'Triplet' or 'Both').")
        end select
      end if

      call getChildValue(child, "NrOfExcitations", ctrl%lrespini%nexc)

      call getChild(child, "StateOfInterest", child2, requested=.false.)
      if (.not. associated(child2)) then
        ctrl%lrespini%nstat = 0
        call setChildValue(child, "StateOfInterest", 0)
      else
        call getChildValue(child2, "", buffer)
        if (tolower(unquote(buffer)) == "brightest") then
          if (ctrl%lrespini%sym /= "S" .or. ctrl%tSpin) then
            call dftbp_error(child2, "Brightest mode only allowed for spin unpolarised singlet&
                & excitations.")
          end if
          ctrl%lrespini%nstat = -1
        else
          call getChildValue(child2, "", ctrl%lrespini%nstat)
          if (ctrl%lrespini%nstat > ctrl%lrespini%nexc) then
            call dftbp_error(child2, "Invalid value, must be within range of NrOfExcitations")
          elseif (ctrl%lrespini%sym == "B" .and. ctrl%lrespini%nstat /= 0) then
            call dftbp_error(child2, "You cannot specify a particular excited state if symmetry&
                & is 'B'")
          end if
        end if
      end if

      call getChildValue(child, "EnergyWindow", ctrl%lrespini%energyWindow, 0.0_dp, &
          & modifier=modifier, child=child2)
      ctrl%lrespini%tEnergyWindow = ctrl%lrespini%energyWindow /= 0.0_dp
      call convertUnitHsd(modifier, energyUnits, child2, ctrl%lrespini%energyWindow)
      call getChildValue(child, "OscillatorWindow", ctrl%lrespini%oscillatorWindow, 0.0_dp, &
          & modifier=modifier,  child=child2)
      ctrl%lrespini%tOscillatorWindow = ctrl%lrespini%oscillatorWindow /= 0.0_dp
      call convertUnitHsd(modifier, dipoleUnits, child2, ctrl%lrespini%oscillatorWindow)
      call getChildValue(child, "CacheCharges", ctrl%lrespini%tCacheCharges, default=.true.)
      call getChildValue(child, "WriteMulliken", ctrl%lrespini%tMulliken, default=.false.)
      call getChildValue(child, "WriteCoefficients", ctrl%lrespini%tCoeffs, default=.false.)
      ctrl%lrespini%tGrndState = .false.
      if (ctrl%lrespini%tCoeffs) then
        call getChildValue(child, "TotalStateCoeffs", ctrl%lrespini%tGrndState, .false.)
      end if
      call getChildValue(child, "WriteEigenvectors", ctrl%lrespini%tPrintEigVecs, .false.)
      call getChildValue(child, "WriteDensityMatrix", ctrl%lrespini%tWriteDensityMatrix, .false.)
      call getChildValue(child, "WriteXplusY", ctrl%lrespini%tXplusY, default=.false.)
      call getChildValue(child, "StateCouplings", ctrl%lrespini%indNACouplings, default=[0, 0])
      if (all(ctrl%lrespini%indNACouplings == 0)) then
        ctrl%lrespini%tNaCoupling = .false.
      else
        ctrl%lrespini%tNaCoupling = .true.
      end if
      call getChildValue(child, "WriteSPTransitions", ctrl%lrespini%tSPTrans, default=.false.)
      call getChildValue(child, "WriteTransitions", ctrl%lrespini%tTrans, default=.false.)
      call getChildValue(child, "WriteTransitionDipole", ctrl%lrespini%tTradip, default=.false.)
      call getChildValue(child, "WriteTransitionCharges", ctrl%lrespini%tTransQ, default=.false.)
      ctrl%lrespini%iLinRespSolver = linRespSolverTypes%None

      call hsd_rename_child(child, "Diagonalizer", "Diagonaliser")
      call getChildValue(child, "Diagonaliser", child2, allowEmptyValue=.true.)
      if (associated(child2)) then
        call getNodeName(child2, buffer)
        select case(buffer)
        case ("arpack")
          if (.not. withArpack) then
            call dftbp_error(child2, 'This DFTB+ binary has been compiled without support for&
                & linear response calculations using the ARPACK/ngARPACK library.')
          end if
          call getChildValue(child2, "WriteStatusArnoldi", ctrl%lrespini%tArnoldi, default=.false.)
          call getChildValue(child2, "TestArnoldi", ctrl%lrespini%tDiagnoseArnoldi, default=.false.)
          ctrl%lrespini%iLinRespSolver = linRespSolverTypes%Arpack
        case ("stratmann")
          ctrl%lrespini%iLinRespSolver = linRespSolverTypes%Stratmann
          call getChildValue(child2, "SubSpaceFactor", ctrl%lrespini%subSpaceFactorStratmann, 20)
        case default
          call dftbp_error(child2, "Invalid diagonaliser method '" // buffer // "'")
        end select
      else
        call dftbp_error(child, "Missing diagonaliser method")
      end if

      call hsd_rename_child(child, "OptimizerCI", "OptimiserCI")
      call getChild(child, "OptimiserCI", child2, requested=.false.)
      if (associated(child2)) then
        call getChildValue(child, "OptimiserCI", child2, child=child3)
        call getNodeName(child2, buffer)
        select case(buffer)
        case ("bearpark")
          ctrl%lrespini%isCIopt = .true.
          call getChildValue(child2, "EnergyShift", ctrl%lrespini%energyShiftCI,&
              & modifier=modifier, default=0.0_dp)
          call convertUnitHsd(modifier, energyUnits, child, ctrl%lrespini%energyShiftCI)
        case ("")
          call dftbp_error(child2, "Missing choice of CI optimiser.")
        case default
          call dftbp_error(child2, "Invalid CI optimiser method '" // buffer // "'")
        end select
      else
        ctrl%lrespini%isCIopt = .false.
      end if

      if (ctrl%tForces .or. ctrl%tPrintForces) then
        call getChildValue(child, "ExcitedStateForces", ctrl%tCasidaForces, default=.true.)
      end if

    end if

    !pp-RPA
    call getChild(node, "PP-RPA", child, requested=.false.)

    if (associated(child)) then

      allocate(ctrl%pprpa)

      if (ctrl%tSpin) then
        ctrl%pprpa%sym = ' '
      else
        call getChildValue(child, "Symmetry", buffer, child=child2)
        select case (unquote(buffer))
        case ("Singlet" , "singlet")
          ctrl%pprpa%sym = 'S'
        case ("Triplet" , "triplet")
          ctrl%pprpa%sym = 'T'
        case ("Both" , "both")
          ctrl%pprpa%sym = 'B'
        case default
          call dftbp_error(child2, "Invalid symmetry value '"  // buffer // &
              & "' (must be 'Singlet', 'Triplet' or 'Both').")
        end select
      end if

      call getChildValue(child, "NrOfExcitations", ctrl%pprpa%nexc)

      call getChildValue(child, "HHubbard", value, child=child2)
      allocate(ctrl%pprpa%hhubbard(geo%nSpecies))
      call readSpeciesList(child2, geo%speciesNames, ctrl%pprpa%hhubbard)

      call getChildValue(child, "TammDancoff", ctrl%pprpa%tTDA, default=.false.)

      call getChild(child, "NrOfVirtualStates", child2, requested=.false.)
      if (.not. associated(child2)) then
        ctrl%pprpa%nvirtual = 0
        ctrl%pprpa%tConstVir = .false.
        call setChildValue(child, "NrOfVirtualStates", 0)
      else
        call getChildValue(child2, "", ctrl%pprpa%nvirtual)
        ctrl%pprpa%tConstVir = .true.
      end if

    end if

  end subroutine readExcited


end module dftbp_dftbplus_parser_excited
