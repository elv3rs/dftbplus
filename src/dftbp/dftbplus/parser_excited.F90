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
  use hsd, only : hsd_rename_child, hsd_get_or_set, hsd_get, hsd_get_table, hsd_set,&
      & hsd_get_attrib, hsd_get_choice, HSD_STAT_OK, hsd_schema_t, hsd_error_t, &
      & schema_init, schema_add_field, schema_validate, schema_destroy, FIELD_REQUIRED, &
      & FIELD_OPTIONAL, FIELD_TYPE_INTEGER, FIELD_TYPE_REAL, FIELD_TYPE_LOGICAL, &
      & FIELD_TYPE_TABLE, FIELD_TYPE_STRING
  use hsd_data, only : hsd_table
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning
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
    integer :: stat

    ! Linear response stuff
    call hsd_get_table(node, "Casida", child, stat, auto_wrap=.true.)

    if (associated(child)) then

      allocate(ctrl%lrespini)
      ctrl%lrespini%tPrintEigVecs = .false.

      if (ctrl%tSpin) then
        ctrl%lrespini%sym = ' '
      else
        call hsd_get(child, "Symmetry", buffer, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing required value: 'Symmetry'")
        call hsd_get_table(child, "Symmetry", child2, stat)
        if (.not. associated(child2)) child2 => child
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

      call hsd_get(child, "NrOfExcitations", ctrl%lrespini%nexc, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing required value: 'NrOfExcitations'")

      call hsd_get_table(child, "StateOfInterest", child2, stat, auto_wrap=.true.)
      if (.not. associated(child2)) then
        ctrl%lrespini%nstat = 0
        call hsd_set(child, "StateOfInterest", 0)
      else
        call hsd_get(child2, "#text", buffer, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value")
        if (tolower(unquote(buffer)) == "brightest") then
          if (ctrl%lrespini%sym /= "S" .or. ctrl%tSpin) then
            call dftbp_error(child2, "Brightest mode only allowed for spin unpolarised singlet&
                & excitations.")
          end if
          ctrl%lrespini%nstat = -1
        else
          call hsd_get(child2, "#text", ctrl%lrespini%nstat, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value")
          if (ctrl%lrespini%nstat > ctrl%lrespini%nexc) then
            call dftbp_error(child2, "Invalid value, must be within range of NrOfExcitations")
          elseif (ctrl%lrespini%sym == "B" .and. ctrl%lrespini%nstat /= 0) then
            call dftbp_error(child2, "You cannot specify a particular excited state if symmetry&
                & is 'B'")
          end if
        end if
      end if

      call hsd_get_or_set(child, "EnergyWindow", ctrl%lrespini%energyWindow, 0.0_dp, child=child2)
      call hsd_get_attrib(child, "EnergyWindow", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      ctrl%lrespini%tEnergyWindow = ctrl%lrespini%energyWindow /= 0.0_dp
      call convertUnitHsd(modifier, energyUnits, child2, ctrl%lrespini%energyWindow)
      call hsd_get_or_set(child, "OscillatorWindow", ctrl%lrespini%oscillatorWindow, 0.0_dp,&
          & child=child2)
      call hsd_get_attrib(child, "OscillatorWindow", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      ctrl%lrespini%tOscillatorWindow = ctrl%lrespini%oscillatorWindow /= 0.0_dp
      call convertUnitHsd(modifier, dipoleUnits, child2, ctrl%lrespini%oscillatorWindow)
      call hsd_get_or_set(child, "CacheCharges", ctrl%lrespini%tCacheCharges, .true.)
      call hsd_get_or_set(child, "WriteMulliken", ctrl%lrespini%tMulliken, .false.)
      call hsd_get_or_set(child, "WriteCoefficients", ctrl%lrespini%tCoeffs, .false.)
      ctrl%lrespini%tGrndState = .false.
      if (ctrl%lrespini%tCoeffs) then
        call hsd_get_or_set(child, "TotalStateCoeffs", ctrl%lrespini%tGrndState, .false.)
      end if
      call hsd_get_or_set(child, "WriteEigenvectors", ctrl%lrespini%tPrintEigVecs, .false.)
      call hsd_get_or_set(child, "WriteDensityMatrix", ctrl%lrespini%tWriteDensityMatrix, .false.)
      call hsd_get_or_set(child, "WriteXplusY", ctrl%lrespini%tXplusY, .false.)
      block
        integer, allocatable :: tmpI(:)
        call hsd_get_or_set(child, "StateCouplings", tmpI, [0, 0])
        ctrl%lrespini%indNACouplings(:min(size(tmpI),2)) = tmpI(:min(size(tmpI),2))
      end block
      if (all(ctrl%lrespini%indNACouplings == 0)) then
        ctrl%lrespini%tNaCoupling = .false.
      else
        ctrl%lrespini%tNaCoupling = .true.
      end if
      call hsd_get_or_set(child, "WriteSPTransitions", ctrl%lrespini%tSPTrans, .false.)
      call hsd_get_or_set(child, "WriteTransitions", ctrl%lrespini%tTrans, .false.)
      call hsd_get_or_set(child, "WriteTransitionDipole", ctrl%lrespini%tTradip, .false.)
      call hsd_get_or_set(child, "WriteTransitionCharges", ctrl%lrespini%tTransQ, .false.)
      ctrl%lrespini%iLinRespSolver = linRespSolverTypes%None

      call hsd_rename_child(child, "Diagonalizer", "Diagonaliser")
      call hsd_get_table(child, "Diagonaliser", child3, stat, auto_wrap=.true.)
      if (.not. associated(child3)) call dftbp_error(child, "Missing required block: 'Diagonaliser'")
      call hsd_get_choice(child3, "", buffer, child2, stat)
      if (associated(child2)) then
        select case(buffer)
        case ("arpack")
          if (.not. withArpack) then
            call dftbp_error(child2, 'This DFTB+ binary has been compiled without support for&
                & linear response calculations using the ARPACK/ngARPACK library.')
          end if
          call hsd_get_or_set(child2, "WriteStatusArnoldi", ctrl%lrespini%tArnoldi, .false.)
          call hsd_get_or_set(child2, "TestArnoldi", ctrl%lrespini%tDiagnoseArnoldi, .false.)
          ctrl%lrespini%iLinRespSolver = linRespSolverTypes%Arpack
        case ("stratmann")
          ctrl%lrespini%iLinRespSolver = linRespSolverTypes%Stratmann
          call hsd_get_or_set(child2, "SubSpaceFactor", ctrl%lrespini%subSpaceFactorStratmann, 20)
        case default
          call dftbp_error(child2, "Invalid diagonaliser method '" // buffer // "'")
        end select
      else
        call dftbp_error(child, "Missing diagonaliser method")
      end if

      call hsd_rename_child(child, "OptimizerCI", "OptimiserCI")
      call hsd_get_table(child, "OptimiserCI", child3, stat, auto_wrap=.true.)
      if (associated(child3)) then
        call hsd_get_choice(child3, "", buffer, child2, stat)
        select case(buffer)
        case ("bearpark")
          ctrl%lrespini%isCIopt = .true.
          call hsd_get_or_set(child2, "EnergyShift", ctrl%lrespini%energyShiftCI, 0.0_dp)
          call hsd_get_attrib(child2, "EnergyShift", modifier, stat)
          if (stat /= HSD_STAT_OK) modifier = ""
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
        call hsd_get_or_set(child, "ExcitedStateForces", ctrl%tCasidaForces, .true.)
      end if

      ! -- Schema validation for Casida (warnings only) --
      block
        type(hsd_schema_t) :: schema
        type(hsd_error_t), allocatable :: schemaErrors(:)
        integer :: iErr

        call schema_init(schema, name="Casida")
        call schema_add_field(schema, "Symmetry", FIELD_REQUIRED, FIELD_TYPE_STRING)
        call schema_add_field(schema, "NrOfExcitations", FIELD_REQUIRED, FIELD_TYPE_INTEGER)
        call schema_add_field(schema, "StateOfInterest", FIELD_OPTIONAL, FIELD_TYPE_INTEGER)
        call schema_add_field(schema, "EnergyWindow", FIELD_OPTIONAL, FIELD_TYPE_REAL)
        call schema_add_field(schema, "OscillatorWindow", FIELD_OPTIONAL, FIELD_TYPE_REAL)
        call schema_add_field(schema, "CacheCharges", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "WriteMulliken", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "WriteCoefficients", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "TotalStateCoeffs", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "WriteEigenvectors", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "WriteDensityMatrix", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "WriteXplusY", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "StateCouplings", FIELD_OPTIONAL, FIELD_TYPE_INTEGER)
        call schema_add_field(schema, "WriteSPTransitions", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "WriteTransitions", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "WriteTransitionDipole", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "WriteTransitionCharges", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "Diagonaliser", FIELD_REQUIRED, FIELD_TYPE_TABLE)
        call schema_add_field(schema, "OptimiserCI", FIELD_OPTIONAL, FIELD_TYPE_TABLE)
        call schema_add_field(schema, "ExcitedStateForces", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_validate(schema, child, schemaErrors)
        if (size(schemaErrors) > 0) then
          do iErr = 1, size(schemaErrors)
            call dftbp_warning(child, "[schema] " // schemaErrors(iErr)%message)
          end do
        end if
        call schema_destroy(schema)
      end block

    end if

    !pp-RPA
    call hsd_get_table(node, "PP-RPA", child, stat, auto_wrap=.true.)

    if (associated(child)) then

      allocate(ctrl%pprpa)

      if (ctrl%tSpin) then
        ctrl%pprpa%sym = ' '
      else
        call hsd_get(child, "Symmetry", buffer, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing required value: 'Symmetry'")
        call hsd_get_table(child, "Symmetry", child2, stat)
        if (.not. associated(child2)) child2 => child
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

      call hsd_get(child, "NrOfExcitations", ctrl%pprpa%nexc, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing required value: 'NrOfExcitations'")

      call hsd_get_table(child, "HHubbard", child2, stat, auto_wrap=.true.)
      if (.not. associated(child2)) call dftbp_error(child, "Missing required block: 'HHubbard'")
      allocate(ctrl%pprpa%hhubbard(geo%nSpecies))
      call readSpeciesList(child2, geo%speciesNames, ctrl%pprpa%hhubbard)

      call hsd_get_or_set(child, "TammDancoff", ctrl%pprpa%tTDA, .false.)

      call hsd_get_table(child, "NrOfVirtualStates", child2, stat, auto_wrap=.true.)
      if (.not. associated(child2)) then
        ctrl%pprpa%nvirtual = 0
        ctrl%pprpa%tConstVir = .false.
        call hsd_set(child, "NrOfVirtualStates", 0)
      else
        call hsd_get(child2, "#text", ctrl%pprpa%nvirtual, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value")
        ctrl%pprpa%tConstVir = .true.
      end if

      ! -- Schema validation for PP-RPA (warnings only) --
      block
        type(hsd_schema_t) :: schema
        type(hsd_error_t), allocatable :: schemaErrors(:)
        integer :: iErr

        call schema_init(schema, name="PP-RPA")
        call schema_add_field(schema, "Symmetry", FIELD_REQUIRED, FIELD_TYPE_STRING)
        call schema_add_field(schema, "NrOfExcitations", FIELD_REQUIRED, FIELD_TYPE_INTEGER)
        call schema_add_field(schema, "HHubbard", FIELD_REQUIRED, FIELD_TYPE_TABLE)
        call schema_add_field(schema, "TammDancoff", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "NrOfVirtualStates", FIELD_OPTIONAL, FIELD_TYPE_INTEGER)
        call schema_validate(schema, child, schemaErrors)
        if (size(schemaErrors) > 0) then
          do iErr = 1, size(schemaErrors)
            call dftbp_warning(child, "[schema] " // schemaErrors(iErr)%message)
          end do
        end if
        call schema_destroy(schema)
      end block

    end if

  end subroutine readExcited


end module dftbp_dftbplus_parser_excited
