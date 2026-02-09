!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Reads the electronic solver settings from the input tree.
module dftbp_dftbplus_parser_solver
  use dftbp_common_accuracy, only : dp
  use dftbp_common_unitconversion, only : energyUnits
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_elecsolvers_elecsolvers, only : electronicSolverTypes
  use dftbp_extlibs_elsiiface, only : withELSI, withPEXSI
  use dftbp_extlibs_poisson, only : TPoissonInfo
  use dftbp_io_hsdutils, only : getChildValue, getNodeName, dftbp_error
  use hsd, only : hsd_get_or_set
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_io_message, only : error
  use hsd_data, only : hsd_table
  use dftbp_type_typegeometry, only : TGeometry
  use dftbp_dftbplus_parser_filling, only : readElectronicFilling
#:if WITH_TRANSPORT
  use dftbp_transport_negfvars, only : TTransPar, TNEGFGreenDensInfo
  use dftbp_dftbplus_parser_transport, only : readGreensFunction
#:endif

  implicit none

  private
  public :: readSolver

contains


  !> Electronic Solver
#:if WITH_TRANSPORT
  subroutine readSolver(node, ctrl, geo, tp, greendens, poisson)
#:else
  subroutine readSolver(node, ctrl, geo, poisson)
#:endif

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer, modifier

    integer :: iTmp

    ! Electronic solver
    call getChildValue(node, "Solver", value1, "RelativelyRobust")
    call getNodeName(value1, buffer)

    select case(buffer)

    case ("qr")
      ctrl%solver%isolver = electronicSolverTypes%qr

    case ("divideandconquer")
      ctrl%solver%isolver = electronicSolverTypes%divideandconquer

    case ("relativelyrobust")
      ctrl%solver%isolver = electronicSolverTypes%relativelyrobust

    case ("magma")
  #:if WITH_MAGMA
      ctrl%solver%isolver = electronicSolverTypes%magmaGvd
      call hsd_get_or_set(value1, "DensityMatrixGPU", ctrl%isDmOnGpu, .true.)
  #:else
      call dftbp_error(node, "DFTB+ must be compiled with MAGMA support in order to enable&
          & this solver")
  #:endif

    case ("elpa")
      allocate(ctrl%solver%elsi)
      call hsd_get_or_set(value1, "Sparse", ctrl%solver%elsi%elsiCsr, .false.)
      if (ctrl%solver%elsi%elsiCsr) then
        ctrl%solver%isolver = electronicSolverTypes%elpadm
      else
        ctrl%solver%isolver = electronicSolverTypes%elpa
      end if
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      call hsd_get_or_set(value1, "Mode", ctrl%solver%elsi%elpaSolver, 2)
      call hsd_get_or_set(value1, "Autotune", ctrl%solver%elsi%elpaAutotune, .false.)
      call getChildValue(value1, "Gpu", ctrl%solver%elsi%elpaGpu, .false., child=child)
      #:if not WITH_GPU
        if (ctrl%solver%elsi%elpaGpu) then
          call dftbp_error(child, "DFTB+ must be compiled with GPU support in order to enable&
              & the GPU acceleration for the ELPA solver")
        end if
      #:endif

    case ("omm")
      ctrl%solver%isolver = electronicSolverTypes%omm
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      call hsd_get_or_set(value1, "nIterationsELPA", ctrl%solver%elsi%ommIterationsElpa, 5)
      call hsd_get_or_set(value1, "Tolerance", ctrl%solver%elsi%ommTolerance, 1.0E-10_dp)
      call hsd_get_or_set(value1, "Choleskii", ctrl%solver%elsi%ommCholesky, .true.)

    case ("pexsi")
      ctrl%solver%isolver = electronicSolverTypes%pexsi
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
    #:if ELSI_VERSION > 2.5
      call hsd_get_or_set(value1, "Method", ctrl%solver%elsi%pexsiMethod, 3)
    #:else
      call hsd_get_or_set(value1, "Method", ctrl%solver%elsi%pexsiMethod, 2)
    #:endif
      select case(ctrl%solver%elsi%pexsiMethod)
      case(1)
        iTmp = 60
      case(2)
        iTmp = 20
      case(3)
        iTmp = 30
      end select
      call hsd_get_or_set(value1, "Poles", ctrl%solver%elsi%pexsiNPole, iTmp)
      if (ctrl%solver%elsi%pexsiNPole < 10) then
        call dftbp_error(value1, "Too few PEXSI poles")
      end if
      select case(ctrl%solver%elsi%pexsiMethod)
      case(1)
        if (mod(ctrl%solver%elsi%pexsiNPole,10) /= 0 .or. ctrl%solver%elsi%pexsiNPole > 120) then
          call dftbp_error(value1, "Unsupported number of PEXSI poles for method 1")
        end if
      case(2,3)
        if (mod(ctrl%solver%elsi%pexsiNPole,5) /= 0 .or. ctrl%solver%elsi%pexsiNPole > 40) then
          call dftbp_error(value1, "Unsupported number of PEXSI poles for this method")
        end if
      end select
      call hsd_get_or_set(value1, "ProcsPerPole", ctrl%solver%elsi%pexsiNpPerPole, 1)
      call hsd_get_or_set(value1, "muPoints", ctrl%solver%elsi%pexsiNMu, 2)
      call hsd_get_or_set(value1, "SymbolicFactorProcs", ctrl%solver%elsi%pexsiNpSymbo, 1)
      call getChildValue(value1, "SpectralRadius", ctrl%solver%elsi%pexsiDeltaE, 10.0_dp,&
          & modifier=modifier, child=child)
      call convertUnitHsd(modifier, energyUnits, child, ctrl%solver%elsi%pexsiDeltaE)

    case ("ntpoly")
      ctrl%solver%isolver = electronicSolverTypes%ntpoly
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      if (ctrl%tSpin) then
        call dftbp_error(value1, "Solver does not currently support spin polarisation")
      end if
      call hsd_get_or_set(value1, "PurificationMethod", ctrl%solver%elsi%ntpolyMethod, 2)
      call hsd_get_or_set(value1, "Tolerance", ctrl%solver%elsi%ntpolyTolerance, 1.0E-5_dp)
      call hsd_get_or_set(value1, "Truncation", ctrl%solver%elsi%ntpolyTruncation, 1.0E-10_dp)

  #:if WITH_TRANSPORT
    case ("greensfunction")
      ctrl%solver%isolver = electronicSolverTypes%GF
      ! need electronic temperature to be read for this solver:
      call readElectronicFilling(node, ctrl, geo)
      if (tp%defined .and. .not.tp%taskUpload) then
        call dftbp_error(node, "greensfunction solver cannot be used "// &
            &  "when task = contactHamiltonian")
      end if
      call readGreensFunction(value1, greendens, tp, ctrl%tempElec)
      ! fixEf also avoids checks of total charge later on in the run
      ctrl%tFixEf = .true.
    case ("transportonly")
      if (tp%defined .and. .not.tp%taskUpload) then
        call dftbp_error(node, "transportonly cannot be used when task = contactHamiltonian")
      end if
      call readGreensFunction(value1, greendens, tp, ctrl%tempElec)
      ctrl%solver%isolver = electronicSolverTypes%OnlyTransport
      ctrl%tFixEf = .true.
  #:endif

    case default
      call dftbp_error(value1, "Unknown electronic solver")

    end select

    if ((ctrl%solver%isolver == electronicSolverTypes%omm .or.&
        & ctrl%solver%isolver == electronicSolverTypes%pexsi ) .and. .not.ctrl%tSpinSharedEf&
        & .and. ctrl%tSpin .and. .not. ctrl%t2Component) then
      call dftbp_error(value1, "This solver currently requires spin values to be relaxed")
    end if
    if (ctrl%solver%isolver == electronicSolverTypes%pexsi .and. .not.withPEXSI) then
      call error("Not compiled with PEXSI support via ELSI")
    end if
    if (any(ctrl%solver%isolver == [electronicSolverTypes%elpa, electronicSolverTypes%omm,&
        & electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly])) then
      if (.not.withELSI) then
        call error("Not compiled with ELSI supported solvers")
      end if
    end if

    if (any(ctrl%solver%isolver == [electronicSolverTypes%omm, electronicSolverTypes%pexsi,&
        & electronicSolverTypes%ntpoly])) then
      call hsd_get_or_set(value1, "Sparse", ctrl%solver%elsi%elsiCsr, .true.)
      if (.not.ctrl%solver%elsi%elsiCsr) then
        if (any(ctrl%solver%isolver == [electronicSolverTypes%pexsi,electronicSolverTypes%ntpoly]))&
            & then
          call hsd_get_or_set(value1, "Threshold", ctrl%solver%elsi%elsi_zero_def, 1.0E-15_dp)
        end if
      end if
    end if

  #:if WITH_TRANSPORT
    if (all(ctrl%solver%isolver /= [electronicSolverTypes%GF,electronicSolverTypes%OnlyTransport])&
        & .and. tp%taskUpload) then
      call dftbp_error(value1, "Eigensolver incompatible with transport calculation&
          & (GreensFunction or TransportOnly required)")
    end if
  #:endif

  end subroutine readSolver


end module dftbp_dftbplus_parser_solver
