!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Reads SCC options, force evaluation, SK truncations, H correction and differentiation settings.
module dftbp_dftbplus_parser_sccoptions
  use dftbp_common_accuracy, only : distFudge, distFudgeOld, dp
  use dftbp_common_unitconversion, only : lengthUnits
  use dftbp_dftb_nonscc, only : diffTypes
  use dftbp_dftb_slakoeqgrid, only : skEqGridNew, skEqGridOld
  use dftbp_dftbplus_forcetypes, only : forceTypes
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_io_charmanip, only : tolower, unquote
  use dftbp_io_hsdutils, only : getChild, getChildValue, getNodeName, getNodeHSDName
  use dftbp_io_hsdutils, only : dftbp_error
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_type_typegeometry, only : TGeometry
  use hsd_data, only : hsd_table
  use dftbp_dftbplus_parser_spin, only : getInitialCharges
  implicit none

  private
  public :: readSccOptions, readForceOptions, SKTruncations, readHCorrection, readDifferentiation

contains


  subroutine readSccOptions(node, ctrl, geo)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    ctrl%tMulliken = .true.

    call getChildValue(node, "ReadInitialCharges", ctrl%tReadChrg, .false.)
    if (.not. ctrl%tReadChrg) then
      call getInitialCharges(node, geo, ctrl%initialCharges)
    end if

    call getChildValue(node, "SCCTolerance", ctrl%sccTol, 1.0e-5_dp)

    ! temporarily removed until debugged
    ! call getChildValue(node, "WriteShifts", ctrl%tWriteShifts, .false.)
    ctrl%tWriteShifts = .false.

    if (geo%tPeriodic) then
      call getChildValue(node, "EwaldParameter", ctrl%ewaldAlpha, 0.0_dp)
      call getChildValue(node, "EwaldTolerance", ctrl%tolEwald, 1.0e-9_dp)
    end if

    if (geo%tHelical) then
      ! Tolerance for k-points being commensurate with C_n rotation
      call getChildValue(node, "HelicalSymmetryTol", ctrl%helicalSymTol, 1.0E-6_dp)
    end if

    ! self consistency required or not to proceed
    call getChildValue(node, "ConvergentSCCOnly", ctrl%isSccConvRequired, .true.)

  end subroutine readSccOptions


  !> Force evaluation options that are need for different hamiltonian choices
  subroutine readForceOptions(node, ctrl)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    type(hsd_table), pointer :: child
    character(len=:), allocatable :: buffer

    call getChildValue(node, "ForceEvaluation", buffer, "Traditional", child=child)
    select case (tolower(unquote(buffer)))
    case("traditional")
      ctrl%forceType = forceTypes%orig
    case("dynamicst0")
      ctrl%forceType = forceTypes%dynamicT0
    case("dynamics")
      ctrl%forceType = forceTypes%dynamicTFinite
    case default
      call dftbp_error(child, "Invalid force evaluation method.")
    end select

  end subroutine readForceOptions


  !> Options for truncation of the SK data sets at a fixed distance
  subroutine SKTruncations(node, truncationCutOff, skInterMeth)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> This is the resulting cutoff distance
    real(dp), intent(out) :: truncationCutOff

    !> Method of the sk interpolation
    integer, intent(in) :: skInterMeth

    logical :: tHardCutOff
    type(hsd_table), pointer :: field
    character(len=:), allocatable :: modifier

    ! Artificially truncate the SK table
    call getChildValue(node, "SKMaxDistance", truncationCutOff, modifier=modifier, child=field)
    call convertUnitHsd(modifier, lengthUnits, field, truncationCutOff)

    call getChildValue(node, "HardCutOff", tHardCutOff, .true.)
    if (tHardCutOff) then
      ! Adjust by the length of the tail appended to the cutoff
      select case(skInterMeth)
      case(skEqGridOld)
        truncationCutOff = truncationCutOff - distFudgeOld
      case(skEqGridNew)
        truncationCutOff = truncationCutOff - distFudge
      end select
    end if
    if (truncationCutOff < epsilon(0.0_dp)) then
      call dftbp_error(field, "Truncation is shorter than the minimum distance over which SK data&
          & goes to 0")
    end if

  end subroutine SKTruncations


  !> Reads numerical differentiation method to be used
  subroutine readDifferentiation(node, ctrl)

    !> relevant node in input tree
    type(hsd_table), pointer, intent(in) :: node

    !> control structure to fill
    type(TControl), intent(inout) :: ctrl


    !> default of a reasonable choice for round off when using a second order finite difference
    !> formula
    real(dp), parameter :: defDelta = epsilon(1.0_dp)**0.25_dp

    character(len=:), allocatable :: buffer, modifier
    type(hsd_table), pointer :: val, child

    call getChildValue(node, "Differentiation", val, "FiniteDiff",&
        & child=child)
    call getNodeName(val, buffer)
    select case (buffer)
    case ("finitediff")
      ctrl%iDerivMethod = diffTypes%finiteDiff
      call getChildValue(val, "Delta", ctrl%deriv1stDelta, defDelta,&
          & modifier=modifier, child=child)
      call convertUnitHsd(modifier, lengthUnits, child,&
          & ctrl%deriv1stDelta)
    case ("richardson")
      ctrl%iDerivMethod = diffTypes%richardson
    case default
      call getNodeHSDName(val, buffer)
      call dftbp_error(child, "Invalid derivative calculation '" &
          & // buffer // "'")
    end select

  end subroutine readDifferentiation


  !> Reads the H corrections (H5, Damp)
  subroutine readHCorrection(node, geo, ctrl)

    !> Node containing the h-bond correction sub-block.
    type(hsd_table), pointer, intent(in) :: node

    !> Geometry.
    type(TGeometry), intent(in) :: geo

    !> Control structure
    type(TControl), intent(inout) :: ctrl

    type(hsd_table), pointer :: value1, child, child2
    character(len=:), allocatable :: buffer
    real(dp) :: h5ScalingDef
    integer :: iSp

    ! X-H interaction corrections including H5 and damping
    ctrl%tDampH = .false.
    call getChildValue(node, "HCorrection", value1, "None", child=child)
    call getNodeName(value1, buffer)

    select case (buffer)

    case ("none")
      ! nothing to do

    case ("damping")
      ! Switch the correction on
      ctrl%tDampH = .true.
      call getChildValue(value1, "Exponent", ctrl%dampExp)

    case ("h5")
      allocate(ctrl%h5Input)
      associate (h5Input => ctrl%h5Input)
        call getChildValue(value1, "RScaling", h5Input%rScale, 0.714_dp)
        call getChildValue(value1, "WScaling", h5Input%wScale, 0.25_dp)
        allocate(h5Input%elementParams(geo%nSpecies))
        call getChild(value1, "H5Scaling", child2, requested=.false., emptyIfMissing=.true.)
        do iSp = 1, geo%nSpecies
          select case (geo%speciesNames(iSp))
          case ("O")
            h5ScalingDef = 0.06_dp
          case ("N")
            h5ScalingDef = 0.18_dp
          case ("S")
            h5ScalingDef = 0.21_dp
          case default
            ! Default value is -1, this indicates that the element should be ignored
            h5ScalingDef = -1.0_dp
          end select
          call getChildValue(child2, geo%speciesNames(iSp), h5Input%elementParams(iSp),&
              & h5ScalingDef)
        end do
        h5Input%speciesNames = geo%speciesNames
      end associate

    case default
      call getNodeHSDName(value1, buffer)
      call dftbp_error(child, "Invalid HCorrection '" // buffer // "'")
    end select

  end subroutine readHCorrection

end module dftbp_dftbplus_parser_sccoptions
