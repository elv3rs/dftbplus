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
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning
  use hsd, only : hsd_get_or_set, hsd_get, hsd_get_table, hsd_get_choice, hsd_get_attrib, &
      & HSD_STAT_OK, hsd_schema_t, hsd_error_t, schema_init, schema_add_field, &
      & schema_validate, schema_destroy, FIELD_OPTIONAL, FIELD_REQUIRED, FIELD_TYPE_INTEGER, &
      & FIELD_TYPE_REAL, FIELD_TYPE_LOGICAL, FIELD_TYPE_TABLE, FIELD_TYPE_STRING, hsd_get_name
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_type_typegeometry, only : TGeometry
  use hsd_data, only : hsd_table, new_table
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

    call hsd_get_or_set(node, "ReadInitialCharges", ctrl%tReadChrg, .false.)
    if (.not. ctrl%tReadChrg) then
      call getInitialCharges(node, geo, ctrl%initialCharges)
    end if

    call hsd_get_or_set(node, "SCCTolerance", ctrl%sccTol, 1.0e-5_dp)

    ctrl%tWriteShifts = .false.

    if (geo%tPeriodic) then
      call hsd_get_or_set(node, "EwaldParameter", ctrl%ewaldAlpha, 0.0_dp)
      call hsd_get_or_set(node, "EwaldTolerance", ctrl%tolEwald, 1.0e-9_dp)
    end if

    if (geo%tHelical) then
      ! Tolerance for k-points being commensurate with C_n rotation
      call hsd_get_or_set(node, "HelicalSymmetryTol", ctrl%helicalSymTol, 1.0E-6_dp)
    end if

    ! self consistency required or not to proceed
    call hsd_get_or_set(node, "ConvergentSCCOnly", ctrl%isSccConvRequired, .true.)

    ! -- Schema validation (warnings only) --
    block
      type(hsd_schema_t) :: schema
      type(hsd_error_t), allocatable :: schemaErrors(:)
      integer :: iErr

      call schema_init(schema, name="SccOptions")
      call schema_add_field(schema, "ReadInitialCharges", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
      call schema_add_field(schema, "InitialCharges", FIELD_OPTIONAL, FIELD_TYPE_TABLE)
      call schema_add_field(schema, "SCCTolerance", FIELD_OPTIONAL, FIELD_TYPE_REAL)
      call schema_add_field(schema, "EwaldParameter", FIELD_OPTIONAL, FIELD_TYPE_REAL)
      call schema_add_field(schema, "EwaldTolerance", FIELD_OPTIONAL, FIELD_TYPE_REAL)
      call schema_add_field(schema, "HelicalSymmetryTol", FIELD_OPTIONAL, FIELD_TYPE_REAL)
      call schema_add_field(schema, "ConvergentSCCOnly", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
      call schema_validate(schema, node, schemaErrors)
      if (size(schemaErrors) > 0) then
        do iErr = 1, size(schemaErrors)
          call dftbp_warning(node, "[schema] " // schemaErrors(iErr)%message)
        end do
      end if
      call schema_destroy(schema)
    end block

  end subroutine readSccOptions


  !> Force evaluation options that are need for different hamiltonian choices
  subroutine readForceOptions(node, ctrl)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    type(hsd_table), pointer :: child
    character(len=:), allocatable :: buffer

    call hsd_get_or_set(node, "ForceEvaluation", buffer, "Traditional", child=child)
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

    ! -- Schema validation (warnings only) --
    block
      type(hsd_schema_t) :: schema
      type(hsd_error_t), allocatable :: schemaErrors(:)
      integer :: iErr

      call schema_init(schema, name="ForceOptions")
      call schema_add_field(schema, "ForceEvaluation", FIELD_OPTIONAL, FIELD_TYPE_STRING)
      call schema_validate(schema, node, schemaErrors)
      if (size(schemaErrors) > 0) then
        do iErr = 1, size(schemaErrors)
          call dftbp_warning(node, "[schema] " // schemaErrors(iErr)%message)
        end do
      end if
      call schema_destroy(schema)
    end block

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
    integer :: stat

    ! Artificially truncate the SK table
    call hsd_get(node, "SKMaxDistance", truncationCutOff, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'SKMaxDistance'")
    call hsd_get_attrib(node, "SKMaxDistance", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    call hsd_get_table(node, "SKMaxDistance", field, stat, auto_wrap=.true.)
    call convertUnitHsd(modifier, lengthUnits, field, truncationCutOff)

    call hsd_get_or_set(node, "HardCutOff", tHardCutOff, .true.)
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

    ! -- Schema validation (warnings only) --
    block
      type(hsd_schema_t) :: schema
      type(hsd_error_t), allocatable :: schemaErrors(:)
      integer :: iErr

      call schema_init(schema, name="SKTruncations")
      call schema_add_field(schema, "SKMaxDistance", FIELD_REQUIRED, FIELD_TYPE_REAL)
      call schema_add_field(schema, "HardCutOff", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
      call schema_validate(schema, node, schemaErrors)
      if (size(schemaErrors) > 0) then
        do iErr = 1, size(schemaErrors)
          call dftbp_warning(node, "[schema] " // schemaErrors(iErr)%message)
        end do
      end if
      call schema_destroy(schema)
    end block

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
    integer :: stat

    call hsd_get_table(node, "Differentiation", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: defTbl, defChild
        call new_table(defTbl, name="differentiation")
        call new_table(defChild, name="finitediff")
        call defTbl%add_child(defChild)
        call node%add_child(defTbl)
      end block
      call hsd_get_table(node, "Differentiation", child, stat, auto_wrap=.true.)
    end if
    call hsd_get_choice(child, "", buffer, val, stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(child, "Invalid or missing choice in 'Differentiation'")
    select case (buffer)
    case ("finitediff")
      ctrl%iDerivMethod = diffTypes%finiteDiff
      call hsd_get_or_set(val, "Delta", ctrl%deriv1stDelta, defDelta)
      call hsd_get_attrib(val, "Delta", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call hsd_get_table(val, "Delta", child, stat, auto_wrap=.true.)
      call convertUnitHsd(modifier, lengthUnits, child,&
          & ctrl%deriv1stDelta)
    case ("richardson")
      ctrl%iDerivMethod = diffTypes%richardson
    case default
      call hsd_get_name(val, buffer)
      call dftbp_error(child, "Invalid derivative calculation '" &
          & // buffer // "'")
    end select

    ! -- Schema validation for differentiation choice (warnings only) --
    if (associated(val)) then
      block
        type(hsd_schema_t) :: schema
        type(hsd_error_t), allocatable :: schemaErrors(:)
        integer :: iErr

        call schema_init(schema, name="Differentiation")
        call schema_add_field(schema, "Delta", FIELD_OPTIONAL, FIELD_TYPE_REAL)
        call schema_validate(schema, val, schemaErrors)
        if (size(schemaErrors) > 0) then
          do iErr = 1, size(schemaErrors)
            call dftbp_warning(val, "[schema] " // schemaErrors(iErr)%message)
          end do
        end if
        call schema_destroy(schema)
      end block
    end if

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
    integer :: iSp, stat

    ! X-H interaction corrections including H5 and damping
    ctrl%tDampH = .false.
    call hsd_get_table(node, "HCorrection", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: defTbl, defChild
        call new_table(defTbl, name="hcorrection")
        call new_table(defChild, name="none")
        call defTbl%add_child(defChild)
        call node%add_child(defTbl)
      end block
      call hsd_get_table(node, "HCorrection", child, stat, auto_wrap=.true.)
    end if
    call hsd_get_choice(child, "", buffer, value1, stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(child, "Invalid or missing choice in 'HCorrection'")

    select case (buffer)

    case ("none")
      ! nothing to do

    case ("damping")
      ! Switch the correction on
      ctrl%tDampH = .true.
      call hsd_get(value1, "Exponent", ctrl%dampExp, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'Exponent'")

    case ("h5")
      allocate(ctrl%h5Input)
      associate (h5Input => ctrl%h5Input)
        call hsd_get_or_set(value1, "RScaling", h5Input%rScale, 0.714_dp)
        call hsd_get_or_set(value1, "WScaling", h5Input%wScale, 0.25_dp)
        allocate(h5Input%elementParams(geo%nSpecies))
        call hsd_get_table(value1, "H5Scaling", child2, stat, auto_wrap=.true.)
        if (.not. associated(child2)) then
          block
            type(hsd_table) :: emptyTbl
            call new_table(emptyTbl, name="h5scaling")
            call value1%add_child(emptyTbl)
          end block
          call hsd_get_table(value1, "H5Scaling", child2, stat, auto_wrap=.true.)
        end if
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
          call hsd_get_or_set(child2, geo%speciesNames(iSp), h5Input%elementParams(iSp),&
              & h5ScalingDef)
        end do
        h5Input%speciesNames = geo%speciesNames
      end associate

    case default
      call hsd_get_name(value1, buffer)
      call dftbp_error(child, "Invalid HCorrection '" // buffer // "'")
    end select

    ! -- Schema validation for HCorrection choice (warnings only) --
    if (associated(value1)) then
      block
        type(hsd_schema_t) :: schema
        type(hsd_error_t), allocatable :: schemaErrors(:)
        integer :: iErr

        call schema_init(schema, name="HCorrection")
        call schema_add_field(schema, "Exponent", FIELD_OPTIONAL, FIELD_TYPE_REAL)
        call schema_add_field(schema, "RScaling", FIELD_OPTIONAL, FIELD_TYPE_REAL)
        call schema_add_field(schema, "WScaling", FIELD_OPTIONAL, FIELD_TYPE_REAL)
        call schema_add_field(schema, "H5Scaling", FIELD_OPTIONAL, FIELD_TYPE_TABLE)
        call schema_validate(schema, value1, schemaErrors)
        if (size(schemaErrors) > 0) then
          do iErr = 1, size(schemaErrors)
            call dftbp_warning(value1, "[schema] " // schemaErrors(iErr)%message)
          end do
        end if
        call schema_destroy(schema)
      end block
    end if

  end subroutine readHCorrection

end module dftbp_dftbplus_parser_sccoptions
