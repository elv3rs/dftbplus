!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Compatibility shim mapping legacy HSD API patterns to hsd-fortran/hsd-data.
!>
!> This module provides call signatures similar to the legacy hsdutils/hsdutils2 modules
!> but operating on hsd_table (from hsd-fortran) instead of fnode (from xmlf90). It enables
!> incremental conversion of parser.F90 from fnode to hsd_table without a big-bang rewrite.
!>
!> Key differences from legacy API:
!>   - type(fnode), pointer  → type(hsd_table), intent(inout)/pointer
!>   - type(string)          → character(len=:), allocatable
!>   - getChildValue(node,name,val,default) → getChildValue(table,name,val,default)
!>   - getChild(node,name,child,requested)  → getChild(table,name,child,requested)
!>   - convertUnitHsd(mod,units,child,val)  → convertUnitHsd(mod,units,child,val)
!>
!> This is a TRANSITIONAL module — it will be removed once all call sites are converted
!> to direct hsd-fortran API calls (Phase 7).
module dftbp_io_hsdcompat
  use dftbp_common_accuracy, only : dp, lc, mc
  use dftbp_common_unitconversion, only : TUnit, convertUnit, statusCodes
  use dftbp_extlibs_hsddata, only : hsd_table, hsd_value, hsd_node, hsd_node_ptr, hsd_iterator, &
      & hsd_get, hsd_get_or, hsd_get_or_set, hsd_get_matrix, &
      & hsd_set, hsd_get_child, hsd_get_table, hsd_has_child, hsd_remove_child, &
      & hsd_child_count, hsd_get_keys, hsd_get_attrib, hsd_has_attrib, hsd_set_attrib, &
      & hsd_rename_child, hsd_get_choice, &
      & new_table, new_value, &
      & HSD_STAT_OK, HSD_STAT_NOT_FOUND, HSD_STAT_TYPE_ERROR
  use dftbp_io_message, only : error, warning
  implicit none

  private

  !> Public compatibility API — generic interfaces matching legacy calling patterns
  public :: getChildValue, getChild, setChildValue
  public :: detailedError, detailedWarning
  public :: convertUnitHsd
  public :: getNodeName2
  public :: warnUnprocessedNodes

  !> Error message prefix for invalid modifiers (matches legacy MSG_INVALID_MODIFIER)
  character(len=*), parameter :: MSG_INVALID_MODIFIER = "Unknown modifier '"

  !> Generic interface for reading child values.
  !>
  !> Maps legacy getChildValue(fnode, name, val, ...) patterns to hsd_get / hsd_get_or_set.
  !> Supports scalar and rank-1 reads for integer, real(dp), logical, character.
  !> For child block dispatch, use the table-returning variant.
  interface getChildValue
    ! Scalar without default (required)
    module procedure :: getChVal_int
    module procedure :: getChVal_real
    module procedure :: getChVal_logical
    module procedure :: getChVal_string
    ! Scalar with default
    module procedure :: getChVal_int_def
    module procedure :: getChVal_real_def
    module procedure :: getChVal_logical_def
    ! Rank-1 arrays without default (required)
    module procedure :: getChVal_intR1
    module procedure :: getChVal_realR1
    module procedure :: getChVal_logicalR1
    ! Rank-2 arrays without default
    module procedure :: getChVal_realR2
    module procedure :: getChVal_intR2
    ! Child table retrieval (for dispatch blocks)
    module procedure :: getChVal_table
  end interface getChildValue

  !> Generic interface for writing child values.
  !>
  !> Maps legacy setChildValue(fnode, name, val) patterns to hsd_set.
  interface setChildValue
    module procedure :: setChVal_int
    module procedure :: setChVal_real
    module procedure :: setChVal_logical
    module procedure :: setChVal_char
    module procedure :: setChVal_intR1
    module procedure :: setChVal_realR1
    module procedure :: setChVal_intR2
    module procedure :: setChVal_realR2
  end interface setChildValue

  !> Generic interface for unit conversion.
  !>
  !> Maps legacy convertUnitHsd(modifier, units, child, val) to
  !> convertUnit + detailedError. Same signature pattern as legacy but
  !> takes hsd_table instead of fnode for the child argument.
  interface convertUnitHsd
    module procedure :: convertUnitHsd_R0
    module procedure :: convertUnitHsd_R1
    module procedure :: convertUnitHsd_R2
  end interface convertUnitHsd

contains

  ! ============================================================
  !  getChildValue — scalar required (no default)
  ! ============================================================

  !> Get required integer child value
  subroutine getChVal_int(node, name, val, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(out) :: val
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    integer :: stat

    call hsd_get(node, name, val, stat=stat)
    if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Missing required integer value: '" // name // "'"))
    end if
    if (present(modifier)) call getModifier_(node, name, modifier)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine getChVal_int

  !> Get required real(dp) child value
  subroutine getChVal_real(node, name, val, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: val
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    integer :: stat

    call hsd_get(node, name, val, stat=stat)
    if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Missing required real value: '" // name // "'"))
    end if
    if (present(modifier)) call getModifier_(node, name, modifier)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine getChVal_real

  !> Get required logical child value
  subroutine getChVal_logical(node, name, val, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    logical, intent(out) :: val
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    integer :: stat

    call hsd_get(node, name, val, stat=stat)
    if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Missing required logical value: '" // name // "'"))
    end if
    if (present(modifier)) call getModifier_(node, name, modifier)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine getChVal_logical

  !> Get string child value, with optional default.
  !>
  !> If default is provided and the child is not found, the default value is used
  !> and written back to the tree (for dftb_pin.hsd generation).
  !> If default is not provided, the child is required and an error is raised if absent.
  subroutine getChVal_string(node, name, val, default, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(out) :: val
    character(len=*), intent(in), optional :: default
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    integer :: stat

    if (present(default)) then
      call hsd_get_or_set(node, name, val, default)
    else
      call hsd_get(node, name, val, stat=stat)
      if (stat /= HSD_STAT_OK) then
        call error(formatErrMsg_(node, "Missing required string value: '" // name // "'"))
      end if
    end if
    if (present(modifier)) call getModifier_(node, name, modifier)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine getChVal_string

  ! ============================================================
  !  getChildValue — scalar with default
  ! ============================================================

  !> Get integer child value with default (writes default back if absent)
  subroutine getChVal_int_def(node, name, val, default, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(out) :: val
    integer, intent(in) :: default
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    call hsd_get_or_set(node, name, val, default)
    if (present(modifier)) call getModifier_(node, name, modifier)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine getChVal_int_def

  !> Get real(dp) child value with default
  subroutine getChVal_real_def(node, name, val, default, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: val
    real(dp), intent(in) :: default
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    call hsd_get_or_set(node, name, val, default)
    if (present(modifier)) call getModifier_(node, name, modifier)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine getChVal_real_def

  !> Get logical child value with default
  subroutine getChVal_logical_def(node, name, val, default, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    logical, intent(out) :: val
    logical, intent(in) :: default
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    call hsd_get_or_set(node, name, val, default)
    if (present(modifier)) call getModifier_(node, name, modifier)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine getChVal_logical_def

  ! ============================================================
  !  getChildValue — rank-1 arrays (required, no default)
  ! ============================================================

  !> Get required integer array child value.
  !> The val array must be pre-allocated to the expected size.
  subroutine getChVal_intR1(node, name, val, nItem, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(out) :: val(:)
    integer, intent(out), optional :: nItem
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    integer, allocatable :: tmp(:)
    integer :: stat, nn

    call hsd_get(node, name, tmp, stat=stat)
    if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Missing required integer array: '" // name // "'"))
    end if
    nn = min(size(tmp), size(val))
    val(:nn) = tmp(:nn)
    if (present(nItem)) nItem = size(tmp)
    if (present(modifier)) call getModifier_(node, name, modifier)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine getChVal_intR1

  !> Get required real(dp) array child value.
  subroutine getChVal_realR1(node, name, val, nItem, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: val(:)
    integer, intent(out), optional :: nItem
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    real(dp), allocatable :: tmp(:)
    integer :: stat, nn

    call hsd_get(node, name, tmp, stat=stat)
    if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Missing required real array: '" // name // "'"))
    end if
    nn = min(size(tmp), size(val))
    val(:nn) = tmp(:nn)
    if (present(nItem)) nItem = size(tmp)
    if (present(modifier)) call getModifier_(node, name, modifier)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine getChVal_realR1

  !> Get required logical array child value.
  subroutine getChVal_logicalR1(node, name, val, nItem, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    logical, intent(out) :: val(:)
    integer, intent(out), optional :: nItem
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    logical, allocatable :: tmp(:)
    integer :: stat, nn

    call hsd_get(node, name, tmp, stat=stat)
    if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Missing required logical array: '" // name // "'"))
    end if
    nn = min(size(tmp), size(val))
    val(:nn) = tmp(:nn)
    if (present(nItem)) nItem = size(tmp)
    if (present(modifier)) call getModifier_(node, name, modifier)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine getChVal_logicalR1

  ! ============================================================
  !  getChildValue — rank-2 arrays (required, no default)
  ! ============================================================

  !> Get required real(dp) matrix child value.
  subroutine getChVal_realR2(node, name, val, nItem, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: val(:,:)
    integer, intent(out), optional :: nItem
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    real(dp), allocatable :: tmp(:,:)
    integer :: stat, nrows, ncols, r, c
    integer :: nr, nc

    call hsd_get_matrix(node, name, tmp, nrows, ncols, stat=stat)
    if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Missing required real matrix: '" // name // "'"))
    end if
    nr = min(nrows, size(val, 1))
    nc = min(ncols, size(val, 2))
    val(:nr, :nc) = tmp(:nr, :nc)
    if (present(nItem)) nItem = nrows * ncols
    if (present(modifier)) call getModifier_(node, name, modifier)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine getChVal_realR2

  !> Get required integer matrix child value.
  subroutine getChVal_intR2(node, name, val, nItem, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(out) :: val(:,:)
    integer, intent(out), optional :: nItem
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    integer, allocatable :: tmp(:,:)
    integer :: stat, nrows, ncols
    integer :: nr, nc

    call hsd_get_matrix(node, name, tmp, nrows, ncols, stat=stat)
    if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Missing required integer matrix: '" // name // "'"))
    end if
    nr = min(nrows, size(val, 1))
    nc = min(ncols, size(val, 2))
    val(:nr, :nc) = tmp(:nr, :nc)
    if (present(nItem)) nItem = nrows * ncols
    if (present(modifier)) call getModifier_(node, name, modifier)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine getChVal_intR2

  ! ============================================================
  !  getChildValue — child table retrieval (for dispatch blocks)
  ! ============================================================

  !> Get a child table for dispatch or block reading.
  !>
  !> This maps the legacy pattern:
  !>   call getChildValue(node, "", child)
  !>   call getNodeName(child, buffer)
  !>   select case (char(buffer))
  !>
  !> To:
  !>   call getChildValue(table, "", child, requested=.true.)
  !>   ! then use child%name for dispatch
  subroutine getChVal_table(node, name, child, requested, modifier, allowEmptyValue, dummyValue)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(hsd_table), pointer, intent(out) :: child
    logical, intent(in), optional :: requested
    character(len=:), allocatable, intent(out), optional :: modifier
    logical, intent(in), optional :: allowEmptyValue
    logical, intent(in), optional :: dummyValue

    integer :: stat
    logical :: isRequired

    isRequired = .true.
    if (present(requested)) isRequired = requested

    if (len_trim(name) == 0) then
      ! Empty name means "get the first table child" (dispatch pattern)
      call getFirstTableChild_(node, child, stat)
    else
      call hsd_get_table(node, name, child, stat)
    end if

    if (stat /= HSD_STAT_OK) then
      child => null()
      if (isRequired) then
        if (len_trim(name) == 0) then
          call error(formatErrMsg_(node, "Missing required child block"))
        else
          call error(formatErrMsg_(node, "Missing required block: '" // name // "'"))
        end if
      end if
      return
    end if

    if (present(modifier) .and. associated(child)) then
      if (len_trim(name) == 0) then
        ! For dispatch blocks, the modifier is on the child itself
        if (allocated(child%attrib)) then
          modifier = child%attrib
        else
          modifier = ""
        end if
      else
        call getModifier_(node, name, modifier)
      end if
    end if

  end subroutine getChVal_table

  ! ============================================================
  !  getChild — get a child table by name
  ! ============================================================

  !> Get a child table by name (maps legacy getChild pattern).
  !>
  !> @param node        Parent table.
  !> @param name        Name of the child to retrieve.
  !> @param child       Pointer to the child table on return (null if not found).
  !> @param requested   If .true. (default), error if not found.
  !> @param modifier    Optional: returns the child's attrib (modifier string).
  subroutine getChild(node, name, child, requested, modifier)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(hsd_table), pointer, intent(out) :: child
    logical, intent(in), optional :: requested
    character(len=:), allocatable, intent(out), optional :: modifier

    integer :: stat
    logical :: isRequired

    isRequired = .true.
    if (present(requested)) isRequired = requested

    call hsd_get_table(node, name, child, stat)

    if (stat /= HSD_STAT_OK) then
      child => null()
      if (isRequired) then
        call error(formatErrMsg_(node, "Missing required block: '" // name // "'"))
      end if
    end if

    if (present(modifier) .and. associated(child)) then
      call getModifier_(node, name, modifier)
    end if

  end subroutine getChild

  ! ============================================================
  !  setChildValue — write values to tree
  ! ============================================================

  !> Set integer child value
  subroutine setChVal_int(node, name, val, replace, child)
    type(hsd_table), intent(inout) :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: val
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    ! In hsd-fortran, hsd_set always replaces (upsert semantics)
    call hsd_set(node, name, val)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine setChVal_int

  !> Set real(dp) child value
  subroutine setChVal_real(node, name, val, replace, child)
    type(hsd_table), intent(inout) :: node
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: val
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    call hsd_set(node, name, val)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine setChVal_real

  !> Set logical child value
  subroutine setChVal_logical(node, name, val, replace, child)
    type(hsd_table), intent(inout) :: node
    character(len=*), intent(in) :: name
    logical, intent(in) :: val
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    call hsd_set(node, name, val)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine setChVal_logical

  !> Set character child value
  subroutine setChVal_char(node, name, val, replace, child)
    type(hsd_table), intent(inout) :: node
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    call hsd_set(node, name, val)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine setChVal_char

  !> Set integer array child value
  subroutine setChVal_intR1(node, name, val, replace, child)
    type(hsd_table), intent(inout) :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: val(:)
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    call hsd_set(node, name, val)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine setChVal_intR1

  !> Set real(dp) array child value
  subroutine setChVal_realR1(node, name, val, replace, child)
    type(hsd_table), intent(inout) :: node
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: val(:)
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    call hsd_set(node, name, val)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine setChVal_realR1

  !> Set integer matrix child value
  subroutine setChVal_intR2(node, name, val, replace, child)
    type(hsd_table), intent(inout) :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: val(:,:)
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    call hsd_set(node, name, val)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine setChVal_intR2

  !> Set real(dp) matrix child value
  subroutine setChVal_realR2(node, name, val, replace, child)
    type(hsd_table), intent(inout) :: node
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: val(:,:)
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    call hsd_set(node, name, val)
    if (present(child)) call hsd_get_table(node, name, child)

  end subroutine setChVal_realR2

  ! ============================================================
  !  Error reporting
  ! ============================================================

  !> Report a fatal error with node context (replaces legacy detailedError)
  subroutine detailedError(node, msg)
    type(hsd_table), intent(in) :: node
    character(len=*), intent(in) :: msg

    call error(formatErrMsg_(node, msg))

  end subroutine detailedError

  !> Report a warning with node context (replaces legacy detailedWarning)
  subroutine detailedWarning(node, msg)
    type(hsd_table), intent(in) :: node
    character(len=*), intent(in) :: msg

    call warning(formatWarnMsg_(node, msg))

  end subroutine detailedWarning

  ! ============================================================
  !  convertUnitHsd — unit conversion with modifier
  ! ============================================================

  !> Convert a scalar value using the unit modifier (rank 0).
  !>
  !> Equivalent to legacy convertUnitHsd(modifier, units, child, val).
  !> If modifier is empty, no conversion is done.
  subroutine convertUnitHsd_R0(modifier, units, child, convertValue, replace, changed)
    character(len=*), intent(in) :: modifier
    type(TUnit), intent(in) :: units(:)
    type(hsd_table), intent(inout) :: child
    real(dp), intent(inout) :: convertValue
    logical, intent(in), optional :: replace
    logical, intent(out), optional :: changed

    logical :: changed_
    integer :: status

    changed_ = len_trim(modifier) > 0
    if (changed_) then
      call convertUnit(units, modifier, convertValue, status)
      if (status /= statusCodes%ok) then
        call detailedError(child, MSG_INVALID_MODIFIER // modifier // "'")
      end if
      if (present(replace)) then
        if (replace) call hsd_set(child, "", convertValue)
      end if
    end if
    if (present(changed)) changed = changed_

  end subroutine convertUnitHsd_R0

  !> Convert a rank-1 array using the unit modifier.
  subroutine convertUnitHsd_R1(modifier, units, child, convertValue, replace, changed)
    character(len=*), intent(in) :: modifier
    type(TUnit), intent(in) :: units(:)
    type(hsd_table), intent(inout) :: child
    real(dp), intent(inout) :: convertValue(:)
    logical, intent(in), optional :: replace
    logical, intent(out), optional :: changed

    logical :: changed_
    integer :: status

    changed_ = len_trim(modifier) > 0
    if (changed_) then
      call convertUnit(units, modifier, convertValue, status)
      if (status /= statusCodes%ok) then
        call detailedError(child, MSG_INVALID_MODIFIER // modifier // "'")
      end if
      if (present(replace)) then
        if (replace) call hsd_set(child, "", convertValue)
      end if
    end if
    if (present(changed)) changed = changed_

  end subroutine convertUnitHsd_R1

  !> Convert a rank-2 matrix using the unit modifier.
  subroutine convertUnitHsd_R2(modifier, units, child, convertValue, replace, changed)
    character(len=*), intent(in) :: modifier
    type(TUnit), intent(in) :: units(:)
    type(hsd_table), intent(inout) :: child
    real(dp), intent(inout) :: convertValue(:,:)
    logical, intent(in), optional :: replace
    logical, intent(out), optional :: changed

    logical :: changed_
    integer :: status

    changed_ = len_trim(modifier) > 0
    if (changed_) then
      call convertUnit(units, modifier, convertValue, status)
      if (status /= statusCodes%ok) then
        call detailedError(child, MSG_INVALID_MODIFIER // modifier // "'")
      end if
      if (present(replace)) then
        if (replace) call hsd_set(child, "", convertValue)
      end if
    end if
    if (present(changed)) changed = changed_

  end subroutine convertUnitHsd_R2

  ! ============================================================
  !  getNodeName2 — get node name safely
  ! ============================================================

  !> Get the name of an HSD node (replaces legacy getNodeName / getNodeName2).
  !>
  !> Returns the node's name as an allocatable character string.
  subroutine getNodeName2(node, nodeName)
    type(hsd_table), intent(in) :: node
    character(len=:), allocatable, intent(out) :: nodeName

    if (allocated(node%name)) then
      nodeName = node%name
    else
      nodeName = ""
    end if

  end subroutine getNodeName2

  ! ============================================================
  !  warnUnprocessedNodes — stub for transition
  ! ============================================================

  !> Warn about unprocessed nodes (stub during transition).
  !>
  !> In the new architecture, unprocessed node detection is handled by
  !> schema_validate_strict. This is a no-op stub during transition.
  subroutine warnUnprocessedNodes(node, tIgnoreUnprocessed)
    type(hsd_table), intent(in) :: node
    logical, intent(in), optional :: tIgnoreUnprocessed

    ! No-op: schema_validate_strict will replace this
    ! once schemas are defined for each parser section.

  end subroutine warnUnprocessedNodes

  ! ============================================================
  !  Internal helpers
  ! ============================================================

  !> Get the modifier (attrib) of a child node by path.
  subroutine getModifier_(parent, name, modifier)
    type(hsd_table), intent(in), target :: parent
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(out) :: modifier

    integer :: stat

    call hsd_get_attrib(parent, name, modifier, stat)
    if (stat /= HSD_STAT_OK) then
      modifier = ""
    end if

  end subroutine getModifier_

  !> Get the first table child of a node (for dispatch patterns with empty name).
  subroutine getFirstTableChild_(node, child, stat)
    type(hsd_table), intent(in), target :: node
    type(hsd_table), pointer, intent(out) :: child
    integer, intent(out) :: stat

    type(hsd_iterator) :: iter
    class(hsd_node), pointer :: cur

    child => null()
    stat = HSD_STAT_NOT_FOUND

    call iter%init(node)
    do while (iter%next(cur))
      select type (t => cur)
      type is (hsd_table)
        child => t
        stat = HSD_STAT_OK
        return
      end select
    end do

  end subroutine getFirstTableChild_

  !> Format an error message with node context.
  function formatErrMsg_(node, msg) result(full_msg)
    type(hsd_table), intent(in) :: node
    character(len=*), intent(in) :: msg
    character(len=:), allocatable :: full_msg

    character(len=20) :: linebuf

    if (allocated(node%name) .and. node%line > 0) then
      write(linebuf, '(i0)') node%line
      full_msg = "Error in '" // node%name // "' (line " // trim(linebuf) // "): " // msg
    else if (allocated(node%name)) then
      full_msg = "Error in '" // node%name // "': " // msg
    else
      full_msg = msg
    end if

  end function formatErrMsg_

  !> Format a warning message with node context.
  function formatWarnMsg_(node, msg) result(full_msg)
    type(hsd_table), intent(in) :: node
    character(len=*), intent(in) :: msg
    character(len=:), allocatable :: full_msg

    character(len=20) :: linebuf

    if (allocated(node%name) .and. node%line > 0) then
      write(linebuf, '(i0)') node%line
      full_msg = "Warning in '" // node%name // "' (line " // trim(linebuf) // "): " // msg
    else if (allocated(node%name)) then
      full_msg = "Warning in '" // node%name // "': " // msg
    else
      full_msg = msg
    end if

  end function formatWarnMsg_

end module dftbp_io_hsdcompat
