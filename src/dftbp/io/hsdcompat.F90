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
  use dftbp_common_status, only : TStatus
  use dftbp_common_unitconversion, only : TUnit, convertUnit, statusCodes
  use dftbp_extlibs_hsddata, only : hsd_table, hsd_value, hsd_node, hsd_node_ptr, hsd_iterator, &
      & hsd_get, hsd_get_or, hsd_get_or_set, hsd_get_matrix, &
      & hsd_set, hsd_get_child, hsd_get_table, hsd_has_child, hsd_remove_child, &
      & hsd_child_count, hsd_get_keys, hsd_get_attrib, hsd_has_attrib, hsd_set_attrib, &
      & hsd_rename_child, hsd_get_choice, &
      & new_table, new_value, &
      & HSD_STAT_OK, HSD_STAT_NOT_FOUND, HSD_STAT_TYPE_ERROR
  use dftbp_io_charmanip, only : i2c, tolower
  use dftbp_io_indexselection, only : getIndexSelection
  use dftbp_io_message, only : error, warning
  implicit none

  private

  !> Public compatibility API — generic interfaces matching legacy calling patterns
  public :: getChildValue, getChild, setChildValue, setChild
  public :: detailedError, detailedWarning
  public :: convertUnitHsd
  public :: getNodeName, getNodeName2, getNodeHSDName
  public :: getChildren, getLength, getItem1, destroyNodeList
  public :: getSelectedAtomIndices, getSelectedIndices
  public :: splitModifier
  public :: warnUnprocessedNodes, setUnprocessed

  !> Constant for text node name (matches xmlf90 textNodeName)
  character(len=*), parameter, public :: textNodeName = "#text"

  !> Error message prefix for invalid modifiers (matches legacy MSG_INVALID_MODIFIER)
  character(len=*), parameter :: MSG_INVALID_MODIFIER = "Unknown modifier '"

  !> Separator character for splitting modifiers
  character, parameter :: sepModifier = ","

  !> Wrapper for a pointer to hsd_table (for arrays of pointers).
  type :: hsd_table_ptr
    type(hsd_table), pointer :: ptr => null()
  end type hsd_table_ptr

  !> List of child tables (replaces fnodeList).
  !>
  !> Created by getChildren(), indexed by getItem1(), sized by getLength(),
  !> freed by destroyNodeList().
  type, public :: hsd_child_list
    integer :: count = 0
    type(hsd_table_ptr), allocatable :: items(:)
  end type hsd_child_list

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
  subroutine getChVal_string(node, name, val, default, modifier, child, multiple)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(out) :: val
    character(len=*), intent(in), optional :: default
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child
    logical, intent(in), optional :: multiple

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
  !>   call getChildValue(node, "Name", value1, "", child=child)
  !>   call getNodeName(value1, buffer)
  !>   select case (buffer)
  !>
  !> The variableValue (child) output is a pointer to the first table child of the named block.
  !> If the first child is a value node (inline data), variableValue is set to null,
  !> and getNodeName will return "#text" for it.
  !>
  !> When called with 3 positional args and child= keyword:
  !>   call getChildValue(node, "Name", value1, child=child)
  !> variableValue = first table child, child = named container
  !>
  !> When called with 4 positional args (with default ""):
  !>   call getChildValue(node, "Name", value1, "", child=child, ...)
  !> Same, but creates block if missing.
  subroutine getChVal_table(node, name, variableValue, default, modifier, child, &
      & list, allowEmptyValue, dummyValue, multiple, requested)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(hsd_table), pointer, intent(out) :: variableValue
    character(len=*), intent(in), optional :: default
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child
    logical, intent(in), optional :: list
    logical, intent(in), optional :: allowEmptyValue
    logical, intent(in), optional :: dummyValue
    logical, intent(in), optional :: multiple
    logical, intent(in), optional :: requested

    type(hsd_table), pointer :: container
    integer :: stat
    logical :: isRequired

    variableValue => null()

    isRequired = .true.
    if (present(requested)) isRequired = requested
    if (present(default)) isRequired = .false.

    if (len_trim(name) == 0) then
      ! Empty name means "get the first child" (dispatch pattern from root)
      container => node
    else
      ! Find the named child table
      call hsd_get_table(node, name, container, stat)
      if (stat /= HSD_STAT_OK) then
        if (present(default)) then
          ! Create the child with default
          allocate(container)
          call new_table(container, name=tolower(name))
          call node%add_child(container)
        else if (isRequired) then
          call error(formatErrMsg_(node, "Missing required block: '" // name // "'"))
          return
        else
          container => null()
          if (present(child)) child => null()
          if (present(modifier)) modifier = ""
          return
        end if
      end if
    end if

    ! Get the first table child of the container for dispatch
    if (associated(container)) then
      call getFirstTableChild_(container, variableValue, stat)
    end if
    ! If no table child found, variableValue stays null (= "#text" / inline data)

    if (present(child)) then
      child => container
    end if

    if (present(modifier)) then
      if (associated(container) .and. allocated(container%attrib)) then
        modifier = container%attrib
      else
        modifier = ""
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
  !  getNodeName — get raw node name (legacy getNodeName from xmlf90)
  ! ============================================================

  !> Get the raw name of a node.
  !>
  !> For hsd_table nodes, returns the node name. For the dispatch pattern where
  !> the table pointer is null (indicating inline value data), returns "#text".
  !> This matches the legacy xmlf90 getNodeName behavior.
  subroutine getNodeName(node, nodeName)
    type(hsd_table), pointer, intent(in) :: node
    character(len=:), allocatable, intent(out) :: nodeName

    if (.not. associated(node)) then
      nodeName = textNodeName
    else if (allocated(node%name)) then
      nodeName = node%name
    else
      nodeName = ""
    end if

  end subroutine getNodeName

  ! ============================================================
  !  getNodeHSDName — get pretty HSD node name (for error messages)
  ! ============================================================

  !> Get the HSD-style name of a node (for error messages).
  !>
  !> Returns the node name in its original case. In the legacy code, this would
  !> look up the attrName attribute. In hsd-fortran, we just return the name.
  subroutine getNodeHSDName(node, nodeName)
    type(hsd_table), intent(in), target :: node
    character(len=:), allocatable, intent(out) :: nodeName

    if (allocated(node%name)) then
      nodeName = node%name
    else
      nodeName = ""
    end if

  end subroutine getNodeHSDName

  ! ============================================================
  !  setChild — create a new child block
  ! ============================================================

  !> Create a new named child block under a parent.
  !>
  !> This replaces the legacy setChild(node, name, child, replace, list, modifier)
  !> pattern that creates a new element node and appends it.
  subroutine setChild(node, name, child, replace, list, modifier)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(hsd_table), pointer, intent(out) :: child
    logical, intent(in), optional :: replace
    logical, intent(in), optional :: list
    character(len=*), intent(in), optional :: modifier

    type(hsd_table), pointer :: newChild

    ! Create and add a new table child
    allocate(newChild)
    call new_table(newChild, name=name)
    call node%add_child(newChild)
    child => newChild

    if (present(modifier)) then
      if (len_trim(modifier) > 0) then
        call hsd_set_attrib(node, name, modifier)
      end if
    end if

  end subroutine setChild

  ! ============================================================
  !  setUnprocessed — no-op stub
  ! ============================================================

  !> Mark a node as unprocessed (no-op in hsd-data).
  !>
  !> In the legacy code, this removed the "processed" attribute.
  !> In hsd-data, processing tracking is not used.
  subroutine setUnprocessed(node)
    type(hsd_table), intent(in) :: node

    ! No-op

  end subroutine setUnprocessed

  ! ============================================================
  !  getChildren / getLength / getItem1 / destroyNodeList
  !  — fnodeList compatibility layer
  ! ============================================================

  !> Get all child tables matching a given name.
  !>
  !> This replaces the legacy getChildren(node, name, children) pattern.
  !> Returns an hsd_child_list containing pointers to matching children.
  subroutine getChildren(node, name, children)
    type(hsd_table), intent(in), target :: node
    character(len=*), intent(in) :: name
    type(hsd_child_list), pointer, intent(out) :: children

    type(hsd_iterator) :: iter
    class(hsd_node), pointer :: cur
    integer :: cnt, idx
    character(len=:), allocatable :: childName

    ! First pass: count matches
    cnt = 0
    call iter%init(node)
    do while (iter%next(cur))
      select type (t => cur)
      type is (hsd_table)
        if (allocated(t%name)) then
          childName = tolower(t%name)
          if (childName == tolower(name)) cnt = cnt + 1
        end if
      end select
    end do

    ! Allocate the list
    allocate(children)
    children%count = cnt
    if (cnt > 0) then
      allocate(children%items(cnt))
    end if

    ! Second pass: fill pointers
    idx = 0
    call iter%init(node)
    do while (iter%next(cur))
      select type (t => cur)
      type is (hsd_table)
        if (allocated(t%name)) then
          childName = tolower(t%name)
          if (childName == tolower(name)) then
            idx = idx + 1
            children%items(idx)%ptr => t
          end if
        end if
      end select
    end do

  end subroutine getChildren

  !> Get the number of items in a child list.
  pure function getLength(children) result(length)
    type(hsd_child_list), pointer, intent(in) :: children
    integer :: length

    if (associated(children)) then
      length = children%count
    else
      length = 0
    end if

  end function getLength

  !> Get the ii-th item from a child list (1-based indexing).
  subroutine getItem1(children, ii, child)
    type(hsd_child_list), pointer, intent(in) :: children
    integer, intent(in) :: ii
    type(hsd_table), pointer, intent(out) :: child

    if (.not. associated(children) .or. ii < 1 .or. ii > children%count) then
      child => null()
      return
    end if
    child => children%items(ii)%ptr

  end subroutine getItem1

  !> Free a child list (replaces legacy destroyNodeList).
  !>
  !> Note: This only frees the list container, not the child nodes themselves
  !> (which are owned by the parent table).
  subroutine destroyNodeList(children)
    type(hsd_child_list), pointer, intent(inout) :: children

    if (associated(children)) then
      if (allocated(children%items)) deallocate(children%items)
      deallocate(children)
    end if
    children => null()

  end subroutine destroyNodeList

  ! ============================================================
  !  splitModifier — split comma-separated modifier string
  ! ============================================================

  !> Split a modifier containing comma-separated list of modifiers into components.
  !>
  !> Replaces the legacy splitModifier that used type(string) arrays.
  !> This version uses allocatable character arrays.
  subroutine splitModifier(modifier, child, modifiers)
    character(len=*), intent(in) :: modifier
    type(hsd_table), intent(in) :: child
    character(len=mc), intent(out) :: modifiers(:)

    integer :: nModif, ii, iStart, iEnd

    nModif = size(modifiers)
    iStart = 1
    do ii = 1, nModif - 1
      iEnd = index(modifier(iStart:), sepModifier)
      if (iEnd == 0) then
        call detailedError(child, "Invalid number of specified modifiers (" &
            & // i2c(ii) // " instead of " // i2c(nModif) // ").")
      end if
      iEnd = iStart + iEnd - 1
      modifiers(ii) = trim(adjustl(modifier(iStart:iEnd-1)))
      iStart = iEnd + 1
    end do
    if (index(modifier(iStart:), sepModifier) /= 0) then
      call detailedError(child, "Invalid number of specified modifiers (" &
          & // "more than " // i2c(nModif) // ").")
    end if
    modifiers(nModif) = trim(adjustl(modifier(iStart:)))

  end subroutine splitModifier

  ! ============================================================
  !  getSelectedAtomIndices / getSelectedIndices
  ! ============================================================

  !> Convert a string containing atom indices, ranges and species names to a list of atom indices.
  !>
  !> This is a copy of the legacy hsdutils version but takes hsd_table instead of fnode.
  subroutine getSelectedAtomIndices(node, selectionExpr, speciesNames, species, selectedIndices,&
        & selectionRange, indexRange)

    !> Top node for detailed errors.
    type(hsd_table), intent(in) :: node

    !> String to convert
    character(len=*), intent(in) :: selectionExpr

    !> Contains the valid species names.
    character(len=*), intent(in) :: speciesNames(:)

    !> Contains for every atom its species index
    integer, intent(in) :: species(:)

    !> Integer list of atom indices on return.
    integer, allocatable, intent(out) :: selectedIndices(:)

    !> The range of indices [from, to] available for selection. Default: [1, size(species)]
    integer, optional, intent(in) :: selectionRange(:)

    !> The range of indices [from, to] available in general.
    integer, optional, intent(in) :: indexRange(:)

    type(TStatus) :: errStatus
    logical, allocatable :: selected(:)
    integer :: selectionRange_(2)
    integer :: ii

    if (present(selectionRange)) then
      selectionRange_(:) = selectionRange
    else
      selectionRange_(:) = [1, size(species)]
    end if

    allocate(selected(selectionRange_(2) - selectionRange_(1) + 1))
    call getIndexSelection(selectionExpr, selectionRange_, selected, errStatus,&
        & indexRange=indexRange, speciesNames=speciesNames, species=species)
    if (errStatus%hasError()) then
      call detailedError(node, "Invalid atom selection expression '" // trim(selectionExpr) &
          & // "': " // errStatus%message)
    end if
    selectedIndices = pack([(ii, ii = selectionRange_(1), selectionRange_(2))], selected)
    if (size(selectedIndices) == 0) then
      call detailedWarning(node, "Atom index selection expression selected no atoms")
    end if

  end subroutine getSelectedAtomIndices

  !> Convert a string containing indices and ranges to a list of indices.
  subroutine getSelectedIndices(node, selectionExpr, selectionRange, selectedIndices, indexRange)

    !> Top node for detailed errors.
    type(hsd_table), intent(in) :: node

    !> String to convert
    character(len=*), intent(in) :: selectionExpr

    !> Range of indices [from, to] available for selection.
    integer, intent(in) :: selectionRange(:)

    !> Integer list of selected indices on return.
    integer, allocatable, intent(out) :: selectedIndices(:)

    !> The range of indices [from, to] available in general.
    integer, optional, intent(in) :: indexRange(:)

    type(TStatus) :: errStatus
    logical, allocatable :: selected(:)
    integer :: ii

    allocate(selected(selectionRange(2) - selectionRange(1) + 1))
    call getIndexSelection(selectionExpr, selectionRange, selected, errStatus,&
        & indexRange=indexRange)
    if (errStatus%hasError()) then
      call detailedError(node, "Invalid index selection expression '" // trim(selectionExpr) &
          & // "': " // errStatus%message)
    end if
    selectedIndices = pack([(ii, ii = selectionRange(1), selectionRange(2))], selected)

  end subroutine getSelectedIndices

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
