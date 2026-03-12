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
!> to direct hsd-fortran API calls.
module dftbp_io_hsdcompat
  use dftbp_common_accuracy, only : dp, lc, mc
  use dftbp_common_status, only : TStatus
  use dftbp_common_unitconversion, only : TUnit, convertUnit, statusCodes
  use dftbp_extlibs_hsd, only : hsd_table, hsd_value, hsd_node, hsd_iterator, &
      & hsd_error_t, &
      & hsd_get, hsd_get_or, hsd_get_or_set, hsd_get_matrix, &
      & hsd_set, hsd_clear_children, &
      & hsd_get_child, hsd_get_table, hsd_has_child, hsd_remove_child, &
      & hsd_child_count, hsd_get_keys, hsd_get_attrib, hsd_has_attrib, hsd_set_attrib, &
      & hsd_rename_child, &
      & new_table, new_value, &
      & hsd_load, hsd_load_string, hsd_dump, hsd_dump_to_string, &
      & hsd_warn_unprocessed, MAX_WARNING_LEN, &
      & HSD_STAT_OK, HSD_STAT_NOT_FOUND, HSD_STAT_TYPE_ERROR
  use dftbp_io_charmanip, only : i2c, tolower, unquote
  use dftbp_io_indexselection, only : getIndexSelection
  use dftbp_io_message, only : error, warning
  use dftbp_io_tokenreader, only : getNextToken, TOKEN_EOS, TOKEN_ERROR, TOKEN_OK
  use dftbp_type_linkedlist, only : append, len, TListComplex, TListComplexR1, TListInt, &
      & TListIntR1, TListReal, TListRealR1, TListString
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
  public :: warnUnprocessedNodes, setUnprocessed, setProcessed
  public :: hasInlineData
  public :: getDescendant, setNodeName, removeChildNodes, destroyNode, removeChild
  public :: renameDescendant, localiseName
  public :: checkError
  public :: getFirstTextChild
  public :: dumpHsd

  !> Re-export core types from hsd-data for convenience
  public :: hsd_table, hsd_value, hsd_node, hsd_iterator, new_table, new_value
  public :: hsd_error_t
  public :: hsd_set, hsd_get, hsd_get_child, hsd_get_table, hsd_has_child
  public :: hsd_remove_child, hsd_rename_child, hsd_set_attrib, hsd_get_attrib
  public :: hsd_child_count, hsd_get_keys, hsd_has_attrib
  public :: hsd_load, hsd_load_string, hsd_dump, hsd_dump_to_string
  public :: HSD_STAT_OK, HSD_STAT_NOT_FOUND, HSD_STAT_TYPE_ERROR

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
    ! Rank-1 arrays with default
    module procedure :: getChVal_intR1_def
    module procedure :: getChVal_realR1_def
    module procedure :: getChVal_logicalR1_def
    ! Rank-2 arrays without default
    module procedure :: getChVal_realR2
    module procedure :: getChVal_intR2
    ! Rank-2 arrays with default
    module procedure :: getChVal_realR2_def
    ! Child table retrieval (for dispatch blocks)
    module procedure :: getChVal_table
    ! Linked list readers (variable-length)
    module procedure :: getChVal_lString
    module procedure :: getChVal_lReal
    module procedure :: getChVal_lRealR1
    module procedure :: getChVal_lInt
    module procedure :: getChVal_lIntR1
    module procedure :: getChVal_lComplex
    module procedure :: getChVal_lComplexR1
    module procedure :: getChVal_lIntR1RealR1
    module procedure :: getChVal_lStringIntR1RealR1
  end interface getChildValue

  !> Generic interface for writing child values.
  !>
  !> Maps legacy setChildValue(fnode, name, val) patterns to hsd_set.
  interface setChildValue
    module procedure :: setChVal_int
    module procedure :: setChVal_real
    module procedure :: setChVal_logical
    module procedure :: setChVal_char
    module procedure :: setChVal_charR1
    module procedure :: setChVal_intR1
    module procedure :: setChVal_realR1
    module procedure :: setChVal_intR2
    module procedure :: setChVal_realR2
    module procedure :: setChVal_intR2RealR2
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

  !> Generic interface for dumping HSD tree (replaces legacy dumpHSD)
  interface dumpHsd
    module procedure :: dumpHsd_file
    module procedure :: dumpHsd_unit
  end interface dumpHsd

contains

  ! ============================================================
  !  Private helpers
  ! ============================================================

  !> Resolve empty name to "#text" for inline-text access.
  !>
  !> The legacy xmlf90 API used getChildValue(node, "", val) to read the text
  !> content of the current element.  In hsd-fortran, inline text is stored as
  !> a child named "#text".  This helper maps "" → "#text" so that
  !> hsd_get / hsd_get_or_set / hsd_get_matrix find the right child.
  pure function textName_(name) result(resolved)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: resolved

    if (len_trim(name) == 0) then
      resolved = "#text"
    else
      resolved = name
    end if

  end function textName_

  !> Add a table child to a parent and return a pointer to the STORED copy.
  !>
  !> add_child() copies via source=, so the local table and the stored table are
  !> distinct objects.  This helper ensures callers get a pointer to the stored
  !> copy so that subsequent modifications (adding sub-children, etc.) affect
  !> the tree.
  subroutine addChildGetPtr_(parent, child, stored)
    type(hsd_table), intent(inout), target :: parent
    type(hsd_table), intent(in) :: child
    type(hsd_table), pointer, intent(out) :: stored

    class(hsd_node), pointer :: storedNode

    call parent%add_child(child)
    call parent%get_child(parent%num_children, storedNode)
    select type (t => storedNode)
    type is (hsd_table)
      stored => t
    class default
      stored => null()
    end select

  end subroutine addChildGetPtr_


  !> Return a pointer to a module-level static pseudo-table named "#text".
  !>
  !> Used by getChVal_table so that getNodeName2(variableValue) returns
  !> textNodeName ("#text") for the dispatch pattern, without modifying
  !> the actual HSD tree.
  function getTextPseudoTable_() result(ptr)
    type(hsd_table), pointer :: ptr

    type(hsd_table), target, save :: singleton
    logical, save :: initialized = .false.

    !$omp critical (getTextPseudoTable_init)
    if (.not. initialized) then
      call new_table(singleton, name=textNodeName)
      initialized = .true.
    end if
    !$omp end critical (getTextPseudoTable_init)
    ptr => singleton

  end function getTextPseudoTable_


  !> Remove any existing child with the given name, then create a new TABLE
  !> child and return a pointer to the stored copy.
  !>
  !> Used by setChVal_* when the caller requests a child pointer. By creating
  !> a TABLE + #text VALUE instead of a bare VALUE, we preserve typed data for
  !> later hsd_get() calls while also providing a valid hsd_table pointer.
  subroutine replaceOrAddTable_(node, name, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(hsd_table), pointer, intent(out) :: child

    integer :: stat
    type(hsd_table) :: newTable

    ! Remove any existing child with this name (VALUE or TABLE)
    call hsd_remove_child(node, name, stat, case_insensitive=.true.)

    ! Create and add a new TABLE
    call new_table(newTable, name=name)
    call addChildGetPtr_(node, newTable, child)

  end subroutine replaceOrAddTable_


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

    call hsd_get(node, textName_(name), val, stat=stat)
    if (stat == HSD_STAT_NOT_FOUND) then
      call error(formatErrMsg_(node, "Missing required integer value: '" // name // "'"))
    else if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Invalid value for '" // name // "': type conversion error"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_int

  !> Get required real(dp) child value
  subroutine getChVal_real(node, name, val, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: val
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    integer :: stat

    call hsd_get(node, textName_(name), val, stat=stat)
    if (stat == HSD_STAT_NOT_FOUND) then
      call error(formatErrMsg_(node, "Missing required real value: '" // name // "'"))
    else if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Invalid value for '" // name // "': type conversion error"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_real

  !> Get required logical child value
  subroutine getChVal_logical(node, name, val, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    logical, intent(out) :: val
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    integer :: stat

    call hsd_get(node, textName_(name), val, stat=stat)
    if (stat == HSD_STAT_NOT_FOUND) then
      call error(formatErrMsg_(node, "Missing required logical value: '" // name // "'"))
    else if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Invalid value for '" // name // "': type conversion error"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

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
      ! First attempt a direct get to handle empty table nodes correctly.
      ! hsd_get_or_set would fall back to the default for an empty table like
      ! "MovedAtoms = {}", because hsd-fortran reports TYPE_ERROR for a table
      ! with no inline text.  The old xmlf90 parser returned "" for such nodes,
      ! so we replicate that behaviour here.
      call hsd_get(node, textName_(name), val, stat=stat)
      if (stat /= HSD_STAT_OK) then
        ! Check whether the child is an EMPTY table (e.g. "MovedAtoms = {}").
        ! hsd_get_or_set would fall back to the default for such nodes because
        ! hsd-fortran reports TYPE_ERROR for a table with no inline text.
        ! The old xmlf90 parser returned "" for empty blocks, so replicate that.
        ! Only truly empty tables (no children) get this treatment — tables WITH
        ! content (e.g. "Prefix = {path/}") should NOT be overridden.
        block
          type(hsd_table), pointer :: tbl
          call hsd_get_table(node, name, tbl)
          if (associated(tbl)) then
            if (tbl%num_children == 0) then
              ! Truly empty table block – return empty string
              val = ""
            else
              ! Table with content (e.g. "Prefix = {path/}") – try extracting text
              call getTextContent_(node, name, val)
              if (len(val) == 0) then
                val = default
                call hsd_set(node, textName_(name), default)
              end if
            end if
          else
            ! Child does not exist at all – use default and write it back
            val = default
            call hsd_set(node, textName_(name), default)
          end if
        end block
      end if
    else
      call hsd_get(node, textName_(name), val, stat=stat)
      if (stat /= HSD_STAT_OK) then
        block
          type(hsd_table), pointer :: tbl
          call hsd_get_table(node, name, tbl)
          if (associated(tbl)) then
            call getTextContent_(node, name, val)
            if (len(val) == 0) then
              call error(formatErrMsg_(node, "Missing required string value: '" // name // "'"))
            end if
          else
            call error(formatErrMsg_(node, "Missing required string value: '" // name // "'"))
          end if
        end block
      end if
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_string

  ! ============================================================
  !  getChildValue — scalar with default
  ! ============================================================

  !> Get integer child value with default (writes default back if absent)
  subroutine getChVal_int_def(node, name, val, default, modifier, child, isDefaultExported)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(out) :: val
    integer, intent(in) :: default
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child
    logical, intent(in), optional :: isDefaultExported

    integer :: stat

    if (present(isDefaultExported)) then
      if (.not. isDefaultExported) then
        call hsd_get_or(node, textName_(name), val, default, stat=stat)
        if (stat == HSD_STAT_TYPE_ERROR) then
          call error(formatErrMsg_(node, "Invalid value for '" // name // "': type conversion error"))
        end if
        if (present(modifier)) then
          call getModifier_(node, name, modifier)
        else
          call checkNoModifier_(node, name)
        end if
        if (present(child)) then
          call hsd_get_table(node, name, child)
          if (.not. associated(child)) child => node
        end if
        return
      end if
    end if
    call hsd_get_or_set(node, textName_(name), val, default, stat=stat)
    if (stat == HSD_STAT_TYPE_ERROR) then
      call error(formatErrMsg_(node, "Invalid value for '" // name // "': type conversion error"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_int_def

  !> Get real(dp) child value with default
  subroutine getChVal_real_def(node, name, val, default, modifier, child, isDefaultExported)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: val
    real(dp), intent(in) :: default
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child
    logical, intent(in), optional :: isDefaultExported

    integer :: stat

    if (present(isDefaultExported)) then
      if (.not. isDefaultExported) then
        call hsd_get_or(node, textName_(name), val, default, stat=stat)
        if (stat == HSD_STAT_TYPE_ERROR) then
          call error(formatErrMsg_(node, "Invalid value for '" // name // "': type conversion error"))
        end if
        if (present(modifier)) then
          call getModifier_(node, name, modifier)
        else
          call checkNoModifier_(node, name)
        end if
        if (present(child)) then
          call hsd_get_table(node, name, child)
          if (.not. associated(child)) child => node
        end if
        return
      end if
    end if
    call hsd_get_or_set(node, textName_(name), val, default, stat=stat)
    if (stat == HSD_STAT_TYPE_ERROR) then
      call error(formatErrMsg_(node, "Invalid value for '" // name // "': type conversion error"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_real_def

  !> Get logical child value with default
  subroutine getChVal_logical_def(node, name, val, default, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    logical, intent(out) :: val
    logical, intent(in) :: default
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    integer :: stat

    call hsd_get_or_set(node, textName_(name), val, default, stat=stat)
    if (stat == HSD_STAT_TYPE_ERROR) then
      call error(formatErrMsg_(node, "Invalid value for '" // name // "': type conversion error"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

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

    call hsd_get(node, textName_(name), tmp, stat=stat)
    if (stat == HSD_STAT_NOT_FOUND) then
      call error(formatErrMsg_(node, "Missing required integer array: '" // name // "'"))
    else if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Invalid value for integer array '" // name // "': type conversion error"))
    end if
    nn = min(size(tmp), size(val))
    val(:nn) = tmp(:nn)
    if (present(nItem)) nItem = size(tmp)
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

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

    call hsd_get(node, textName_(name), tmp, stat=stat)
    if (stat == HSD_STAT_NOT_FOUND) then
      call error(formatErrMsg_(node, "Missing required real array: '" // name // "'"))
    else if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Invalid value for real array '" // name // "': type conversion error"))
    end if
    nn = min(size(tmp), size(val))
    val(:nn) = tmp(:nn)
    if (present(nItem)) nItem = size(tmp)
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

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

    call hsd_get(node, textName_(name), tmp, stat=stat)
    if (stat == HSD_STAT_NOT_FOUND) then
      call error(formatErrMsg_(node, "Missing required logical array: '" // name // "'"))
    else if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Invalid value for logical array '" // name // "': type conversion error"))
    end if
    nn = min(size(tmp), size(val))
    val(:nn) = tmp(:nn)
    if (present(nItem)) nItem = size(tmp)
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

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
    integer :: stat, nrows, ncols

    call hsd_get_matrix(node, textName_(name), tmp, nrows, ncols, stat=stat)
    if (stat == HSD_STAT_NOT_FOUND) then
      call error(formatErrMsg_(node, "Missing required real matrix: '" // name // "'"))
    else if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Invalid value for real matrix '" // name // "': type conversion error"))
    end if
    ! Legacy compatibility: the old xmlf90 code read matrix values sequentially
    ! (row by row from text) and filled Fortran arrays in column-major order.
    ! hsd_get_matrix returns tmp(nrows, ncols) in text layout (rows = lines).
    ! We must flatten in row-major (text) order and refill in column-major.
    call fillColumnMajorFromTextMatrix_real_(tmp, nrows, ncols, val)
    if (present(nItem)) nItem = nrows * ncols
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

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

    call hsd_get_matrix(node, textName_(name), tmp, nrows, ncols, stat=stat)
    if (stat == HSD_STAT_NOT_FOUND) then
      call error(formatErrMsg_(node, "Missing required integer matrix: '" // name // "'"))
    else if (stat /= HSD_STAT_OK) then
      call error(formatErrMsg_(node, "Invalid value for integer matrix '" // name // "': type conversion error"))
    end if
    call fillColumnMajorFromTextMatrix_int_(tmp, nrows, ncols, val)
    if (present(nItem)) nItem = nrows * ncols
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_intR2

  !> Get real(dp) matrix child value with default.
  subroutine getChVal_realR2_def(node, name, val, default, nItem, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: val(:,:)
    real(dp), intent(in) :: default(:,:)
    integer, intent(out), optional :: nItem
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    real(dp), allocatable :: tmp(:,:)
    integer :: stat, nrows, ncols
    integer :: nr, nc

    call hsd_get_matrix(node, textName_(name), tmp, nrows, ncols, stat=stat)
    if (stat == HSD_STAT_TYPE_ERROR) then
      call error(formatErrMsg_(node, "Invalid value for '" // name // "': type conversion error"))
    else if (stat /= HSD_STAT_OK) then
      nr = min(size(default, 1), size(val, 1))
      nc = min(size(default, 2), size(val, 2))
      val(:nr, :nc) = default(:nr, :nc)
      call hsd_set(node, textName_(name), default)
      if (present(nItem)) nItem = size(default)
    else
      call fillColumnMajorFromTextMatrix_real_(tmp, nrows, ncols, val)
      if (present(nItem)) nItem = nrows * ncols
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_realR2_def

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
        ! Not found as a table — check if it exists as a value node.
        ! (Same value/table mismatch handling as getChild.)
        block
          class(hsd_node), pointer :: anyChild
          call hsd_get_child(node, name, anyChild, stat)
          if (stat == HSD_STAT_OK .and. associated(anyChild)) then
            select type (anyChild)
            type is (hsd_value)
              ! Wrap value in a table with #text child (like getChild does)
              block
                character(len=:), allocatable :: valStr
                integer :: getStat
                type(hsd_table) :: wrapTbl
                call anyChild%get_string(valStr, getStat)
                if (getStat /= 0) valStr = ""
                call new_table(wrapTbl, name=anyChild%name, attrib=anyChild%attrib, &
                    & line=anyChild%line)
                block
                  type(hsd_value) :: txt
                  call new_value(txt, name="#text", line=anyChild%line)
                  call txt%set_string(valStr)
                  call wrapTbl%add_child(txt)
                end block
                call hsd_remove_child(node, name, stat, case_insensitive=.true.)
                call addChildGetPtr_(node, wrapTbl, container)
                if (associated(container)) container%processed = .true.
              end block
            class default
              container => null()
              stat = HSD_STAT_NOT_FOUND
            end select
          else
            stat = HSD_STAT_NOT_FOUND
          end if
        end block
      end if

      if (.not. associated(container)) then
        if (present(default)) then
          ! Create the child with default.
          ! In the legacy dispatch pattern, the default string becomes the name
          ! of a child element inside the container (e.g. default="None" creates
          ! a child table named "None" so getNodeName returns "none").
          block
            type(hsd_table) :: defContainer
            call new_table(defContainer, name=tolower(name))
            if (len_trim(default) > 0) then
              block
                type(hsd_table) :: default_child
                call new_table(default_child, name=default)
                call defContainer%add_child(default_child)
              end block
            end if
            call addChildGetPtr_(node, defContainer, container)
          end block
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
      ! Mark the dispatch child as processed so that warnUnprocessedNodes does
      ! not report it as unread.  The caller's select-case dispatch will read
      ! its children (marking them individually); the dispatch child itself is
      ! only used to determine the method name via getNodeName / getNodeName2.
      if (associated(variableValue)) variableValue%processed = .true.
      ! If no table child found but text content exists, return a pointer to a
      ! static pseudo-table named "#text" so that the dispatch pattern
      ! (getNodeName2 + select case) can detect inline text data — matching the
      ! old xmlf90 behavior where TEXT_NODE children were returned as fnode ptrs.
      ! The pseudo-table is NOT added to the tree to avoid polluting it.
      if (.not. associated(variableValue) .and. container%num_children > 0) then
        block
          type(hsd_iterator) :: it
          class(hsd_node), pointer :: cur
          logical :: hasText
          hasText = .false.
          call it%init(container)
          do while (it%next(cur))
            select type (cur)
            type is (hsd_value)
              hasText = .true.
              exit
            end select
          end do
          if (hasText) then
            variableValue => getTextPseudoTable_()
          end if
        end block
      end if
    end if

    if (present(child)) then
      child => container
    end if

    if (present(modifier)) then
      if (associated(container) .and. allocated(container%attrib)) then
        modifier = container%attrib
      else
        modifier = ""
      end if
    else
      if (associated(container) .and. allocated(container%attrib) &
          & .and. len_trim(container%attrib) > 0) then
        call error(formatErrMsg_(container, "Entity is not allowed to have a modifier"))
      end if
    end if

  end subroutine getChVal_table

  ! ============================================================
  !  getChildValue — rank-1 arrays with default
  ! ============================================================

  !> Get integer array child value with default.
  subroutine getChVal_intR1_def(node, name, val, default, nItem, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(out) :: val(:)
    integer, intent(in) :: default(:)
    integer, intent(out), optional :: nItem
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    integer, allocatable :: tmp(:)
    integer :: stat, nn

    call hsd_get(node, textName_(name), tmp, stat=stat)
    if (stat == HSD_STAT_TYPE_ERROR) then
      call error(formatErrMsg_(node, "Invalid value for '" // name // "': type conversion error"))
    else if (stat /= HSD_STAT_OK) then
      nn = min(size(default), size(val))
      val(:nn) = default(:nn)
      call hsd_set(node, textName_(name), default)
      if (present(nItem)) nItem = size(default)
    else
      nn = min(size(tmp), size(val))
      val(:nn) = tmp(:nn)
      if (present(nItem)) nItem = size(tmp)
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_intR1_def

  !> Get real(dp) array child value with default.
  subroutine getChVal_realR1_def(node, name, val, default, nItem, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: val(:)
    real(dp), intent(in) :: default(:)
    integer, intent(out), optional :: nItem
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    real(dp), allocatable :: tmp(:)
    integer :: stat, nn

    call hsd_get(node, textName_(name), tmp, stat=stat)
    if (stat == HSD_STAT_TYPE_ERROR) then
      call error(formatErrMsg_(node, "Invalid value for '" // name // "': type conversion error"))
    else if (stat /= HSD_STAT_OK) then
      nn = min(size(default), size(val))
      val(:nn) = default(:nn)
      call hsd_set(node, textName_(name), default)
      if (present(nItem)) nItem = size(default)
    else
      nn = min(size(tmp), size(val))
      val(:nn) = tmp(:nn)
      if (present(nItem)) nItem = size(tmp)
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_realR1_def

  !> Get logical array child value with default.
  subroutine getChVal_logicalR1_def(node, name, val, default, nItem, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    logical, intent(out) :: val(:)
    logical, intent(in) :: default(:)
    integer, intent(out), optional :: nItem
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    logical, allocatable :: tmp(:)
    integer :: stat, nn

    call hsd_get(node, textName_(name), tmp, stat=stat)
    if (stat == HSD_STAT_TYPE_ERROR) then
      call error(formatErrMsg_(node, "Invalid value for '" // name // "': type conversion error"))
    else if (stat /= HSD_STAT_OK) then
      nn = min(size(default), size(val))
      val(:nn) = default(:nn)
      call hsd_set(node, textName_(name), default)
      if (present(nItem)) nItem = size(default)
    else
      nn = min(size(tmp), size(val))
      val(:nn) = tmp(:nn)
      if (present(nItem)) nItem = size(tmp)
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_logicalR1_def

  ! ============================================================
  !  getChildValue — linked-list readers (variable-length data)
  ! ============================================================

  !> Get a list of strings from a child's text content.
  !>
  !> Reads all whitespace-separated tokens from the named child's text.
  !> Each token is unquoted and appended to variableValue.
  subroutine getChVal_lString(node, name, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(TListString), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text, token
    integer :: iStart, iErr

    call getTextContent_(node, name, text)
    iStart = 1
    call getNextToken(text, token, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(variableValue, trim(unquote(token)))
      call getNextToken(text, token, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call error(formatErrMsg_(node, "Invalid string value in '" // name // "'"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_lString

  !> Get a list of real values from a child's text content.
  subroutine getChVal_lReal(node, name, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(TListReal), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    real(dp) :: buffer
    integer :: iStart, iErr

    call getTextContent_(node, name, text)
    iStart = 1
    call getNextToken(text, buffer, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call error(formatErrMsg_(node, "Invalid real value in '" // name // "'"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_lReal

  !> Get a list of real arrays (rank-1 of given dimension) from a child's text content.
  !>
  !> Reads dim reals at a time as groups, appending each group to the list.
  subroutine getChVal_lRealR1(node, name, dim, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim
    type(TListRealR1), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    real(dp), allocatable :: buffer(:)
    integer :: iStart, iErr, nItem

    allocate(buffer(dim))
    call getTextContent_(node, name, text)
    iStart = 1
    call getNextToken(text, buffer, iStart, iErr, nItem)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr, nItem)
    end do
    if (iErr == TOKEN_ERROR) then
      call error(formatErrMsg_(node, "Invalid real value in '" // name // "'"))
    else if (iErr == TOKEN_EOS .and. nItem /= 0) then
      call error(formatErrMsg_(node, "Unexpected end of data in '" // name // "'"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_lRealR1

  !> Get a list of integer values from a child's text content.
  subroutine getChVal_lInt(node, name, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(TListInt), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    integer :: buffer
    integer :: iStart, iErr

    call getTextContent_(node, name, text)
    iStart = 1
    call getNextToken(text, buffer, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call error(formatErrMsg_(node, "Invalid integer value in '" // name // "'"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_lInt

  !> Get a list of integer arrays (rank-1 of given dimension) from a child's text content.
  subroutine getChVal_lIntR1(node, name, dim, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim
    type(TListIntR1), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    integer, allocatable :: buffer(:)
    integer :: iStart, iErr, nItem

    allocate(buffer(dim))
    call getTextContent_(node, name, text)
    iStart = 1
    call getNextToken(text, buffer, iStart, iErr, nItem)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr, nItem)
    end do
    if (iErr == TOKEN_ERROR) then
      call error(formatErrMsg_(node, "Invalid integer value in '" // name // "'"))
    else if (iErr == TOKEN_EOS .and. nItem /= 0) then
      call error(formatErrMsg_(node, "Unexpected end of data in '" // name // "'"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_lIntR1

  !> Get a list of complex values from a child's text content.
  subroutine getChVal_lComplex(node, name, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(TListComplex), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    complex(dp) :: buffer
    integer :: iStart, iErr

    call getTextContent_(node, name, text)
    iStart = 1
    call getNextToken(text, buffer, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call error(formatErrMsg_(node, "Invalid complex value in '" // name // "'"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_lComplex

  !> Get a list of complex arrays (rank-1 of given dimension) from a child's text content.
  subroutine getChVal_lComplexR1(node, name, dim, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim
    type(TListComplexR1), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    complex(dp), allocatable :: buffer(:)
    integer :: iStart, iErr, nItem

    allocate(buffer(dim))
    call getTextContent_(node, name, text)
    iStart = 1
    call getNextToken(text, buffer, iStart, iErr, nItem)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr, nItem)
    end do
    if (iErr == TOKEN_ERROR) then
      call error(formatErrMsg_(node, "Invalid complex value in '" // name // "'"))
    else if (iErr == TOKEN_EOS .and. nItem /= 0) then
      call error(formatErrMsg_(node, "Unexpected end of data in '" // name // "'"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_lComplexR1

  !> Get paired lists of integer and real arrays from a child's text content.
  !>
  !> Reads alternating groups: dimInt integers then dimReal reals per record.
  subroutine getChVal_lIntR1RealR1(node, name, dimInt, valueInt, dimReal, valueReal, &
      & modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: dimInt
    type(TListIntR1), intent(inout) :: valueInt
    integer, intent(in) :: dimReal
    type(TListRealR1), intent(inout) :: valueReal
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    integer, allocatable :: bufferInt(:)
    real(dp), allocatable :: bufferReal(:)
    integer :: iStart, iErr, nItem

    allocate(bufferInt(dimInt))
    allocate(bufferReal(dimReal))
    call getTextContent_(node, name, text)
    iStart = 1
    iErr = TOKEN_OK
    do while (iErr == TOKEN_OK)
      call getNextToken(text, bufferInt, iStart, iErr, nItem)
      if (iErr == TOKEN_ERROR) then
        call error(formatErrMsg_(node, "Invalid integer value in '" // name // "'"))
      end if
      if (iErr == TOKEN_EOS .and. nItem /= 0) then
        call error(formatErrMsg_(node, "Unexpected end of data in '" // name // "'"))
      end if
      if (iErr == TOKEN_OK) then
        call append(valueInt, bufferInt)
        call getNextToken(text, bufferReal, iStart, iErr, nItem)
        if (iErr == TOKEN_ERROR) then
          call error(formatErrMsg_(node, "Invalid real value in '" // name // "'"))
        end if
        if (iErr == TOKEN_EOS .and. nItem /= 0) then
          call error(formatErrMsg_(node, "Unexpected end of data in '" // name // "'"))
        end if
        if (iErr == TOKEN_OK) then
          call append(valueReal, bufferReal)
        end if
      end if
    end do
    if (len(valueInt) /= len(valueReal)) then
      call error(formatErrMsg_(node, "Unexpected end of data in '" // name // "'"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_lIntR1RealR1

  !> Get combined string, integer-array and real-array lists from a child's text content.
  !>
  !> Reads records consisting of: one string, dimInt integers, dimReal reals.
  subroutine getChVal_lStringIntR1RealR1(node, name, valueStr, dimInt, valueInt, dimReal, &
      & valueReal, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(TListString), intent(inout) :: valueStr
    integer, intent(in) :: dimInt
    type(TListIntR1), intent(inout) :: valueInt
    integer, intent(in) :: dimReal
    type(TListRealR1), intent(inout) :: valueReal
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text, bufferStr
    integer, allocatable :: bufferInt(:)
    real(dp), allocatable :: bufferReal(:)
    integer :: iStart, iErr, nItem

    allocate(bufferInt(dimInt))
    allocate(bufferReal(dimReal))
    call getTextContent_(node, name, text)
    iStart = 1
    iErr = TOKEN_OK
    do while (iErr == TOKEN_OK)
      call getNextToken(text, bufferStr, iStart, iErr)
      if (iErr == TOKEN_ERROR) then
        call error(formatErrMsg_(node, "Invalid string value in '" // name // "'"))
      end if
      if (iErr == TOKEN_EOS) exit
      call append(valueStr, bufferStr)

      call getNextToken(text, bufferInt, iStart, iErr, nItem)
      if (iErr /= TOKEN_OK) then
        call error(formatErrMsg_(node, "Invalid integer value in '" // name // "'"))
      end if
      call append(valueInt, bufferInt)

      call getNextToken(text, bufferReal, iStart, iErr, nItem)
      if (iErr /= TOKEN_OK) then
        call error(formatErrMsg_(node, "Invalid real value in '" // name // "'"))
      end if
      call append(valueReal, bufferReal)
    end do
    if (len(valueStr) /= len(valueInt) .or. len(valueInt) /= len(valueReal)) then
      call error(formatErrMsg_(node, "Unexpected end of data in '" // name // "'"))
    end if
    if (present(modifier)) then
      call getModifier_(node, name, modifier)
    else
      call checkNoModifier_(node, name)
    end if
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if

  end subroutine getChVal_lStringIntR1RealR1

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
  !> @param emptyIfMissing  If .true., create an empty child if not found (instead of error).
  subroutine getChild(node, name, child, requested, modifier, emptyIfMissing)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(hsd_table), pointer, intent(out) :: child
    logical, intent(in), optional :: requested
    character(len=:), allocatable, intent(out), optional :: modifier
    logical, intent(in), optional :: emptyIfMissing

    integer :: stat
    logical :: isRequired
    logical :: emptyIfMissing_
    class(hsd_node), pointer :: anyChild

    isRequired = .true.
    if (present(requested)) isRequired = requested
    emptyIfMissing_ = .false.
    if (present(emptyIfMissing)) emptyIfMissing_ = emptyIfMissing

    call hsd_get_table(node, name, child, stat)

    if (stat /= HSD_STAT_OK) then
      ! Not found as a table — check if it exists as a value node.
      ! The legacy xmlf90 DOM used uniform fnode elements for both blocks and
      ! scalar values.  hsd-fortran distinguishes hsd_table from hsd_value.
      ! To preserve the legacy behaviour where getChild() could locate value
      ! nodes (e.g. `ParserVersion = 4`), we wrap value children in a new
      ! table and splice them into the tree so that subsequent calls for the
      ! same name also succeed.
      call hsd_get_child(node, name, anyChild, stat)
      if (stat == HSD_STAT_OK .and. associated(anyChild)) then
        select type (anyChild)
        type is (hsd_value)
          ! Wrap the value node in a table child so the caller sees an
          ! hsd_table that contains the original value as inline text.
          block
            character(len=:), allocatable :: valStr
            integer :: getStat
            type(hsd_table) :: wrapTbl
            call anyChild%get_string(valStr, getStat)
            if (getStat /= 0) valStr = ""
            call new_table(wrapTbl, name=anyChild%name, attrib=anyChild%attrib, line=anyChild%line)
            block
              type(hsd_value) :: txt
              call new_value(txt, name="#text", line=anyChild%line)
              call txt%set_string(valStr)
              call wrapTbl%add_child(txt)
            end block
            ! Replace the value node with the new table in the parent
            call hsd_remove_child(node, name, stat, case_insensitive=.true.)
            call addChildGetPtr_(node, wrapTbl, child)
            if (associated(child)) child%processed = .true.
          end block
        class default
          child => null()
        end select
      else
        child => null()
      end if
    end if

    if (.not. associated(child)) then
      if (emptyIfMissing_ .and. .not. isRequired) then
        ! Create an empty child block
        block
          type(hsd_table) :: emptyChild
          call new_table(emptyChild, name=tolower(name))
          call addChildGetPtr_(node, emptyChild, child)
        end block
      else if (isRequired) then
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
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: val
    !> Intentionally unused: hsd_set is inherently upsert (always replaces if the key exists),
    !> so no separate replace logic is needed. Kept for API compatibility with the old xmlf90 layer.
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    if (present(child)) then
      ! Create TABLE + typed #text VALUE so that both the child pointer
      ! (for detailedWarning/setUnprocessed) and hsd_get() work correctly.
      call replaceOrAddTable_(node, name, child)
      block
        type(hsd_value) :: txt
        call new_value(txt, name="#text")
        call txt%set_integer(val)
        call child%add_child(txt)
      end block
    else
      call hsd_set(node, textName_(name), val)
    end if

  end subroutine setChVal_int

  !> Set real(dp) child value
  subroutine setChVal_real(node, name, val, replace, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: val
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    if (present(child)) then
      call replaceOrAddTable_(node, name, child)
      block
        type(hsd_value) :: txt
        call new_value(txt, name="#text")
        call txt%set_real(val)
        call child%add_child(txt)
      end block
    else
      call hsd_set(node, textName_(name), val)
    end if

  end subroutine setChVal_real

  !> Set logical child value
  subroutine setChVal_logical(node, name, val, replace, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    logical, intent(in) :: val
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    if (present(child)) then
      call replaceOrAddTable_(node, name, child)
      block
        type(hsd_value) :: txt
        call new_value(txt, name="#text")
        call txt%set_logical(val)
        call child%add_child(txt)
      end block
    else
      call hsd_set(node, textName_(name), val)
    end if

  end subroutine setChVal_logical

  !> Set character child value
  subroutine setChVal_char(node, name, val, replace, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    if (present(child)) then
      call replaceOrAddTable_(node, name, child)
      block
        type(hsd_value) :: txt
        call new_value(txt, name="#text")
        call txt%set_string(val)
        call child%add_child(txt)
      end block
    else
      ! Use textName_ to convert empty name ("") to "#text" for inline text,
      ! since hsd_set with an empty path silently fails.
      call hsd_set(node, textName_(name), val)
    end if

  end subroutine setChVal_char

  !> Set character array child value
  subroutine setChVal_charR1(node, name, val, replace, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val(:)
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: combined
    integer :: ii

    combined = ""
    do ii = 1, size(val)
      if (ii > 1) combined = combined // " "
      combined = combined // trim(val(ii))
    end do

    if (present(child)) then
      call replaceOrAddTable_(node, name, child)
      block
        type(hsd_value) :: txt
        call new_value(txt, name="#text")
        call txt%set_raw(combined)
        call child%add_child(txt)
      end block
    else
      call hsd_set(node, textName_(name), combined)
    end if

  end subroutine setChVal_charR1

  !> Set integer array child value
  subroutine setChVal_intR1(node, name, val, replace, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: val(:)
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    if (present(child)) then
      call replaceOrAddTable_(node, name, child)
      block
        type(hsd_value) :: txt
        character(len=32) :: buffer
        character(len=:), allocatable :: raw
        integer :: ii
        raw = ""
        do ii = 1, size(val)
          write(buffer, '(I0)') val(ii)
          if (ii > 1) raw = raw // " "
          raw = raw // trim(adjustl(buffer))
        end do
        call new_value(txt, name="#text")
        call txt%set_raw(raw)
        call child%add_child(txt)
      end block
    else
      call hsd_set(node, textName_(name), val)
    end if

  end subroutine setChVal_intR1

  !> Set real(dp) array child value
  subroutine setChVal_realR1(node, name, val, replace, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: val(:)
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    if (present(child)) then
      call replaceOrAddTable_(node, name, child)
      block
        type(hsd_value) :: txt
        character(len=32) :: buffer
        character(len=:), allocatable :: raw
        integer :: ii
        raw = ""
        do ii = 1, size(val)
          write(buffer, '(G0)') val(ii)
          if (ii > 1) raw = raw // " "
          raw = raw // trim(adjustl(buffer))
        end do
        call new_value(txt, name="#text")
        call txt%set_raw(raw)
        call child%add_child(txt)
      end block
    else
      call hsd_set(node, textName_(name), val)
    end if

  end subroutine setChVal_realR1

  !> Set integer matrix child value
  subroutine setChVal_intR2(node, name, val, replace, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: val(:,:)
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    if (present(child)) then
      call replaceOrAddTable_(node, name, child)
      block
        type(hsd_value) :: txt
        character(len=32) :: buffer
        character(len=:), allocatable :: raw
        integer :: ir, ic
        raw = ""
        do ic = 1, size(val, 2)
          if (ic > 1) raw = raw // new_line('a')
          do ir = 1, size(val, 1)
            write(buffer, '(I0)') val(ir, ic)
            if (ir > 1) raw = raw // " "
            raw = raw // trim(adjustl(buffer))
          end do
        end do
        call new_value(txt, name="#text")
        call txt%set_raw(raw)
        call child%add_child(txt)
      end block
    else
      ! hsd_set writes rows as text lines, but DFTB+ convention expects
      ! columns as text lines, so transpose before passing to hsd_set.
      call hsd_set(node, textName_(name), transpose(val))
    end if

  end subroutine setChVal_intR2

  !> Set real(dp) matrix child value
  subroutine setChVal_realR2(node, name, val, replace, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: val(:,:)
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    if (present(child)) then
      call replaceOrAddTable_(node, name, child)
      block
        type(hsd_value) :: txt
        character(len=32) :: buffer
        character(len=:), allocatable :: raw
        integer :: ir, ic
        raw = ""
        do ic = 1, size(val, 2)
          if (ic > 1) raw = raw // new_line('a')
          do ir = 1, size(val, 1)
            write(buffer, '(G0)') val(ir, ic)
            if (ir > 1) raw = raw // " "
            raw = raw // trim(adjustl(buffer))
          end do
        end do
        call new_value(txt, name="#text")
        call txt%set_raw(raw)
        call child%add_child(txt)
      end block
    else
      ! hsd_set writes rows as text lines, but DFTB+ convention expects
      ! columns as text lines, so transpose before passing to hsd_set.
      call hsd_set(node, textName_(name), transpose(val))
    end if

  end subroutine setChVal_realR2


  !> Set mixed integer + real rank-2 child value (TypesAndCoordinates pattern)
  subroutine setChVal_intR2RealR2(node, name, intValue, realValue, replace, child, modifier)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: intValue(:,:)
    real(dp), intent(in) :: realValue(:,:)
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child
    character(len=*), intent(in), optional :: modifier

    character(len=100) :: buffer
    character(len=:), allocatable :: strBuffer
    integer :: nRow, nCol1, nCol2, ii, jj

    nRow = size(intValue, dim=2)
    nCol1 = size(intValue, dim=1)
    nCol2 = size(realValue, dim=1)
    strBuffer = ""
    do ii = 1, nRow
      do jj = 1, nCol1
        write(buffer, *) intValue(jj, ii)
        strBuffer = strBuffer // " " // trim(adjustl(buffer))
      end do
      do jj = 1, nCol2
        write(buffer, *) realValue(jj, ii)
        strBuffer = strBuffer // " " // trim(adjustl(buffer))
      end do
      if (ii < nRow) strBuffer = strBuffer // new_line('a')
    end do

    if (present(child)) then
      call replaceOrAddTable_(node, name, child)
      block
        type(hsd_value) :: txt
        call new_value(txt, name="#text")
        call txt%set_raw(trim(adjustl(strBuffer)))
        call child%add_child(txt)
      end block
      if (present(modifier)) child%attrib = modifier
    else
      call hsd_set(node, textName_(name), trim(adjustl(strBuffer)))
      if (present(modifier)) call hsd_set_attrib(node, name, modifier)
    end if

  end subroutine setChVal_intR2RealR2


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
    integer :: status, setStat

    changed_ = len_trim(modifier) > 0
    if (changed_) then
      call convertUnit(units, modifier, convertValue, status)
      if (status /= statusCodes%ok) then
        call detailedError(child, MSG_INVALID_MODIFIER // modifier // "'")
      end if
      if (present(replace)) then
        if (replace) then
          call hsd_set(child, "#text", convertValue, stat=setStat)
          if (setStat /= 0) call detailedError(child, "Failed to write converted value back to HSD tree")
        end if
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
    integer :: status, setStat

    changed_ = len_trim(modifier) > 0
    if (changed_) then
      call convertUnit(units, modifier, convertValue, status)
      if (status /= statusCodes%ok) then
        call detailedError(child, MSG_INVALID_MODIFIER // modifier // "'")
      end if
      if (present(replace)) then
        if (replace) then
          call hsd_set(child, "#text", convertValue, stat=setStat)
          if (setStat /= 0) call detailedError(child, "Failed to write converted value back to HSD tree")
        end if
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
    integer :: status, setStat

    changed_ = len_trim(modifier) > 0
    if (changed_) then
      call convertUnit(units, modifier, convertValue, status)
      if (status /= statusCodes%ok) then
        call detailedError(child, MSG_INVALID_MODIFIER // modifier // "'")
      end if
      if (present(replace)) then
        if (replace) then
          call hsd_set(child, "#text", transpose(convertValue), stat=setStat)
          if (setStat /= 0) call detailedError(child, "Failed to write converted value back to HSD tree")
        end if
      end if
    end if
    if (present(changed)) changed = changed_

  end subroutine convertUnitHsd_R2

  ! ============================================================
  !  getNodeName2 — get node name safely
  ! ============================================================

  !> Get the name of an HSD node (replaces legacy getNodeName / getNodeName2).
  !>
  !> Returns the node's name as an allocatable character string, lowercased
  !> to match the legacy xmlf90 HSD parser behavior.
  !> If node pointer is null, returns "".
  subroutine getNodeName2(node, nodeName)
    type(hsd_table), pointer, intent(in) :: node
    character(len=:), allocatable, intent(out) :: nodeName

    if (.not. associated(node)) then
      nodeName = ""
      return
    end if

    if (allocated(node%name)) then
      nodeName = tolower(node%name)
    else
      nodeName = ""
    end if

  end subroutine getNodeName2


  !> Check whether inline data was provided in a getChildValue table result.
  !>
  !> After calling getChildValue(node, name, variableValue, "", child=child),
  !> getNodeName2(variableValue) returns "" for BOTH empty defaults AND inline
  !> text data (because variableValue is null in both cases). This function
  !> distinguishes the two by checking if the container (child) has text content.
  !>
  !> Usage:
  !>   call getChildValue(node, "Velocities", value1, "", child=child, ...)
  !>   if (.not. associated(value1) .and. .not. hasInlineData(child)) then
  !>     ! Truly empty — no data provided
  !>   else
  !>     ! Data provided (inline text or dispatch method)
  !>   end if
  function hasInlineData(container) result(has)
    type(hsd_table), pointer, intent(in) :: container
    logical :: has

    if (.not. associated(container)) then
      has = .false.
    else
      has = hasTextChildren_(container)
    end if

  end function hasInlineData

  ! ============================================================
  !  warnUnprocessedNodes — detect unused input nodes
  ! ============================================================

  !> Warn about unprocessed nodes in the HSD tree.
  !>
  !> Walks the tree and collects any nodes that were not accessed during parsing.
  !> If tIgnoreUnprocessed is .true., warnings are suppressed; if .false. (the
  !> default) and unprocessed nodes exist, the program stops with an error.
  subroutine warnUnprocessedNodes(node, tIgnoreUnprocessed)
    type(hsd_table), intent(in) :: node
    logical, intent(in), optional :: tIgnoreUnprocessed

    character(len=MAX_WARNING_LEN), allocatable :: warnings(:)
    logical :: tIgnore
    integer :: ii

    tIgnore = .false.
    if (present(tIgnoreUnprocessed)) tIgnore = tIgnoreUnprocessed

    call hsd_warn_unprocessed(node, warnings)

    if (size(warnings) > 0) then
      if (tIgnore) then
        do ii = 1, size(warnings)
          call warning(trim(warnings(ii)))
        end do
      else
        do ii = 1, size(warnings)
          call warning(trim(warnings(ii)))
        end do
        call error("The following " // i2c(size(warnings)) &
            & // " node(s) have been ignored by the parser.")
      end if
    end if

  end subroutine warnUnprocessedNodes

  ! ============================================================
  !  getNodeName — get raw node name (legacy getNodeName from xmlf90)
  ! ============================================================

  !> Get the raw name of a node.
  !>
  !> For hsd_table nodes, returns the node name. For the dispatch pattern where
  !> the table pointer is null (indicating inline value data), returns "#text".
  !> This matches the legacy xmlf90 getNodeName behavior where names were lowercased.
  subroutine getNodeName(node, nodeName)
    type(hsd_table), pointer, intent(in) :: node
    character(len=:), allocatable, intent(out) :: nodeName

    if (.not. associated(node)) then
      nodeName = textNodeName
    else if (allocated(node%name)) then
      nodeName = tolower(node%name)
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
  !> If node pointer is null, returns "".
  subroutine getNodeHSDName(node, nodeName)
    type(hsd_table), pointer, intent(in) :: node
    character(len=:), allocatable, intent(out) :: nodeName

    if (.not. associated(node)) then
      nodeName = ""
      return
    end if

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

    type(hsd_table) :: newTable
    class(hsd_node), pointer :: storedChild

    ! If replace requested, remove any existing child with this name first
    if (present(replace)) then
      if (replace) then
        call hsd_remove_child(node, name, case_insensitive=.true.)
      end if
    end if

    ! Create a table and add it (add_child makes a copy via source=)
    call new_table(newTable, name=name)
    call node%add_child(newTable)

    ! Get pointer to the STORED copy (not the local original)
    call node%get_child(node%num_children, storedChild)
    select type (t => storedChild)
    type is (hsd_table)
      child => t
    class default
      child => null()
    end select

    if (present(modifier)) then
      if (len_trim(modifier) > 0) then
        child%attrib = modifier
      end if
    end if

  end subroutine setChild

  ! ============================================================
  !  setUnprocessed — clear processed flag
  ! ============================================================

  !> Mark a node as unprocessed.
  !>
  !> Clears the processed flag so that warnUnprocessedNodes will report the node
  !> as unread. Used by oldcompat when restructuring nodes for re-parsing.
  subroutine setUnprocessed(node)
    type(hsd_table), intent(inout) :: node

    node%processed = .false.

  end subroutine setUnprocessed

  !> Mark a node (and optionally all descendants) as processed.
  !>
  !> Sets the processed flag so that warnUnprocessedNodes will not report the
  !> node as unread. When recursive is .true., all table children are also
  !> marked recursively.
  recursive subroutine setProcessed(node, recursive)
    type(hsd_table), pointer, intent(in) :: node
    logical, intent(in), optional :: recursive

    integer :: ii
    class(hsd_node), pointer :: ch
    logical :: doRecurse

    if (.not. associated(node)) return
    node%processed = .true.

    doRecurse = .false.
    if (present(recursive)) doRecurse = recursive
    if (doRecurse) then
      do ii = 1, node%num_children
        call node%get_child(ii, ch)
        if (.not. associated(ch)) cycle
        ch%processed = .true.
        select type (t => ch)
        type is (hsd_table)
          block
            type(hsd_table), pointer :: tptr
            tptr => t
            call setProcessed(tptr, recursive=.true.)
          end block
        end select
      end do
    end if

  end subroutine setProcessed

  !> Check tokenization error flag and report detailed error if needed.
  !>
  !> Replaces the legacy checkError from hsdutils.
  subroutine checkError(node, iErr, msg)
    type(hsd_table), intent(in) :: node
    integer, intent(in) :: iErr
    character(len=*), intent(in) :: msg

    if (iErr == TOKEN_ERROR) then
      call detailedError(node, msg)
    else if (iErr == TOKEN_EOS) then
      call detailedError(node, "Unexpected end of data")
    end if

  end subroutine checkError

  ! ============================================================
  !  getDescendant — path-based child lookup (oldcompat support)
  ! ============================================================

  !> Navigate to a descendant node along a "/" separated path.
  !>
  !> Replaces the legacy getDescendant(root, path, child, requested, parent).
  subroutine getDescendant(root, path, child, requested, parent)
    type(hsd_table), intent(inout), target :: root
    character(len=*), intent(in) :: path
    type(hsd_table), pointer, intent(out) :: child
    logical, intent(in), optional :: requested
    type(hsd_table), pointer, intent(out), optional :: parent

    integer :: iSlash
    logical :: isRequired
    type(hsd_table), pointer :: cur
    character(len=:), allocatable :: remaining, segment

    isRequired = .false.
    if (present(requested)) isRequired = requested

    ! Walk down the path "/" by "/"
    cur => root
    remaining = path
    child => null()

    do while (len_trim(remaining) > 0)
      iSlash = index(remaining, "/")
      if (iSlash > 0) then
        segment = remaining(:iSlash - 1)
        remaining = remaining(iSlash + 1:)
      else
        segment = remaining
        remaining = ""
      end if

      if (present(parent)) parent => cur

      ! Use getChild which handles VALUE→TABLE wrapping, so that
      ! descendant paths like "Options/CalculateForces" work even when
      ! the leaf is a VALUE node (the legacy xmlf90 DOM made no such
      ! distinction).
      call getChild(cur, segment, child, requested=.false.)
      if (.not. associated(child)) then
        if (isRequired) then
          call error(formatErrMsg_(root, "Required path not found: '" // path // "'"))
        end if
        return
      end if
      cur => child
    end do

  end subroutine getDescendant

  ! ============================================================
  !  setNodeName — rename a node directly
  ! ============================================================

  !> Rename a node (replaces legacy setNodeName from hsdutils2).
  !>
  !> Directly modifies the node's name field. If parent is provided,
  !> invalidates the parent's hash index so lookups use the new name.
  subroutine setNodeName(node, name, updateHsdName, parent)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: updateHsdName
    type(hsd_table), intent(inout), optional :: parent

    node%name = name
    if (present(parent)) call parent%invalidate_index()

  end subroutine setNodeName


  !> Rename a descendant node by path, handling both table and value nodes.
  !>
  !> This wraps hsd_rename_child with path-based navigation and issues a warning.
  !> Unlike getDescendant + setNodeName, this works for value nodes too.
  subroutine renameDescendant(root, path, newName, warningMsg)
    type(hsd_table), intent(inout), target :: root
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: newName
    character(len=*), intent(in), optional :: warningMsg

    integer :: stat

    call hsd_rename_child(root, path, newName, stat, case_insensitive=.true.)
    if (stat == HSD_STAT_OK) then
      if (present(warningMsg)) then
        call detailedWarning(root, warningMsg)
      end if
    end if

  end subroutine renameDescendant


  ! ============================================================
  !  removeChildNodes — remove all children from a node
  ! ============================================================

  !> Remove all children from a node.
  !>
  !> Replaces the legacy removeChildNodes from xmlutils.
  subroutine removeChildNodes(node)
    type(hsd_table), intent(inout), target :: node

    call hsd_clear_children(node)

  end subroutine removeChildNodes

  ! ============================================================
  !  destroyNode — remove a child from its parent
  ! ============================================================

  !> Remove/destroy a child node. This is a no-op stub — hsd_table nodes
  !> are managed by their parent. The caller should use hsd_remove_child
  !> on the parent if actual removal is needed.
  subroutine destroyNode(node)
    type(hsd_table), pointer, intent(inout) :: node

    ! In hsd-fortran, nodes are owned by their parent table.
    ! We just nullify the pointer. The actual memory is managed by the parent.
    node => null()

  end subroutine destroyNode

  ! ============================================================
  !  removeChild — remove a specific child from its parent (compat)
  ! ============================================================

  !> Removes a specific child node from its parent, compatible with legacy
  !> `dummy => removeChild(parent, child)` pattern.  Returns a null pointer
  !> (the old xmlf90 returned the removed node, but callers always discard it).
  function removeChild(parent, child) result(removed)
    type(hsd_table), pointer, intent(inout) :: parent
    type(hsd_table), pointer, intent(inout) :: child
    type(hsd_table), pointer :: removed

    if (associated(child)) then
      call hsd_remove_child(parent, child%name)
      child => null()
    end if
    removed => null()

  end function removeChild

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
  !>
  !> When name is empty, returns the attrib of the parent node itself (legacy
  !> behavior for inline-value patterns like ``getChildValue(node, "", val,
  !> modifier=modifier)``).  When name is non-empty, returns the attrib of the
  !> named child.
  subroutine getModifier_(parent, name, modifier)
    type(hsd_table), intent(in), target :: parent
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(out) :: modifier

    integer :: stat

    if (len_trim(name) == 0) then
      ! Empty name: return the attrib of the node itself.
      ! This matches the legacy xmlf90 behavior where the modifier lived on the
      ! enclosing element (e.g. Temperature [Kelvin] = 1273.15).
      if (allocated(parent%attrib)) then
        modifier = parent%attrib
      else
        modifier = ""
      end if
    else
      call hsd_get_attrib(parent, name, modifier, stat)
      if (stat /= HSD_STAT_OK) then
        modifier = ""
      end if
    end if

  end subroutine getModifier_


  !> Check that a node does not carry a modifier (attrib).
  !>
  !> If the named child (or the node itself when name is empty) has a non-empty
  !> attrib, the user supplied a modifier (e.g. [Kelvin]) where none is accepted.
  !> This matches the legacy MSG_NOMODIFIER check.
  subroutine checkNoModifier_(parent, name)
    type(hsd_table), intent(in), target :: parent
    character(len=*), intent(in) :: name

    character(len=:), allocatable :: attrib
    integer :: stat

    if (len_trim(name) == 0) then
      if (allocated(parent%attrib) .and. len_trim(parent%attrib) > 0) then
        call error(formatErrMsg_(parent, "Entity is not allowed to have a modifier"))
      end if
    else
      call hsd_get_attrib(parent, name, attrib, stat)
      if (stat == HSD_STAT_OK .and. len_trim(attrib) > 0) then
        call error(formatErrMsg_(parent, "Entity '" // name // "' is not allowed to have a modifier"))
      end if
    end if

  end subroutine checkNoModifier_


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


  !> Check whether a table node has any text (value) children.
  !>
  !> Returns .true. if the node contains at least one hsd_value child with
  !> non-empty text content. Used to distinguish empty default containers
  !> from containers with inline text data.
  function hasTextChildren_(node) result(has)
    type(hsd_table), intent(in) :: node
    logical :: has

    type(hsd_iterator) :: iter
    class(hsd_node), pointer :: cur

    has = .false.
    call iter%init(node)
    do while (iter%next(cur))
      select type (cur)
      type is (hsd_value)
        has = .true.
        return
      end select
    end do

  end function hasTextChildren_


  !> Get the text content of a node (public wrapper around getTextContent_).
  !>
  !> Replaces the legacy xmlf90 getFirstTextChild + getNodeValue pattern.
  !> Returns the concatenated text content of all unnamed value children.
  subroutine getFirstTextChild(node, text)
    type(hsd_table), intent(inout), target :: node
    character(len=:), allocatable, intent(out) :: text

    call getTextContent_(node, "", text)

  end subroutine getFirstTextChild


  !> Get the text content of a named child from an hsd_table.
  !>
  !> This extracts the raw text from either a named child value node or,
  !> if name is empty, from the first unnamed value child of the table.
  !> It collects all unnamed value children's text content (for multi-line data).
  subroutine getTextContent_(node, name, text)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(out) :: text

    type(hsd_table), pointer :: container
    type(hsd_iterator) :: iter
    class(hsd_node), pointer :: cur
    character(len=:), allocatable :: piece
    integer :: stat

    text = ""

    if (len_trim(name) == 0) then
      container => node
    else
      call hsd_get_table(node, name, container, stat)
      if (stat /= HSD_STAT_OK) then
        ! Try getting as a simple value
        call hsd_get(node, name, text, stat=stat)
        if (stat /= HSD_STAT_OK) then
          call error(formatErrMsg_(node, "Missing required data: '" // name // "'"))
        end if
        return
      end if
    end if

    ! Collect all unnamed or #text value children's text
    call iter%init(container)
    do while (iter%next(cur))
      select type (v => cur)
      type is (hsd_value)
        if (.not. allocated(v%name) .or. len_trim(v%name) == 0 &
            & .or. v%name == "#text") then
          call v%get_string(piece, stat)
          if (stat == HSD_STAT_OK .and. allocated(piece)) then
            if (len_trim(piece) > 0) then
              if (len(text) > 0) then
                text = text // " " // piece
              else
                text = piece
              end if
            end if
          end if
        end if
      end select
    end do

    if (len(text) == 0) then
      ! Fallback: try hsd_get as string (handles simple named values)
      call hsd_get(container, "#text", text, stat=stat)
      if (stat /= HSD_STAT_OK) text = ""
    end if

  end subroutine getTextContent_

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


  ! ============================================================
  !  dumpHsd — dump HSD tree (replaces legacy dumpHSD)
  ! ============================================================

  !> Dump HSD tree to a named file
  subroutine dumpHsd_file(table, fileName)
    type(hsd_table), intent(in) :: table
    character(len=*), intent(in) :: fileName

    type(hsd_error_t), allocatable :: err

    call hsd_dump(table, fileName, err)
    if (allocated(err)) then
      call error("Error writing HSD file '" // fileName // "'")
    end if

  end subroutine dumpHsd_file


  !> Dump HSD tree to an open file unit
  subroutine dumpHsd_unit(table, unit)
    type(hsd_table), intent(in) :: table
    integer, intent(in) :: unit

    character(len=:), allocatable :: str

    call hsd_dump_to_string(table, str)
    write(unit, '(A)') str

  end subroutine dumpHsd_unit


  ! ============================================================
  !  fillColumnMajorFromTextMatrix — legacy matrix-fill helpers
  ! ============================================================

  !> Fill a Fortran rank-2 real array in column-major order from a text-layout matrix.
  !>
  !> hsd_get_matrix returns tmp(nrows, ncols) where rows correspond to text lines
  !> and columns to values within each line. The legacy xmlf90 code read all values
  !> sequentially (row by row) and filled the Fortran array in column-major order.
  !> This helper reproduces that behavior.
  subroutine fillColumnMajorFromTextMatrix_real_(tmp, nrows, ncols, val)
    real(dp), intent(in) :: tmp(:,:)
    integer, intent(in) :: nrows, ncols
    real(dp), intent(out) :: val(:,:)

    integer :: totalText, totalVal, fillSize, k, i, j

    totalText = nrows * ncols
    totalVal = size(val, 1) * size(val, 2)
    fillSize = min(totalText, totalVal)

    ! Fill val in column-major order from text values in row-major order
    k = 0
    outer: do j = 1, size(val, 2)
      do i = 1, size(val, 1)
        k = k + 1
        if (k > fillSize) exit outer
        ! k-th value in text order is at tmp row (k-1)/ncols+1, col mod(k-1,ncols)+1
        val(i, j) = tmp((k - 1) / ncols + 1, mod(k - 1, ncols) + 1)
      end do
    end do outer

  end subroutine fillColumnMajorFromTextMatrix_real_

  !> Fill a Fortran rank-2 integer array in column-major order from a text-layout matrix.
  subroutine fillColumnMajorFromTextMatrix_int_(tmp, nrows, ncols, val)
    integer, intent(in) :: tmp(:,:)
    integer, intent(in) :: nrows, ncols
    integer, intent(out) :: val(:,:)

    integer :: totalText, totalVal, fillSize, k, i, j

    totalText = nrows * ncols
    totalVal = size(val, 1) * size(val, 2)
    fillSize = min(totalText, totalVal)

    k = 0
    outer: do j = 1, size(val, 2)
      do i = 1, size(val, 1)
        k = k + 1
        if (k > fillSize) exit outer
        val(i, j) = tmp((k - 1) / ncols + 1, mod(k - 1, ncols) + 1)
      end do
    end do outer

  end subroutine fillColumnMajorFromTextMatrix_int_


  !> Legacy alias for renameChildren/hsd_rename_child.
  !>
  !> Renames a child node (e.g. American → British spelling). If the child
  !> does not exist the call is silently ignored.
  subroutine localiseName(node, localName, anglicisedName)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: localName
    character(len=*), intent(in) :: anglicisedName

    integer :: stat

    call hsd_rename_child(node, localName, anglicisedName, stat, case_insensitive=.true.)

  end subroutine localiseName


end module dftbp_io_hsdcompat
