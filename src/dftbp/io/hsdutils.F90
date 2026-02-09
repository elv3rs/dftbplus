!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> DFTB+-specific IO utilities for HSD input parsing.
!>
!> This module provides all DFTB+-specific convenience wrappers around
!> hsd-fortran's low-level tree API:
!>
!>   - dftbp_error / dftbp_warning — formatted error/warning with node context
!>   - getChildValue / getChild / setChildValue / setChild — high-level accessors
!>     with error reporting, modifier extraction, and default write-back
!>   - getDescendant / renameDescendant / setNodeName / removeChild — tree manipulation
!>   - Atom/index selection helpers
!>   - Node name and modifier extraction helpers
!>   - Processed-flag management (setProcessed / setUnprocessed)
!>
!> This module is a permanent part of the DFTB+ IO layer.  It has no legacy-API
!> surface — callers use hsd-fortran types natively.
module dftbp_io_hsdutils
  use dftbp_common_accuracy, only : dp, lc, mc
  use dftbp_common_status, only : TStatus
  use hsd, only : hsd_table, hsd_value, hsd_node, hsd_node_ptr, hsd_iterator, &
      & hsd_format_error, hsd_format_warning, hsd_get, hsd_get_or, hsd_get_or_set, &
      & hsd_get_matrix, hsd_get_child, hsd_get_table, hsd_get_attrib, hsd_has_child, &
      & hsd_set, hsd_set_attrib, hsd_remove_child, hsd_rename_child, hsd_clear_children, &
      & HSD_STAT_OK, HSD_STAT_NOT_FOUND, HSD_STAT_TYPE_ERROR
  use hsd_data, only : new_table, new_value
  use dftbp_io_charmanip, only : i2c, tolower
  use dftbp_io_indexselection, only : getIndexSelection
  use dftbp_io_message, only : error, warning
  use dftbp_io_tokenreader, only : TOKEN_EOS, TOKEN_ERROR
  implicit none

  private

  ! --- Public: error / warning helpers ---
  public :: dftbp_error, dftbp_warning

  ! --- Public: getChildValue / getChild / setChildValue / setChild ---
  public :: getChildValue, getChild, setChildValue, setChild

  ! --- Public: tree manipulation ---
  public :: getDescendant, renameDescendant, setNodeName, removeChild

  ! --- Public: atom / index selection ---
  public :: getSelectedAtomIndices, getSelectedIndices

  ! --- Public: modifier / name helpers ---
  public :: splitModifier, getModifier
  public :: getNodeName, getNodeName2, getNodeHSDName

  ! --- Public: processed flags ---
  public :: setProcessed, setUnprocessed

  ! --- Public: misc ---
  public :: checkError, hasInlineData, getFirstTextChild
  public :: textNodeName

  ! ============================================================
  !  Constants
  ! ============================================================

  !> Name used for inline text data nodes
  character(len=*), parameter :: textNodeName = "#text"

  !> Separator character for splitting modifiers
  character, parameter :: sepModifier = ","

  ! ============================================================
  !  Types
  ! ============================================================

  ! ============================================================
  !  Generic interfaces
  ! ============================================================

  !> Generic interface for reading child values.
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
  end interface getChildValue

  !> Generic interface for writing child values.
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

contains

  ! ============================================================
  !  Private helpers
  ! ============================================================

  !> Resolve empty name to "#text" for inline-text access.
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


  !> Remove any existing child with the given name, then create a new TABLE
  !> child and return a pointer to the stored copy.
  subroutine replaceOrAddTable_(node, name, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(hsd_table), pointer, intent(out) :: child

    integer :: stat
    type(hsd_table) :: newTable

    call hsd_remove_child(node, name, stat, case_insensitive=.true.)
    call new_table(newTable, name=name)
    call addChildGetPtr_(node, newTable, child)

  end subroutine replaceOrAddTable_


  !> Mark a named child of a parent table as processed.
  subroutine markChildProcessed_(node, name)
    type(hsd_table), intent(in) :: node
    character(len=*), intent(in) :: name

    class(hsd_node), pointer :: childNode
    integer :: stat

    call hsd_get_child(node, name, childNode, stat)
    if (stat == HSD_STAT_OK .and. associated(childNode)) then
      childNode%processed = .true.
    end if

  end subroutine markChildProcessed_


  !> Get the first table child of a node (for dispatch patterns).
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


  !> Get the text content of a named child from an hsd_table.
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
        call hsd_get(node, name, text, stat=stat)
        if (stat /= HSD_STAT_OK) then
          call dftbp_error(node, "Missing required data: '" // name // "'")
        end if
        return
      end if
    end if

    call iter%init(container)
    do while (iter%next(cur))
      select type (v => cur)
      type is (hsd_value)
        if (.not. allocated(v%name) .or. len_trim(v%name) == 0 &
            & .or. v%name == "#text") then
          call v%get_string(piece, stat)
          if (stat == HSD_STAT_OK .and. allocated(piece)) then
            if (len(text) > 0) then
              text = text // " " // piece
            else
              text = piece
            end if
          end if
        end if
      end select
    end do

    if (len(text) == 0) then
      call hsd_get(container, "#text", text, stat=stat)
      if (stat /= HSD_STAT_OK) text = ""
    end if

  end subroutine getTextContent_


  !> Fill a Fortran rank-2 real array in column-major order from a text-layout matrix.
  subroutine fillColumnMajorFromTextMatrix_real_(tmp, nrows, ncols, val)
    real(dp), intent(in) :: tmp(:,:)
    integer, intent(in) :: nrows, ncols
    real(dp), intent(out) :: val(:,:)

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


  ! ============================================================
  !  Error / warning helpers
  ! ============================================================

  !> Report a fatal error with node context.
  !>
  !> Thin wrapper: hsd_format_error → DFTB+ error().
  !> Replaces the legacy detailedError(node, msg) pattern.
  subroutine dftbp_error(node, msg)
    type(hsd_table), intent(in) :: node
    character(len=*), intent(in) :: msg

    character(len=:), allocatable :: formatted

    call hsd_format_error(node, msg, formatted)
    call error(formatted)

  end subroutine dftbp_error


  !> Report a warning with node context.
  !>
  !> Thin wrapper: hsd_format_warning → DFTB+ warning().
  !> Replaces the legacy detailedWarning(node, msg) pattern.
  subroutine dftbp_warning(node, msg)
    type(hsd_table), intent(in) :: node
    character(len=*), intent(in) :: msg

    character(len=:), allocatable :: formatted

    call hsd_format_warning(node, msg, formatted)
    call warning(formatted)

  end subroutine dftbp_warning


  ! ============================================================
  !  splitModifier
  ! ============================================================

  !> Split a modifier containing comma-separated list of modifiers into components.
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
        call dftbp_error(child, "Invalid number of specified modifiers (" &
            & // i2c(ii) // " instead of " // i2c(nModif) // ").")
      end if
      iEnd = iStart + iEnd - 1
      modifiers(ii) = trim(adjustl(modifier(iStart:iEnd-1)))
      iStart = iEnd + 1
    end do
    if (index(modifier(iStart:), sepModifier) /= 0) then
      call dftbp_error(child, "Invalid number of specified modifiers (" &
          & // "more than " // i2c(nModif) // ").")
    end if
    modifiers(nModif) = trim(adjustl(modifier(iStart:)))

  end subroutine splitModifier


  ! ============================================================
  !  Atom / index selection helpers
  ! ============================================================

  !> Convert a string containing atom indices, ranges and species names to a list of atom indices.
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

    !> The range of indices [from, to] available for selection.
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
      call dftbp_error(node, "Invalid atom selection expression '" // trim(selectionExpr) &
          & // "': " // errStatus%message)
    end if
    selectedIndices = pack([(ii, ii = selectionRange_(1), selectionRange_(2))], selected)
    if (size(selectedIndices) == 0) then
      call dftbp_warning(node, "Atom index selection expression selected no atoms")
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
      call dftbp_error(node, "Invalid index selection expression '" // trim(selectionExpr) &
          & // "': " // errStatus%message)
    end if
    selectedIndices = pack([(ii, ii = selectionRange(1), selectionRange(2))], selected)

  end subroutine getSelectedIndices


  ! ============================================================
  !  checkError — tokenizer error bridge
  ! ============================================================

  !> Check tokenization error flag and report detailed error if needed.
  subroutine checkError(node, iErr, msg)
    type(hsd_table), intent(in) :: node
    integer, intent(in) :: iErr
    character(len=*), intent(in) :: msg

    if (iErr == TOKEN_ERROR) then
      call dftbp_error(node, msg)
    else if (iErr == TOKEN_EOS) then
      call dftbp_error(node, "Unexpected end of data")
    end if

  end subroutine checkError


  ! ============================================================
  !  hasInlineData — check for text value children
  ! ============================================================

  !> Check whether a container table has any value children (inline data).
  function hasInlineData(container) result(has)
    type(hsd_table), pointer, intent(in) :: container
    logical :: has

    type(hsd_iterator) :: iter
    class(hsd_node), pointer :: cur

    if (.not. associated(container)) then
      has = .false.
      return
    end if

    has = .false.
    call iter%init(container)
    do while (iter%next(cur))
      select type (cur)
      type is (hsd_value)
        has = .true.
        return
      end select
    end do

  end function hasInlineData


  ! ============================================================
  !  getFirstTextChild — concatenated text content
  ! ============================================================

  !> Get the concatenated text content of all unnamed/text value children.
  subroutine getFirstTextChild(node, text)
    type(hsd_table), intent(inout), target :: node
    character(len=:), allocatable, intent(out) :: text

    type(hsd_iterator) :: iter
    class(hsd_node), pointer :: cur
    character(len=:), allocatable :: piece
    integer :: stat

    text = ""
    call iter%init(node)
    do while (iter%next(cur))
      select type (v => cur)
      type is (hsd_value)
        if (.not. allocated(v%name) .or. len_trim(v%name) == 0 &
            & .or. v%name == "#text") then
          call v%get_string(piece, stat)
          if (stat == HSD_STAT_OK .and. allocated(piece)) then
            if (len(text) > 0) then
              text = text // " " // piece
            else
              text = piece
            end if
          end if
        end if
      end select
    end do

    if (len(text) == 0) then
      call hsd_get(node, "#text", text, stat=stat)
      if (stat /= HSD_STAT_OK) text = ""
    end if

  end subroutine getFirstTextChild


  ! ============================================================
  !  Node name helpers
  ! ============================================================

  !> Get the lowercased name of a node. Null pointer → "#text".
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


  !> Get the lowercased name of a node. Null pointer → "".
  subroutine getNodeName2(node, nodeName)
    type(hsd_table), pointer, intent(in) :: node
    character(len=:), allocatable, intent(out) :: nodeName

    if (.not. associated(node)) then
      nodeName = ""
    else if (allocated(node%name)) then
      nodeName = tolower(node%name)
    else
      nodeName = ""
    end if

  end subroutine getNodeName2


  !> Get the original-case name of a node. Null pointer → "".
  subroutine getNodeHSDName(node, nodeName)
    type(hsd_table), pointer, intent(in) :: node
    character(len=:), allocatable, intent(out) :: nodeName

    if (.not. associated(node)) then
      nodeName = ""
    else if (allocated(node%name)) then
      nodeName = node%name
    else
      nodeName = ""
    end if

  end subroutine getNodeHSDName


  ! ============================================================
  !  Processed flag helpers
  ! ============================================================

  !> Mark a node as unprocessed. Null-safe.
  subroutine setUnprocessed(node)
    type(hsd_table), pointer, intent(in) :: node

    if (.not. associated(node)) return
    node%processed = .false.

  end subroutine setUnprocessed


  !> Mark a node (and optionally all descendants) as processed. Null-safe.
  recursive subroutine setProcessed(node, recursive)
    type(hsd_table), pointer, intent(in) :: node
    logical, intent(in), optional :: recursive

    logical :: doRecurse
    integer :: ii
    class(hsd_node), pointer :: childNode

    if (.not. associated(node)) return
    node%processed = .true.

    doRecurse = .false.
    if (present(recursive)) doRecurse = recursive
    if (.not. doRecurse) return

    do ii = 1, node%num_children
      call node%get_child(ii, childNode)
      if (.not. associated(childNode)) cycle
      childNode%processed = .true.
      select type (t => childNode)
      type is (hsd_table)
        block
          type(hsd_table), pointer :: tPtr
          tPtr => t
          call setProcessed(tPtr, .true.)
        end block
      end select
    end do

  end subroutine setProcessed


  ! ============================================================
  !  Modifier extraction helper
  ! ============================================================

  !> Extract the modifier (attrib) string for a child node.
  !>
  !> If name is empty, returns the attrib of the parent node itself.
  !> Otherwise, uses hsd_get_attrib to find the named child's attrib.
  subroutine getModifier(parent, name, modifier)
    type(hsd_table), intent(in), target :: parent
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(out) :: modifier

    integer :: stat

    if (len_trim(name) == 0) then
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

  end subroutine getModifier


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
    if (stat /= HSD_STAT_OK) then
      call dftbp_error(node, "Missing required integer value: '" // name // "'")
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

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
    if (stat /= HSD_STAT_OK) then
      call dftbp_error(node, "Missing required real value: '" // name // "'")
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

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
    if (stat /= HSD_STAT_OK) then
      call dftbp_error(node, "Missing required logical value: '" // name // "'")
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

  end subroutine getChVal_logical


  !> Get string child value, with optional default.
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
      call hsd_get_or_set(node, textName_(name), val, default)
    else
      call hsd_get(node, textName_(name), val, stat=stat)
      if (stat /= HSD_STAT_OK) then
        call dftbp_error(node, "Missing required string value: '" // name // "'")
      end if
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

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

    call hsd_get_or_set(node, textName_(name), val, default)
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

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

    call hsd_get_or_set(node, textName_(name), val, default)
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

  end subroutine getChVal_real_def


  !> Get logical child value with default
  subroutine getChVal_logical_def(node, name, val, default, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    logical, intent(out) :: val
    logical, intent(in) :: default
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    call hsd_get_or_set(node, textName_(name), val, default)
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

  end subroutine getChVal_logical_def


  ! ============================================================
  !  getChildValue — rank-1 arrays (required, no default)
  ! ============================================================

  !> Get required integer array child value.
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
    if (stat /= HSD_STAT_OK) then
      call dftbp_error(node, "Missing required integer array: '" // name // "'")
    end if
    nn = min(size(tmp), size(val))
    val(:nn) = tmp(:nn)
    if (present(nItem)) nItem = size(tmp)
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

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
    if (stat /= HSD_STAT_OK) then
      call dftbp_error(node, "Missing required real array: '" // name // "'")
    end if
    nn = min(size(tmp), size(val))
    val(:nn) = tmp(:nn)
    if (present(nItem)) nItem = size(tmp)
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

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
    if (stat /= HSD_STAT_OK) then
      call dftbp_error(node, "Missing required logical array: '" // name // "'")
    end if
    nn = min(size(tmp), size(val))
    val(:nn) = tmp(:nn)
    if (present(nItem)) nItem = size(tmp)
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

  end subroutine getChVal_logicalR1


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
    if (stat /= HSD_STAT_OK) then
      nn = min(size(default), size(val))
      val(:nn) = default(:nn)
      call hsd_set(node, name, default)
      if (present(nItem)) nItem = size(default)
    else
      nn = min(size(tmp), size(val))
      val(:nn) = tmp(:nn)
      if (present(nItem)) nItem = size(tmp)
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

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
    if (stat /= HSD_STAT_OK) then
      nn = min(size(default), size(val))
      val(:nn) = default(:nn)
      call hsd_set(node, name, default)
      if (present(nItem)) nItem = size(default)
    else
      nn = min(size(tmp), size(val))
      val(:nn) = tmp(:nn)
      if (present(nItem)) nItem = size(tmp)
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

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
    if (stat /= HSD_STAT_OK) then
      nn = min(size(default), size(val))
      val(:nn) = default(:nn)
      if (present(nItem)) nItem = size(default)
    else
      nn = min(size(tmp), size(val))
      val(:nn) = tmp(:nn)
      if (present(nItem)) nItem = size(tmp)
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

  end subroutine getChVal_logicalR1_def


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
    if (stat /= HSD_STAT_OK) then
      call dftbp_error(node, "Missing required real matrix: '" // name // "'")
    end if
    call fillColumnMajorFromTextMatrix_real_(tmp, nrows, ncols, val)
    if (present(nItem)) nItem = nrows * ncols
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

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
    if (stat /= HSD_STAT_OK) then
      call dftbp_error(node, "Missing required integer matrix: '" // name // "'")
    end if
    call fillColumnMajorFromTextMatrix_int_(tmp, nrows, ncols, val)
    if (present(nItem)) nItem = nrows * ncols
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

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
    if (stat /= HSD_STAT_OK) then
      nr = min(size(default, 1), size(val, 1))
      nc = min(size(default, 2), size(val, 2))
      val(:nr, :nc) = default(:nr, :nc)
      call hsd_set(node, name, reshape(default, [size(default)]))
      if (present(nItem)) nItem = size(default)
    else
      call fillColumnMajorFromTextMatrix_real_(tmp, nrows, ncols, val)
      if (present(nItem)) nItem = nrows * ncols
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) then
      call hsd_get_table(node, name, child)
      if (.not. associated(child)) child => node
    end if
    call markChildProcessed_(node, name)

  end subroutine getChVal_realR2_def


  ! ============================================================
  !  getChildValue — child table retrieval (for dispatch blocks)
  ! ============================================================

  !> Get a child table for dispatch or block reading.
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
      container => node
    else
      call hsd_get_table(node, name, container, stat)
      if (stat /= HSD_STAT_OK) then
        block
          class(hsd_node), pointer :: anyChild
          call hsd_get_child(node, name, anyChild, stat)
          if (stat == HSD_STAT_OK .and. associated(anyChild)) then
            select type (anyChild)
            type is (hsd_value)
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
          call dftbp_error(node, "Missing required block: '" // name // "'")
          return
        else
          container => null()
          if (present(child)) child => null()
          if (present(modifier)) modifier = ""
          return
        end if
      end if
    end if

    if (associated(container)) then
      call getFirstTableChild_(container, variableValue, stat)
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
    end if
    call markChildProcessed_(node, name)

    if (associated(variableValue)) then
      variableValue%processed = .true.
    end if

  end subroutine getChVal_table


  ! ============================================================
  !  getChild — get a child table by name
  ! ============================================================

  !> Get a child table by name.
  !!
  !! Uses hsd_get_table with auto_wrap to automatically promote value nodes
  !! to table nodes with a #text child (handles bare "Key = value" patterns).
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

    isRequired = .true.
    if (present(requested)) isRequired = requested
    emptyIfMissing_ = .false.
    if (present(emptyIfMissing)) emptyIfMissing_ = emptyIfMissing

    call hsd_get_table(node, name, child, stat, auto_wrap=.true.)

    if (.not. associated(child)) then
      if (emptyIfMissing_ .and. .not. isRequired) then
        block
          type(hsd_table) :: emptyChild
          call new_table(emptyChild, name=tolower(name))
          call addChildGetPtr_(node, emptyChild, child)
        end block
      else if (isRequired) then
        call dftbp_error(node, "Missing required block: '" // name // "'")
      end if
    end if
    if (present(modifier) .and. associated(child)) then
      call getModifier(node, name, modifier)
    end if
    if (associated(child)) call markChildProcessed_(node, name)

  end subroutine getChild


  ! ============================================================
  !  setChildValue — write values to tree
  ! ============================================================

  !> Set integer child value
  subroutine setChVal_int(node, name, val, replace, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: val
    logical, intent(in), optional :: replace
    type(hsd_table), pointer, intent(out), optional :: child

    if (present(child)) then
      call replaceOrAddTable_(node, name, child)
      block
        type(hsd_value) :: txt
        call new_value(txt, name="#text")
        call txt%set_integer(val)
        call child%add_child(txt)
      end block
    else
      call hsd_set(node, name, val)
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
      call hsd_set(node, name, val)
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
      call hsd_set(node, name, val)
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
      call hsd_set(node, name, val)
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
      call hsd_set(node, name, combined)
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
      call hsd_set(node, name, val)
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
      call hsd_set(node, name, val)
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
        do ir = 1, size(val, 1)
          if (ir > 1) raw = raw // new_line('a')
          do ic = 1, size(val, 2)
            write(buffer, '(I0)') val(ir, ic)
            if (ic > 1) raw = raw // " "
            raw = raw // trim(adjustl(buffer))
          end do
        end do
        call new_value(txt, name="#text")
        call txt%set_raw(raw)
        call child%add_child(txt)
      end block
    else
      call hsd_set(node, name, val)
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
        do ir = 1, size(val, 1)
          if (ir > 1) raw = raw // new_line('a')
          do ic = 1, size(val, 2)
            write(buffer, '(G0)') val(ir, ic)
            if (ic > 1) raw = raw // " "
            raw = raw // trim(adjustl(buffer))
          end do
        end do
        call new_value(txt, name="#text")
        call txt%set_raw(raw)
        call child%add_child(txt)
      end block
    else
      call hsd_set(node, name, val)
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
    else
      call hsd_set(node, name, trim(adjustl(strBuffer)))
    end if
    if (present(modifier)) call hsd_set_attrib(node, name, modifier)

  end subroutine setChVal_intR2RealR2


  ! ============================================================
  !  setChild — create a new child block
  ! ============================================================

  !> Create a new named child block under a parent.
  subroutine setChild(node, name, child, replace, list, modifier)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(hsd_table), pointer, intent(out) :: child
    logical, intent(in), optional :: replace
    logical, intent(in), optional :: list
    character(len=*), intent(in), optional :: modifier

    type(hsd_table) :: newTable
    class(hsd_node), pointer :: storedChild

    call new_table(newTable, name=name)
    call node%add_child(newTable)

    call node%get_child(node%num_children, storedChild)
    select type (t => storedChild)
    type is (hsd_table)
      child => t
    class default
      child => null()
    end select

    if (present(modifier)) then
      if (len_trim(modifier) > 0) then
        call hsd_set_attrib(node, name, modifier)
      end if
    end if

  end subroutine setChild


  ! ============================================================
  !  getDescendant — path-based child lookup
  ! ============================================================

  !> Navigate to a descendant node along a "/" separated path.
  subroutine getDescendant(root, path, child, requested, processed, parent)
    type(hsd_table), intent(inout), target :: root
    character(len=*), intent(in) :: path
    type(hsd_table), pointer, intent(out) :: child
    logical, intent(in), optional :: requested
    logical, intent(in), optional :: processed
    type(hsd_table), pointer, intent(out), optional :: parent

    integer :: iSlash
    logical :: isRequired
    type(hsd_table), pointer :: cur
    character(len=:), allocatable :: remaining, segment

    isRequired = .false.
    if (present(requested)) isRequired = requested

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

      call getChild(cur, segment, child, requested=.false.)
      if (.not. associated(child)) then
        if (isRequired) then
          call dftbp_error(root, "Required path not found: '" // path // "'")
        end if
        return
      end if
      cur => child
    end do

  end subroutine getDescendant


  ! ============================================================
  !  setNodeName — rename a node directly
  ! ============================================================

  !> Rename a node.
  subroutine setNodeName(node, name, updateHsdName, parent)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: updateHsdName
    type(hsd_table), intent(inout), optional :: parent

    node%name = tolower(name)
    if (present(parent)) call parent%invalidate_index()

  end subroutine setNodeName


  ! ============================================================
  !  renameDescendant — rename a child by path
  ! ============================================================

  !> Rename a descendant node by path, handling both table and value nodes.
  subroutine renameDescendant(root, path, newName, warningMsg)
    type(hsd_table), intent(inout), target :: root
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: newName
    character(len=*), intent(in), optional :: warningMsg

    integer :: stat

    call hsd_rename_child(root, path, newName, stat, case_insensitive=.true.)
    if (stat == HSD_STAT_OK) then
      if (present(warningMsg)) then
        call dftbp_warning(root, warningMsg)
      end if
    end if

  end subroutine renameDescendant


  ! ============================================================
  !  removeChild — remove a specific child from its parent
  ! ============================================================

  !> Removes a specific child node from its parent.
  function removeChild(parent, child) result(removed)
    type(hsd_table), pointer, intent(inout) :: parent
    type(hsd_table), pointer, intent(inout) :: child
    type(hsd_table), pointer :: removed

    call hsd_remove_child(parent, child%name)
    child => null()
    removed => null()

  end function removeChild


end module dftbp_io_hsdutils
