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
!>   - getChild — child table lookup with error reporting and modifier extraction
!>   - getDescendant / renameDescendant / setNodeName / removeChild — tree manipulation
!>   - Atom/index selection helpers
!>   - Node name and modifier extraction helpers
!>   - Processed-flag management (setProcessed / setUnprocessed)
!>
!> This module is a permanent part of the DFTB+ IO layer.  It has no legacy-API
!> surface — callers use hsd-fortran types natively.
module dftbp_io_hsdutils
  use dftbp_common_accuracy, only : mc
  use dftbp_common_status, only : TStatus
  use hsd, only : hsd_table, hsd_value, hsd_node, hsd_iterator, &
      & hsd_format_error, hsd_format_warning, hsd_get, &
      & hsd_get_child, hsd_get_table, hsd_get_attrib, &
      & hsd_remove_child, hsd_rename_child, &
      & HSD_STAT_OK
  use hsd_data, only : new_table
  use dftbp_io_charmanip, only : i2c, tolower
  use dftbp_io_indexselection, only : getIndexSelection
  use dftbp_io_message, only : error, warning
  use dftbp_io_tokenreader, only : TOKEN_EOS, TOKEN_ERROR
  implicit none

  private

  ! --- Public: error / warning helpers ---
  public :: dftbp_error, dftbp_warning

  ! --- Public: getChild ---
  public :: getChild

  ! --- Public: tree manipulation ---
  public :: getDescendant, renameDescendant, setNodeName, removeChild

  ! --- Public: atom / index selection ---
  public :: getSelectedAtomIndices, getSelectedIndices

  ! --- Public: modifier / name helpers ---
  public :: splitModifier, getModifier
  public :: getNodeName, getNodeName2, getNodeHSDName

  ! --- Public: processed flags ---
  public :: setProcessed

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

contains

  ! ============================================================
  !  Private helpers
  ! ============================================================


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
