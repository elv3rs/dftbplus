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

!>   - Atom/index selection helpers
!>   - Node name and modifier extraction helpers
!>   - Processed-flag management (setProcessed / setUnprocessed)
!>
!> This module is a permanent part of the DFTB+ IO layer.  It has no legacy-API
!> surface — callers use hsd-fortran types natively.
module dftbp_io_hsdutils
  use dftbp_common_accuracy, only : mc
  use dftbp_common_status, only : TStatus
  use hsd, only : hsd_table, &
      & hsd_format_error, hsd_format_warning
  use dftbp_io_charmanip, only : i2c
  use dftbp_io_indexselection, only : getIndexSelection
  use dftbp_io_message, only : error, warning
  implicit none

  private

  ! --- Public: error / warning helpers ---
  public :: dftbp_error, dftbp_warning


  ! --- Public: atom / index selection ---
  public :: getSelectedAtomIndices, getSelectedIndices

  ! --- Public: modifier / name helpers ---
  public :: splitModifier


  ! ============================================================
  !  Constants
  ! ============================================================

  !> Separator character for splitting modifiers
  character, parameter :: sepModifier = ","

  ! ============================================================
  !  Types
  ! ============================================================

contains

  ! ============================================================
  !  Private helpers
  ! ============================================================




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




end module dftbp_io_hsdutils
