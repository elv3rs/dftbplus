!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Unit conversion helpers bridging HSD modifier strings to DFTB+'s TUnit system.
!>
!> Extracted from the legacy compatibility layer as part of the hsd-fortran direct-API migration.
!> This module provides the convertUnitHsd interface that:
!>   1. Reads the modifier string (e.g. "[Angstrom]") from an HSD node's attrib
!>   2. Looks up the unit in the DFTB+ TUnit table
!>   3. Applies the conversion factor to the value(s)
!>
!> It depends on hsd-fortran for error formatting and hsd_set, and on DFTB+'s
!> own unitconversion module for the TUnit type and convertUnit procedure.
module dftbp_io_unitconv
  use dftbp_common_accuracy, only : dp
  use dftbp_common_unitconversion, only : TUnit, convertUnit, statusCodes
  use hsd, only : hsd_table, hsd_format_error, hsd_set
  use dftbp_io_message, only : error
  implicit none

  private

  !> Generic interface for unit conversion with HSD modifier.
  !>
  !> Maps the modifier string (from hsd_get_attrib) to a unit conversion factor
  !> via DFTB+'s TUnit table and applies it.  Reports a formatted error if the
  !> modifier names an unknown unit.
  interface convertUnitHsd
    module procedure :: convertUnitHsd_R0
    module procedure :: convertUnitHsd_R1
    module procedure :: convertUnitHsd_R2
  end interface convertUnitHsd

  public :: convertUnitHsd

  !> Error message prefix for invalid modifiers
  character(len=*), parameter :: MSG_INVALID_MODIFIER = "Unknown modifier '"

contains

  !> Convert a scalar value using the unit modifier (rank 0).
  subroutine convertUnitHsd_R0(modifier, units, child, convertValue, replace, changed)
    character(len=*), intent(in) :: modifier
    type(TUnit), intent(in) :: units(:)
    type(hsd_table), intent(inout) :: child
    real(dp), intent(inout) :: convertValue
    logical, intent(in), optional :: replace
    logical, intent(out), optional :: changed

    logical :: changed_
    integer :: status
    character(len=:), allocatable :: msg

    changed_ = len_trim(modifier) > 0
    if (changed_) then
      call convertUnit(units, modifier, convertValue, status)
      if (status /= statusCodes%ok) then
        call hsd_format_error(child, MSG_INVALID_MODIFIER // modifier // "'", msg)
        call error(msg)
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
    character(len=:), allocatable :: msg

    changed_ = len_trim(modifier) > 0
    if (changed_) then
      call convertUnit(units, modifier, convertValue, status)
      if (status /= statusCodes%ok) then
        call hsd_format_error(child, MSG_INVALID_MODIFIER // modifier // "'", msg)
        call error(msg)
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
    character(len=:), allocatable :: msg

    changed_ = len_trim(modifier) > 0
    if (changed_) then
      call convertUnit(units, modifier, convertValue, status)
      if (status /= statusCodes%ok) then
        call hsd_format_error(child, MSG_INVALID_MODIFIER // modifier // "'", msg)
        call error(msg)
      end if
      if (present(replace)) then
        if (replace) call hsd_set(child, "", convertValue)
      end if
    end if
    if (present(changed)) changed = changed_

  end subroutine convertUnitHsd_R2

end module dftbp_io_unitconv
