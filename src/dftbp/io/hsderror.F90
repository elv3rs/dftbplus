!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Error and warning reporting bridge for the hsd-fortran/hsd-data migration.
!>
!> Replaces the legacy detailedError/detailedWarning (from hsdutils) which extracted path/line
!> information from xmlf90 DOM attributes. This module constructs error messages using hsd_node
!> metadata (name, line) and delegates to dftbp_io_message for output.
module dftbp_io_hsderror
  use dftbp_extlibs_hsddata, only : hsd_node, hsd_table, hsd_value, hsd_error_t
  use dftbp_io_message, only : error, warning
  implicit none

  private
  public :: hsd_fatal_error, hsd_warning, hsd_error_from_stat

  !> Newline character
  character(len=*), parameter :: nl = new_line('a')

contains


  !> Report a fatal error associated with an HSD node and abort.
  !>
  !> Constructs an error message including the node name and line number,
  !> then delegates to the DFTB+ error() routine which terminates the program.
  !>
  !> @param node  The HSD node where the error occurred.
  !> @param msg   Human-readable error description.
  !> @param path  Optional HSD path context (e.g., "Hamiltonian/DFTB").
  subroutine hsd_fatal_error(node, msg, path)
    class(hsd_node), intent(in) :: node
    character(len=*), intent(in) :: msg
    character(len=*), intent(in), optional :: path

    call error(format_message_(node, msg, path))

  end subroutine hsd_fatal_error


  !> Report a warning associated with an HSD node.
  !>
  !> Constructs a warning message including the node name and line number,
  !> then delegates to the DFTB+ warning() routine.
  !>
  !> @param node  The HSD node where the warning occurred.
  !> @param msg   Human-readable warning description.
  !> @param path  Optional HSD path context (e.g., "Hamiltonian/DFTB").
  subroutine hsd_warning(node, msg, path)
    class(hsd_node), intent(in) :: node
    character(len=*), intent(in) :: msg
    character(len=*), intent(in), optional :: path

    call warning(format_message_(node, msg, path))

  end subroutine hsd_warning


  !> Convert an hsd_error_t to a fatal DFTB+ error.
  !>
  !> If the error is allocated (i.e. an error occurred), report it and abort.
  !> If not allocated, this is a no-op.
  !>
  !> @param err    The hsd_error_t from an hsd-data operation.
  !> @param prefix Optional prefix message for context.
  subroutine hsd_error_from_stat(err, prefix)
    type(hsd_error_t), allocatable, intent(in) :: err
    character(len=*), intent(in), optional :: prefix

    character(len=:), allocatable :: full_msg

    if (.not. allocated(err)) return

    if (present(prefix)) then
      full_msg = trim(prefix) // ": " // trim(err%message)
    else
      full_msg = trim(err%message)
    end if

    call error(full_msg)

  end subroutine hsd_error_from_stat


  ! ---------------------------------------------------------------------------
  !  Private helpers
  ! ---------------------------------------------------------------------------


  !> Format an error/warning message with node location info.
  function format_message_(node, msg, path) result(formatted)
    class(hsd_node), intent(in) :: node
    character(len=*), intent(in) :: msg
    character(len=*), intent(in), optional :: path
    character(len=:), allocatable :: formatted

    character(len=20) :: line_str

    formatted = trim(msg)

    ! Append path context
    if (present(path)) then
      if (len_trim(path) > 0) then
        if (allocated(node%name) .and. len_trim(node%name) > 0) then
          formatted = formatted // nl // "Path: " // trim(path) // "/" // trim(node%name)
        else
          formatted = formatted // nl // "Path: " // trim(path)
        end if
      else if (allocated(node%name) .and. len_trim(node%name) > 0) then
        formatted = formatted // nl // "Path: " // trim(node%name)
      end if
    else if (allocated(node%name) .and. len_trim(node%name) > 0) then
      formatted = formatted // nl // "Path: " // trim(node%name)
    end if

    ! Append line number
    if (node%line > 0) then
      write(line_str, '(I0)') node%line
      formatted = formatted // nl // "Line: " // trim(line_str)
    end if

  end function format_message_

end module dftbp_io_hsderror
