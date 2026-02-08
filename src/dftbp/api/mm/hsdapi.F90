!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides access to HSD manipulation functions
module dftbp_hsdapi
  use hsd, only : hsd_table, hsd_dump, hsd_dump_to_string, hsd_error_t
  implicit none
  private

  public :: hsd_table
  public :: hsd_dump

  !> Backward-compatible dumpHsd interface (file or unit overloads)
  interface dumpHsd
    module procedure :: dumpHsd_file
    module procedure :: dumpHsd_unit
  end interface dumpHsd

  public :: dumpHsd

contains

  !> Dump HSD tree to a named file
  subroutine dumpHsd_file(table, fileName)
    type(hsd_table), intent(in) :: table
    character(len=*), intent(in) :: fileName

    type(hsd_error_t), allocatable :: err

    call hsd_dump(table, fileName, err)

  end subroutine dumpHsd_file

  !> Dump HSD tree to an open file unit
  subroutine dumpHsd_unit(table, unit)
    type(hsd_table), intent(in) :: table
    integer, intent(in) :: unit

    character(len=:), allocatable :: str

    call hsd_dump_to_string(table, str)
    write(unit, '(A)') str

  end subroutine dumpHsd_unit

end module dftbp_hsdapi
