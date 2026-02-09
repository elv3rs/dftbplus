!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module dftbp_dftbplus_input_fileaccess
  use dftbp_common_file, only : fileAccessValues
  use dftbp_io_charmanip, only : tolower, unquote
  use hsd, only : hsd_table, hsd_get
  use dftbp_io_hsdutils, only : getChild, setChildValue
  use dftbp_io_hsdutils, only : dftbp_error
  implicit none

  private
  public :: readBinaryAccessTypes

contains


  !> Reads in the file acess types
  subroutine readBinaryAccessTypes(node, accessTypes)

    !> Parent note which should contain the "BinaryAccessTypes" subnode
    type(hsd_table), pointer, intent(in) :: node

    !> Read and write access types on exit (defaulting to ["stream", "stream"])
    character(*), intent(out) :: accessTypes(:)

    type(hsd_table), pointer :: child
    character(:), allocatable :: stringArr(:)
    integer :: nStr, ii

    @:ASSERT(size(accessTypes) == 2)

    call getChild(node, "BinaryAccessTypes", child, requested=.false.)
    if (.not. associated(child)) then
      call setChildValue(node, "BinaryAccessTypes", ["stream"], child=child)
    end if
    call hsd_get(node, "BinaryAccessTypes", stringArr)
    nStr = size(stringArr)
    if (nStr < 1 .or. nStr > 2) then
      call dftbp_error(child, "BinaryAccessTypes needs one or two arguments")
    end if
    accessTypes(1:nStr) = stringArr(1:nStr)
    if (nStr == 1) then
      accessTypes(2) = accessTypes(1)
    end if
    accessTypes(:) = tolower(unquote(accessTypes))
    do ii = 1, size(accessTypes)
      if (.not. any(accessTypes(ii) == fileAccessValues)) then
        call dftbp_error(child, "Invalid file access type '" // trim(accessTypes(ii)) // "'")
      end if
    end do

  end subroutine readBinaryAccessTypes

end module dftbp_dftbplus_input_fileaccess
