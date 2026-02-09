!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Reads the parallel/BLACS block from the HSD input.
module dftbp_dftbplus_parser_parallel
  use hsd_data, only : hsd_table, new_table
  use dftbp_common_globalenv, only : withMpi, withScalapack
  use dftbp_dftbplus_inputdata, only : TBlacsOpts, TInputData
  use hsd, only : hsd_get_or_set, hsd_get_table
  use dftbp_io_hsdutils, only : dftbp_warning
  implicit none

  private
  public :: readParallel

contains


  !> Reads the parallel block.
  subroutine readParallel(root, input)

    !> Root node eventually containing the current block
    type(hsd_table), pointer, intent(in) :: root

    !> Input structure to be filled
    type(TInputData), intent(inout) :: input

    type(hsd_table), pointer :: node, pTmp
    integer :: stat

    call hsd_get_table(root, "Parallel", node, stat, auto_wrap=.true.)
    if (.not. associated(node) .and. withMpi) then
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="parallel")
        call root%add_child(tmpTbl)
      end block
      call hsd_get_table(root, "Parallel", node, stat)
    end if
    if (associated(node)) then
      if (.not. withMpi) then
        call dftbp_warning(node, "Settings will be read but ignored (compiled without MPI&
            & support)")
      end if
      allocate(input%ctrl%parallelOpts)
      call hsd_get_or_set(node, "Groups", input%ctrl%parallelOpts%nGroup, 1, child=pTmp)
      call hsd_get_or_set(node, "UseOmpThreads", input%ctrl%parallelOpts%tOmpThreads, .not. withMpi)
      call readBlacs(node, input%ctrl%parallelOpts%blacsOpts)
    end if

  end subroutine readParallel


  !> Reads the blacs block
  subroutine readBlacs(root, blacsOpts)

    !> Root node eventually containing the current block
    type(hsd_table), pointer, intent(in) :: root

    !> Blacs settings
    type(TBlacsOpts), intent(inout) :: blacsOpts

    type(hsd_table), pointer :: node
    integer :: stat

    call hsd_get_table(root, "Blacs", node, stat, auto_wrap=.true.)
    if (.not. associated(node) .and. withScalapack) then
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="blacs")
        call root%add_child(tmpTbl)
      end block
      call hsd_get_table(root, "Blacs", node, stat)
    end if
    if (associated(node)) then
      if (.not. withScalapack) then
        call dftbp_warning(node, "Settings will be read but ignored (compiled without SCALAPACK&
            & support)")
      end if
      call hsd_get_or_set(node, "BlockSize", blacsOpts%blockSize, 32)
    end if

  end subroutine readBlacs


end module dftbp_dftbplus_parser_parallel
