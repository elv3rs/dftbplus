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
  use dftbp_common_globalenv, only : withMpi, withScalapack
  use dftbp_dftbplus_inputdata, only : TBlacsOpts, TInputData
  use dftbp_io_hsdcompat, only : hsd_table, detailedWarning, getChild, getChildValue
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

    call getChild(root, "Parallel", child=node, requested=.false., emptyIfMissing=withMpi)
    if (associated(node)) then
      if (.not. withMpi) then
        call detailedWarning(node, "Settings will be read but ignored (compiled without MPI&
            & support)")
      end if
      allocate(input%ctrl%parallelOpts)
      call getChildValue(node, "Groups", input%ctrl%parallelOpts%nGroup, 1, child=pTmp)
      call getChildValue(node, "UseOmpThreads", input%ctrl%parallelOpts%tOmpThreads, .not. withMpi)
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

    call getChild(root, "Blacs", child=node, requested=.false., emptyIfMissing=withScalapack)
    if (associated(node)) then
      if (.not. withScalapack) then
        call detailedWarning(node, "Settings will be read but ignored (compiled without SCALAPACK&
            & support)")
      end if
      call getChildValue(node, "BlockSize", blacsOpts%blockSize, 32)
    end if

  end subroutine readBlacs


end module dftbp_dftbplus_parser_parallel
