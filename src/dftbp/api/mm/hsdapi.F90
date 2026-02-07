!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides access to HSD manipulation functions
module dftbp_hsdapi
  use dftbp_io_hsdcompat, only : hsd_table, hsd_child_list, &
      & dumpHsd, getChild, getChildren, getChildValue, setChild, setChildValue
  implicit none
  private

  public :: hsd_table, hsd_child_list
  public :: getChild, getChildren, setChild, getChildValue, setChildValue
  public :: dumpHsd

end module dftbp_hsdapi
