!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Shared parser utility routines exposed for unit testing.
!> This module re-exports selected routines from dftbp_dftbplus_parser.
module dftbp_dftbplus_parser_shared_utils
  use dftbp_dftbplus_parser, only : readMDInitTemp
  implicit none

  public :: readMDInitTemp

end module dftbp_dftbplus_parser_shared_utils
