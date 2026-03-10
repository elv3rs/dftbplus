!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Parallel input parsing routines exposed for unit testing.
!> This module re-exports selected routines from dftbp_dftbplus_parser.
module dftbp_dftbplus_parser_parallel
  use dftbp_dftbplus_parser, only : readParallel
  implicit none

  public :: readParallel

end module dftbp_dftbplus_parser_parallel
