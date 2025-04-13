!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_waveplot_cache
  use dftbp_common_accuracy, only : dp
  use waveplot_molorb2, only: findCacheSize
  use fortuno_serial, only : suite => serial_suite_item, test_list, all_close
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains
    
  $:TEST("waveplot_findCacheSize")
    !> The grids basis vectors
    real(dp) :: gridVecs(3,3) 

    !> Size of cache from center to cutoff
    integer :: nPointsHalved(3) 
    
    ! Max radial size
    real(dp) :: cutoff 

    call findCacheSize(gridVecs, cutoff, nPointsHalved)
    



  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("wpcache", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_waveplot_cache
