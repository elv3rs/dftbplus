!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_waveplot_cache
  use dftbp_common_accuracy, only : dp
  use waveplot_molorb2, only: gridBoxFromSphereRadius
  use fortuno_serial, only : suite => serial_suite_item, test_list, all_close
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains


  $:TEST("gridBoxFromSphereRadiusOrthogonal")
    real(dp) :: gridVecs(3,3)
    integer :: nPointsHalved(3)
    real(dp) :: cutoff

    ! Friendly cartesian orthonormal (identity) basis 
    gridVecs = 0.0_dp
    gridVecs(1,1) = 1.0_dp
    gridVecs(2,2) = 1.0_dp
    gridVecs(3,3) = 1.0_dp

    ! Simplest case
    cutoff = 5.0_dp
    call gridBoxFromSphereRadius(gridVecs, cutoff, nPointsHalved)
    @:ASSERT(all(nPointsHalved == [5, 5, 5]))

    ! Ensure we dont truncate the wavefunction
    cutoff = 5.1_dp
    call gridBoxFromSphereRadius(gridVecs, cutoff, nPointsHalved)
    @:ASSERT(all(nPointsHalved == [6, 6, 6]))

  $:END_TEST()

  $:TEST("gridBoxFromSphereRadiusSkewed")
    ! Test with non-orthogonal, non-unit vectors
    real(dp) :: gridVecs(3,3)
    integer :: nPointsHalved(3)
    real(dp) :: cutoff

    ! Basis: = [[1, 1, 0], [0, 1, 0], [0, 0, 1]]
    gridVecs = reshape([ 1.0_dp, 0.0_dp, 0.0_dp, &
                         1.0_dp, 1.0_dp, 0.0_dp, & 
                         0.0_dp, 0.0_dp, 1.0_dp], [3,3])
    cutoff = 10.0_dp
    ! Simple verification:
    ! Consider circle in 2d space:
    ! |r| = sqrt(x^2 + y^2)
    ! x = a + b
    ! y = b
    ! |r| = sqrt(a^2 + 2ab + 2b^2)
    ! r^2 = a^2 + 2ab + 2b^2
    ! Set r=10
    ! 100 = a^2 + 2ab + 2b^2
    ! Solve quadratic in terms of a:
    ! a^2 + 2ab  + b^2 = 100 - b^2
    !          (a+b)^2 = 100-b^2
    !                a = sqrt(100-b^2)-b
    ! Find max b:
    !  => 0 = d/db a 
    ! <=> 0 = -b (100-b^2)^{-1/2} -1
    ! <=> -1 = b (100-b^2)^{-1/2}
    ! <=> 100-b^2 = b^2
    ! <=> b = sqrt(50) = 5 * sqrt(2)
    ! Substituting into a:
    ! a = sqrt(100 - 50) - sqrt(50)
    !   = 10 sqrt(2) \approx 14.4


    call gridBoxFromSphereRadius(gridVecs, cutoff, nPointsHalved)
    @:ASSERT(all(nPointsHalved == [15, 10, 10]))

    ! Ensure we dont truncate the wavefunction
    cutoff = 3.0_dp
    call gridBoxFromSphereRadius(gridVecs, cutoff, nPointsHalved)
    @:ASSERT(all(nPointsHalved == [5, 3, 3]))
  $:END_TEST()


  $:TEST("TOrbitalCacheAlignOrthogonal")
    use waveplot_molorb2, only: TOrbitalCache
    type(TOrbitalCache) :: cache
    real(dp) :: gridVecs(3,3)
    integer :: gridDims(4)
    real(dp) :: shiftedPos(3)
    integer :: iOffset(3), iChunk(3), iMain(3,2)

    ! Friendly cartesian orthonormal (identity) basis
    gridVecs = 0.0_dp
    gridVecs(1,1) = 1.0_dp
    gridVecs(2,2) = 1.0_dp
    gridVecs(3,3) = 1.0_dp

    ! Main grid dimensions
    gridDims = [20, 20, 20, 1]

    ! Setup cache instance (supply parameters for align)
    cache%gridVecs = gridVecs
    cache%subdivisionFactor = 10
    cache%nPointsHalved = [5, 5, 5]


    !---- Case 1
    shiftedPos = [0.0_dp, 0.0_dp, 0.0_dp]
    call cache%align(gridDims, shiftedPos, iOffset, iChunk, iMain)
    @:ASSERT(all(iOffset == [-1, -1, -1]))
    @:ASSERT(all(iChunk == [0, 0, 0]))
    @:ASSERT(all(iMain(:,1) == [1, 1, 1]))
    @:ASSERT(all(iMain(:,2) == [6, 6, 6]))
    !---- Case 2
    shiftedPos = [2.5_dp, 3.0_dp, 4.75_dp]
    call cache%align(gridDims, shiftedPos, iOffset, iChunk, iMain)
    @:ASSERT(all(iOffset == [-3, -4, -5])) 
    @:ASSERT(all(iChunk == [5, 0, 7]))
    @:ASSERT(all(iMain(:,1) == [1, 1, 1]))
    @:ASSERT(all(iMain(:,2) == [8, 9, 10]))
    !---- Case 3
    shiftedPos = [18.0_dp, 19.5_dp, 0.25_dp]
    call cache%align(gridDims, shiftedPos, iOffset, iChunk, iMain)
    @:ASSERT(.true.)
    @:ASSERT(all(iOffset == [-19, -20, -1])) 
    @:ASSERT(all(iChunk == [0, 5, 2]))
    @:ASSERT(all(iMain(:,1) == [14, 15, 1]))
    @:ASSERT(all(iMain(:,2) == [20, 20, 6]))

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
