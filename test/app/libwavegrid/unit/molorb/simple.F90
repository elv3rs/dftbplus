!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_libwavegrid_simple
  use dftbp_common_accuracy, only : dp
  use libwavegrid, only : TSlaterOrbital
  use fortuno_serial, only : suite => serial_suite_item, test_list, all_close
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains

  subroutine initHydrogenSto(sto, useRadialLut)
    type(TSlaterOrbital), intent(out) :: sto
    logical, intent(in), optional :: useRadialLut

    integer, parameter :: sto_ll = 0
    real(dp), parameter :: lut_resolution = 0.01_dp
    real(dp), parameter :: sto_cutoff = 6.0_dp
    real(dp), parameter :: sto_alpha(3) = [0.5_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: sto_aa(3,3) = reshape([ &
        -2.2765228685565400_dp, 17.453716731738609_dp, -12.701455953438341_dp, &
         0.26641083435260687_dp, -5.4229751699602433_dp, -6.5568796727250120_dp, &
        -7.9427553748566294E-003_dp, 0.96370929548055750_dp, -0.85307020704514269_dp], [3,3])

    call sto%init(sto_aa, sto_alpha, sto_ll, lut_resolution, sto_cutoff, useRadialLut=useRadialLut)
  end subroutine initHydrogenSto

  subroutine initOxygenSto(sto, useRadialLut)
    type(TSlaterOrbital), intent(out) :: sto
    logical, intent(in), optional :: useRadialLut

    integer, parameter :: sto_ll = 0
    real(dp), parameter :: lut_resolution = 0.01_dp
    real(dp), parameter :: sto_cutoff = 6.0_dp
    real(dp), parameter :: sto_alpha(3) = [2.0_dp, 4.0_dp, 6.0_dp]
    real(dp), parameter :: sto_aa(3,3) = reshape([ &
        -7.8382725000000000E-002_dp, 1.6669999999999999E+000_dp, -1.3059999999999999E+000_dp, &
         1.1930000000000000E-002_dp, -2.7699999999999999E-001_dp, -3.3499999999999999E-001_dp, &
        -3.7870000000000001E-004_dp, 4.8999999999999998E-002_dp, -4.3299999999999998E-002_dp], [3,3])

    call sto%init(sto_aa, sto_alpha, sto_ll, lut_resolution, sto_cutoff, useRadialLut=useRadialLut)
  end subroutine initOxygenSto



  $:TEST("TOrbitalCacheInitialiseSimple")
  
    @:ASSERT(.true.)
  ! Cases to handle:
  ! Both cpu, gpu.
  ! real : default, totChrg, addDensities, lut
  ! cmplx: default, totChrg, lut

    ! TODO: Check that orbital matches expectations#
    ! if WITH_CUDA, duplicate all and set useGPU = .true.
    ! H-chain, H2O
    ! using macros: for every <name> in list
    ! New testcase_<name>
    ! allocate buffer
    ! load input args from file
    ! call molorb
    ! load expected output from file
    ! compare all_close(actual, expected)
    ! compare sum(actual), sum(expected)

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

end module test_libwavegrid_simple
