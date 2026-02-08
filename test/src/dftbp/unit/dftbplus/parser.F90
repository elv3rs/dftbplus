!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

!> Unit tests for parser sub-module subroutines.
module test_dftbplus_parser
  use fortuno_serial, only : suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : Boltzmann
  use hsd, only : hsd_table, hsd_error_t, hsd_get, hsd_set, HSD_STAT_OK, hsd_has_child, &
      & hsd_get_table
  use hsd_data, only : new_table, data_load_string, DATA_FMT_HSD
  use dftbp_dftbplus_parser_shared_utils, only : readMDInitTemp
  use dftbp_dftbplus_parser_parallel, only : readParallel
  use dftbp_dftbplus_inputdata, only : TInputData, init, destruct
  use dftbp_common_globalenv, only : withScalapack
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains


  $:TEST("readMDInitTemp_kelvin", label="parser")
    !! Read initial temperature in Kelvin
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot
    type(hsd_error_t), allocatable :: err
    real(dp) :: tempAtom

    call data_load_string("InitialTemperature [K] = 300.0", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    pRoot => root
    call readMDInitTemp(pRoot, tempAtom, 0.0_dp)
    ! 300 K in atomic units: 300 * Boltzmann
    @:ASSERT(abs(tempAtom - 300.0_dp * Boltzmann) < 1.0e-12_dp)
  $:END_TEST()


  $:TEST("readMDInitTemp_hartree", label="parser")
    !! Read initial temperature in Hartree (default units)
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot
    type(hsd_error_t), allocatable :: err
    real(dp) :: tempAtom

    call data_load_string("InitialTemperature = 0.001", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    pRoot => root
    call readMDInitTemp(pRoot, tempAtom, 0.0_dp)
    @:ASSERT(abs(tempAtom - 0.001_dp) < 1.0e-12_dp)
  $:END_TEST()


  $:TEST("readMDInitTemp_minimum", label="parser")
    !! Temperature below minimum is clamped
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot
    type(hsd_error_t), allocatable :: err
    real(dp) :: tempAtom

    call data_load_string("InitialTemperature = 0.00001", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    pRoot => root
    call readMDInitTemp(pRoot, tempAtom, 0.001_dp)
    ! Result should be clamped to minimum
    @:ASSERT(abs(tempAtom - 0.001_dp) < 1.0e-12_dp)
  $:END_TEST()


  $:TEST("readParallel_defaults", label="parser")
    !! Read Parallel block with defaults
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot
    type(hsd_error_t), allocatable :: err
    type(TInputData) :: input

    call data_load_string("Parallel {}", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call init(input)
    pRoot => root
    call readParallel(pRoot, input)
    ! With default values: nGroup=1, tOmpThreads depends on MPI build
    @:ASSERT(allocated(input%ctrl%parallelOpts))
    @:ASSERT(input%ctrl%parallelOpts%nGroup == 1)
    ! blockSize only meaningful with scalapack; skip check otherwise
    if (withScalapack) then
      @:ASSERT(input%ctrl%parallelOpts%blacsOpts%blockSize == 32)
    end if
    call destruct(input)
  $:END_TEST()


  $:TEST("readParallel_custom", label="parser")
    !! Read Parallel block with custom values
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot
    type(hsd_error_t), allocatable :: err
    type(TInputData) :: input

    call data_load_string(&
        & "Parallel { Groups = 4" // new_line("a") // &
        & "  UseOmpThreads = Yes" // new_line("a") // &
        & "  Blacs { BlockSize = 64 }" // new_line("a") // &
        & "}", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call init(input)
    pRoot => root
    call readParallel(pRoot, input)
    @:ASSERT(allocated(input%ctrl%parallelOpts))
    @:ASSERT(input%ctrl%parallelOpts%nGroup == 4)
    @:ASSERT(input%ctrl%parallelOpts%tOmpThreads .eqv. .true.)
    @:ASSERT(input%ctrl%parallelOpts%blacsOpts%blockSize == 64)
    call destruct(input)
  $:END_TEST()


  $:TEST("readParallel_absent", label="parser")
    !! No Parallel block → parallelOpts not allocated
    type(hsd_table), target :: root
    type(hsd_table), pointer :: pRoot
    type(hsd_error_t), allocatable :: err
    type(TInputData) :: input

    call data_load_string("SomeOther = 1", root, DATA_FMT_HSD, err)
    @:ASSERT(.not. allocated(err))
    call init(input)
    pRoot => root
    call readParallel(pRoot, input)
    @:ASSERT(.not. allocated(input%ctrl%parallelOpts))
    call destruct(input)
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("parser", test_list([&
            $:TEST_ITEMS(label="parser")
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()
  end function tests


end module test_dftbplus_parser
