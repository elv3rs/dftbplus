!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_libwavegrid_simple
  use dftbp_common_accuracy, only : dp
  use libwavegrid, only : TSlaterOrbital, TSpeciesBasis, TMolecularOrbital, &
                       & getValue, getAtomicDensities, getTotalChrg
  use dftbp_type_typegeometry, only : TGeometry
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_boundarycond, only : TBoundaryConds, boundaryCondsEnum, TBoundaryConds_init
  use fortuno_serial, only : suite => serial_suite_item, test_list, all_close
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains

  subroutine initOrbitalHydrogenS(sto, useRadialLut)
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
  end subroutine initOrbitalHydrogenS


  subroutine initOrbitalOxygenS(sto, useRadialLut)
    type(TSlaterOrbital), intent(out) :: sto
    logical, intent(in), optional :: useRadialLut

    integer, parameter :: sto_ll = 0
    real(dp), parameter :: lut_resolution = 0.01_dp
    real(dp), parameter :: sto_cutoff = 6.0_dp
    real(dp), parameter :: sto_alpha(4) = [0.5_dp, 1.26_dp, 3.17_dp, 8.0_dp]
    real(dp), parameter :: sto_aa(3,4) = reshape([ &
        0.21323488915449521_dp, -0.031152441012403959_dp, 0.0011303933349092960_dp, -9.0596686234211106_dp, &
        4.0675254100925002_dp, -0.60689938592758674_dp, 10.444780695232501_dp, -3.1373217750406419_dp, &
        4.4644627824691749_dp, 8.9215282105208260_dp, 7.2210396633361826_dp, 16.146571373535430_dp], [3,4])
    call sto%init(sto_aa, sto_alpha, sto_ll, lut_resolution, sto_cutoff, useRadialLut=useRadialLut)
  end subroutine initOrbitalOxygenS


  subroutine initOrbitalOxygenP(sto, useRadialLut)
    type(TSlaterOrbital), intent(out) :: sto
    logical, intent(in), optional :: useRadialLut

    integer, parameter :: sto_ll = 1
    real(dp), parameter :: lut_resolution = 0.01_dp
    real(dp), parameter :: sto_cutoff = 6.0_dp
    real(dp), parameter :: sto_alpha(4) = [0.5_dp, 1.26_dp, 3.17_dp, 8.0_dp]
    real(dp), parameter :: sto_aa(3,4) = reshape([ &
        -0.021351405651207991_dp, 0.0028544859270132768_dp, -9.4141846289124166E-005_dp, 1.8517392789336220_dp, &
        -0.79114942586812875_dp, 0.10094277989615121_dp, 16.210706533770320_dp, -10.077615056451849_dp, &
        7.7615980276616314_dp, -1.7017045797631820_dp, -10.773616241206961_dp, -35.439076485248712_dp], [3,4])
    call sto%init(sto_aa, sto_alpha, sto_ll, lut_resolution, sto_cutoff, useRadialLut=useRadialLut)
  end subroutine initOrbitalOxygenP

  subroutine initSpeciesBasisH(speciesBasis, useRadialLut)
    type(TSpeciesBasis), intent(out) :: speciesBasis(1)
    logical, intent(in), optional :: useRadialLut

    ! Hydrogen
    speciesBasis(1)%atomicNumber = 1
    speciesBasis(1)%nOrb = 1
    allocate(speciesBasis(1)%stos(1))
    call initOrbitalHydrogenS(speciesBasis(1)%stos(1), useRadialLut=useRadialLut)

  end subroutine initSpeciesBasisH

  subroutine initSpeciesBasisH2O(speciesBasis, useRadialLut)
    type(TSpeciesBasis), intent(out) :: speciesBasis(2)
    logical, intent(in), optional :: useRadialLut
  
    ! Oxygen
    speciesBasis(1)%atomicNumber = 8
    speciesBasis(1)%nOrb = 2
    allocate(speciesBasis(1)%stos(2))
    call initOrbitalOxygenS(speciesBasis(1)%stos(1), useRadialLut=useRadialLut)
    call initOrbitalOxygenP(speciesBasis(1)%stos(2), useRadialLut=useRadialLut)
    
    ! Hydrogen 
    speciesBasis(2)%atomicNumber = 1
    speciesBasis(2)%nOrb = 1
    allocate(speciesBasis(2)%stos(1))
    call initOrbitalHydrogenS(speciesBasis(2)%stos(1), useRadialLut=useRadialLut)

  end subroutine initSpeciesBasisH2O


  subroutine initGeometryH2O(geometry)
    type(TGeometry), intent(out) :: geometry
    geometry%nAtom = 3
    geometry%tPeriodic = .false.
    geometry%tFracCoord = .false.
    allocate(geometry%species(3))
    geometry%species = [1, 2, 2]
    allocate(geometry%coords(3,3))
    geometry%coords(:,1) = [0.0_dp, 0.0_dp, 0.0_dp]
    geometry%coords(:,2) = [0.0_dp, 0.0_dp, 1.4797763915205659_dp]
    geometry%coords(:,3) = [0.0_dp, -1.8897259885789233_dp, -1.4797763915205659_dp]
    geometry%nSpecies = 2
    geometry%speciesNames = ["O", "H"]
    geometry%tHelical = .false.
  end subroutine initGeometryH2O


  subroutine initMolorbH2O(molorb, useRadialLut)
    type(TMolecularOrbital), intent(out) :: molorb
    logical, intent(in), optional :: useRadialLut
    
    real(dp), parameter :: gridOrigin(3) = [-6.0_dp, -7.8897259885789230_dp, -7.4797763915205664_dp]
    real(dp), parameter :: gridVecs(3,3) = reshape([ &
        0.12_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.13889725988578924_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.14959552783041133_dp], [3,3])
    type(TSpeciesBasis) :: speciesBasis(2)
    type(TGeometry) :: geometry
    type(TBoundaryConds) :: bconds
    type(TStatus) :: status

    call initSpeciesBasisH2O(speciesBasis, useRadialLut=useRadialLut)
    call initGeometryH2O(geometry)
    call TBoundaryConds_init(bconds, boundaryCondsEnum%cluster, errStatus=status)
    @:ASSERT(status%isOk())

    call molorb%init(geometry, bconds, speciesBasis, gridOrigin, gridVecs)
  end subroutine initMolorbH2O






  $:TEST("LibwavegridH2O")
    type(TMolecularOrbital) :: molorb
    logical, parameter :: useRadialLut = .false.
  
    call initMolorbH2O(molorb, useRadialLut=useRadialLut)

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
