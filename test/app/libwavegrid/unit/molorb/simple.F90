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
  use dftbp_math_simplealgebra, only : determinant33
  use dftbp_dftb_boundarycond, only : TBoundaryConds, boundaryCondsEnum, TBoundaryConds_init
  use fortuno_serial, only : suite => serial_suite_item, test_list, all_close, is_close
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

  type :: spotCheck
    integer :: idx(4)
    real(dp) :: expected
  end type spotCheck

  type :: spotCheckCmpl
    integer :: idx(4)
    complex(dp) :: expected
  end type spotCheckCmpl

  !> Allow 0.01% relative error (LUT interpolation etc.)
  real(dp), parameter :: rtol = 1.0e-4_dp
  
  !> Real Eigenvectors for H2O molecule
  real(dp), parameter :: eigVecsReal(6,4) = reshape([ &
       0.93075335285360816_dp,       5.4749851762277632E-003_dp,   6.5189438567696949E-019_dp,  0.0000000000000000_dp,       &
      -0.11131233240146322_dp,      -0.11131233240146322_dp,      -0.27821574567737223_dp,     -0.73551457100773243_dp,      &
      -7.8726661848432967E-018_dp,   0.0000000000000000_dp,      -0.33117590117481399_dp,      -0.33117590117481405_dp,      &
       2.7154435594648282E-018_dp,  -9.9931138111063277E-018_dp,  0.73160204739108692_dp,       1.1102230246251565E-016_dp,  &
       0.39813631802212995_dp,      -0.39813631802212995_dp,      3.7150709608707145E-033_dp,  -8.1237292785459783E-033_dp,  &
       9.7504438592516575E-017_dp,   1.0000000000000000_dp,      -6.6858833593547624E-017_dp,   6.6858833593547636E-017_dp], &
       & [6,4])
  !> Atomic density occupations for H2O per Orb(distributed evenly over same angMom)
  real(dp), parameter :: occupationAtomDens(6,1) = reshape([2.0_dp, 4.0_dp / 3.0_dp, 4.0_dp / 3.0_dp, 4.0_dp / 3.0_dp, 1.0_dp,&
  1.0_dp], [6,1])
  !> Occupation per state as computed by DFTB+
  real(dp), parameter :: occupationVecH2O(4) = [2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp]

  !> Complex eigenvectors for H chain
  complex(dp), parameter :: eigVecsCmpl(1,4) = reshape([ &
       (0.62981679983963734_dp, 0.0_dp), (0.68885855972586518_dp, 0.0_dp), &
       (0.81991408755170758_dp, 0.0_dp), (1.0501565770378785_dp, 0.0_dp) &
       ], [1,4])
  !> Occupation per state as computed by DFTB+
  real(dp), parameter :: occupationVecHchain(4) = [0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp]
  !> k-points for H chain
  real(dp), parameter :: kPointsHchain(3,8) = reshape([ &
       0.0_dp, 0.0_dp, 0.19634954084936207_dp, 0.0_dp, 0.0_dp, 0.58904862254808621_dp, 0.0_dp, 0.0_dp, &
       0.98174770424681035_dp, 0.0_dp, 0.0_dp, 1.3744467859455345_dp, 0.0_dp, 0.0_dp, 1.7671458676442586_dp, &
       0.0_dp, 0.0_dp, 2.1598449493429825_dp, 0.0_dp, 0.0_dp, 2.5525440310417071_dp, 0.0_dp, 0.0_dp, &
       2.9452431127404308_dp], [3,8])
  integer, parameter :: kIndexesHchain(4) = [1, 2, 3, 4]


contains

  subroutine initOrbitalHydrogenS(sto, useRadialLut)
    type(TSlaterOrbital), intent(out) :: sto
    logical, intent(in), optional :: useRadialLut

    integer, parameter :: sto_ll = 0
    real(dp), parameter :: lut_resolution = 0.01_dp
    real(dp), parameter :: sto_cutoff = 6.0_dp
    real(dp), parameter :: sto_alpha(3) = [0.5_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: sto_aa(3,3) = reshape([ &
        -2.2765228685565400_dp, 0.26641083435260687_dp, -7.9427553748566294E-003_dp, &
         17.453716731738609_dp, -5.4229751699602433_dp, 0.96370929548055750_dp, &
        -12.701455953438341_dp, -6.5568796727250120_dp, -0.85307020704514269_dp], [3,3])
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
    geometry%coords(:,1) = [0.0_dp, -1.8897259885789233_dp, 0.0_dp]
    geometry%coords(:,2) = [0.0_dp, 0.0_dp, 1.4797763915205659_dp]
    geometry%coords(:,3) = [0.0_dp, 0.0_dp, -1.4797763915205659_dp]

    geometry%nSpecies = 2
    allocate(geometry%speciesNames(2))
    geometry%speciesNames = ["O", "H"]
    geometry%tHelical = .false.
  end subroutine initGeometryH2O


  subroutine initGeometryHchain(geometry)
    type(TGeometry), intent(out) :: geometry
    geometry%nAtom = 1
    geometry%tPeriodic = .true.
    geometry%tFracCoord = .false.
    allocate(geometry%species(1))
    geometry%species = [1]
    allocate(geometry%coords(3,1))
    geometry%coords(:,1) = 0.0_dp
    geometry%nSpecies = 1
    allocate(geometry%origin(3))
    geometry%origin = 0.0_dp
    allocate(geometry%latVecs(3,3))
    geometry%latVecs = 0.0_dp
    geometry%latVecs(1,1) = 188.97259885789234_dp
    geometry%latVecs(2,2) = 188.97259885789234_dp
    geometry%latVecs(3,3) = 1.5117807908631387_dp
    allocate(geometry%recVecs2p(3,3))
    geometry%recVecs2p = 0.0_dp
    geometry%recVecs2p(1,1) = 5.2917724899999999E-003_dp
    geometry%recVecs2p(2,2) = 5.2917724899999999E-003_dp
    geometry%recVecs2p(3,3) = 0.66147156125000006_dp
    allocate(geometry%speciesNames(1))
    geometry%speciesNames = ["H"]
    geometry%tHelical = .false.
  end subroutine initGeometryHchain


  
  subroutine initMolorbHchain(molorb, useRadialLut)
    type(TMolecularOrbital), intent(out) :: molorb
    logical, intent(in), optional :: useRadialLut
    
    real(dp), parameter :: gridOrigin(3) = [-5.0_dp, -5.0_dp, -5.0_dp]
    real(dp), parameter :: gridVecs(3,3) = reshape([ &
        0.1_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.1_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.1_dp], [3,3])
    type(TSpeciesBasis) :: speciesBasis(1)
    type(TGeometry) :: geometry
    type(TBoundaryConds) :: bconds
    type(TStatus) :: status

    call initSpeciesBasisH(speciesBasis, useRadialLut=useRadialLut)
    call initGeometryHchain(geometry)
    call TBoundaryConds_init(bconds, boundaryCondsEnum%pbc3d, errStatus=status)
    @:ASSERT(status%isOk())

    call molorb%init(geometry, bconds, speciesBasis, gridOrigin, gridVecs)
  end subroutine initMolorbHchain



  subroutine initMolorbH2O(molorb, useRadialLut)
    type(TMolecularOrbital), intent(out) :: molorb
    logical, intent(in), optional :: useRadialLut
    
    real(dp), parameter :: gridOrigin(3) = [-5.0_dp, -5.0_dp, -5.0_dp]
    real(dp), parameter :: gridVecs(3,3) = reshape([ &
        0.1_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.1_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.1_dp], [3,3])
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

  ! ### Real (H2O) Test cases ###

  ! -- Real (H2O) : Total charge calculation --
  $:TEST("molorb_real_totChrg")
    type(TMolecularOrbital) :: molorb
    logical, parameter :: useRadialLut = .true.
    logical, parameter :: useGPU = .false.
    real(dp) :: valueOnGrid(100,100,100,1)
    real(dp) :: gridVol, actual, expected
    integer :: i, idx(4)
    type(spotCheck) :: spotChecks(6) = [ &
        spotCheck([50, 31, 50, 1], 1.009617874235818_dp), & ! O at (0,-1.89,0)
        spotCheck([1, 1, 1, 1], 0.0_dp), & ! Grid bounds
        spotCheck([100, 100, 100, 1], 0.0_dp), & ! Grid bounds
        spotCheck([51, 51, 51, 1], 0.6092945991929660E-01_dp), & ! Grid Center
        spotCheck([50, 40, 51, 1], 0.5131651070097035_dp), & ! Symmetry: molecule in YZ plane
        spotCheck([52, 40, 51, 1], 0.5131651070097035_dp) & ! Symmetric: molecule in YZ plane
    ]

    call initMolorbH2O(molorb, useRadialLut)
    gridVol = abs(determinant33(molorb%system%gridVecs))

    ! -- Real (H2O) : Total charge calculation --
    call getTotalChrg(molorb, eigVecsReal, valueOnGrid, occupationVecH2O, useGPU)
    ! Check sum over grid
    expected = 8.0040445629655839_dp
    actual = sum(valueOnGrid) * gridVol
    print *, "Total charge:", actual
    @:CHECK(is_close(actual, expected, rtol=rtol))
    ! Spot check values
    do i = 1, size(spotChecks)
      idx = spotChecks(i)%idx
      expected = spotChecks(i)%expected
      actual = valueOnGrid(idx(1), idx(2), idx(3), idx(4))
      @:CHECK(is_close(actual, expected, rtol=rtol))
    end do
  $:END_TEST()

  ! -- Real (H2O) : Atomic densities calculation --
  $:TEST("molorb_real_atomicDensities")
    type(TMolecularOrbital) :: molorb
    logical, parameter :: useRadialLut = .true.
    logical, parameter :: useGPU = .false.
    real(dp) :: valueOnGrid(100,100,100,1)
    real(dp) :: gridVol, actual, expected
    integer :: i, idx(4)
    type(spotCheck) :: spotChecks(6) = [ &
        spotCheck([50, 31, 50, 1], 1.238914874283217_dp), & ! O at (0,-1.89,0)
        spotCheck([1, 1, 1, 1], 0.0_dp), & ! Grid bounds
        spotCheck([100, 100, 100, 1], 0.0_dp), & ! Grid bounds
        spotCheck([51, 51, 51, 1], 0.5835226183394406E-01_dp), & ! Grid Center
        spotCheck([50, 40, 51, 1], 0.6846187059267557_dp), & ! Symmetry: molecule in YZ plane
        spotCheck([52, 40, 51, 1], 0.6846187059267551_dp) & ! Symmetric: molecule in YZ plane
    ]

    call initMolorbH2O(molorb, useRadialLut)
    gridVol = abs(determinant33(molorb%system%gridVecs))

    ! Real (H2O) : Atomic densities calculation
    call getAtomicDensities(molorb, occupationAtomDens, valueOnGrid, useGPU)
    ! Check sum over grid
    expected = 8.004352_dp
    actual = sum(valueOnGrid) * gridVol
    print *, "Total charge:", actual
    @:CHECK(is_close(actual, expected, rtol=rtol))
    ! Spot check values
    do i = 1, size(spotChecks)
      idx = spotChecks(i)%idx
      expected = spotChecks(i)%expected
      actual = valueOnGrid(idx(1), idx(2), idx(3), idx(4))
      @:CHECK(is_close(actual, expected, rtol=rtol))
    end do
  $:END_TEST()

  ! -- Real (H2O) : All states calculation --
  $:TEST("molorb_real_allStates")
    type(TMolecularOrbital) :: molorb
    logical, parameter :: useRadialLut = .true.
    logical, parameter :: useGPU = .false.
    real(dp) :: valueOnGrid(100,100,100,4)
    real(dp) :: gridVol, actual, expected
    integer :: i, idx(4)
    ! 2 randomly chosen points per state.
    ! [random.randrange(1,101) for _ in "123"] 
    type(spotCheck) :: spotChecks(8) = [ &
        spotCheck([31, 82, 52, 1], -0.4237591897694508E-03_dp), & 
        spotCheck([93, 19, 63, 1], -0.8172698875983751E-04_dp), & 
        spotCheck([80, 82, 83, 2], -0.1274610827226662E-03_dp), & 
        spotCheck([32, 28, 73, 2], -0.2681614272853517E-03_dp), & 
        spotCheck([68, 17, 95, 3],  0.1151590378991745E-03), & 
        spotCheck([39, 88, 09, 3], -0.9397073043782955E-04_dp), & 
        spotCheck([09, 78, 23, 4],  0.4296017603685945E-20_dp), & 
        spotCheck([47, 59, 06, 4], -0.9990029155728657E-06_dp)  & 
    ]

    call initMolorbH2O(molorb, useRadialLut)
    gridVol = abs(determinant33(molorb%system%gridVecs))
    
    ! Real (H2O) : All states calculation
    call getValue(molorb, eigVecsReal, valueOnGrid, useGPU)
    ! Check sum over grid
    expected = -11.36846731051559_dp
    actual = sum(valueOnGrid) * gridVol
    print *, "Total charge:", actual
    @:CHECK(is_close(actual, expected, rtol=rtol))
    ! Spot check values
    do i = 1, size(spotChecks)
      idx = spotChecks(i)%idx
      expected = spotChecks(i)%expected
      actual = valueOnGrid(idx(1), idx(2), idx(3), idx(4))
      @:CHECK(is_close(actual, expected, rtol=rtol))
    end do
  $:END_TEST()


  ! ### Complex (H-chain) Test cases ###
  ! -- Complex (H-chain) : Total charge calculation --
  $:TEST("molorb_cmplx_totChrg")
    type(TMolecularOrbital) :: molorb
    logical, parameter :: useRadialLut = .false.
    logical, parameter :: useGPU = .false.
    real(dp) :: valueOnGrid(100,100,100,1)
    real(dp) :: gridVol, actual, expected
    integer :: i, idx(4)
    ! Randomly chosen points.
    ! [random.randrange(1,101) for _ in "123"] 
    type(spotCheck) :: spotChecks(6) = [ &
        spotCheck([31, 82, 52, 1], 0.1854186607568482E-04_dp), & 
        spotCheck([93, 19, 63, 1], 0.4554475419034117E-08_dp), & 
        spotCheck([80, 82, 83, 1], 0.1400485004792566E-05_dp), & 
        spotCheck([32, 28, 73, 1], 0.2821278119060284E-03_dp), & 
        spotCheck([68, 17, 95, 1], 0.1144197239608134E-04_dp), & 
        spotCheck([39, 88, 09, 1], 0.7705316328219835E-05_dp) & 
    ]
    call initMolorbHchain(molorb, useRadialLut)
    gridVol = abs(determinant33(molorb%system%gridVecs))
    ! -- Complex (H-chain) : Total charge calculation --
    call getTotalChrg(molorb, eigVecsCmpl, kPointsHchain, kIndexesHchain, valueOnGrid, occupationVecHchain, useGPU)
    ! Check sum over grid
    expected = 6.637263751554149_dp
    actual = sum(valueOnGrid) * gridVol
    print *, "Total charge:", actual
    @:CHECK(is_close(actual, expected, rtol=rtol))
    ! Spot check values
    do i = 1, size(spotChecks)
      idx = spotChecks(i)%idx
      expected = spotChecks(i)%expected
      actual = valueOnGrid(idx(1), idx(2), idx(3), idx(4))
      @:CHECK(is_close(actual, expected, rtol=rtol))
    end do
    $:END_TEST()


    ! -- Complex (H-chain) : All states calculation --
    $:TEST("molorb_cmplx_allStates")
      type(TMolecularOrbital) :: molorb
      logical, parameter :: useRadialLut = .false.
      logical, parameter :: useGPU = .false.
      complex(dp) :: valueOnGrid(100,100,100,4)
      real(dp) :: gridVol
      complex(dp) :: actual, expected
      integer :: i, idx(4)
      ! 2 randomly chosen points per state.
      ! [random.randrange(1,101) for _ in "123"] 
      type(spotCheckCmpl) :: spotChecks(8) = [ &
          spotCheckCmpl([31, 82, 52, 1], (0.4618562431412820E-02_dp, 0.5805855442516971E-04_dp)), &
          spotCheckCmpl([93, 19, 63, 1], (0.8330079947653005E-04_dp, 0.8753533677015034E-05_dp)), &
          spotCheckCmpl([80, 82, 83, 2], (0.1222651639518469E-02_dp, 0.8048672887188940E-04_dp)), &
          spotCheckCmpl([32, 28, 73, 2], (0.1690763857502260E-01_dp, 0.4646391064170055E-02_dp)), &
          spotCheckCmpl([68, 17, 95, 3], (0.2091620118572572E-02_dp, 0.2596145058129883E-02_dp)), &
          spotCheckCmpl([39, 88, 09, 3], (0.2671119431633243E-02_dp, 0.5944536018057613E-03_dp)), &
          spotCheckCmpl([09, 78, 23, 4], (0.8227151295698438E-04_dp, 0.1991128564447004E-04_dp)), &
          spotCheckCmpl([47, 59, 06, 4], (0.2967129806717818_dp, 0.8175874030009455E-02_dp)) &
      ]

      call initMolorbHchain(molorb, useRadialLut)
      gridVol = abs(determinant33(molorb%system%gridVecs))
      ! Complex (H-chain) : All states calculation
      call getValue(molorb, eigVecsCmpl, kPointsHchain, kIndexesHchain, valueOnGrid, useGPU)
      ! Check sum over grid
      expected = (113.9831688331955_dp, 45.30115350686194_dp)
      actual = sum(valueOnGrid) * gridVol
      print *, "Total charge:", actual
      @:CHECK(is_close(real(actual), real(expected), rtol=rtol))
      @:CHECK(is_close(aimag(actual), aimag(expected), rtol=rtol))

      ! Spot check values
      do i = 1, size(spotChecks)
        idx = spotChecks(i)%idx
        expected = spotChecks(i)%expected
        actual = valueOnGrid(idx(1), idx(2), idx(3), idx(4))

        @:CHECK(is_close(real(actual), real(expected), rtol=rtol))
        @:CHECK(is_close(aimag(actual), aimag(expected), rtol=rtol))
      end do
    $:END_TEST()

    ! ### Summary / todo ###
    ! Cases to handle:
    ! Parametrisation: Run on both cpu, gpu, lut ond and off (template 4x).
    ! Idea: Generate reference data with old code, spot check values & check sum
    !
    ! Additional tests:
    ! Init an sto, and spot check values
    ! Init sto using lut, check if lut access works
    !
    ! todo: Allow mixed lut/direct evaluation in gpu kernel by resampling to uniform luts in fortran
    ! Tolerance considerations:
    ! OLD   
    ! CPU     8.004 0445629655839
    ! GPU     8.004 0445629655839
    ! CPU LUT 8.004 1122133449658
    ! GPU LUT 8.004 1087880444035


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("libwavegrid", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_libwavegrid_simple
