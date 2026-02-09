!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

program test_extpot
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbplus, only : dumpHsd, hsd_table, getDftbPlusApi, getDftbPlusBuild, &
      & TDftbPlus, TDftbPlus_init, TDftbPlusInput
  use hsd, only : hsd_set, new_table
  use hsd_data, only : hsd_node
  use extchargepot, only : getPointChargeGradients, getPointChargePotential
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: nAtom = 3

  integer, parameter :: nExtChrg = 2

  ! H2O coordinates (atomic units)
  real(dp), parameter :: initialCoords(3, nAtom) = reshape([&
      & 0.000000000000000E+00_dp, -0.188972598857892E+01_dp,  0.000000000000000E+00_dp,&
      & 0.000000000000000E+00_dp,  0.000000000000000E+00_dp,  0.147977639152057E+01_dp,&
      & 0.000000000000000E+00_dp,  0.000000000000000E+00_dp, -0.147977639152057E+01_dp], [3, nAtom])

  ! H2O atom types
  integer, parameter :: species(nAtom) = [1, 2, 2]

  ! External charges (positions and charges, again atomic units)
  real(dp), parameter :: extCharges(4, nExtChrg) = reshape([&
      &-0.94486343888717805E+00_dp,-0.94486343888717794E+01_dp, 0.17007541899969201E+01_dp, 2.5_dp,&
      & 0.43463718188810203E+01_dp,-0.58581533211004997E+01_dp, 0.26456176288841000E+01_dp, -1.9_dp&
      &], [4, nExtChrg])

  call main_()

contains

  !! Main test routine
  !!
  !! All non-constant variables must be defined here to ensure that they are all explicitely
  !! deallocated before the program finishes (avoiding residual memory that tools like valgrind
  !! notice).
  !!
  subroutine main_()

    type(TDftbPlus) :: dftbp
    type(TDftbPlusInput) :: input

    real(dp) :: merminEnergy
    real(dp) :: coords(3, nAtom), gradients(3, nAtom), extPot(nAtom), extPotGrad(3, nAtom)
    real(dp) :: atomCharges(nAtom), cm5Charges(nAtom), extChargeGrads(3, nExtChrg)
    real(dp) :: esps(2)
    real(dp), parameter :: espLocations(3,2) = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp,0.1_dp,0.0_dp],&
        & [3,2])
    type(hsd_table), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pType2Files, pAnalysis, pCm5
    type(hsd_table), pointer :: pParserOpts

    character(:), allocatable :: DftbVersion
    integer :: major, minor, patch

    !integer :: devNull

    call getDftbPlusBuild(DftbVersion)
    write(*,*)'DFTB+ build: ' // "'" // trim(DftbVersion) // "'"
    call getDftbPlusApi(major, minor, patch)
    write(*,"(1X,A,1X,I0,'.',I0,'.',I0)")'API version:', major, minor, patch

    ! Note: setting the global standard output to /dev/null will also suppress run-time error
    ! messages
    !open(newunit=devNull, file="/dev/null", action="write")
    !call TDftbPlus_init(dftbp, outputUnit=devNull)
    call TDftbPlus_init(dftbp)

    call dftbp%getEmptyInput(input)
    call input%getRootNode(pRoot)
    call create_child_table(pRoot, "geometry", pGeo)
    call hsd_set(pGeo, "periodic", .false.)
    call hsd_set(pGeo, "typenames", ["O", "H"])
    coords(:,:) = 0.0_dp
    call setTypesAndCoords(pGeo, "typesandcoordinates", reshape(species, [1, size(species)]), coords)
    call create_child_table(pRoot, "hamiltonian", pHam)
    call create_child_table(pHam, "dftb", pDftb)
    call hsd_set(pDftb, "scc", .true.)
    call hsd_set(pDftb, "scctolerance", 1e-12_dp)

    ! sub-block inside hamiltonian for the maximum angular momenta
    call create_child_table(pDftb, "maxangularmomentum", pMaxAng)
    ! explicitly set the maximum angular momenta for the species
    call hsd_set(pMaxAng, "o", "p")
    call hsd_set(pMaxAng, "h", "s")

    ! get the SK data
    ! You should provide the skfiles as found in the external/slakos/origin/mio-1-1/ folder. These
    ! can be downloaded with the utils/get_opt_externals script
    call create_child_table(pDftb, "slaterkosterfiles", pSlakos)
    call create_child_table(pSlakos, "type2filenames", pType2Files)
    call hsd_set(pType2Files, "prefix", "./")
    call hsd_set(pType2Files, "separator", "-")
    call hsd_set(pType2Files, "suffix", ".skf")

    !  set up analysis options
    call create_child_table(pRoot, "analysis", pAnalysis)
    call hsd_set(pAnalysis, "calculateforces", .true.)
    call create_child_table(pAnalysis, "cm5", pCm5)

    call create_child_table(pRoot, "parseroptions", pParserOpts)
    call hsd_set(pParserOpts, "parserversion", 5)

    print "(A)", 'Input tree in HSD format:'
    call dumpHsd(input%hsdTree, output_unit)

    ! initialise the DFTB+ calculator
    call dftbp%setupCalculator(input)

    ! Replace coordinates
    coords(:,:) = initialCoords
    call dftbp%setGeometry(coords)

    ! add external point charges
    call getPointChargePotential(extCharges(1:3,:), extCharges(4,:), coords, extPot, extPotGrad)
    call dftbp%setExternalPotential(atomPot=extPot, potGrad=extPotGrad)

    ! get results
    call dftbp%getEnergy(merminEnergy)
    call dftbp%getGradients(gradients)
    call dftbp%getGrossCharges(atomCharges)
    call dftbp%getCM5Charges(cm5Charges)
    call dftbp%getElStatPotential(esps, espLocations)
    call getPointChargeGradients(coords, atomCharges, extCharges(1:3,:), extCharges(4,:),&
        & extChargeGrads)

    print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
    print "(A,3F15.10)", 'Obtained gross charges:', atomCharges
    print "(A,3F15.10)", 'Obtained cm5 charges:', cm5Charges
    print "(A,2F15.10)", 'Obtained electrostatic potentials:', esps
    print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)
    print "(A,3F15.10)", 'Obtained gradient of atom 2:', gradients(:,2)
    print "(A,3F15.10)", 'Obtained gradient of atom 3:', gradients(:,3)
    print "(A,3F15.10)", 'Obtained gradient of charge 1:', extChargeGrads(:,1)
    print "(A,3F15.10)", 'Obtained gradient of charge 2:', extChargeGrads(:,2)

    ! Write file for internal test system
    call writeAutotestTag(merminEnergy=merminEnergy, gradients=gradients, grossCharges=atomCharges,&
        & extChargeGradients=extChargeGrads, potential=esps, cm5Charges=cm5Charges)

  end subroutine main_


  subroutine create_child_table(parent, name, child)
    type(hsd_table), intent(inout), target :: parent
    character(len=*), intent(in) :: name
    type(hsd_table), pointer, intent(out) :: child

    type(hsd_table) :: tmp
    class(hsd_node), pointer :: stored

    call new_table(tmp, name=name)
    call parent%add_child(tmp)
    call parent%get_child(parent%num_children, stored)
    select type (t => stored)
    type is (hsd_table)
      child => t
    end select
  end subroutine create_child_table


  subroutine setTypesAndCoords(node, name, types, coords)
    type(hsd_table), intent(inout) :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: types(:,:)
    real(dp), intent(in) :: coords(:,:)

    character(len=100) :: buffer
    character(len=:), allocatable :: strBuffer
    integer :: ii, jj, nRow, nCol1, nCol2

    nRow = size(types, dim=2)
    nCol1 = size(types, dim=1)
    nCol2 = size(coords, dim=1)
    strBuffer = ""
    do ii = 1, nRow
      do jj = 1, nCol1
        write(buffer, *) types(jj, ii)
        strBuffer = strBuffer // " " // trim(adjustl(buffer))
      end do
      do jj = 1, nCol2
        write(buffer, *) coords(jj, ii)
        strBuffer = strBuffer // " " // trim(adjustl(buffer))
      end do
      if (ii < nRow) strBuffer = strBuffer // new_line('a')
    end do
    call hsd_set(node, name, trim(adjustl(strBuffer)))
  end subroutine setTypesAndCoords


end program test_extpot
