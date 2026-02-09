!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> example code for changing internal charges inside calculation
program test_intcharges
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbplus, only : convertAtomTypesToSpecies, dumpHsd, hsd_table, getDftbPlusApi, getDftbPlusBuild,&
      & getMaxAngFromSlakoFile, TDftbPlus, TDftbPlus_init, TDftbPlusInput
  use hsd, only : hsd_set, new_table
  use hsd_data, only : hsd_node
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: nAtom = 3

  ! H2O coordinates, atomic units
  real(dp), parameter :: initialCoords(3, nAtom) = reshape([&
      & 0.000000000000000E+00_dp, -0.188972598857892E+01_dp,  0.000000000000000E+00_dp,&
      & 0.000000000000000E+00_dp,  0.000000000000000E+00_dp,  0.147977639152057E+01_dp,&
      & 0.000000000000000E+00_dp,  0.000000000000000E+00_dp, -0.147977639152057E+01_dp], [3, nAtom])

  ! Atomic number of each atom
  integer, parameter :: atomTypes(nAtom) = [8, 1, 1]

  ! list of atoms by atomic number
  character(2), parameter :: atomTypeNames(10) = [character(2) ::&
      & "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"]

  character(100), parameter :: slakoFiles(2, 2) = reshape([character(100) :: &
      & "./O-O.skf", "./H-O.skf", "./O-H.skf", "./H-H.skf"], [2, 2])

  character(1), parameter :: maxAngNames(4) = ["s", "p", "d", "f"]

  call main_()

contains


  !! Main test routine
  !!
  !! All non-constant variables must be defined here to ensure that they are all explicitely
  !! deallocated before the program finishes  (avoiding residual memory that tools like valgrind notice).
  !!
  subroutine main_()

    type(TDftbPlus) :: dftbp
    type(TDftbPlusInput) :: input

    integer, allocatable :: species(:)
    character(2), allocatable :: speciesNames(:)
    real(dp) :: merminEnergy
    real(dp) :: coords(3, nAtom), gradients(3, nAtom)
    real(dp) :: atomCharges(nAtom), cm5Charges(nAtom), atomMasses(nAtom)
    real(dp) :: q0(nAtom)
    type(hsd_table), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pAnalysis, pCm5
    type(hsd_table), pointer :: pParserOpts

    character(:), allocatable :: DftbVersion
    integer :: major, minor, patch
    integer :: ii

    call getDftbPlusBuild(DftbVersion)
    write(*,*)'DFTB+ build: ' // "'" // trim(DftbVersion) // "'"
    call getDftbPlusApi(major, minor, patch)
    write(*,"(1X,A,1X,I0,'.',I0,'.',I0)")'API version:', major, minor, patch

    ! Note: could pass an optional file unit to use for standard output as the argument outputUnit
    call TDftbPlus_init(dftbp)

    call dftbp%getEmptyInput(input)
    call input%getRootNode(pRoot)
    call create_child_table(pRoot, "geometry", pGeo)
    call hsd_set(pGeo, "periodic", .false.)

    ! Demonstrates how to convert the atom types if they are not numbered from 1 but use atomic
    ! numbers instead. The atomTypeNames array is optional, if not present, the resulting type names
    ! (which will have to be used at other places) will be X1, X2, etc.
    print "(A)", "Converting atom types"
    call convertAtomTypesToSpecies(atomTypes, species, speciesNames, atomTypeNames)

    call hsd_set(pGeo, "typenames", speciesNames)
    coords(:,:) = 0.0_dp
    call setTypesAndCoords(pGeo, "typesandcoordinates", reshape(species, [1, size(species)]), coords)
    call create_child_table(pRoot, "hamiltonian", pHam)
    call create_child_table(pHam, "dftb", pDftb)
    call hsd_set(pDftb, "scc", .true.)
    call hsd_set(pDftb, "scctolerance", 1e-12_dp)
    call create_child_table(pDftb, "maxangularmomentum", pMaxAng)

    ! read angular momenta from SK data
    call hsd_set(pMaxAng, to_lower(speciesNames(1)),&
        & maxAngNames(getMaxAngFromSlakoFile(slakoFiles(1, 1)) + 1))
    call hsd_set(pMaxAng, to_lower(speciesNames(2)),&
        & maxAngNames(getMaxAngFromSlakoFile(slakoFiles(2, 2)) + 1))

    ! set up locations for SK file data
    call create_child_table(pDftb, "slaterkosterfiles", pSlakos)
    call hsd_set(pSlakos, "o-o", trim(slakoFiles(1, 1)))
    call hsd_set(pSlakos, "h-o", trim(slakoFiles(2, 1)))
    call hsd_set(pSlakos, "o-h", trim(slakoFiles(1, 2)))
    call hsd_set(pSlakos, "h-h", trim(slakoFiles(2, 2)))
    call create_child_table(pRoot, "analysis", pAnalysis)
    call hsd_set(pAnalysis, "calculateforces", .true.)
    call create_child_table(pAnalysis, "cm5", pCm5)
    call create_child_table(pRoot, "parseroptions", pParserOpts)
    call hsd_set(pParserOpts, "parserversion", 13)

    print "(A)", 'Input tree in HSD format:'
    call dumpHsd(input%hsdTree, output_unit)

    ! convert input into settings for the DFTB+ calculator
    call dftbp%setupCalculator(input)

    ! replace QM atom coordinates
    coords(:,:) = initialCoords
    call dftbp%setGeometry(coords)

    ! get energy, charges and forces
    call dftbp%getEnergy(merminEnergy)
    call dftbp%getGradients(gradients)
    call dftbp%getGrossCharges(atomCharges)
    call dftbp%getCM5Charges(cm5Charges)
    call dftbp%getAtomicMasses(atomMasses)

    print "(A)", "Neutral system"
    print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
    print "(A,3F15.10)", 'Obtained total charge:', sum(atomCharges)
    print "(A,3F15.10)", 'Obtained individual gross charges:', atomCharges
    print "(A,3F15.10)", 'Obtained CM5 charges:', cm5Charges
    print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)
    print "(A,3F15.10)", 'Obtained gradient of atom 2:', gradients(:,2)
    print "(A,3F15.10)", 'Obtained gradient of atom 3:', gradients(:,3)

    call dftbp%getRefCharges(q0)
    print "(A)", "Original reference atomic charges"
    do ii = 1, 3
      print "(I2,F8.4)", ii, q0(ii)
    end do

    print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy

    q0(1) = q0(1) + 0.5_dp
    call dftbp%setRefCharges(q0)

    call dftbp%getRefCharges(q0)
    print "(A)", "New reference charges"
    do ii = 1, 3
      print "(I2,F8.4)", ii, q0(ii)
    end do

    ! get energy, charges and forces
    call dftbp%getEnergy(merminEnergy)
    call dftbp%getGradients(gradients)
    call dftbp%getGrossCharges(atomCharges)
    call dftbp%getCM5Charges(cm5Charges)
    call dftbp%getAtomicMasses(atomMasses)

    print "(A)", "Neutral system"
    print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
    print "(A,3F15.10)", 'Obtained total charge:', sum(atomCharges)
    print "(A,3F15.10)", 'Obtained individual gross charges:', atomCharges
    print "(A,3F15.10)", 'Obtained CM5 charges:', cm5Charges
    print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)
    print "(A,3F15.10)", 'Obtained gradient of atom 2:', gradients(:,2)
    print "(A,3F15.10)", 'Obtained gradient of atom 3:', gradients(:,3)

    ! Write file for internal test system
    call writeAutotestTag(merminEnergy=merminEnergy, gradients=gradients, grossCharges=atomCharges,&
        & atomMasses=atomMasses, cm5Charges=cm5Charges)

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


  function to_lower(str) result(lower)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: lower
    integer :: ii, ic

    lower = str
    do ii = 1, len(lower)
      ic = iachar(lower(ii:ii))
      if (ic >= iachar('A') .and. ic <= iachar('Z')) lower(ii:ii) = achar(ic + 32)
    end do
  end function to_lower


end program test_intcharges
