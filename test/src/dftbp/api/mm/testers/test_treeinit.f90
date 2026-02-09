!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Test code which builds the DFTB input tree and then evaluates energy and forces. Example is a
!> periodic geometry with k-points.
program test_treeinit
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbplus, only : dumpHsd, hsd_table, getDftbPlusApi, getDftbPlusBuild, &
      & TDftbPlus, TDftbPlus_init, TDftbPlusInput
  use hsd, only : hsd_set, new_table
  use hsd_data, only : hsd_node
  use dftbp_common_constants, only : AA__Bohr
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! reference coordinates, (xyz,:nAtom) in atomic units
  real(dp), parameter :: initialCoords(3, 2) = reshape([&
      & 0.0000000000000000_dp, 0.0000000000000000_dp, 0.0000000000000000_dp,&
      & 2.5639291987021915_dp, 2.5639291987021915_dp, 2.5639291987021915_dp], [3, 2])

  ! lattice vectors in atomic units
  real(dp), parameter :: initialLatVecs(3, 3) = reshape([&
      & 5.1278583974043830_dp, 5.1278583974043830_dp, 0.0000000000000000_dp,&
      & 0.0000000000000000_dp, 5.1278583974043830_dp, 5.1278583974043830_dp,&
      & 5.1278583974043830_dp, 0.0000000000000000_dp, 5.1278583974043830_dp], [3, 3])

  call main_()

contains


  !! Main test routine
  !!
  !! All non-constant variables must be defined here to ensure that they are all explicitely
  !! deallocated before the program finishes (avoiding residual memory that tools like valgrind notice).
  !!
  subroutine main_()

    ! DFTB+ calculation itself
    type(TDftbPlus) :: dftbp
    ! input settings
    type(TDftbPlusInput) :: input

    real(dp) :: merminEnergy
    real(dp) :: coords(3, 2), latVecs(3, 3), gradients(3, 2), stressTensor(3,3)

    ! pointers to the parts of the input tree that will be set
    type(hsd_table), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pOptions, pParserOpts

    character(:), allocatable :: DftbVersion
    integer :: major, minor, patch

    !integer :: devNull

    call getDftbPlusBuild(DftbVersion)
    write(*,*)'DFTB+ build: ' // "'" // trim(DftbVersion) // "'"
    call getDftbPlusApi(major, minor, patch)
    write(*,"(1X,A,1X,I0,'.',I0,'.',I0)")'API version:', major, minor, patch

    ! Note: setting the global standard output to /dev/null will also suppress run-time error messages
    !open(newunit=devNull, file="/dev/null", action="write")
    !call TDftbPlus_init(dftbp, outputUnit=devNull)
    call TDftbPlus_init(dftbp)

    ! You should provide the skfiles found in the external/slakos/origin/pbc-0-3/ folder. These can be
    ! downloaded with the utils/get_opt_externals script.

    ! initialise a DFTB input tree and populate entries which do not have relevant or appropriate
    ! default values
    call dftbp%getEmptyInput(input)
    call input%getRootNode(pRoot)
    call create_child_table(pRoot, "geometry", pGeo)
    call hsd_set(pGeo, "periodic", .true.)
    call hsd_set(pGeo, "latticevectors", initialLatVecs)
    call hsd_set(pGeo, "typenames", ["Si"])
    coords(:,:) = 0.0_dp
    call setTypesAndCoords(pGeo, "typesandcoordinates", reshape([1, 1], [1, 2]), coords)
    call create_child_table(pRoot, "hamiltonian", pHam)
    call create_child_table(pHam, "dftb", pDftb)
    call hsd_set(pDftb, "scc", .false.)
    call create_child_table(pDftb, "maxangularmomentum", pMaxAng)
    call hsd_set(pMaxAng, "si", "p")
    call create_child_table(pDftb, "slaterkosterfiles", pSlakos)
    call hsd_set(pSlakos, "si-si", "./Si-Si.skf")
    call hsd_set(pDftb, "kpointsandweights", reshape([&
        &  0.25_dp,  0.25_dp, 0.25_dp, 1.00_dp,&
        & -0.25_dp,  0.25_dp, 0.25_dp, 1.00_dp,&
        &  0.25_dp, -0.25_dp, 0.25_dp, 1.00_dp,&
        & -0.25_dp, -0.25_dp, 0.25_dp, 1.00_dp], [4, 4], order=[2, 1]))
    call create_child_table(pRoot, "options", pOptions)
    call hsd_set(pOptions, "calculateforces", .true.)
    call create_child_table(pRoot, "parseroptions", pParserOpts)
    call hsd_set(pParserOpts, "parserversion", 3)

    ! print resulting input file, including defaults
    print *, 'Input tree in HSD format:'
    call dumpHsd(input%hsdTree, output_unit)
    print *

    ! parse the input for the DFTB+ instance
    call dftbp%setupCalculator(input)

    ! set the lattice vectors and coordinates in the document tree
    latVecs(:,:) = initialLatVecs
    coords(:,:) = initialCoords
    call dftbp%setGeometry(coords, latVecs)
    call dftbp%getEnergy(merminEnergy)
    call dftbp%getGradients(gradients)
    print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
    print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)

    ! make a small displacement in the lattice vectors and coordinates
    latVecs(1, 1) = latVecs(1, 1) + 0.1_dp * AA__Bohr
    coords(1, 1) = coords(1, 1) + 0.1_dp * AA__Bohr
    call dftbp%setGeometry(coords, latVecs)

    ! re-calculate energy and forces
    call dftbp%getEnergy(merminEnergy)
    call dftbp%getGradients(gradients)
    call dftbp%getStressTensor(stressTensor)
    print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
    print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)

    ! Write file for internal test system
    call writeAutotestTag(merminEnergy=merminEnergy, gradients=gradients, stressTensor=stressTensor)

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


end program test_treeinit
