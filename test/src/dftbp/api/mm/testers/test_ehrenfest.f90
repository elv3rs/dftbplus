!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

program test_ehrenfest
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbplus, only : dumpHsd, hsd_table, getDftbPlusApi, getDftbPlusBuild, &
      & TDftbPlus, TDftbPlus_init, TDftbPlusInput
  use hsd, only : hsd_set, new_table
  use hsd_data, only : hsd_node
  use dftbp_common_constants, only : AA__Bohr, eV__Hartree, fs__au, imag, pi, V_m__au
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: nAtom = 12

  ! benzene coordinates (atomic units)
  real(dp), parameter :: initialCoords(3, nAtom) = reshape([&
       &    0.00000000_dp,      0.41727209_dp,      1.34035331_dp,&
       &    0.00000000_dp,      1.36277581_dp,      0.31264346_dp,&
       &    0.00000000_dp,      0.94549248_dp,     -1.02003049_dp,&
       &    0.00000000_dp,     -0.41727209_dp,     -1.32501204_dp,&
       &    0.00000000_dp,     -1.36277581_dp,     -0.29730219_dp,&
       &    0.00000000_dp,     -0.94549248_dp,      1.03537176_dp,&
       &    0.00000000_dp,      0.74550740_dp,      2.38862741_dp,&
       &    0.00000000_dp,      2.43471932_dp,      0.55254238_dp,&
       &    0.00000000_dp,      1.68922144_dp,     -1.82842183_dp,&
       &    0.00000000_dp,     -0.74550739_dp,     -2.37328614_dp,&
       &    0.00000000_dp,     -2.43471932_dp,     -0.53720110_dp,&
       &    0.00000000_dp,     -1.68922144_dp,      1.84376309_dp], [3, nAtom])

  ! benzene atom types
  integer, parameter :: species(nAtom) = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]

  integer, parameter :: nsteps = 2000
  real(dp), parameter :: timestep = 0.2_dp ! au
  real(dp), parameter :: fstrength = 0.01_dp ! V/angstrom
  real(dp), parameter :: poldir(3) = [ 0.0_dp, 1.0_dp, 1.0_dp ]
  real(dp), parameter :: omega = 6.795 ! eV

  call main_()

contains

  !! Main test routine
  !!
  !! All non-constant variables must be defined here to ensure that they are all explicitly
  !! deallocated before the program finishes (avoiding residual memory that tools like valgrind
  !! notice).
  !!
  subroutine main_()

    type(TDftbPlus) :: dftbp
    type(TDftbPlusInput) :: input

    real(dp) :: coords(3, nAtom), merminEnergy, dipole(3, 1), energy, atomNetCharges(nAtom, 1)
    real(dp) :: force(3, nAtom)
    real(dp) :: norm, fielddir(3), angFreq, envelope, field(3), time
    real(dp) :: time0 = 0.0_dp, time1 = 6.0_dp ! fs
    type(hsd_table), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pType2Files, pElecDyn
    type(hsd_table), pointer :: pPerturb, pLaser

    character(:), allocatable :: DftbVersion
    integer :: major, minor, patch, istep, ii


    call getDftbPlusBuild(DftbVersion)
    write(*,*)'DFTB+ build: ' // "'" // trim(DftbVersion) // "'"
    call getDftbPlusApi(major, minor, patch)
    write(*,"(1X,A,1X,I0,'.',I0,'.',I0)")'API version:', major, minor, patch

    call TDftbPlus_init(dftbp)

    call dftbp%getEmptyInput(input)
    call input%getRootNode(pRoot)
    call create_child_table(pRoot, "geometry", pGeo)
    call hsd_set(pGeo, "periodic", .false.)
    call hsd_set(pGeo, "typenames", ["C", "H"])
    coords(:,:) = 0.0_dp
    call setTypesAndCoords(pGeo, "typesandcoordinates", reshape(species, [1, size(species)]), coords)
    call create_child_table(pRoot, "hamiltonian", pHam)
    call create_child_table(pHam, "dftb", pDftb)
    call hsd_set(pDftb, "scc", .true.)
    call hsd_set(pDftb, "scctolerance", 1e-10_dp)

    call create_child_table(pDftb, "maxangularmomentum", pMaxAng)
    call hsd_set(pMaxAng, "c", "p")
    call hsd_set(pMaxAng, "h", "s")

    call create_child_table(pDftb, "slaterkosterfiles", pSlakos)
    call create_child_table(pSlakos, "type2filenames", pType2Files)
    call hsd_set(pType2Files, "prefix", "slakos/origin/mio-1-1/")
    call hsd_set(pType2Files, "separator", "-")
    call hsd_set(pType2Files, "suffix", ".skf")

    !  set up electron dynamics options
    call create_child_table(pRoot, "electrondynamics", pElecDyn)
    call hsd_set(pElecDyn, "steps", nsteps)
    call hsd_set(pElecDyn, "timestep", timestep)
    call hsd_set(pElecDyn, "fieldstrength", fstrength*1.0e10_dp*V_m__au)
    call hsd_set(pElecDyn, "iondynamics", .true.)
    call hsd_set(pElecDyn, "initialtemperature", 0.0_dp)

    call create_child_table(pElecDyn, "perturbation", pPerturb)
    call create_child_table(pPerturb, "laser", pLaser)
    ! these twovalues will be overriden
    call hsd_set(pLaser, "polarisationdirection", [0.0_dp, 1.0_dp, 1.0_dp])
    call hsd_set(pLaser, "laserenergy", 1.0_dp)

    norm = sqrt(dot_product(poldir, poldir))
    fielddir = poldir / norm
    angFreq = omega * eV__Hartree
    time0 = time0 * fs__au
    time1 = time1 * fs__au

    print "(A)", 'Input tree in HSD format:'
    call dumpHsd(input%hsdTree, output_unit)

    ! initialise the DFTB+ calculator
    call dftbp%setupCalculator(input)

    ! Replace coordinates
    coords(:,:) = initialCoords*AA__Bohr
    call dftbp%setGeometry(coords)

    ! get ground state
    call dftbp%getEnergy(merminEnergy)

    call dftbp%initializeTimeProp(timestep, .true., .false.)

    do istep = 1, nsteps
      ! calculate field at present timestep
      time = istep * timestep
      if (time >= time0 .and. time <= time1) then
        envelope = sin(pi*(time - time0)/(time1 - time0))**2
      else
        envelope = 0.0_dp
      end if
      field = fstrength * 1.0e10_dp * V_m__au * envelope * aimag(exp(imag*time*angFreq) * fielddir)
      call dftbp%setTdElectricField(field)
      call dftbp%doOneTdStep(istep, dipole=dipole, energy=energy, atomNetCharges=atomNetCharges,&
          & coord=coords, force=force)
    end do

    print "(A,F15.10)", 'Final total Energy:', energy
    print "(A,3F15.10)", 'Final dipole:', (dipole(ii,1), ii=1,3)
    print "(A,100F15.10)", 'Final net atomic charges:', (atomNetCharges(ii,1), ii=1,nAtom)
    print "(A,100F15.10)", 'Final coordinates:', (coords(:,ii), ii=1,nAtom)

    ! Write file for internal test system
    call writeAutotestTag(tdEnergy=energy, tdDipole=dipole, tdCharges=atomNetCharges,&
        & tdCoords=coords, tdForces=force)

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


end program test_ehrenfest
