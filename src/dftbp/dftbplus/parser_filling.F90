!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Reads the electronic filling settings from the HSD input.
module dftbp_dftbplus_parser_filling
  use hsd_data, only : hsd_table, new_table
  use dftbp_common_accuracy, only : dp, lc, minTemp
  use dftbp_common_constants, only : Boltzmann
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_unitconversion, only : energyUnits
  use dftbp_dftb_etemp, only : fillingTypes
  use dftbp_dftbplus_inputdata, only : TControl
  use hsd, only : hsd_get_or_set, hsd_get, hsd_get_table, hsd_get_choice, hsd_get_attrib, &
      & HSD_STAT_OK
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, getNodeName2
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_type_typegeometry, only : TGeometry
  implicit none

  private
  public :: readFilling, readElectronicFilling

contains


  !> Filling of electronic levels
  subroutine readFilling(node, ctrl, geo, temperatureDefault)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to test for periodicity
    type(TGeometry), intent(in) :: geo

    !> Default temperature for filling
    real(dp), intent(in) :: temperatureDefault

    type(hsd_table), pointer :: value1, child, child2, child3, field
    character(len=:), allocatable :: buffer, modifier
    character(lc) :: errorStr
    integer :: stat

    call hsd_get_table(node, "Filling", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: defTbl, defChild
        call new_table(defTbl, name="filling")
        call new_table(defChild, name="fermi")
        call defTbl%add_child(defChild)
        call node%add_child(defTbl)
      end block
      call hsd_get_table(node, "Filling", child, stat, auto_wrap=.true.)
    end if
    call hsd_get_choice(child, "", buffer, value1, stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(child, "Invalid or missing choice in 'Filling'")

    select case (buffer)
    case ("fermi")
      ctrl%iDistribFn = fillingTypes%Fermi ! Fermi function
    case ("gaussian")
      ctrl%iDistribFn = fillingTypes%Methfessel ! Gauss function broadening of levels (0th order MP)
    case ("methfesselpaxton")
      ! Set the order of the Methfessel-Paxton step function approximation, defaulting to 1st order
      call hsd_get_or_set(value1, "Order", ctrl%iDistribFn, 1)
      if (ctrl%iDistribFn < 1) then
        call getNodeName2(value1, buffer)
        select case(ctrl%iDistribFn)
        case (0)
          write(errorStr, "(A)")"Methfessel-Paxton filling order 0 is equivalent to gaussian&
              & smearing"
          call dftbp_warning(child, errorStr)
        case default
          write(errorStr, "(A,A,A,I4)")"Filling order must be above zero '", buffer,"' :",&
              &ctrl%iDistribFn
          call dftbp_error(child, errorStr)
        end select
      end if
      ctrl%iDistribFn = ctrl%iDistribFn + fillingTypes%Methfessel
    case default
      call getNodeName2(value1, buffer)
      call dftbp_error(child, "Invalid filling method '" //buffer// "'")
    end select

    if (.not. ctrl%tSetFillingTemp) then
      call hsd_get_or_set(value1, "Temperature", ctrl%tempElec, temperatureDefault)
      call hsd_get_attrib(value1, "Temperature", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call hsd_get_table(value1, "Temperature", field, stat, auto_wrap=.true.)
      call convertUnitHsd(modifier, energyUnits, field, ctrl%tempElec)
      if (ctrl%tempElec < minTemp) then
        ctrl%tempElec = minTemp
      end if
    end if

    call hsd_get_table(value1, "FixedFermiLevel", child2, stat, auto_wrap=.true.)
    ctrl%tFixEf = associated(child2)
    if (ctrl%tFixEf) then
      if (ctrl%tSpin .and. .not.ctrl%t2Component) then
        allocate(ctrl%Ef(2))
      else
        allocate(ctrl%Ef(1))
      end if
      call hsd_get(child2, "#text", ctrl%Ef, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value in 'FixedFermiLevel'")
      if (allocated(child2%attrib)) then
        modifier = child2%attrib
      else
        modifier = ""
      end if
      child3 => child2
      call convertUnitHsd(modifier, energyUnits, child3, ctrl%Ef)
    end if

    if (geo%tPeriodic .and. .not.ctrl%tFixEf) then
      call hsd_get_or_set(value1, "IndependentKFilling", ctrl%tFillKSep, .false.)
    end if

  end subroutine readFilling


  !> Parses for electronic filling temperature (should only read if not either REKS or electron
  !> dynamics from a supplied density matrix)
  subroutine readElectronicFilling(hamNode, ctrl, geo)

    !> Relevant node in input tree
    type(hsd_table), pointer :: hamNode

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to test for periodicity
    type(TGeometry), intent(in) :: geo

    select case(ctrl%hamiltonian)
    case(hamiltonianTypes%xtb)
      call readFilling(hamNode, ctrl, geo, 300.0_dp*Boltzmann)
    case(hamiltonianTypes%dftb)
      call readFilling(hamNode, ctrl, geo, 0.0_dp)
    end select

  end subroutine readElectronicFilling

end module dftbp_dftbplus_parser_filling
