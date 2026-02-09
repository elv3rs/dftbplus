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
  use hsd_data, only : hsd_table
  use dftbp_common_accuracy, only : dp, lc, minTemp
  use dftbp_common_constants, only : Boltzmann
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_unitconversion, only : energyUnits
  use dftbp_dftb_etemp, only : fillingTypes
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_io_hsdutils, only : getChild, getChildValue
  use hsd, only : hsd_get_or_set
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, getNodeName, getNodeHSDName
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

    call getChildValue(node, "Filling", value1, "Fermi", child=child)
    call getNodeName(value1, buffer)

    select case (buffer)
    case ("fermi")
      ctrl%iDistribFn = fillingTypes%Fermi ! Fermi function
    case ("gaussian")
      ctrl%iDistribFn = fillingTypes%Methfessel ! Gauss function broadening of levels (0th order MP)
    case ("methfesselpaxton")
      ! Set the order of the Methfessel-Paxton step function approximation, defaulting to 1st order
      call hsd_get_or_set(value1, "Order", ctrl%iDistribFn, 1)
      if (ctrl%iDistribFn < 1) then
        call getNodeHSDName(value1, buffer)
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
      call getNodeHSDName(value1, buffer)
      call dftbp_error(child, "Invalid filling method '" //buffer// "'")
    end select

    if (.not. ctrl%tSetFillingTemp) then
      call getChildValue(value1, "Temperature", ctrl%tempElec, temperatureDefault, &
          &modifier=modifier, child=field)
      call convertUnitHsd(modifier, energyUnits, field, ctrl%tempElec)
      if (ctrl%tempElec < minTemp) then
        ctrl%tempElec = minTemp
      end if
    end if

    call getChild(value1, "FixedFermiLevel", child=child2, modifier=modifier, requested=.false.)
    ctrl%tFixEf = associated(child2)
    if (ctrl%tFixEf) then
      if (ctrl%tSpin .and. .not.ctrl%t2Component) then
        allocate(ctrl%Ef(2))
      else
        allocate(ctrl%Ef(1))
      end if
      call getChildValue(child2, "", ctrl%Ef, modifier=modifier, child=child3)
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
