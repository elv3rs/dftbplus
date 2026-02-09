!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Reads electrostatics and medium-range DFTB parameters from HSD input.
module dftbp_dftbplus_parser_electrostatics
  use dftbp_common_accuracy, only : dp
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_extlibs_poisson, only : TPoissonInfo, withPoisson
  use dftbp_io_hsdutils, only : getNodeName2, dftbp_error
  use hsd, only : hsd_get_or_set, hsd_get_table, hsd_get_choice, HSD_STAT_OK
  use dftbp_type_typegeometry, only : TGeometry
  use hsd_data, only : hsd_table, new_table
#:if WITH_TRANSPORT
  use dftbp_transport_negfvars, only : TTransPar
#:endif
#:if WITH_POISSON
  use dftbp_dftbplus_parser_transport, only : readPoisson
#:endif
  implicit none

  private
  public :: readElectrostatics, readMdftb

contains


#:if WITH_TRANSPORT
  subroutine readElectrostatics(node, ctrl, geo, tp, poisson)
#:else
  subroutine readElectrostatics(node, ctrl, geo, poisson)
#:endif

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp
  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer
    integer :: stat

    ctrl%tPoisson = .false.

    ! Read in which kind of electrostatics method to use.
    call hsd_get_table(node, "Electrostatics", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: defTbl, defChild
        call new_table(defTbl, name="electrostatics")
        call new_table(defChild, name="gammafunctional")
        call defTbl%add_child(defChild)
        call node%add_child(defTbl)
      end block
      call hsd_get_table(node, "Electrostatics", child, stat, auto_wrap=.true.)
    end if
    call hsd_get_choice(child, "", buffer, value1, stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(child, "Invalid or missing choice in 'Electrostatics'")

    select case (buffer)

    case ("gammafunctional")
    #:if WITH_TRANSPORT
      if (tp%taskUpload .and. ctrl%tSCC) then
        call dftbp_error(value1, "GammaFunctional not available, if you upload contacts in an SCC&
            & calculation.")
      end if
    #:endif

    case ("poisson")
      if (.not. withPoisson) then
        call dftbp_error(value1, "Poisson not available as binary was built without the Poisson&
            &-solver")
      end if
      #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
        ctrl%tPoisson = .true.
        #:if WITH_TRANSPORT
          call readPoisson(value1, poisson, geo%tPeriodic, tp, geo%latVecs, ctrl%updateSccAfterDiag)
        #:else
          call readPoisson(value1, poisson, geo%tPeriodic, geo%latVecs, ctrl%updateSccAfterDiag)
        #:endif
      #:endblock

    case default
      call getNodeName2(value1, buffer)
      call dftbp_error(child, "Unknown electrostatics '" // buffer // "'")
    end select

  end subroutine readElectrostatics


  !> Read in the mdftb parameters
  subroutine readMdftb(node, ctrl, geo)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    type(hsd_table), pointer :: value1, child, child2
    character(len=:), allocatable :: buffer
    integer :: iSp1, stat

    ctrl%isMdftb = .false.
    if (ctrl%tSCC) then
      call hsd_get_table(node, "Mdftb", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) then
        block
          type(hsd_table) :: defTbl, defChild
          call new_table(defTbl, name="mdftb")
          call new_table(defChild, name="none")
          call defTbl%add_child(defChild)
          call node%add_child(defTbl)
        end block
        call hsd_get_table(node, "Mdftb", child, stat, auto_wrap=.true.)
      end if
      call hsd_get_choice(child, "", buffer, value1, stat)
      if (associated(value1)) then
        select case(buffer)
        case("onecenterapproximation")
          ctrl%isMdftb = .true.
          allocate(ctrl%mdftbAtomicIntegrals)
          allocate(ctrl%mdftbAtomicIntegrals%DScaling(geo%nSpecies), source=1.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%QScaling(geo%nSpecies), source=1.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%SXPx(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%PxXDxxyy(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%PxXDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%PyYDxxyy(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%PzZDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%SXXS(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%PxXXPx(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%PyXXPy(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%SXXDxxyy(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%SXXDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%SYYDxxyy(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%SZZDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DxyXXDxy(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DyzXXDyz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DxxyyXXDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DzzXXDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DxxyyYYDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DzzZZDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DxzXZDzz(geo%nSpecies), source=0.0_dp)
          allocate(ctrl%mdftbAtomicIntegrals%DyzYZDxxyy(geo%nSpecies), source=0.0_dp)

          call hsd_get_table(value1, "AtomDIntegralScalings", child2, stat, auto_wrap=.true.)
          if (associated(child2)) then
            do iSp1 = 1, geo%nSpecies
              call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1)),&
                  & ctrl%mdftbAtomicIntegrals%DScaling(iSp1), 1.0_dp)
            end do
          end if

          call hsd_get_table(value1, "AtomQIntegralScalings", child2, stat, auto_wrap=.true.)
          if (associated(child2)) then
            do iSp1 = 1, geo%nSpecies
              call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1)),&
                  & ctrl%mdftbAtomicIntegrals%QScaling(iSp1), 1.0_dp)
            end do
          end if

          call hsd_get_table(value1, "OneCenterAtomIntegrals", child2, stat, auto_wrap=.true.)
          if (.not. associated(child2)) call dftbp_error(value1, &
              & "Missing required block: 'OneCenterAtomIntegrals'")
          do iSp1 = 1, geo%nSpecies
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":S|X|Px",&
                & ctrl%mdftbAtomicIntegrals%SXPx(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Px|X|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%PxXDxxyy (iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Px|X|Dzz",&
                & ctrl%mdftbAtomicIntegrals%PxXDzz(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Py|Y|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%PyYDxxyy(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Pz|Z|Dzz",&
                & ctrl%mdftbAtomicIntegrals%PzZDzz(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":S|XX|S",&
                & ctrl%mdftbAtomicIntegrals%SXXS(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Px|XX|Px",&
                & ctrl%mdftbAtomicIntegrals%PxXXPx(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Py|XX|Py",&
                & ctrl%mdftbAtomicIntegrals%PyXXPy(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":S|XX|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%SXXDxxyy(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":S|XX|Dzz",&
                & ctrl%mdftbAtomicIntegrals%SXXDzz(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":S|YY|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%SYYDxxyy(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":S|ZZ|Dzz",&
                & ctrl%mdftbAtomicIntegrals%SZZDzz(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Dxy|XX|Dxy",&
                & ctrl%mdftbAtomicIntegrals%DxyXXDxy(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Dyz|XX|Dyz",&
                & ctrl%mdftbAtomicIntegrals%DyzXXDyz(iSp1), 0.0_dp)
            !call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Dxx-yy|XX|Dzz",&
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Dzz|XX|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%DxxyyXXDzz(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Dzz|XX|Dzz",&
                & ctrl%mdftbAtomicIntegrals%DzzXXDzz(iSp1), 0.0_dp)
            !call getChildValue(child2, trim(geo%speciesNames(iSp1))//":Dxx-yy|YY|Dzz",&
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Dzz|YY|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%DxxyyYYDzz(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Dzz|ZZ|Dzz",&
                & ctrl%mdftbAtomicIntegrals%DzzZZDzz(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Dxz|XZ|Dzz",&
                & ctrl%mdftbAtomicIntegrals%DxzXZDzz(iSp1), 0.0_dp)
            call hsd_get_or_set(child2, trim(geo%speciesNames(iSp1))//":Dyz|YZ|Dxx-yy",&
                & ctrl%mdftbAtomicIntegrals%DyzYZDxxyy(iSp1), 0.0_dp)
          end do
        case("none")
          ctrl%isMdftb = .false.
        case default
          call dftbp_error(child,"Unknown functions :"// buffer)
        end select
      end if
    end if

  end subroutine readMdftb

end module dftbp_dftbplus_parser_electrostatics
