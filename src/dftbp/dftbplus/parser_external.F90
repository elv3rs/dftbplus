!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Reads external field and potential settings from the input tree.
module dftbp_dftbplus_parser_external
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_common_file, only : closeFile, openFile, TFileDescr
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_unitconversion, only : EFieldUnits, energyUnits, freqUnits, lengthUnits
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_io_charmanip, only : unquote
  use dftbp_io_hsdutils, only : hsd_child_list, &
      & getChild, getChildren, getChildValue, &
      & getLength, getItem1, destroyNodeList
  use hsd, only : hsd_get_or_set
  use dftbp_io_hsdutils, only : dftbp_error, textNodeName, getNodeName
  use dftbp_io_message, only : error
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_type_linkedlist, only : append, asArray, destruct, init, intoArray, len, &
      & TListInt, TListReal, TListRealR1, TListRealR2
  use dftbp_type_typegeometry, only : TGeometry
  use hsd_data, only : hsd_table
  implicit none

  private
  public :: readExternal

contains

  ! External field(s) and potential(s)
  subroutine readExternal(node, ctrl, geo)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    type(hsd_table), pointer :: value1, child, child2, child3
    type(hsd_child_list), pointer :: children
    character(len=:), allocatable :: modifier, buffer, buffer2
    real(dp) :: rTmp
    type(TFileDescr) :: file
    integer :: ind, ii, iErr, nElem
    real(dp), allocatable :: tmpR1(:), tmpR2(:,:)
    type(TListRealR2) :: lCharges
    type(TListRealR1) :: lBlurs, lr1
    type(TListReal) :: lr
    type(TListInt) :: li

    call getChildValue(node, "ElectricField", value1, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! external applied field
    call getChild(child, "External", child2, requested=.false.)
    if (associated(child2)) then
      allocate(ctrl%electricField)
      ctrl%tMulliken = .true.
      call getChildValue(child2, "Strength", ctrl%electricField%EFieldStrength, modifier=modifier,&
          & child=child3)
      call convertUnitHsd(modifier, EFieldUnits, child3, ctrl%electricField%EFieldStrength)
      call getChildValue(child2, "Direction", ctrl%electricField%EfieldVector)
      if (sum(ctrl%electricField%EfieldVector**2) < 1e-8_dp) then
        call dftbp_error(child2,"Vector too small")
      else
        ctrl%electricField%EfieldVector = ctrl%electricField%EfieldVector&
            & / sqrt(sum(ctrl%electricField%EfieldVector**2))
      end if
      call getChildValue(child2, "Frequency", ctrl%electricField%EFieldOmega, 0.0_dp, &
          & modifier=modifier, child=child3)
      call convertUnitHsd(modifier, freqUnits, child3, ctrl%electricField%EFieldOmega)
      if (ctrl%electricField%EFieldOmega > 0.0) then
        ! angular frequency
        ctrl%electricField%EFieldOmega = 2.0_dp * pi * ctrl%electricField%EFieldOmega
        ctrl%electricField%isTDEfield = .true.
      else
        ctrl%electricField%isTDEfield = .false.
        ctrl%electricField%EFieldOmega = 0.0_dp
      end if
      ctrl%electricField%EfieldPhase = 0
      if (ctrl%electricField%isTDEfield) then
        call hsd_get_or_set(child2, "Phase", ctrl%electricField%EfieldPhase, 0)
      end if
    end if

    ctrl%nExtChrg = 0
    if (ctrl%hamiltonian == hamiltonianTypes%dftb) then

      call getChildren(child, "PointCharges", children)
      if (getLength(children) > 0) then
        ! Point charges present
        if (.not.ctrl%tSCC) then
          call error("External charges can only be used in an SCC calculation")
        end if
        call init(lCharges)
        call init(lBlurs)
        ctrl%nExtChrg = 0
        do ii = 1, getLength(children)
          call getItem1(children, ii, child2)
          call getChildValue(child2, "CoordsAndCharges", value1, modifier=modifier, child=child3)
          call getNodeName(value1, buffer)
          select case(buffer)
          case (textNodeName)
            call init(lr1)
            call getChildValue(child3, "", 4, lr1, modifier=modifier)
            allocate(tmpR2(4, len(lr1)))
            call asArray(lr1, tmpR2)
            ctrl%nExtChrg = ctrl%nExtChrg + len(lr1)
            call destruct(lr1)
          case ("directread")
            call getChildValue(value1, "Records", ind)
            call getChildValue(value1, "File", buffer2)
            allocate(tmpR2(4, ind))
            call openFile(file, unquote(buffer2), mode="r", iostat=iErr)
            if (iErr /= 0) then
              call dftbp_error(value1, "Could not open file '"&
                  & // trim(unquote(buffer2)) // "' for direct reading" )
            end if
            read(file%unit, *, iostat=iErr) tmpR2
            if (iErr /= 0) then
              call dftbp_error(value1, "Error during direct reading '"&
                  & // trim(unquote(buffer2)) // "'")
            end if
            call closeFile(file)
            ctrl%nExtChrg = ctrl%nExtChrg + ind
          case default
            call dftbp_error(value1, "Invalid block name")
          end select
          call convertUnitHsd(modifier, lengthUnits, child3, tmpR2(1:3,:))
          call append(lCharges, tmpR2)
          call getChildValue(child2, "GaussianBlurWidth", rTmp, 0.0_dp, modifier=modifier,&
              & child=child3)
          if (rTmp < 0.0_dp) then
            call dftbp_error(child3, "Gaussian blur width may not be negative")
          end if
          call convertUnitHsd(modifier, lengthUnits, child3, rTmp)
          allocate(tmpR1(size(tmpR2, dim=2)))
          tmpR1(:) = rTmp
          call append(lBlurs, tmpR1)
          deallocate(tmpR1)
          deallocate(tmpR2)
        end do

        allocate(ctrl%extChrg(4, ctrl%nExtChrg))
        ind = 1
        do ii = 1, len(lCharges)
          call intoArray(lCharges, ctrl%extChrg(:, ind:), nElem, ii)
          ind = ind + nElem
        end do
        call destruct(lCharges)

        allocate(ctrl%extChrgBlurWidth(ctrl%nExtChrg))
        ind = 1
        do ii = 1, len(lBlurs)
          call intoArray(lBlurs, ctrl%extChrgBlurWidth(ind:), nElem, ii)
          ind = ind + nElem
        end do
        call destruct(lBlurs)
        call destroyNodeList(children)
      end if

    else

      call getChildren(child, "PointCharges", children)
      if (getLength(children) > 0) then
        call dftbp_error(child, "External charges are not currently supported for this model")
      end if

    end if

    call getChild(node, "AtomSitePotential", child, requested=.false.)
    if (associated(child)) then
      allocate(ctrl%atomicExtPotential)

      call getChild(child, "Net", child2, requested=.false.)
      if (associated(child2)) then
        ! onsites
        ctrl%tNetAtomCharges = .true.
        call init(li)
        call init(lr)
        call getChildValue(child2, "Atoms", li)
        call getChildValue(child2, "Vext", lr, modifier=modifier, child=child3)
        if (len(li) /= len(lr)) then
          call dftbp_error(child2, "Mismatch in number of sites and potentials")
        end if
        allocate(ctrl%atomicExtPotential%iAtOnSite(len(li)))
        call asArray(li, ctrl%atomicExtPotential%iAtOnSite)
        allocate(ctrl%atomicExtPotential%VextOnSite(len(lr)))
        call asArray(lr, ctrl%atomicExtPotential%VextOnSite)
        call convertUnitHsd(modifier, energyUnits, child3, ctrl%atomicExtPotential%VextOnSite)
        call destruct(li)
        call destruct(lr)
      end if

      call getChild(child, "Gross", child2, requested=.false.)
      if (associated(child2)) then
        ! atomic
        call init(li)
        call init(lr)
        call getChildValue(child2, "Atoms", li)
        call getChildValue(child2, "Vext", lr, modifier=modifier, child=child3)
        if (len(li) /= len(lr)) then
          call dftbp_error(child2, "Mismatch in number of sites and potentials")
        end if
        allocate(ctrl%atomicExtPotential%iAt(len(li)))
        call asArray(li, ctrl%atomicExtPotential%iAt)
        allocate(ctrl%atomicExtPotential%Vext(len(lr)))
        call asArray(lr, ctrl%atomicExtPotential%Vext)
        call convertUnitHsd(modifier, energyUnits, child3, ctrl%atomicExtPotential%Vext)
        call destruct(li)
        call destruct(lr)
      end if

      if (.not.allocated(ctrl%atomicExtPotential%iAt)&
          & .and. .not.allocated(ctrl%atomicExtPotential%iAtOnSite)) then
        call dftbp_error(child, "No atomic potentials specified")
      end if

    end if

  end subroutine readExternal

end module dftbp_dftbplus_parser_external
