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
  use dftbp_io_hsdutils, only : getChild, getChildValue
  use hsd, only : hsd_get_or_set, hsd_get, hsd_get_matrix, hsd_get_attrib, hsd_table_ptr, &
      & hsd_get_child_tables
  use dftbp_io_hsdutils, only : dftbp_error, textNodeName, getNodeName
  use dftbp_io_message, only : error
  use dftbp_io_unitconv, only : convertUnitHsd

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
    type(hsd_table_ptr), allocatable :: children(:)
    character(len=:), allocatable :: modifier, buffer, buffer2
    real(dp) :: rTmp
    type(TFileDescr) :: file
    integer :: ind, ii, iErr, stat, nrows, ncols, prevN
    real(dp), allocatable :: tmpR2(:,:)
    real(dp), allocatable :: allCharges(:,:), allBlurs(:)

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

      call hsd_get_child_tables(child, "PointCharges", children)
      if (size(children) > 0) then
        ! Point charges present
        if (.not.ctrl%tSCC) then
          call error("External charges can only be used in an SCC calculation")
        end if
        allocate(allCharges(4, 0))
        allocate(allBlurs(0))
        ctrl%nExtChrg = 0
        do ii = 1, size(children)
          child2 => children(ii)%ptr
          call getChildValue(child2, "CoordsAndCharges", value1, modifier=modifier, child=child3)
          call getNodeName(value1, buffer)
          select case(buffer)
          case (textNodeName)
            call hsd_get_matrix(child2, "CoordsAndCharges", tmpR2, nrows, ncols, stat=stat, &
                & order="column-major")
            if (stat /= 0) call dftbp_error(child2, "Error reading CoordsAndCharges")
            ctrl%nExtChrg = ctrl%nExtChrg + ncols
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
          prevN = size(allCharges, 2)
          block
            real(dp), allocatable :: tmpChrg(:,:)
            call move_alloc(allCharges, tmpChrg)
            allocate(allCharges(4, prevN + size(tmpR2, 2)))
            if (prevN > 0) allCharges(:, 1:prevN) = tmpChrg
            allCharges(:, prevN+1:) = tmpR2
          end block
          call getChildValue(child2, "GaussianBlurWidth", rTmp, 0.0_dp, modifier=modifier,&
              & child=child3)
          if (rTmp < 0.0_dp) then
            call dftbp_error(child3, "Gaussian blur width may not be negative")
          end if
          call convertUnitHsd(modifier, lengthUnits, child3, rTmp)
          prevN = size(allBlurs)
          block
            real(dp), allocatable :: tmpBlr(:)
            call move_alloc(allBlurs, tmpBlr)
            allocate(allBlurs(prevN + size(tmpR2, 2)))
            if (prevN > 0) allBlurs(1:prevN) = tmpBlr
            allBlurs(prevN+1:) = rTmp
          end block
          deallocate(tmpR2)
        end do

        call move_alloc(allCharges, ctrl%extChrg)
        call move_alloc(allBlurs, ctrl%extChrgBlurWidth)
      end if

    else

      call hsd_get_child_tables(child, "PointCharges", children)
      if (size(children) > 0) then
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
        call hsd_get(child2, "Atoms", ctrl%atomicExtPotential%iAtOnSite)
        call hsd_get(child2, "Vext", ctrl%atomicExtPotential%VextOnSite)
        if (size(ctrl%atomicExtPotential%iAtOnSite) /= &
            & size(ctrl%atomicExtPotential%VextOnSite)) then
          call dftbp_error(child2, "Mismatch in number of sites and potentials")
        end if
        call hsd_get_attrib(child2, "Vext", modifier)
        call getChild(child2, "Vext", child3)
        call convertUnitHsd(modifier, energyUnits, child3, ctrl%atomicExtPotential%VextOnSite)
      end if

      call getChild(child, "Gross", child2, requested=.false.)
      if (associated(child2)) then
        ! atomic
        call hsd_get(child2, "Atoms", ctrl%atomicExtPotential%iAt)
        call hsd_get(child2, "Vext", ctrl%atomicExtPotential%Vext)
        if (size(ctrl%atomicExtPotential%iAt) /= size(ctrl%atomicExtPotential%Vext)) then
          call dftbp_error(child2, "Mismatch in number of sites and potentials")
        end if
        call hsd_get_attrib(child2, "Vext", modifier)
        call getChild(child2, "Vext", child3)
        call convertUnitHsd(modifier, energyUnits, child3, ctrl%atomicExtPotential%Vext)
      end if

      if (.not.allocated(ctrl%atomicExtPotential%iAt)&
          & .and. .not.allocated(ctrl%atomicExtPotential%iAtOnSite)) then
        call dftbp_error(child, "No atomic potentials specified")
      end if

    end if

  end subroutine readExternal

end module dftbp_dftbplus_parser_external
