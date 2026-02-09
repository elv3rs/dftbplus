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
  use hsd, only : hsd_get_or_set, hsd_get, hsd_get_matrix, hsd_get_attrib, hsd_table_ptr, &
      & hsd_get_child_tables, hsd_get_table, hsd_get_choice, HSD_STAT_OK
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

    call hsd_get_table(node, "ElectricField", child, stat, auto_wrap=.true.)

    ! external applied field
    child2 => null()
    if (associated(child)) then
      call hsd_get_table(child, "External", child2, stat, auto_wrap=.true.)
    end if
    if (associated(child2)) then
      allocate(ctrl%electricField)
      ctrl%tMulliken = .true.
      call hsd_get(child2, "Strength", ctrl%electricField%EFieldStrength, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value: 'Strength'")
      call hsd_get_attrib(child2, "Strength", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call hsd_get_table(child2, "Strength", child3, stat, auto_wrap=.true.)
      if (.not. associated(child3)) child3 => child2
      call convertUnitHsd(modifier, EFieldUnits, child3, ctrl%electricField%EFieldStrength)
      block
        real(dp), allocatable :: tmpArr(:)
        call hsd_get(child2, "Direction", tmpArr, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value: 'Direction'")
        ctrl%electricField%EfieldVector = tmpArr(:3)
      end block
      if (sum(ctrl%electricField%EfieldVector**2) < 1e-8_dp) then
        call dftbp_error(child2,"Vector too small")
      else
        ctrl%electricField%EfieldVector = ctrl%electricField%EfieldVector&
            & / sqrt(sum(ctrl%electricField%EfieldVector**2))
      end if
      call hsd_get_or_set(child2, "Frequency", ctrl%electricField%EFieldOmega, 0.0_dp)
      call hsd_get_attrib(child2, "Frequency", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call hsd_get_table(child2, "Frequency", child3, stat, auto_wrap=.true.)
      if (.not. associated(child3)) child3 => child2
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

      if (associated(child)) then
        call hsd_get_child_tables(child, "PointCharges", children)
      else
        allocate(children(0))
      end if
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
          call hsd_get_table(child2, "CoordsAndCharges", child3, stat, auto_wrap=.true.)
          if (.not. associated(child3)) call dftbp_error(child2, &
              & "Missing required block: 'CoordsAndCharges'")
          call hsd_get_attrib(child2, "CoordsAndCharges", modifier, stat)
          if (stat /= HSD_STAT_OK) modifier = ""
          call hsd_get_choice(child3, "", buffer, value1, stat)
          if (.not. associated(value1)) buffer = textNodeName
          select case(buffer)
          case (textNodeName)
            call hsd_get_matrix(child2, "CoordsAndCharges", tmpR2, nrows, ncols, stat=stat, &
                & order="column-major")
            if (stat /= 0) call dftbp_error(child2, "Error reading CoordsAndCharges")
            ctrl%nExtChrg = ctrl%nExtChrg + ncols
          case ("directread")
            call hsd_get(value1, "Records", ind, stat=stat)
            if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'Records'")
            call hsd_get(value1, "File", buffer2, stat=stat)
            if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'File'")
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
          call hsd_get_or_set(child2, "GaussianBlurWidth", rTmp, 0.0_dp)
          call hsd_get_attrib(child2, "GaussianBlurWidth", modifier, stat)
          if (stat /= HSD_STAT_OK) modifier = ""
          call hsd_get_table(child2, "GaussianBlurWidth", child3, stat, auto_wrap=.true.)
          if (.not. associated(child3)) child3 => child2
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

      if (associated(child)) then
        call hsd_get_child_tables(child, "PointCharges", children)
      else
        allocate(children(0))
      end if
      if (size(children) > 0) then
        call dftbp_error(child, "External charges are not currently supported for this model")
      end if

    end if

    call hsd_get_table(node, "AtomSitePotential", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      allocate(ctrl%atomicExtPotential)

      call hsd_get_table(child, "Net", child2, stat, auto_wrap=.true.)
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
        call hsd_get_table(child2, "Vext", child3, stat, auto_wrap=.true.)
        if (.not. associated(child3)) call dftbp_error(child2, "Missing required block: 'Vext'")
        call convertUnitHsd(modifier, energyUnits, child3, ctrl%atomicExtPotential%VextOnSite)
      end if

      call hsd_get_table(child, "Gross", child2, stat, auto_wrap=.true.)
      if (associated(child2)) then
        ! atomic
        call hsd_get(child2, "Atoms", ctrl%atomicExtPotential%iAt)
        call hsd_get(child2, "Vext", ctrl%atomicExtPotential%Vext)
        if (size(ctrl%atomicExtPotential%iAt) /= size(ctrl%atomicExtPotential%Vext)) then
          call dftbp_error(child2, "Mismatch in number of sites and potentials")
        end if
        call hsd_get_attrib(child2, "Vext", modifier)
        call hsd_get_table(child2, "Vext", child3, stat, auto_wrap=.true.)
        if (.not. associated(child3)) call dftbp_error(child2, "Missing required block: 'Vext'")
        call convertUnitHsd(modifier, energyUnits, child3, ctrl%atomicExtPotential%Vext)
      end if

      if (.not.allocated(ctrl%atomicExtPotential%iAt)&
          & .and. .not.allocated(ctrl%atomicExtPotential%iAtOnSite)) then
        call dftbp_error(child, "No atomic potentials specified")
      end if

    end if

  end subroutine readExternal

end module dftbp_dftbplus_parser_external
