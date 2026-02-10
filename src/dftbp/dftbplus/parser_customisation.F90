!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Subroutines for reading customised Hubbard and reference occupation parameters.
module dftbp_dftbplus_parser_customisation
  use dftbp_common_accuracy, only : dp, sc
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, getSelectedAtomIndices
  use hsd, only : hsd_rename_child, hsd_table_ptr, hsd_get_child_tables, hsd_get, hsd_get_or_set,&
      & hsd_get_table, HSD_STAT_OK, hsd_schema_t, hsd_error_t, schema_init, schema_add_field, &
      & schema_validate, schema_destroy, FIELD_OPTIONAL, FIELD_REQUIRED, FIELD_TYPE_REAL, &
      & FIELD_TYPE_TABLE, FIELD_TYPE_STRING
  use hsd_data, only : hsd_table
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_orbitals, only : getShellnames
  use dftbp_type_typegeometry, only : TGeometry
  use dftbp_type_wrappedintr, only : TWrappedInt1

  implicit none

  private
  public :: readCustomisedHubbards, readCustomReferenceOcc, is_numeric

contains


  !> Reads customised Hubbard U values that over-ride the SK file values
  subroutine readCustomisedHubbards(node, geo, orb, tShellResolvedScc, hubbU)

    !> input data to parse
    type(hsd_table), pointer, intent(in) :: node

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> is this a shell resolved calculation, or only one U value per atom
    logical, intent(in) :: tShellResolvedScc

    !> hubbard U values on exit
    real(dp), allocatable, intent(out) :: hubbU(:,:)

    type(hsd_table), pointer :: child, child2
    integer :: iSp1
    integer :: stat

    call hsd_rename_child(node, "CustomizedHubbards", "CustomisedHubbards")
    call hsd_get_table(node, "CustomisedHubbards", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      allocate(hubbU(orb%mShell, geo%nSpecies))
      hubbU(:,:) = 0.0_dp
      do iSp1 = 1, geo%nSpecies
        call hsd_get_table(child, geo%speciesNames(iSp1), child2, stat, auto_wrap=.true.)
        if (.not. associated(child2)) then
          cycle
        end if
        if (tShellResolvedScc) then
          block
            real(dp), allocatable :: tmpR(:)
            call hsd_get(child2, "#text", tmpR, stat=stat)
            if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required values")
            hubbU(:min(size(tmpR),orb%nShell(iSp1)), iSp1) = &
                & tmpR(:min(size(tmpR),orb%nShell(iSp1)))
          end block
        else
          call hsd_get(child2, "#text", hubbU(1, iSp1), stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value")
          hubbU(:orb%nShell(iSp1), iSp1) = hubbU(1, iSp1)
        end if
      end do

      ! -- Schema validation (warnings only) --
      block
        type(hsd_schema_t) :: schema
        type(hsd_error_t), allocatable :: schemaErrors(:)
        integer :: iErr

        call schema_init(schema, name="CustomisedHubbards")
        call schema_validate(schema, child, schemaErrors)
        if (size(schemaErrors) > 0) then
          do iErr = 1, size(schemaErrors)
            call dftbp_warning(child, "[schema] " // schemaErrors(iErr)%message)
          end do
        end if
        call schema_destroy(schema)
      end block
    end if

  end subroutine readCustomisedHubbards


  !> This subroutine overrides the neutral (reference) atom electronic occupation
  subroutine readCustomReferenceOcc(root, orb, referenceOcc, geo, iAtInRegion, customOcc)

    !> Node to be parsed
    type(hsd_table), pointer, intent(in) :: root

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Default reference occupations
    real(dp), intent(in) :: referenceOcc(:,:)

    !> Geometry information
    type(TGeometry), intent(in) :: geo

    !> Atom indices corresponding to user defined reference atomic charges
    type(TWrappedInt1), allocatable, intent(out) :: iAtInRegion(:)

    !> User-defined reference atomic charges
    real(dp), allocatable, intent(out) :: customOcc(:,:)

    type(hsd_table), pointer :: node, container, child
    type(hsd_table_ptr), allocatable :: nodes(:)
    character(len=:), allocatable :: buffer
    integer :: nCustomOcc, iCustomOcc, iShell, iSpecies, nAtom
    character(sc), allocatable :: shellNamesTmp(:)
    logical, allocatable :: atomOverriden(:)
    integer :: stat

    call hsd_rename_child(root, "CustomizedOccupations", "CustomisedOccupations")
    call hsd_get_table(root, "CustomisedOccupations", container, stat, auto_wrap=.true.)
    if (.not. associated(container)) then
      return
    end if

    call hsd_get_child_tables(container, "ReferenceOccupation", nodes)
    nCustomOcc = size(nodes)
    nAtom = size(geo%species)
    allocate(iAtInRegion(nCustomOcc))
    allocate(customOcc(orb%mShell, nCustomOcc))
    allocate(atomOverriden(nAtom))
    atomOverriden(:) = .false.
    customOcc(:,:) = 0.0_dp

    do iCustomOcc = 1, nCustomOcc
      node => nodes(iCustomOcc)%ptr
      call hsd_get(node, "Atoms", buffer, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'Atoms'")
      call hsd_get_table(node, "Atoms", child, stat)
      if (.not. associated(child)) child => node
      call getSelectedAtomIndices(child, buffer, geo%speciesNames, geo%species,&
          & iAtInRegion(iCustomOcc)%data)
      if (any(atomOverriden(iAtInRegion(iCustomOcc)%data))) then
        call dftbp_error(child, "Atom region contains atom(s) which have already been overridden")
      end if
      atomOverriden(iAtInRegion(iCustomOcc)%data) = .true.
      iSpecies = geo%species(iAtInRegion(iCustomOcc)%data(1))
      if (any(geo%species(iAtInRegion(iCustomOcc)%data) /= iSpecies)) then
        call dftbp_error(child, "All atoms in a ReferenceOccupation declaration must have the&
            & same type.")
      end if
      call getShellNames(iSpecies, orb, shellNamesTmp)
      do iShell = 1, orb%nShell(iSpecies)
          call hsd_get_or_set(node, shellNamesTmp(iShell), customOcc(iShell, iCustomOcc), &
            & referenceOcc(iShell, iSpecies))
      end do
      deallocate(shellNamesTmp)

      ! -- Schema validation for ReferenceOccupation (warnings only) --
      block
        type(hsd_schema_t) :: schema
        type(hsd_error_t), allocatable :: schemaErrors(:)
        integer :: iErr

        call schema_init(schema, name="ReferenceOccupation")
        call schema_add_field(schema, "Atoms", FIELD_REQUIRED, FIELD_TYPE_STRING)
        call schema_validate(schema, node, schemaErrors)
        if (size(schemaErrors) > 0) then
          do iErr = 1, size(schemaErrors)
            call dftbp_warning(node, "[schema] " // schemaErrors(iErr)%message)
          end do
        end if
        call schema_destroy(schema)
      end block
    end do

  end subroutine readCustomReferenceOcc

  function is_numeric(string) result(is)
    character(len=*), intent(in) :: string
    logical :: is

    real :: x
    integer :: err

    read(string,*,iostat=err) x
    is = (err == 0)
  end function is_numeric

end module dftbp_dftbplus_parser_customisation
