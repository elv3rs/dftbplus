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
  use dftbp_io_hsdutils, only : &
      & getChild, getChildValue
  use dftbp_io_hsdutils, only : dftbp_error, getSelectedAtomIndices
  use hsd, only : hsd_rename_child, hsd_table_ptr, hsd_get_child_tables
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

    call hsd_rename_child(node, "CustomizedHubbards", "CustomisedHubbards")
    call getChild(node, "CustomisedHubbards", child, requested=.false.)
    if (associated(child)) then
      allocate(hubbU(orb%mShell, geo%nSpecies))
      hubbU(:,:) = 0.0_dp
      do iSp1 = 1, geo%nSpecies
        call getChild(child, geo%speciesNames(iSp1), child2, requested=.false.)
        if (.not. associated(child2)) then
          cycle
        end if
        if (tShellResolvedScc) then
          call getChildValue(child2, "", hubbU(:orb%nShell(iSp1), iSp1))
        else
          call getChildValue(child2, "", hubbU(1, iSp1))
          hubbU(:orb%nShell(iSp1), iSp1) = hubbU(1, iSp1)
        end if
      end do
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

    call hsd_rename_child(root, "CustomizedOccupations", "CustomisedOccupations")
    call getChild(root, "CustomisedOccupations", container, requested=.false.)
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
      call getChildValue(node, "Atoms", buffer, child=child, multiple=.true.)
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
          call getChildValue(node, shellNamesTmp(iShell), customOcc(iShell, iCustomOcc), &
            & default=referenceOcc(iShell, iSpecies))
      end do
      deallocate(shellNamesTmp)
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
