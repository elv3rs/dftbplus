!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines to deal with HSD species lists
module dftbp_dftbplus_specieslist
  use dftbp_common_accuracy, only : dp
  use dftbp_common_unitconversion, only : TUnit
  use hsd, only : hsd_get, hsd_get_or_set, hsd_get_attrib, hsd_get_table, HSD_STAT_OK
  use dftbp_io_hsdutils, only : dftbp_error
  use dftbp_io_unitconv, only : convertUnitHsd
  use hsd_data, only : hsd_table, hsd_node
  implicit none

  private
  public :: readSpeciesList


  !> Generic wrapper for species specific data
  interface readSpeciesList
    module procedure :: readSpeciesListReal
    module procedure :: readSpeciesListInt
  end interface readSpeciesList


contains


  !> Mark a node (and optionally all descendants) as processed. Null-safe.
  recursive subroutine setProcessed(node, recursive)
    type(hsd_table), pointer, intent(in) :: node
    logical, intent(in), optional :: recursive

    logical :: doRecurse
    integer :: ii
    class(hsd_node), pointer :: childNode

    if (.not. associated(node)) return
    node%processed = .true.

    doRecurse = .false.
    if (present(recursive)) doRecurse = recursive
    if (.not. doRecurse) return

    do ii = 1, node%num_children
      call node%get_child(ii, childNode)
      if (.not. associated(childNode)) cycle
      childNode%processed = .true.
      select type (t => childNode)
      type is (hsd_table)
        block
          type(hsd_table), pointer :: tPtr
          tPtr => t
          call setProcessed(tPtr, .true.)
        end block
      end select
    end do

  end subroutine setProcessed


  !> Read a list of real valued species data
  subroutine readSpeciesListReal(node, speciesNames, array, default, units, markAllProcessed)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Names of all species
    character(len=*), intent(in) :: speciesNames(:)

    !> Data array to read
    real(dp), intent(inout) :: array(:)

    !> Optional default values of data array to be read
    real(dp), intent(in), optional :: default(:)

    !> Conversion factor
    type(TUnit), intent(in), optional :: units(:)

    !> Whether all unread subnodes should also be marked as processed (default: .true.)
    logical, optional, intent(in) :: markAllProcessed

    type(hsd_table), pointer :: child
    character(len=:), allocatable :: modifier
    integer :: iSp, stat
    logical :: markAllProcessed_

    if (present(default)) then
      if (present(units)) then
        do iSp = 1, size(speciesNames)
          call hsd_get_or_set(node, speciesNames(iSp), array(iSp), default(iSp))
          call hsd_get_attrib(node, speciesNames(iSp), modifier, stat)
          if (stat /= HSD_STAT_OK) modifier = ""
          call hsd_get_table(node, speciesNames(iSp), child, stat, auto_wrap=.true.)
          call convertUnitHsd(modifier, units, child, array(iSp))
        end do
      else
        do iSp = 1, size(speciesNames)
          call hsd_get_or_set(node, speciesNames(iSp), array(iSp), default(iSp))
        end do
      end if
    else
      if (present(units)) then
        do iSp = 1, size(speciesNames)
          call hsd_get(node, speciesNames(iSp), array(iSp), stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: '" &
              & // trim(speciesNames(iSp)) // "'")
          call hsd_get_attrib(node, speciesNames(iSp), modifier, stat)
          if (stat /= HSD_STAT_OK) modifier = ""
          call hsd_get_table(node, speciesNames(iSp), child, stat, auto_wrap=.true.)
          call convertUnitHsd(modifier, units, child, array(iSp))
        end do
      else
        do iSp = 1, size(speciesNames)
          call hsd_get(node, speciesNames(iSp), array(iSp), stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: '" &
              & // trim(speciesNames(iSp)) // "'")
        end do
      end if
    end if

    markAllProcessed_ = .true.
    if (present(markAllProcessed)) markAllProcessed_ = markAllProcessed
    if (markAllProcessed_) call setProcessed(node, recursive=.true.)

  end subroutine readSpeciesListReal


  !> Read a list of integer valued species data
  subroutine readSpeciesListInt(node, speciesNames, array, default, markAllProcessed)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Names of all species
    character(len=*), intent(in) :: speciesNames(:)

    !> Data array to read
    integer, intent(out) :: array(:)

    !> Data array to read
    integer, optional, intent(in) :: default(:)

    !> Whether all unread subnodes should also be marked as processed (default: .true.)
    logical, optional, intent(in) :: markAllProcessed

    integer :: iSp, stat
    logical :: markAllProcessed_

    if (present(default)) then
      do iSp = 1, size(speciesNames)
        call hsd_get_or_set(node, speciesNames(iSp), array(iSp), default(iSp))
      end do
    else
      do iSp = 1, size(speciesNames)
        call hsd_get(node, speciesNames(iSp), array(iSp), stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: '" &
            & // trim(speciesNames(iSp)) // "'")
      end do
    end if

    markAllProcessed_ = .true.
    if (present(markAllProcessed)) markAllProcessed_ = markAllProcessed
    if (markAllProcessed_) call setProcessed(node, recursive=.true.)

  end subroutine readSpeciesListInt

end module dftbp_dftbplus_specieslist
