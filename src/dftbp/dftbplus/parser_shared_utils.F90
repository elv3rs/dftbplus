#:include 'common.fypp'
#:include 'error.fypp'

!> Shared utility routines used by multiple parser sub-modules.
module dftbp_dftbplus_parser_shared_utils
  use dftbp_common_accuracy, only : dp
  use dftbp_common_unitconversion, only : energyUnits, massUnits
  use dftbp_io_charmanip, only : i2c
  use dftbp_io_hsdcompat, only : hsd_table, hsd_child_list, detailedError, detailedWarning, &
      & getChild, getChildren, getChildValue, getSelectedAtomIndices, getLength, getItem1, &
      & destroyNodeList, convertUnitHsd
  use dftbp_type_typegeometry, only : TGeometry

  implicit none

  private
  public :: getInputMasses, readMDInitTemp

contains


  subroutine getInputMasses(node, geo, masses)

    !> relevant node of input data
    type(hsd_table), pointer :: node

    !> geometry object, which contains atomic species information
    type(TGeometry), intent(in) :: geo

    !> masses to be returned
    real(dp), allocatable, intent(out) :: masses(:)

    type(hsd_table), pointer :: child, child2, child3, val
    type(hsd_child_list), pointer :: children
    integer, allocatable :: pTmpI1(:)
    character(len=:), allocatable :: buffer, modifier
    real(dp) :: rTmp
    integer :: ii, jj, iAt

    call getChildValue(node, "Masses", val, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! Read individual atom specifications
    call getChildren(child, "Mass", children)
    if (getLength(children) == 0) then
      call destroyNodeList(children)
      return
    end if

    allocate(masses(geo%nAtom))
    masses(:) = -1.0_dp
    do ii = 1, getLength(children)
      call getItem1(children, ii, child2)
      call getChildValue(child2, "Atoms", buffer, child=child3, multiple=.true.)
      call getSelectedAtomIndices(child3, buffer, geo%speciesNames, geo%species, pTmpI1)
      call getChildValue(child2, "MassPerAtom", rTmp, modifier=modifier, child=child)
      call convertUnitHsd(modifier, massUnits, child, rTmp)
      do jj = 1, size(pTmpI1)
        iAt = pTmpI1(jj)
        if (masses(iAt) >= 0.0_dp) then
          call detailedWarning(child3, "Previous setting for the mass  of atom" // i2c(iAt) //&
              & " overwritten")
        end if
        masses(iAt) = rTmp
      end do
      deallocate(pTmpI1)
    end do
    call destroyNodeList(children)

  end subroutine getInputMasses


  subroutine readMDInitTemp(node, tempAtom, minimumTemp)

    !> input data to parse
    type(hsd_table), pointer :: node

    !> Ionic temperature
    real(dp), intent(out) :: tempAtom

    !> Lowest possible ion temperature
    real(dp), intent(in) :: minimumTemp

    type(hsd_table), pointer :: child
    character(len=:), allocatable :: modifier

    call getChildValue(node, "InitialTemperature", tempAtom, modifier=modifier, child=child)
    if (tempAtom < 0.0_dp) then
      call detailedError(node, "Negative temperature")
    end if
    call convertUnitHsd(modifier, energyUnits, node, tempAtom)
    tempAtom = max(tempAtom, minimumTemp)

  end subroutine readMDInitTemp


end module dftbp_dftbplus_parser_shared_utils
