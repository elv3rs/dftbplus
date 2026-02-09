#:include 'common.fypp'
#:include 'error.fypp'

!> Shared utility routines used by multiple parser sub-modules.
module dftbp_dftbplus_parser_shared_utils
  use dftbp_common_accuracy, only : dp
  use dftbp_common_unitconversion, only : energyUnits, massUnits
  use dftbp_io_charmanip, only : i2c
  use dftbp_io_hsdutils, only : getChildValue
  use hsd, only : hsd_table_ptr, hsd_get_child_tables
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, getSelectedAtomIndices
  use dftbp_io_unitconv, only : convertUnitHsd
  use hsd_data, only : hsd_table
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
    type(hsd_table_ptr), allocatable :: children(:)
    integer, allocatable :: pTmpI1(:)
    character(len=:), allocatable :: buffer, modifier
    real(dp) :: rTmp
    integer :: ii, jj, iAt

    call getChildValue(node, "Masses", val, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! Read individual atom specifications
    call hsd_get_child_tables(child, "Mass", children)
    if (size(children) == 0) then
      return
    end if

    allocate(masses(geo%nAtom))
    masses(:) = -1.0_dp
    do ii = 1, size(children)
      child2 => children(ii)%ptr
      call getChildValue(child2, "Atoms", buffer, child=child3, multiple=.true.)
      call getSelectedAtomIndices(child3, buffer, geo%speciesNames, geo%species, pTmpI1)
      call getChildValue(child2, "MassPerAtom", rTmp, modifier=modifier, child=child)
      call convertUnitHsd(modifier, massUnits, child, rTmp)
      do jj = 1, size(pTmpI1)
        iAt = pTmpI1(jj)
        if (masses(iAt) >= 0.0_dp) then
          call dftbp_warning(child3, "Previous setting for the mass  of atom" // i2c(iAt) //&
              & " overwritten")
        end if
        masses(iAt) = rTmp
      end do
      deallocate(pTmpI1)
    end do

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
      call dftbp_error(node, "Negative temperature")
    end if
    call convertUnitHsd(modifier, energyUnits, node, tempAtom)
    tempAtom = max(tempAtom, minimumTemp)

  end subroutine readMDInitTemp


end module dftbp_dftbplus_parser_shared_utils
