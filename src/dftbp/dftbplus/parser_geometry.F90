#:include 'common.fypp'
#:include 'error.fypp'

!> Parser routine for the Geometry block.
module dftbp_dftbplus_parser_geometry
  use dftbp_dftbplus_inputdata, only : TInputData
  use hsd, only : hsd_get_choice
  use hsd_data, only : hsd_table
  use dftbp_type_typegeometryhsd, only : readTGeometryGen, readTGeometryHsd, readTGeometryLammps,&
      & readTGeometryVasp, readTGeometryXyz
  implicit none

  private
  public :: readGeometry

contains


  !> Read in Geometry
  subroutine readGeometry(node, input)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Input structure to be filled
    type(TInputData), intent(inout) :: input

    type(hsd_table), pointer :: value1
    character(len=:), allocatable :: buffer
    integer :: stat

    call hsd_get_choice(node, "", buffer, value1, stat)
    input%geom%tPeriodic = .false.
    input%geom%tHelical = .false.
    select case (buffer)
    case ("genformat")
      call readTGeometryGen(value1, input%geom)
    case ("xyzformat")
      call readTGeometryXyz(value1, input%geom)
    case ("vaspformat")
      call readTGeometryVasp(value1, input%geom)
    case ("lammpsformat")
      call readTGeometryLammps(value1, input%geom)
    case default
      if (associated(value1)) value1%processed = .false.
      call readTGeometryHSD(node, input%geom)
    end select

  end subroutine readGeometry

end module dftbp_dftbplus_parser_geometry
