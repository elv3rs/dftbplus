!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Module to read input from HSD tree
module dftbp_dftbplus_input_geoopt
  use dftbp_common_accuracy, only : dp
  use dftbp_common_unitconversion, only : energyUnits, forceUnits, lengthUnits, timeUnits
  use dftbp_geoopt_package, only : TFilterInput, TFireInput, TLbfgsInput, TOptimizerInput,&
      & TOptTolerance, TRationalFuncInput, TSteepdescInput
  use dftbp_io_charmanip, only : unquote
  use hsd, only : hsd_rename_child, hsd_get_or_set, hsd_get_table, hsd_get_attrib, hsd_get_choice,&
      & HSD_STAT_OK, hsd_get_name
  use dftbp_io_hsdutils, only : getSelectedAtomIndices, dftbp_error
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_type_typegeometry, only : TGeometry
  use hsd_data, only : hsd_table, new_table
  implicit none

  private
  public :: readGeoOptInput, readOptimizerInput, TGeoOptInput


  !> General input wrapper for optimisers in this package
  type :: TGeoOptInput

    !> Input for coordinate transformation and filter step
    type(TFilterInput) :: filter

    !> Optimiser input choice
    class(TOptimizerInput), allocatable :: optimiser

    !> Tolerances for optimization
    type(TOptTolerance) :: tolerance

    !> Number of allowed geometry optimization steps
    integer :: nGeoSteps

    !> Prefix of the output file name
    character(len=:), allocatable :: outFile

  end type TGeoOptInput


contains

  !> General entry point to read geometry optimization
  subroutine readGeoOptInput(node, geom, input, atomsRange)

    !> Node to get the information from
    type(hsd_table), pointer, intent(in) :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    !> Control structure to be filled
    type(TGeoOptInput), intent(out) :: input

    !> Default range of moving atoms (may be restricted for example by contacts in transport
    !> calculations)
    character(len=*), intent(in) :: atomsRange

    type(hsd_table), pointer :: child, value1
    character(len=:), allocatable :: buffer
    integer :: stat

    call hsd_rename_child(node, "Optimizer", "Optimiser")
    call hsd_get_choice(node, "Optimiser", buffer, child, stat)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: defContainer, defChild
        call new_table(defContainer, name="optimiser")
        call new_table(defChild, name="Rational")
        call defContainer%add_child(defChild)
        call node%add_child(defContainer)
      end block
      call hsd_get_choice(node, "Optimiser", buffer, child, stat)
    end if
    call readOptimizerInput(child, input%optimiser)

    call readFilterInput(node, geom, input%filter, atomsRange)

    call hsd_get_table(node, "Convergence", child, stat, auto_wrap=.true.)
    if (.not.associated(child)) then
      block
        type(hsd_table) :: tmpTbl
        call new_table(tmpTbl, name="convergence")
        call node%add_child(tmpTbl)
      end block
      call hsd_get_table(node, "Convergence", child, stat)
    end if
    call readOptTolerance(child, input%tolerance)

    call hsd_get_or_set(node, "MaxSteps", input%nGeoSteps, 20*geom%nAtom)

    call hsd_get_or_set(node, "OutputPrefix", buffer, "geo_end")
    input%outFile = trim(unquote(buffer))

  end subroutine readGeoOptInput


  !> Reads the optimiser
  subroutine readOptimizerInput(node, input)

    !> Optimiser node
    type(hsd_table), pointer, intent(in) :: node

    !> Control structure to be filled
    class(TOptimizerInput), allocatable, intent(out) :: input

    type(TSteepDescInput), allocatable :: steepDescInput
    type(TFireInput), allocatable :: fireInput
    type(TLbfgsInput), allocatable :: lbfgsInput
    type(TRationalFuncInput), allocatable :: rationalFuncInput
    character(len=:), allocatable :: buffer

    call hsd_get_name(node, buffer, "#text")
    select case (buffer)
    case default
      call dftbp_error(node, "Invalid optimiser name.")
    case("steepestdescent")
        allocate(steepDescInput)
        call readSteepDescInput(node, steepDescInput)
        call move_alloc(steepDescInput, input)
    case("fire")
        allocate(fireInput)
        call readFireInput(node, fireInput)
        call move_alloc(fireInput, input)
    case("lbfgs")
      allocate(lbfgsInput)
      call readLbfgsInput(node, lbfgsInput)
      call move_alloc(lbfgsInput, input)
    case("rational")
      allocate(rationalFuncInput)
      call readRationalFuncInput(node, rationalFuncInput)
      call move_alloc(rationalFuncInput, input)
    end select

  end subroutine readOptimizerInput


  !> Entry point for reading input for cartesian geometry transformation filter
  subroutine readFilterInput(node, geom, input, atomsRange)

    !> Node to get the information from
    type(hsd_table), pointer, intent(in) :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    !> Control structure to be filled
    type(TFilterInput), intent(out) :: input

    !> Default range of moving atoms (may be restricted for example by contacts in transport
    !> calculations)
    character(len=*), intent(in) :: atomsRange

    type(hsd_table), pointer :: child
    character(len=:), allocatable :: buffer
    integer :: stat

    call hsd_get_or_set(node, "LatticeOpt", input%lattice, .false.)
    if (input%lattice) then
      call hsd_get_or_set(node, "FixAngles", input%fixAngles, .false.)
      block
        logical, allocatable :: tmpFixLen(:)
        call hsd_get_or_set(node, "FixLengths", tmpFixLen, [.false., .false., .false.])
        input%fixLength(:) = tmpFixLen(:min(size(tmpFixLen), size(input%fixLength)))
      end block
      if (input%fixAngles .and. all(input%fixLength)) then
        call dftbp_error(node, "LatticeOpt with all lattice vectors fixed is not possible")
      end if
      call hsd_get_or_set(node, "Isotropic", input%isotropic, .false.)
    end if
    call hsd_get_or_set(node, "MovedAtoms", buffer, trim(atomsRange))
    call hsd_get_table(node, "MovedAtoms", child, stat, auto_wrap=.true.)
    call getSelectedAtomIndices(child, buffer, geom%speciesNames, geom%species, &
        & input%indMovedAtom)

  end subroutine readFilterInput


  !> Entry point for reading convergence thresholds
  subroutine readOptTolerance(node, input)

    !> Node to get the information from
    type(hsd_table), pointer, intent(in) :: node

    !> Control structure to be filled
    type(TOptTolerance), intent(out) :: input

    type(hsd_table), pointer :: field
    character(len=:), allocatable :: modifier
    integer :: stat

    call hsd_get_or_set(node, "Energy", input%energy, huge(1.0_dp))
    call hsd_get_attrib(node, "Energy", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    call hsd_get_table(node, "Energy", field, stat, auto_wrap=.true.)
    call convertUnitHsd(modifier, energyUnits, field, input%energy)

    call hsd_get_or_set(node, "GradNorm", input%gradNorm, huge(1.0_dp))
    call hsd_get_attrib(node, "GradNorm", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    call hsd_get_table(node, "GradNorm", field, stat, auto_wrap=.true.)
    call convertUnitHsd(modifier, forceUnits, field, input%gradNorm)

    call hsd_get_or_set(node, "GradElem", input%gradElem, 1.0e-4_dp)
    call hsd_get_attrib(node, "GradElem", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    call hsd_get_table(node, "GradElem", field, stat, auto_wrap=.true.)
    call convertUnitHsd(modifier, forceUnits, field, input%gradElem)

    call hsd_get_or_set(node, "DispNorm", input%dispNorm, huge(1.0_dp))
    call hsd_get_attrib(node, "DispNorm", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    call hsd_get_table(node, "DispNorm", field, stat, auto_wrap=.true.)
    call convertUnitHsd(modifier, lengthUnits, field, input%dispNorm)

    call hsd_get_or_set(node, "DispElem", input%dispElem, huge(1.0_dp))
    call hsd_get_attrib(node, "DispElem", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    call hsd_get_table(node, "DispElem", field, stat, auto_wrap=.true.)
    call convertUnitHsd(modifier, lengthUnits, field, input%dispElem)

  end subroutine readOptTolerance


  !> Entry point for reading input for SteepestDescent
  subroutine readSteepDescInput(node, input)

    !> Node to get the information from
    type(hsd_table), intent(in), pointer :: node

    !> Control structure to be filled
    type(TSteepDescInput), intent(out) :: input

    call hsd_get_or_set(node, "ScalingFactor", input%scalingFactor, 1.0_dp)

  end subroutine readSteepDescInput


  !> Entry point for reading input for FIRE
  subroutine readFireInput(node, input)

    !> Node to get the information from
    type(hsd_table), pointer, intent(in) :: node

    !> Control structure to be filled
    type(TFireInput), intent(out) :: input

    type(hsd_table), pointer :: field
    character(len=:), allocatable :: modifier
    integer :: stat

    call hsd_get_or_set(node, "nMin", input%nMin, 5)
    call hsd_get_or_set(node, "aPar", input%a_start, 0.1_dp)
    call hsd_get_or_set(node, "fInc", input%f_inc, 1.1_dp)
    call hsd_get_or_set(node, "fDec", input%f_dec, 0.5_dp)
    call hsd_get_or_set(node, "fAlpha", input%f_alpha, 0.99_dp)
    call hsd_get_or_set(node, "StepSize", input%dt_max, 1.0_dp)
    call hsd_get_attrib(node, "StepSize", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    call hsd_get_table(node, "StepSize", field, stat, auto_wrap=.true.)
    call convertUnitHsd(modifier, timeUnits, field, input%dt_max)

  end subroutine readFireInput


  !> Entry point for reading input for LBFGS
  subroutine readLbfgsInput(node, input)

    !> Node to get the information from
    type(hsd_table), pointer, intent(in) :: node

    !> Control structure to be filled
    type(TLbfgsInput), intent(out) :: input

    call hsd_get_or_set(node, "Memory", input%memory, 20)

  end subroutine readLbfgsInput


  !> Entry point for reading input for rational function optimiser
  subroutine readRationalFuncInput(node, input)

    !> Node to get the information from
    type(hsd_table), pointer, intent(in) :: node

    !> Control structure to be filled
    type(TRationalFuncInput), intent(out) :: input

    call hsd_get_or_set(node, "DiagLimit", input%diagLimit, 1.0e-2_dp)

  end subroutine readRationalFuncInput


end module dftbp_dftbplus_input_geoopt
