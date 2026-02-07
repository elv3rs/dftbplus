#:include 'common.fypp'
#:include 'error.fypp'

!> Parser routines for Driver/MD and electron dynamics input blocks.
module dftbp_dftbplus_parser_driver
  use dftbp_common_accuracy, only : dp, lc, mc, minTemp, sc
  use dftbp_common_constants, only : Bohr__AA
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_unitconversion, only : angularUnits, EFieldUnits, energyUnits, forceUnits, &
      & freqUnits, lengthUnits, massUnits, pressureUnits, timeUnits, VelocityUnits
  use dftbp_dftbplus_input_geoopt, only : readGeoOptInput
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_dftbplus_parser_shared_utils, only : getInputMasses, readMDInitTemp
  use dftbp_elecsolvers_elecsolvers, only : electronicSolverTypes
  use dftbp_extlibs_plumed, only : withPlumed
  use dftbp_geoopt_geoopt, only : geoOptTypes
  use dftbp_io_charmanip, only : i2c, unquote
  use dftbp_io_hsdcompat, only : hsd_table, hsd_child_list, textNodeName, &
      & detailedError, detailedWarning, getChild, getChildren, getChildValue, &
      & getSelectedAtomIndices, getNodeName, getNodeHSDName, getNodeName2, &
      & convertUnitHsd, hsd_rename_child
  use dftbp_io_message, only : error
  use dftbp_md_tempprofile, only : identifyTempProfile, tempProfileTypes, TTempProfileInput
  use dftbp_md_thermostats, only : thermostatTypes, TThermostatInput
  use dftbp_md_xlbomd, only : TXlbomdInp
  use dftbp_timedep_timeprop, only : envTypes, pertTypes, tdSpinTypes, TElecDynamicsInp
  use dftbp_transport_negfvars, only : TTransPar
  use dftbp_type_linkedlist, only : asArray, asVector, destruct, init, len, &
      & TListIntR1, TListRealR1, TListString
  use dftbp_type_typegeometry, only : TGeometry
#:if WITH_SOCKETS
  use dftbp_io_ipisocket, only : IPI_PROTOCOLS
#:endif

  implicit none

  private
  public :: readDriver, readElecDynamics

contains

  !> Read in driver properties
#:if WITH_TRANSPORT
  subroutine readDriver(node, parent, geom, ctrl, transpar)
#:else
  subroutine readDriver(node, parent, geom, ctrl)
#:endif

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Parent of node (for error messages)
    type(hsd_table), pointer :: parent

    !> geometry of the system
    type(TGeometry), intent(in) :: geom

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(in) :: transpar
  #:endif

    type(hsd_table), pointer :: child, child2, child3, value1, value2, field

    character(len=:), allocatable :: buffer, buffer2, modifier
  #:if WITH_SOCKETS
    character(lc) :: sTmp
  #:endif

    ! Default range of atoms to move (may be adjusted if contacts present)
    character(mc) :: atomsRange

    character(mc) :: modeName
    logical :: isMaxStepNeeded

    ctrl%isGeoOpt = .false.
    ctrl%tCoordOpt = .false.
    ctrl%tLatOpt = .false.

    ctrl%iGeoOpt = geoOptTypes%none
    ctrl%tMD = .false.
    ctrl%tForces = .false.
    ctrl%tSetFillingTemp = .false.

    atomsRange = "1:-1"
  #:if WITH_TRANSPORT
    if (transpar%defined) then
      ! only those atoms in the device region
      write(atomsRange,"(I0,':',I0)")transpar%idxdevice
    end if
  #:endif

    call hsd_rename_child(parent, "GeometryOptimization", "GeometryOptimisation")
    call getNodeName2(node, buffer)
    driver: select case (buffer)
    case ("")
      modeName = ""
      continue
    case ("none")
      modeName = ""
      continue
    case ("geometryoptimisation")
      modeName = "geometry optimisation"

      if (geom%tHelical) then
        call detailedError(node, "GeometryOptimisation driver currently does not support helical&
            & geometries")
      end if

      allocate(ctrl%geoOpt)

      call readGeoOptInput(node, geom, ctrl%geoOpt, atomsRange)

      call getChildValue(node, "AppendGeometries", ctrl%tAppendGeo, .false.)

      ! Geometry optimisation drivers
      ctrl%iGeoOpt = geoOptTypes%geometryoptimisation
      ctrl%tForces = .true.
      ctrl%restartFreq = 1

    case ("steepestdescent")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ! Steepest downhill optimisation
      ctrl%iGeoOpt = geoOptTypes%steepestDesc
      call commonGeoOptions(node, ctrl, geom, atomsRange)

    case ("conjugategradient")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')// "Please use the GeometryOptimisation driver instead.")

      ! Conjugate gradient location optimisation
      ctrl%iGeoOpt = geoOptTypes%conjugateGrad
      call commonGeoOptions(node, ctrl, geom, atomsRange)

    case("gdiis")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ! Gradient DIIS optimisation, only stable in the quadratic region
      ctrl%iGeoOpt = geoOptTypes%diis
      call getChildValue(node, "alpha", ctrl%deltaGeoOpt, 1.0E-1_dp)
      call getChildValue(node, "Generations", ctrl%iGenGeoOpt, 8)
      call commonGeoOptions(node, ctrl, geom, atomsRange)

    case ("lbfgs")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ctrl%iGeoOpt = geoOptTypes%lbfgs

      allocate(ctrl%lbfgsInp)
      call getChildValue(node, "Memory", ctrl%lbfgsInp%memory, 20)

      call getChildValue(node, "LineSearch", ctrl%lbfgsInp%isLineSearch, .false.)

      isMaxStepNeeded = .not. ctrl%lbfgsInp%isLineSearch
      if (isMaxStepNeeded) then
        call getChildValue(node, "setMaxStep", ctrl%lbfgsInp%isLineSearch, isMaxStepNeeded)
        ctrl%lbfgsInp%MaxQNStep = isMaxStepNeeded
      else
        call getChildValue(node, "oldLineSearch", ctrl%lbfgsInp%isOldLS, .false.)
      end if

      call commonGeoOptions(node, ctrl, geom, atomsRange, ctrl%lbfgsInp%isLineSearch)

    case ("fire")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ctrl%iGeoOpt = geoOptTypes%fire
      call commonGeoOptions(node, ctrl, geom, atomsRange, .false.)
      call getChildValue(node, "TimeStep", ctrl%deltaT, 1.0_dp, modifier=modifier, child=field)
      call convertUnitHsd(modifier, timeUnits, field, ctrl%deltaT)

    case("secondderivatives")
      ! currently only numerical derivatives of forces is implemented

      modeName = "second derivatives"

      ctrl%tDerivs = .true.
      ctrl%tForces = .true.

      call getChildValue(node, "Atoms", buffer2, trim(atomsRange), child=child,&
          & multiple=.true.)
      call getSelectedAtomIndices(child, buffer2, geom%speciesNames, geom%species,&
          & ctrl%indDerivAtom)
      if (size(ctrl%indDerivAtom) == 0) then
        call error("No atoms specified for derivatives calculation.")
      end if

      call getChild(node, "MovedAtoms", child, requested=.false.)
      if (associated(child)) then
        if (.not. isContiguousRange(ctrl%indDerivAtom)) then
          call detailedError(child,&
            & "Atoms for calculation of partial Hessian must be a contiguous range.")
        end if
        call getChildValue(child, "", buffer2, child=child2, multiple=.true.)
        call getSelectedAtomIndices(child2, buffer2, geom%speciesNames, geom%species, &
           & ctrl%indMovedAtom)
        if (.not. isContiguousRange(ctrl%indMovedAtom)) then
          call detailedError(child2, "MovedAtoms for calculation of partial Hessian must be a &
              & contiguous range.")
        end if
        if (.not. containsAll(ctrl%indDerivAtom, ctrl%indMovedAtom)) then
          call detailedError(child2, "MovedAtoms has indices not contained in Atoms.")
        end if
      else
        ctrl%indMovedAtom = ctrl%indDerivAtom
      end if
      ctrl%nrMoved = size(ctrl%indMovedAtom)

      call getChildValue(node, "Delta", ctrl%deriv2ndDelta, 1.0E-4_dp, &
          & modifier=modifier, child=field)
      call convertUnitHsd(modifier, lengthUnits, field, ctrl%deriv2ndDelta)

    case ("velocityverlet")
      ! molecular dynamics

      modeName = "molecular dynamics"

      ctrl%tForces = .true.
      ctrl%tMD = .true.

      call getChildValue(node, "MDRestartFrequency", ctrl%restartFreq, 1)
      call getChildValue(node, "MovedAtoms", buffer2, trim(atomsRange), child=child, &
          &multiple=.true.)
      call getSelectedAtomIndices(child, buffer2, geom%speciesNames, geom%species, &
          & ctrl%indMovedAtom)
      ctrl%nrMoved = size(ctrl%indMovedAtom)
      if (ctrl%nrMoved == 0) then
        call error("No atoms specified for molecular dynamics.")
      end if
      call readInitialVelocities(node, ctrl, geom%nAtom)

      call getChildValue(node, "KeepStationary", ctrl%tMDstill,.true.)
      if (ctrl%tMDstill .and. geom%nAtom == 1) then
        call error("Removing translational freedom with only one atom not&
            & possible.")
      end if

      call getChildValue(node, "TimeStep", ctrl%deltaT, modifier=modifier, &
          &child=field)
      call convertUnitHsd(modifier, timeUnits, field, ctrl%deltaT)

      call parseThermostat(node, ctrl%deltaT, ctrl%tReadMDVelocities, ctrl%maxRun,&
          & ctrl%thermostatInp, ctrl%tempProfileInp)

      if (ctrl%maxRun < -1) then
        call getChildValue(node, "Steps", ctrl%maxRun)
      end if

      call getChildValue(node, "OutputPrefix", buffer2, "geo_end")
      ctrl%outFile = unquote(buffer2)

      call getChildValue(node, "Plumed", ctrl%tPlumed, default=.false., child=child)
      if (ctrl%tPlumed .and. .not. withPlumed) then
        call detailedError(child, "Metadynamics can not be used since code has been compiled&
            & without PLUMED support")
      end if

      if (geom%tPeriodic) then

        ctrl%tBarostat = .false.
        ctrl%pressure = 0.0_dp
        ctrl%BarostatStrength = 0.0_dp
        call getChild(node, "Barostat", child, requested=.false.)
        if (associated(child)) then
          if (allocated(ctrl%hybridXcInp)) then
            call detailedError(node, "Barostating not currently implemented for hybrid functionals")
          end if
          if (ctrl%nrMoved /= geom%nAtom) then
            call error("Dynamics for a subset of atoms is not currently&
                & possible when using a barostat")
          end if
          call getChildValue(child, "Pressure", ctrl%pressure, &
              & modifier=modifier, child=child2)
          call convertUnitHsd(modifier, pressureUnits, child2, &
              & ctrl%pressure)
          call getChild(child, "Coupling", child=child2, requested=.false.)
          if (associated(child2)) then
            call getChildValue(child2, "", ctrl%BarostatStrength)
            call getChild(child, "Timescale",child=child2,modifier=modifier,&
                &requested=.false.)
            if (associated(child2)) call error("Only Coupling strength OR &
                &Timescale can be set for Barostatting.")
          else
            call getChild(child, "Timescale",child=child2,modifier=modifier,&
                &requested=.false.)
            if (associated(child2)) then
              call getChildValue(child2, "", ctrl%BarostatStrength, &
                  & modifier=modifier, child=child3)
              call convertUnitHsd(modifier, timeUnits, child3, &
                  & ctrl%BarostatStrength)
              ctrl%BarostatStrength = ctrl%deltaT / ctrl%BarostatStrength
            else
              call error("Either Coupling strength or Timescale must be set&
                  & for Barostatting.")
            end if
          end if
          call getChildValue(child, "Isotropic", ctrl%tIsotropic, .true.)
          ctrl%tBarostat = .true.
        end if
      end if

      if (ctrl%hamiltonian == hamiltonianTypes%dftb) then
        call readXlbomdOptions(node, ctrl%xlbomd)
      end if

      call getInputMasses(node, geom, ctrl%masses)

    case ("socket")
      ! external socket control of the run (once initialised from input)

      modeName = "socket control"

    #:if WITH_SOCKETS
      ctrl%tForces = .true.
      allocate(ctrl%socketInput)
      call getChild(node, 'File', child=child2, requested=.false.)
      call getChild(node, 'Host', child=child3, requested=.false.)
      if (associated(child2) .eqv. associated(child3)) then
        call error('Either Host or File (but not both) must be set for socket&
            & communication')
      end if

      ! File communication
      if (associated(child2)) then
        call getChildValue(child2, "", buffer2)
        ctrl%socketInput%host = unquote(buffer2)
        ! zero it to signal to initprogram to use a unix file
        ctrl%socketInput%port = 0
      else
        call getChildValue(child3, "", buffer2)
        ctrl%socketInput%host = unquote(buffer2)
        call getChildValue(node, "Port", ctrl%socketInput%port, child=field)
        if (ctrl%socketInput%port <= 0) then
          call detailedError(field, "Invalid port number")
        end if
      end if

      call getChildValue(node, "Protocol", value1, "i-PI", child=child)
      call getNodeName(value1, buffer)
      select case(buffer)
      case("i-pi")
        ctrl%socketInput%protocol = IPI_PROTOCOLS%IPI_1
        ! want a file path
        if (ctrl%socketInput%port == 0) then
          call getChildValue(node, "Prefix", buffer2, "/tmp/ipi_")
          sTmp = unquote(buffer2)
          ctrl%socketInput%host = trim(sTmp) // trim(ctrl%socketInput%host)
        end if

      case default
        call detailedError(child, "Invalid protocol '" // buffer // "'")
      end select
      call getChildValue(node, "Verbosity", ctrl%socketInput%verbosity, 0)
      call getChildValue(node, "MaxSteps", ctrl%maxRun, 200)

    #:else
      call detailedError(node, "Program had been compiled without socket support")
    #:endif

    case default

      call getNodeHSDName(node, buffer)
      call detailedError(parent, "Invalid driver '" // buffer // "'")

    end select driver

  #:if WITH_TRANSPORT
    if (ctrl%solver%isolver == electronicSolverTypes%OnlyTransport .and. trim(modeName) /= "") then
      call detailederror(node, "transportOnly solver cannot be used with "//trim(modeName))
    end if
  #:endif

  end subroutine readDriver

  !> Simple function to check that an array of indices is a contigous range
  function isContiguousRange(indices) result(isContiguous)

    !> Array of atomic indices
    integer, intent(in) :: indices(:)

    !> whether indices are contigous
    logical :: isContiguous

    isContiguous = all(indices(: size(indices) - 1) + 1 == indices(2:))

  end function isContiguousRange


  !> checks that the array subindices is contained in indices
  function containsAll(indices, subindices)

    !> Array of atomic indices to check against
    integer, intent(in) :: indices(:)

    !> Array of atomic indices to check
    integer, intent(in) :: subindices(:)

    !> whether indices are contigous
    logical :: containsAll

    integer :: kk

    containsAll = .false.
    do kk = 1, size(subindices)
      if (.not. any(indices == subindices(kk))) return
    end do
    containsAll = .true.

  end function containsAll


  !> Common geometry optimisation settings for various drivers
  subroutine commonGeoOptions(node, ctrl, geom, atomsRange, isMaxStepNeeded)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> geometry of the system
    type(TGeometry), intent(in) :: geom

    !> Default range of moving atoms (may be restricted for example by contacts in transport
    !> calculations)
    character(len=*), intent(in) :: atomsRange

    !> Is the maximum step size relevant for this driver
    logical, intent(in), optional :: isMaxStepNeeded

    type(hsd_table), pointer :: child, field
    character(len=:), allocatable :: buffer2, modifier
    logical :: isMaxStep

    if (present(isMaxStepNeeded)) then
      isMaxStep = isMaxStepNeeded
    else
      isMaxStep = .true.
    end if

    ctrl%tForces = .true.
    ctrl%restartFreq = 1

    call getChildValue(node, "LatticeOpt", ctrl%tLatOpt, .false.)
    if (ctrl%tLatOpt) then
      if (allocated(ctrl%hybridXcInp)) then
        call detailedError(node, "Lattice optimisation not currently implemented for hybrid&
            & functionals")
      end if
      call getChildValue(node, "Pressure", ctrl%pressure, 0.0_dp, modifier=modifier, child=child)
      call convertUnitHsd(modifier, pressureUnits, child, ctrl%pressure)
      call getChildValue(node, "FixAngles", ctrl%tLatOptFixAng, .false.)
      if (ctrl%tLatOptFixAng) then
        call getChildValue(node, "FixLengths", ctrl%tLatOptFixLen, [.false.,.false.,.false.])
      else
        call getChildValue(node, "Isotropic", ctrl%tLatOptIsotropic, .false.)
      end if
      if (isMaxStep) then
        call getChildValue(node, "MaxLatticeStep", ctrl%maxLatDisp, 0.2_dp)
      end if
    end if
    call getChildValue(node, "MovedAtoms", buffer2, trim(atomsRange), child=child, multiple=.true.)
    call getSelectedAtomIndices(child, buffer2, geom%speciesNames, geom%species,&
        & ctrl%indMovedAtom)

    ctrl%nrMoved = size(ctrl%indMovedAtom)
    ctrl%tCoordOpt = (ctrl%nrMoved /= 0)
    if (ctrl%tCoordOpt) then
      if (isMaxStep) then
        call getChildValue(node, "MaxAtomStep", ctrl%maxAtomDisp, 0.2_dp)
      end if
    end if
    call getChildValue(node, "MaxForceComponent", ctrl%maxForce, 1e-4_dp, modifier=modifier,&
        & child=field)
    call convertUnitHsd(modifier, forceUnits, field, ctrl%maxForce)
    call getChildValue(node, "MaxSteps", ctrl%maxRun, 200)
    call getChildValue(node, "StepSize", ctrl%deltaT, 100.0_dp, modifier=modifier, child=field)
    call convertUnitHsd(modifier, timeUnits, field, ctrl%deltaT)
    call getChildValue(node, "OutputPrefix", buffer2, "geo_end")
    ctrl%outFile = unquote(buffer2)
    call getChildValue(node, "AppendGeometries", ctrl%tAppendGeo, .false.)
    call readGeoConstraints(node, ctrl, geom%nAtom)
    if (ctrl%tLatOpt) then
      if (ctrl%nrConstr/=0) then
        call error("Lattice optimisation and constraints currently incompatible.")
      end if
      if (ctrl%nrMoved/=0.and.ctrl%nrMoved<geom%nAtom) then
        call error("Subset of optimising atoms not currently possible with lattice optimisation.")
      end if
    end if
    ctrl%isGeoOpt = ctrl%tLatOpt .or. ctrl%tCoordOpt

  end subroutine commonGeoOptions


  !> Extended lagrangian options for XLBOMD
  subroutine readXlbomdOptions(node, input)

    !> node in the input tree
    type(hsd_table), pointer :: node

    !> extracted settings on exit
    type(TXLBOMDInp), allocatable, intent(out) :: input

    type(hsd_table), pointer :: pXlbomd, pXlbomdFast, pRoot, pChild
    logical :: tXlbomdFast

    call getChild(node, 'Xlbomd', pXlbomd, requested=.false.)
    call getChild(node, 'XlbomdFast', pXlbomdFast, requested=.false.)
    if (.not. (associated(pXlbomd) .or. associated(pXlbomdFast))) then
      return
    end if
    if (associated(pXlbomd) .and. associated(pXlbomdFast)) then
      call detailedError(pXlbomdFast, "Blocks 'Xlbomd' and 'XlbomdFast' are&
          & mutually exclusive")
    end if
    if (associated(pXlbomdFast)) then
      tXlbomdFast = .true.
      pRoot => pXlbomdFast
    else
      tXlbomdFast = .false.
      pRoot => pXlbomd
    end if
    allocate(input)
    call getChildValue(pRoot, 'IntegrationSteps', input%nKappa, 5, child=pChild)
    if (all([5, 6, 7] /= input%nKappa)) then
      call detailedError(pChild, 'Invalid number of integration steps (must be&
          & 5, 6 or 7)')
    end if
    call getChildValue(pRoot, 'PreSteps', input%nPreSteps, 0)

    ! Since support for inverse Jacobian has been removed, we can set FullSccSteps
    ! to its minimal value (no averaging of inverse Jacobians is done anymore)
    input%nFullSccSteps = input%nKappa + 1

    if (tXlbomdFast) then
      call getChildValue(pRoot, 'TransientSteps', input%nTransientSteps, 10)
      input%minSccIter = 1
      input%maxSccIter = 1
      ! Dummy value as minSccIter and maxSccIter have been set to 1.
      input%sccTol = 1e-5_dp
      call getChildValue(pRoot, 'Scale', input%scale, 1.0_dp, child=pChild)
      if (input%scale <= 0.0_dp .or. input%scale > 1.0_dp) then
        call detailedError(pChild, 'Scaling value must be in the interval&
            & (0.0, 1.0]')
      end if

    else
      input%nTransientSteps = 0
      call getChildValue(pRoot, 'MinSccIterations', input%minSCCIter, 1)
      call getChildValue(pRoot, 'MaxSccIterations', input%maxSCCIter, 200)
      if (input%maxSCCIter <= 0) then
        call detailedError(pRoot,"MaxSccIterations must be >= 1");
      end if
      call getChildValue(pRoot, 'SccTolerance', input%sccTol, 1e-5_dp)
      input%scale = 1.0_dp
    end if

  end subroutine readXlbomdOptions


  !> Reads geometry constraints
  subroutine readGeoConstraints(node, ctrl, nAtom)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer
    type(TListIntR1) :: intBuffer
    type(TListRealR1) :: realBuffer

    call getChildValue(node, "Constraints", value1, "", child=child, allowEmptyValue=.true.)
    call getNodeName2(value1, buffer)
    if (buffer == "") then
      ctrl%nrConstr = 0
    else
      call init(intBuffer)
      call init(realBuffer)
      call getChildValue(child, "", 1, intBuffer, 3, realBuffer)
      ctrl%nrConstr = len(intBuffer)
      allocate(ctrl%conAtom(ctrl%nrConstr))
      allocate(ctrl%conVec(3, ctrl%nrConstr))
      call asVector(intBuffer, ctrl%conAtom)
      if (.not.all(ctrl%conAtom<=nAtom)) then
        call detailedError(node,"Non-existent atom specified in constraint")
      end if
      call asArray(realBuffer, ctrl%conVec)
      call destruct(intBuffer)
      call destruct(realBuffer)
    end if

  end subroutine readGeoConstraints


  !> Reads MD velocities
  subroutine readInitialVelocities(node, ctrl, nAtom)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Total number of all atoms
    integer, intent(in) :: nAtom

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer, modifier
    type(TListRealR1) :: realBuffer
    integer :: nVelocities
    real(dp), allocatable :: tmpVelocities(:,:)

    call getChildValue(node, "Velocities", value1, "", child=child, modifier=modifier,&
        & allowEmptyValue=.true.)
    call getNodeName2(value1, buffer)
    if (buffer == "") then
      ctrl%tReadMDVelocities = .false.
    else
      call init(realBuffer)
      call getChildValue(child, "", 3, realBuffer, modifier=modifier)
      nVelocities = len(realBuffer)
      if (nVelocities /= nAtom) then
        call detailedError(node, "Incorrect number of specified velocities: " &
            & // i2c(3*nVelocities) // " supplied, " &
            & // i2c(3*nAtom) // " required.")
      end if
      allocate(tmpVelocities(3, nVelocities))
      call asArray(realBuffer, tmpVelocities)
      if (len(modifier) > 0) then
        call convertUnitHsd(modifier, VelocityUnits, child, &
            & tmpVelocities)
      end if
      call destruct(realBuffer)
      allocate(ctrl%initialVelocities(3, ctrl%nrMoved))
      ctrl%initialVelocities(:,:) = tmpVelocities(:,ctrl%indMovedAtom(:))
      ctrl%tReadMDVelocities = .true.
    end if

  end subroutine readInitialVelocities


  subroutine readTemperature(node, tempProfInp)

    !> Temperature node
    type(hsd_table), pointer :: node

    !> Temperature profile input data on exit
    type(TTempProfileInput), intent(out) :: tempProfInp

    character(len=:), allocatable :: modifier
    real(dp) :: temp

    call getChildValue(node, "", temp, modifier=modifier)
    call convertUnitHsd(modifier, energyUnits, node, temp)
    if (temp < 0.0_dp) call detailedError(node, "Negative temperature.")
    temp = max(minTemp, temp)
    tempProfInp%tempInts = [huge(1)]
    tempProfInp%tempValues = [temp]
    tempProfInp%tempMethods = [tempProfileTypes%constant]

  end subroutine readTemperature


  !> reads a temperature profile for MD with correctness checking of the input
  subroutine readTemperatureProfile(node, modifier, tempProfInp)

    !> Temperature profile node
    type(hsd_table), pointer :: node

    !> unit modifier of the node
    character(len=*), intent(in) :: modifier

    !> Temperature profile input data on exit
    type(TTempProfileInput), intent(out) :: tempProfInp

    type(TListString) :: ls
    type(TListIntR1) :: li1
    type(TListRealR1) :: lr1
    character(len=20), allocatable :: tmpC1(:)
    integer :: ii
    logical :: success

    call init(ls)
    call init(li1)
    call init(lr1)
    call getChildValue(node, "", ls, 1, li1, 1, lr1)
    if (len(ls) < 1) then
      call detailedError(node, "At least one annealing step must be specified.")
    end if
    allocate(tmpC1(len(ls)))
    allocate(tempProfInp%tempInts(len(li1)))
    allocate(tempProfInp%tempValues(len(lr1)))
    call asArray(ls, tmpC1)
    call asVector(li1, tempProfInp%tempInts)
    call asVector(lr1, tempProfInp%tempValues)
    call destruct(ls)
    call destruct(li1)
    call destruct(lr1)
    allocate(tempProfInp%tempMethods(size(tmpC1)))
    do ii = 1, size(tmpC1)
      call identifyTempProfile(tempProfInp%tempMethods(ii), tmpC1(ii), success)
      if (success) then
        cycle
      end if
      call detailedError(node, "Invalid annealing method name '" // trim(tmpC1(ii)) // "'.")
    end do

    if (any(tempProfInp%tempInts < 0)) then
      call detailedError(node, "Step values must not be negative.")
    end if

    if (sum(tempProfInp%tempInts) == 0) then
      call detailedError(node, "Sum of steps in the profile must be greater than zero.")
    end if

    if (any(tempProfInp%tempValues < 0.0_dp)) then
      call detailedError(node, "Negative temperature.")
    end if

    call convertUnitHsd(modifier, energyUnits, node, tempProfInp%tempValues)
    if (any(tempProfInp%tempValues < minTemp)) then
      tempProfInp%tempValues = max(tempProfInp%tempValues, minTemp)
    end if

  end subroutine readTemperatureProfile

  subroutine readElecDynamics(node, input, geom, masses)

    !> input data to parse
    type(hsd_table), pointer :: node

    !> ElecDynamicsInp instance
    type(TElecDynamicsInp), intent(inout) :: input

    !> geometry of the system
    type(TGeometry), intent(in) :: geom

    !> masses to be returned
    real(dp), allocatable, intent(inout) :: masses(:)

    type(hsd_table), pointer :: child, value1
    character(len=:), allocatable :: buffer, buffer2, modifier
    logical :: ppRangeInvalid, tNeedFieldStrength
    real (dp) :: defPpRange(2)
    logical :: defaultWrite

    call getChildValue(node, "Steps", input%steps)
    call getChildValue(node, "TimeStep", input%dt, modifier=modifier, child=child)
    call convertUnitHsd(modifier, timeUnits, child, input%dt)

    call getChildValue(node, "Populations", input%tPopulations, .false.)
    call getChildValue(node, "WriteFrequency", input%writeFreq, 50)
    call getChildValue(node, "Restart", input%tReadRestart, .false.)
    if (input%tReadRestart) then
      call getChildValue(node, "RestartFromAscii", input%tReadRestartAscii, .false.)
    end if
    call getChildValue(node, "WriteRestart", input%tWriteRestart, .true.)
    if (input%tWriteRestart) then
      call getChildValue(node, "WriteAsciiRestart", input%tWriteRestartAscii, .false.)
    end if
    call getChildValue(node, "RestartFrequency", input%restartFreq, max(input%Steps / 10, 1))
    call getChildValue(node, "Forces", input%tForces, .false.)
    call getChildValue(node, "WriteBondEnergy", input%tBondE, .false.)
    call getChildValue(node, "WriteBondPopulation", input%tBondP, .false.)
    call getChildValue(node, "WriteAtomicEnergies", input%tWriteAtomEnergies, .false.)
    call getChildValue(node, "Pump", input%tPump, .false.)
    call getChildValue(node, "FillingsFromFile", input%tFillingsFromFile, .false.)

    if (input%tPump) then
      call getChildValue(node, "PumpProbeFrames", input%tdPPFrames)
      defPpRange = [0.0_dp, input%steps * input%dt]
      call getChildValue(node, "PumpProbeRange", input%tdPpRange, defPprange, modifier=modifier,&
          & child=child)
      call convertUnitHsd(modifier, timeUnits, child, input%tdPpRange)

      ppRangeInvalid = (input%tdPpRange(2) <= input%tdPpRange(1))&
          & .or. (input%tdPprange(1) < defPpRange(1))&
          & .or. (input%tdPpRange(2) > defPpRange(2))
      if (ppRangeInvalid) then
        call detailederror(child, "Wrong definition of PumpProbeRange, either incorrect order&
            & or outside of simulation time range")
      end if
    end if

    call getChildValue(node, "Probe", input%tProbe, .false.)
    if (input%tPump .and. input%tProbe) then
      call detailedError(child, "Pump and probe cannot be simultaneously true.")
    end if

    call getChildValue(node, "EulerFrequency", input%eulerFreq, 0)

    call getChildValue(node, "VerboseDynamics", input%tVerboseDyn, .true.)

    if ((input%eulerFreq < 50) .and. (input%eulerFreq > 0)) then
      call detailedError(child, "Wrong number of Euler steps, should be above 50")
    end if
    if (input%eulerFreq >= 50) then
      input%tEulers = .true.
    else
      input%tEulers = .false.
    end if

    ! assume this is required (needed for most perturbations, but not none)
    tNeedFieldStrength = .true.

    defaultWrite = .true.

    !! Different perturbation types
    call getChildValue(node, "Perturbation", value1, "None", child=child)
    call getNodeName(value1, buffer)
    select case(buffer)

    case ("kick")
      input%pertType = pertTypes%kick
      call hsd_rename_child(value1, "PolarizationDirection", "PolarisationDirection")
      call getChildValue(value1, "PolarisationDirection", buffer2)
      input%polDir = directionConversion(unquote(buffer2), value1)

      call getChildValue(value1, "SpinType", buffer2, "Singlet")
      select case(unquote(buffer2))
      case ("singlet", "Singlet")
        input%spType = tdSpinTypes%singlet
      case ("triplet", "Triplet")
        input%spType = tdSpinTypes%triplet
      case default
        call detailedError(value1, "Unknown spectrum spin type " // buffer2)
      end select

      defaultWrite = .false.

    case ("laser")
      input%pertType = pertTypes%laser
      call hsd_rename_child(value1, "PolarizationDirection", "PolarisationDirection")
      call getChildValue(value1, "PolarisationDirection", input%reFieldPolVec)
      call hsd_rename_child(value1, "ImagPolarizationDirection", "ImagPolarisationDirection")
      call getChildValue(value1, "ImagPolarisationDirection", input%imFieldPolVec, &
          & [0.0_dp, 0.0_dp, 0.0_dp])
      call getChildValue(value1, "LaserEnergy", input%omega, modifier=modifier, child=child)
      call convertUnitHsd(modifier, energyUnits, child, input%omega)
      call getChildValue(value1, "Phase", input%phase, 0.0_dp, modifier=modifier, child=child)
      call convertUnitHsd(modifier, angularUnits, child, input%phase)
      call getChildValue(value1, "ExcitedAtoms", buffer, "1:-1", child=child, multiple=.true.)
      call getSelectedAtomIndices(child, buffer, geom%speciesNames, geom%species,&
          & input%indExcitedAtom)

      input%nExcitedAtom = size(input%indExcitedAtom)
      if (input%nExcitedAtom == 0) then
        call error("No atoms specified for laser excitation.")
      end if

      defaultWrite = .true.

    case ("kickandlaser")
      input%pertType = pertTypes%kickAndLaser
      call getChildValue(value1, "KickPolDir", buffer2)
      input%polDir = directionConversion(unquote(buffer2), value1)
      call getChildValue(value1, "SpinType", input%spType, tdSpinTypes%singlet)
      call getChildValue(value1, "LaserPolDir", input%reFieldPolVec)
      call getChildValue(value1, "LaserImagPolDir", input%imFieldPolVec, [0.0_dp, 0.0_dp, 0.0_dp])
      call getChildValue(value1, "LaserEnergy", input%omega, modifier=modifier, child=child)
      call convertUnitHsd(modifier, energyUnits, child, input%omega)
      call getChildValue(value1, "Phase", input%phase, 0.0_dp, modifier=modifier, child=child)
      call convertUnitHsd(modifier, angularUnits, child, input%phase)
      call getChildValue(value1, "LaserStrength", input%tdLaserField, modifier=modifier,&
          & child=child)
      call convertUnitHsd(modifier, EFieldUnits, child, input%tdLaserField)

      call getChildValue(value1, "ExcitedAtoms", buffer, "1:-1", child=child, multiple=.true.)
      call getSelectedAtomIndices(child, buffer, geom%speciesNames, geom%species,&
          & input%indExcitedAtom)
      input%nExcitedAtom = size(input%indExcitedAtom)
      if (input%nExcitedAtom == 0) then
        call error("No atoms specified for laser excitation.")
      end if

      defaultWrite = .false.

    case ("none")
      input%pertType = pertTypes%noTDPert
      tNeedFieldStrength = .false.

      defaultWrite = .true.

    case default
      call detailedError(child, "Unknown perturbation type " // buffer)
    end select

    if (tNeedFieldStrength) then
      call getChildValue(node, "FieldStrength", input%tdfield, modifier=modifier, child=child)
      call convertUnitHsd(modifier, EFieldUnits, child, input%tdfield)
    end if

    call getChildValue(node, "WriteEnergyAndCharges", input%tdWriteExtras, defaultWrite)

    !! Different envelope functions
    call getChildValue(node, "EnvelopeShape", value1, "Constant")
    call getNodeName(value1, buffer)
    select case(buffer)

    case("constant")
      input%envType = envTypes%constant

    case("gaussian")
      input%envType = envTypes%gaussian
      call getChildValue(value1, "Time0", input%time0, 0.0_dp, modifier=modifier, child=child)
      call convertUnitHsd(modifier, timeUnits, child, input%Time0)

      call getChildValue(value1, "Time1", input%time1, modifier=modifier, child=child)
      call convertUnitHsd(modifier, timeUnits, child, input%Time1)

    case("sin2")
      input%envType = envTypes%sin2
      call getChildValue(value1, "Time0", input%time0, 0.0_dp, modifier=modifier, child=child)
      call convertUnitHsd(modifier, timeUnits, child, input%Time0)

      call getChildValue(value1, "Time1", input%time1, modifier=modifier, child=child)
      call convertUnitHsd(modifier, timeUnits, child, input%Time1)

    case("fromfile")
      input%envType = envTypes%fromFile
      call getChildValue(value1, "Time0", input%time0, 0.0_dp, modifier=modifier, child=child)
      call convertUnitHsd(modifier, timeUnits, child, input%Time0)

    case default
      call detailedError(value1, "Unknown envelope shape " // buffer)
    end select

    !! Non-adiabatic molecular dynamics
    call getChildValue(node, "IonDynamics", input%tIons, .false.)
    if (input%tIons) then
      call getChildValue(node, "MovedAtoms", buffer, "1:-1", child=child, multiple=.true.)
      call getSelectedAtomIndices(child, buffer, geom%speciesNames, geom%species,&
          & input%indMovedAtom)

      input%nMovedAtom = size(input%indMovedAtom)
      call readInitialVelocitiesNAMD(node, input, geom%nAtom)
      if (input%tReadMDVelocities) then
        ! without a thermostat, if we know the initial velocities, we do not need a temperature, so
        ! just set it to something 'safe'
        input%tempAtom = minTemp
      else
        if (.not. input%tReadRestart) then
          ! previously lower limit was minTemp:
          call readMDInitTemp(node, input%tempAtom, 0.0_dp)
        end if
        call getInputMasses(node, geom, masses)
      end if
    end if

  end subroutine readElecDynamics


  !> Converts direction label text string into corresponding numerical value
  function directionConversion(direction, node) result(iX)

    !> Direction label
    character(*), intent(in) :: direction

    !> input tree for error return
    type(hsd_table), pointer :: node

    !> direction indicator (1 - 4) for (x,y,z,all)
    integer :: iX

    select case(trim(direction))
    case ("x", "X")
      iX = 1
    case ("y", "Y")
      iX = 2
    case ("z", "Z")
      iX = 3
    case ("all", "All", "ALL")
      iX = 4
    case default
      call detailedError(node, "Wrongly specified polarisation direction " // trim(direction)&
          & // ". Must be x, y, z or all.")
    end select

  end function directionConversion


  !> Reads MD velocities
  subroutine readInitialVelocitiesNAMD(node, input, nAtom)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> ElecDynamicsInp object structure to be filled
    type(TElecDynamicsInp), intent(inout) :: input

    !> Total number of all atoms
    integer, intent(in) :: nAtom

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer, modifier
    type(TListRealR1) :: realBuffer
    integer :: nVelocities
    real(dp), pointer :: tmpVelocities(:,:)

    call getChildValue(node, "Velocities", value1, "", child=child, modifier=modifier,&
        & allowEmptyValue=.true.)
    call getNodeName2(value1, buffer)
    if (buffer == "") then
       input%tReadMDVelocities = .false.
    else
       call init(realBuffer)
       call getChildValue(child, "", 3, realBuffer, modifier=modifier)
       nVelocities = len(realBuffer)
       if (nVelocities /= nAtom) then
          call detailedError(node, "Incorrect number of specified velocities: " &
               & // i2c(3*nVelocities) // " supplied, " &
               & // i2c(3*nAtom) // " required.")
       end if
       allocate(tmpVelocities(3, nVelocities))
       call asArray(realBuffer, tmpVelocities)
       if (len(modifier) > 0) then
          call convertUnitHsd(modifier, VelocityUnits, child, &
               & tmpVelocities)
       end if
       call destruct(realBuffer)
       allocate(input%initialVelocities(3, input%nMovedAtom))
       input%initialVelocities(:,:) = tmpVelocities(:, input%indMovedAtom(:))
       input%tReadMDVelocities = .true.
    end if

  end subroutine readInitialVelocitiesNAMD

  subroutine parseThermostat(node, deltaT, hasInitVelocities, maxRun, thermostatInp, tempProfileInp)

    !> Parent node of the thermostat node
    type(hsd_table), pointer, intent(in) :: node
    
    !> Time step
    real(dp), intent(in) :: deltaT
    
    !> Whether initial velocities had been specified for the MD run
    logical, intent(in) :: hasInitVelocities
    
    !> Number of MD timesteps, will be updated by adding up the steps in the temperature profile
    integer, intent(inout) :: maxRun
    
    !> Thermostat input filled up from the HSD data
    type(TThermostatInput), allocatable, intent(out) :: thermostatInp
    
    !> Temperature profile input filled up from the HSD data
    type(TTempProfileInput), allocatable, intent(out) :: tempProfileInp

    type(hsd_table), pointer :: thermNode, child, child2, child3
    character(len=:), allocatable :: thermName, modifier

    allocate(thermostatInp, tempProfileInp)
    call getChildValue(node, "Thermostat", thermNode, child=child)
    call getNodeName(thermNode, thermName)

    select case(thermName)

    case ("berendsen")

      thermostatInp%thermostatType = thermostatTypes%berendsen
      allocate(thermostatInp%berendsen)
      associate (inp => thermostatInp%berendsen)
        call readTempOrTempProfile_(thermNode, maxRun, tempProfileInp)
        call getChild(thermNode, "CouplingStrength", child=child2, requested=.false.)
        if (associated(child2)) then
          call getChildValue(child2, "", inp%coupling)
          call getChild(thermNode, "Timescale", child=child2, modifier=modifier, requested=.false.)
          if (associated(child2)) then
            call error("Only Coupling strength OR Timescale can be set for Berendsen thermostats.")
          end if
        else
          call getChild(thermNode, "Timescale", child=child2, modifier=modifier, requested=.false.)
          if (associated(child2)) then
            call getChildValue(child2, "", inp%coupling, modifier=modifier, child=child3)
            call convertUnitHsd(modifier, timeUnits, child3, inp%coupling)
            inp%coupling = deltaT / inp%coupling
          else
            call error("Either CouplingStrength or Timescale must be set for Berendsen thermostats.")
          end if
        end if
      end associate

    case ("nosehoover")

      thermostatInp%thermostatType = thermostatTypes%nhc
      allocate(thermostatInp%nhc)
      associate (inp => thermostatInp%nhc)
        call readTempOrTempProfile_(thermNode, maxRun, tempProfileInp)
        call getChildValue(thermNode, "CouplingStrength", inp%coupling, modifier=modifier,&
            & child=child2)
        call convertUnitHsd(modifier, freqUnits, child2, inp%coupling)

        call getChildValue(thermNode, "ChainLength", inp%chainLength, 3)
        call getChildValue(thermNode, "Order", inp%expOrder, 3, child=child2)
        if (.not. any(inp%expOrder == [3, 5])) then
          call detailedError(child2, "Order of Nose-Hoover thermostat must be either 3 or 5")
        end if
        call getChildValue(thermNode, "IntegratorSteps", inp%nExpSteps, 1)
        call getChild(thermNode, "Restart",  child=child2, requested=.false.)
        if (associated(child2)) then
          allocate(inp%xnose(inp%chainLength))
          allocate(inp%vnose(inp%chainLength))
          allocate(inp%gnose(inp%chainLength))
          call getChildValue(child2,"x",inp%xnose)
          call getChildValue(child2,"v",inp%vnose)
          call getChildValue(child2,"g",inp%gnose)
        end if
      end associate

    case ("andersen")

      thermostatInp%thermostatType = thermostatTypes%andersen
      allocate(thermostatInp%andersen)
      associate (inp => thermostatInp%andersen)
        call readTempOrTempProfile_(thermNode, maxRun, tempProfileInp)
        call getChildValue(thermNode, "ReselectProbability", inp%rescaleProb, child=child2)
        if (inp%rescaleProb <= 0.0_dp .or. inp%rescaleProb > 1.0_dp) then
          call detailedError(child2, "ReselectProbability must be in the range (0,1]!")
        end if
        call getChildValue(thermNode, "ReselectIndividually", inp%rescaleIndiv)
      end associate

    case ("none")

      ! Create a fake thermostat with a single constant temperature value
      ! It will only used to generate the initial velocities for the MD anyway.
      thermostatInp%thermostatType = thermostatTypes%dummy
      tempProfileInp%tempInts = [huge(1)]
      tempProfileInp%tempMethods = [tempProfileTypes%constant]
      tempProfileInp%tempValues = [minTemp]
      if (.not. hasInitVelocities) then
        ! Initial velocities had not been provided, overwrite 'safe' default value (minTemp)
        ! by reading the temperature explicitly (needed for generating the initial velocities)
        call readMDInitTemp(thermNode, tempProfileInp%tempValues(1), minTemp)
      end if

    case default
      call getNodeHSDName(thermNode, thermName)
      call detailedError(child, "Invalid thermostat '" // thermName // "'")

    end select

  contains

    !> Reads the temperature or the temperature profile
    subroutine readTempOrTempProfile_(thermNode, maxRun, tempProfileInp)
      type(hsd_table), pointer, intent(in) :: thermNode
      integer, intent(inout) :: maxRun
      type(TTempProfileInput), intent(out) :: tempProfileInp

      type(hsd_table), pointer :: value, child
      character(len=:), allocatable :: buffer, modifier

      call getChildValue(thermNode, "Temperature", value, modifier=modifier, child=child)
      call getNodeName(value, buffer)
      select case(buffer)
      case (textNodeName)
        call readTemperature(child, tempProfileInp)
      case ("temperatureprofile")
        call readTemperatureProfile(value, modifier, tempProfileInp)
        maxRun = sum(tempProfileInp%tempInts) - 1
      case default
        call detailedError(value, "Invalid method name.")
      end select

    end subroutine readTempOrTempProfile_

  end subroutine parseThermostat


end module dftbp_dftbplus_parser_driver
