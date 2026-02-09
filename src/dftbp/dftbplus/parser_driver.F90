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
  use dftbp_io_charmanip, only : i2c, tolower, unquote
  use hsd, only : hsd_rename_child, hsd_get, hsd_get_or_set, hsd_get_matrix, hsd_get_table, &
      & hsd_get_attrib, hsd_get_choice, HSD_STAT_OK
  use hsd_data, only : hsd_table, hsd_node
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, getSelectedAtomIndices,&
      & getNodeName, getNodeName2, hasInlineData, &
      & getFirstTextChild
  use dftbp_io_tokenreader, only : getNextToken, TOKEN_EOS, TOKEN_ERROR, TOKEN_OK
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_io_message, only : error
  use dftbp_md_tempprofile, only : identifyTempProfile, tempProfileTypes, TTempProfileInput
  use dftbp_md_thermostats, only : thermostatTypes, TThermostatInput
  use dftbp_md_xlbomd, only : TXlbomdInp
  use dftbp_timedep_timeprop, only : envTypes, pertTypes, tdSpinTypes, TElecDynamicsInp
  use dftbp_transport_negfvars, only : TTransPar
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
    integer :: stat

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
        call dftbp_error(node, "GeometryOptimisation driver currently does not support helical&
            & geometries")
      end if

      allocate(ctrl%geoOpt)

      call readGeoOptInput(node, geom, ctrl%geoOpt, atomsRange)

      call hsd_get_or_set(node, "AppendGeometries", ctrl%tAppendGeo, .false.)

      ! Geometry optimisation drivers
      ctrl%iGeoOpt = geoOptTypes%geometryoptimisation
      ctrl%tForces = .true.
      ctrl%restartFreq = 1

    case ("steepestdescent")

      modeName = "geometry relaxation"
      call dftbp_warning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ! Steepest downhill optimisation
      ctrl%iGeoOpt = geoOptTypes%steepestDesc
      call commonGeoOptions(node, ctrl, geom, atomsRange)

    case ("conjugategradient")

      modeName = "geometry relaxation"
      call dftbp_warning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')// "Please use the GeometryOptimisation driver instead.")

      ! Conjugate gradient location optimisation
      ctrl%iGeoOpt = geoOptTypes%conjugateGrad
      call commonGeoOptions(node, ctrl, geom, atomsRange)

    case("gdiis")

      modeName = "geometry relaxation"
      call dftbp_warning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ! Gradient DIIS optimisation, only stable in the quadratic region
      ctrl%iGeoOpt = geoOptTypes%diis
      call hsd_get_or_set(node, "alpha", ctrl%deltaGeoOpt, 1.0E-1_dp)
      call hsd_get_or_set(node, "Generations", ctrl%iGenGeoOpt, 8)
      call commonGeoOptions(node, ctrl, geom, atomsRange)

    case ("lbfgs")

      modeName = "geometry relaxation"
      call dftbp_warning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ctrl%iGeoOpt = geoOptTypes%lbfgs

      allocate(ctrl%lbfgsInp)
      call hsd_get_or_set(node, "Memory", ctrl%lbfgsInp%memory, 20)

      call hsd_get_or_set(node, "LineSearch", ctrl%lbfgsInp%isLineSearch, .false.)

      isMaxStepNeeded = .not. ctrl%lbfgsInp%isLineSearch
      if (isMaxStepNeeded) then
        call hsd_get_or_set(node, "setMaxStep", ctrl%lbfgsInp%isLineSearch, isMaxStepNeeded)
        ctrl%lbfgsInp%MaxQNStep = isMaxStepNeeded
      else
        call hsd_get_or_set(node, "oldLineSearch", ctrl%lbfgsInp%isOldLS, .false.)
      end if

      call commonGeoOptions(node, ctrl, geom, atomsRange, ctrl%lbfgsInp%isLineSearch)

    case ("fire")

      modeName = "geometry relaxation"
      call dftbp_warning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ctrl%iGeoOpt = geoOptTypes%fire
      call commonGeoOptions(node, ctrl, geom, atomsRange, .false.)
      call hsd_get_or_set(node, "TimeStep", ctrl%deltaT, 1.0_dp)
      call hsd_get_attrib(node, "TimeStep", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, timeUnits, node, ctrl%deltaT)

    case("secondderivatives")
      ! currently only numerical derivatives of forces is implemented

      modeName = "second derivatives"

      ctrl%tDerivs = .true.
      ctrl%tForces = .true.

      call hsd_get_or_set(node, "Atoms", buffer2, trim(atomsRange))
      call hsd_get_table(node, "Atoms", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) child => node
      call getSelectedAtomIndices(child, buffer2, geom%speciesNames, geom%species,&
          & ctrl%indDerivAtom)
      if (size(ctrl%indDerivAtom) == 0) then
        call error("No atoms specified for derivatives calculation.")
      end if

      call hsd_get_table(node, "MovedAtoms", child, stat, auto_wrap=.true.)
      if (associated(child)) then
        if (.not. isContiguousRange(ctrl%indDerivAtom)) then
          call dftbp_error(child,&
            & "Atoms for calculation of partial Hessian must be a contiguous range.")
        end if
        call hsd_get(child, "#text", buffer2, stat=stat)
        if (stat /= HSD_STAT_OK) buffer2 = ""
        child2 => child
        call getSelectedAtomIndices(child2, buffer2, geom%speciesNames, geom%species, &
           & ctrl%indMovedAtom)
        if (.not. isContiguousRange(ctrl%indMovedAtom)) then
          call dftbp_error(child2, "MovedAtoms for calculation of partial Hessian must be a &
              & contiguous range.")
        end if
        if (.not. containsAll(ctrl%indDerivAtom, ctrl%indMovedAtom)) then
          call dftbp_error(child2, "MovedAtoms has indices not contained in Atoms.")
        end if
      else
        ctrl%indMovedAtom = ctrl%indDerivAtom
      end if
      ctrl%nrMoved = size(ctrl%indMovedAtom)

      call hsd_get_or_set(node, "Delta", ctrl%deriv2ndDelta, 1.0E-4_dp)
      call hsd_get_attrib(node, "Delta", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, lengthUnits, node, ctrl%deriv2ndDelta)

    case ("velocityverlet")
      ! molecular dynamics

      modeName = "molecular dynamics"

      ctrl%tForces = .true.
      ctrl%tMD = .true.

      call hsd_get_or_set(node, "MDRestartFrequency", ctrl%restartFreq, 1)
      call hsd_get_or_set(node, "MovedAtoms", buffer2, trim(atomsRange))
      call hsd_get_table(node, "MovedAtoms", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) child => node
      call getSelectedAtomIndices(child, buffer2, geom%speciesNames, geom%species, &
          & ctrl%indMovedAtom)
      ctrl%nrMoved = size(ctrl%indMovedAtom)
      if (ctrl%nrMoved == 0) then
        call error("No atoms specified for molecular dynamics.")
      end if
      call readInitialVelocities(node, ctrl, geom%nAtom)

      call hsd_get_or_set(node, "KeepStationary", ctrl%tMDstill,.true.)
      if (ctrl%tMDstill .and. geom%nAtom == 1) then
        call error("Removing translational freedom with only one atom not&
            & possible.")
      end if

      call hsd_get(node, "TimeStep", ctrl%deltaT, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'TimeStep'")
      call hsd_get_attrib(node, "TimeStep", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, timeUnits, node, ctrl%deltaT)

      call parseThermostat(node, ctrl%deltaT, ctrl%tReadMDVelocities, ctrl%maxRun,&
          & ctrl%thermostatInp, ctrl%tempProfileInp)

      if (ctrl%maxRun < -1) then
        call hsd_get(node, "Steps", ctrl%maxRun, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'Steps'")
      end if

      call hsd_get_or_set(node, "OutputPrefix", buffer2, "geo_end")
      ctrl%outFile = unquote(buffer2)

      call hsd_get_or_set(node, "Plumed", ctrl%tPlumed, .false.)
      if (ctrl%tPlumed .and. .not. withPlumed) then
        call dftbp_error(node, "Metadynamics can not be used since code has been compiled&
            & without PLUMED support")
      end if

      if (geom%tPeriodic) then

        ctrl%tBarostat = .false.
        ctrl%pressure = 0.0_dp
        ctrl%BarostatStrength = 0.0_dp
        call hsd_get_table(node, "Barostat", child, stat, auto_wrap=.true.)
        if (associated(child)) then
          if (allocated(ctrl%hybridXcInp)) then
            call dftbp_error(node, "Barostating not currently implemented for hybrid functionals")
          end if
          if (ctrl%nrMoved /= geom%nAtom) then
            call error("Dynamics for a subset of atoms is not currently&
                & possible when using a barostat")
          end if
          call hsd_get(child, "Pressure", ctrl%pressure, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing required value: 'Pressure'")
          call hsd_get_attrib(child, "Pressure", modifier, stat)
          if (stat /= HSD_STAT_OK) modifier = ""
          call convertUnitHsd(modifier, pressureUnits, child, &
              & ctrl%pressure)
          call hsd_get_table(child, "Coupling", child2, stat, auto_wrap=.true.)
          if (associated(child2)) then
            call hsd_get(child2, "#text", ctrl%BarostatStrength, stat=stat)
            if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value")
            call hsd_get_table(child, "Timescale", child2, stat, auto_wrap=.true.)
            if (associated(child2)) call error("Only Coupling strength OR &
                &Timescale can be set for Barostatting.")
          else
            call hsd_get_table(child, "Timescale", child2, stat, auto_wrap=.true.)
            if (associated(child2)) then
              call hsd_get(child2, "#text", ctrl%BarostatStrength, stat=stat)
              if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value")
              if (allocated(child2%attrib)) then
                modifier = child2%attrib
              else
                modifier = ""
              end if
              call convertUnitHsd(modifier, timeUnits, child2, &
                  & ctrl%BarostatStrength)
              ctrl%BarostatStrength = ctrl%deltaT / ctrl%BarostatStrength
            else
              call error("Either Coupling strength or Timescale must be set&
                  & for Barostatting.")
            end if
          end if
          call hsd_get_or_set(child, "Isotropic", ctrl%tIsotropic, .true.)
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
      call hsd_get_table(node, 'File', child2, stat, auto_wrap=.true.)
      call hsd_get_table(node, 'Host', child3, stat, auto_wrap=.true.)
      if (associated(child2) .eqv. associated(child3)) then
        call error('Either Host or File (but not both) must be set for socket&
            & communication')
      end if

      ! File communication
      if (associated(child2)) then
        call hsd_get(child2, "#text", buffer2, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value")
        ctrl%socketInput%host = unquote(buffer2)
        ! zero it to signal to initprogram to use a unix file
        ctrl%socketInput%port = 0
      else
        call hsd_get(child3, "#text", buffer2, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child3, "Missing required value")
        ctrl%socketInput%host = unquote(buffer2)
        call hsd_get(node, "Port", ctrl%socketInput%port, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'Port'")
        if (ctrl%socketInput%port <= 0) then
          call dftbp_error(node, "Invalid port number")
        end if
      end if

      call hsd_get_choice(node, "Protocol", buffer, value1, stat)
      if (stat /= HSD_STAT_OK) buffer = "i-PI"
      buffer = tolower(buffer)
      select case(buffer)
      case("i-pi")
        ctrl%socketInput%protocol = IPI_PROTOCOLS%IPI_1
        ! want a file path
        if (ctrl%socketInput%port == 0) then
          call hsd_get_or_set(node, "Prefix", buffer2, "/tmp/ipi_")
          sTmp = unquote(buffer2)
          ctrl%socketInput%host = trim(sTmp) // trim(ctrl%socketInput%host)
        end if

      case default
        call dftbp_error(node, "Invalid protocol '" // buffer // "'")
      end select
      call hsd_get_or_set(node, "Verbosity", ctrl%socketInput%verbosity, 0)
      call hsd_get_or_set(node, "MaxSteps", ctrl%maxRun, 200)

    #:else
      call dftbp_error(node, "Program had been compiled without socket support")
    #:endif

    case default

      call getNodeName2(node, buffer)
      call dftbp_error(parent, "Invalid driver '" // buffer // "'")

    end select driver

  #:if WITH_TRANSPORT
    if (ctrl%solver%isolver == electronicSolverTypes%OnlyTransport .and. trim(modeName) /= "") then
      call dftbp_error(node, "transportOnly solver cannot be used with "//trim(modeName))
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
    integer :: stat

    if (present(isMaxStepNeeded)) then
      isMaxStep = isMaxStepNeeded
    else
      isMaxStep = .true.
    end if

    ctrl%tForces = .true.
    ctrl%restartFreq = 1

    call hsd_get_or_set(node, "LatticeOpt", ctrl%tLatOpt, .false.)
    if (ctrl%tLatOpt) then
      if (allocated(ctrl%hybridXcInp)) then
        call dftbp_error(node, "Lattice optimisation not currently implemented for hybrid&
            & functionals")
      end if
      call hsd_get_or_set(node, "Pressure", ctrl%pressure, 0.0_dp)
      call hsd_get_attrib(node, "Pressure", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, pressureUnits, node, ctrl%pressure)
      call hsd_get_or_set(node, "FixAngles", ctrl%tLatOptFixAng, .false.)
      if (ctrl%tLatOptFixAng) then
        block
          logical, allocatable :: tmpLogArr(:)
          call hsd_get_or_set(node, "FixLengths", tmpLogArr, [.false.,.false.,.false.])
          ctrl%tLatOptFixLen = tmpLogArr
        end block
      else
        call hsd_get_or_set(node, "Isotropic", ctrl%tLatOptIsotropic, .false.)
      end if
      if (isMaxStep) then
        call hsd_get_or_set(node, "MaxLatticeStep", ctrl%maxLatDisp, 0.2_dp)
      end if
    end if
    call hsd_get_or_set(node, "MovedAtoms", buffer2, trim(atomsRange))
    call hsd_get_table(node, "MovedAtoms", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) child => node
    call getSelectedAtomIndices(child, buffer2, geom%speciesNames, geom%species,&
        & ctrl%indMovedAtom)

    ctrl%nrMoved = size(ctrl%indMovedAtom)
    ctrl%tCoordOpt = (ctrl%nrMoved /= 0)
    if (ctrl%tCoordOpt) then
      if (isMaxStep) then
        call hsd_get_or_set(node, "MaxAtomStep", ctrl%maxAtomDisp, 0.2_dp)
      end if
    end if
    call hsd_get_or_set(node, "MaxForceComponent", ctrl%maxForce, 1e-4_dp)
    call hsd_get_attrib(node, "MaxForceComponent", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    call convertUnitHsd(modifier, forceUnits, node, ctrl%maxForce)
    call hsd_get_or_set(node, "MaxSteps", ctrl%maxRun, 200)
    call hsd_get_or_set(node, "StepSize", ctrl%deltaT, 100.0_dp)
    call hsd_get_attrib(node, "StepSize", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    call convertUnitHsd(modifier, timeUnits, node, ctrl%deltaT)
    call hsd_get_or_set(node, "OutputPrefix", buffer2, "geo_end")
    ctrl%outFile = unquote(buffer2)
    call hsd_get_or_set(node, "AppendGeometries", ctrl%tAppendGeo, .false.)
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
    integer :: stat

    call hsd_get_table(node, 'Xlbomd', pXlbomd, stat, auto_wrap=.true.)
    call hsd_get_table(node, 'XlbomdFast', pXlbomdFast, stat, auto_wrap=.true.)
    if (.not. (associated(pXlbomd) .or. associated(pXlbomdFast))) then
      return
    end if
    if (associated(pXlbomd) .and. associated(pXlbomdFast)) then
      call dftbp_error(pXlbomdFast, "Blocks 'Xlbomd' and 'XlbomdFast' are&
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
    call hsd_get_or_set(pRoot, 'IntegrationSteps', input%nKappa, 5, child=pChild)
    if (all([5, 6, 7] /= input%nKappa)) then
      call dftbp_error(pChild, 'Invalid number of integration steps (must be&
          & 5, 6 or 7)')
    end if
    call hsd_get_or_set(pRoot, 'PreSteps', input%nPreSteps, 0)

    ! Since support for inverse Jacobian has been removed, we can set FullSccSteps
    ! to its minimal value (no averaging of inverse Jacobians is done anymore)
    input%nFullSccSteps = input%nKappa + 1

    if (tXlbomdFast) then
      call hsd_get_or_set(pRoot, 'TransientSteps', input%nTransientSteps, 10)
      input%minSccIter = 1
      input%maxSccIter = 1
      ! Dummy value as minSccIter and maxSccIter have been set to 1.
      input%sccTol = 1e-5_dp
      call hsd_get_or_set(pRoot, 'Scale', input%scale, 1.0_dp, child=pChild)
      if (input%scale <= 0.0_dp .or. input%scale > 1.0_dp) then
        call dftbp_error(pChild, 'Scaling value must be in the interval&
            & (0.0, 1.0]')
      end if

    else
      input%nTransientSteps = 0
      call hsd_get_or_set(pRoot, 'MinSccIterations', input%minSCCIter, 1)
      call hsd_get_or_set(pRoot, 'MaxSccIterations', input%maxSCCIter, 200)
      if (input%maxSCCIter <= 0) then
        call dftbp_error(pRoot,"MaxSccIterations must be >= 1");
      end if
      call hsd_get_or_set(pRoot, 'SccTolerance', input%sccTol, 1e-5_dp)
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
    real(dp), allocatable :: constraintMatrix(:,:)
    integer :: nMatRows, nMatCols, ii, stat

    call hsd_get_table(node, "Constraints", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      ctrl%nrConstr = 0
    else
      call hsd_get_choice(child, "", buffer, value1, stat)
      if (stat /= HSD_STAT_OK) buffer = ""
      if (buffer == "" .and. .not. hasInlineData(child)) then
        ctrl%nrConstr = 0
      else
        call hsd_get_matrix(node, "Constraints", constraintMatrix, nMatRows, nMatCols, &
            & stat=stat, order="column-major")
        if (stat /= 0 .or. nMatRows /= 4) then
          call dftbp_error(child, "Invalid constraint data (expected rows of 1 int + 3 reals)")
        end if
        ctrl%nrConstr = nMatCols
        allocate(ctrl%conAtom(ctrl%nrConstr))
        allocate(ctrl%conVec(3, ctrl%nrConstr))
        do ii = 1, ctrl%nrConstr
          ctrl%conAtom(ii) = nint(constraintMatrix(1, ii))
        end do
        if (.not.all(ctrl%conAtom<=nAtom)) then
          call dftbp_error(node,"Non-existent atom specified in constraint")
        end if
        ctrl%conVec(:,:) = constraintMatrix(2:4, :)
      end if
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
    integer :: nVelocities, nMatRows, nMatCols, stat
    real(dp), allocatable :: tmpVelocities(:,:)

    call hsd_get_table(node, "Velocities", child, stat, auto_wrap=.true.)
    modifier = ""
    buffer = ""
    if (associated(child)) then
      if (allocated(child%attrib)) modifier = child%attrib
      call hsd_get_choice(child, "", buffer, value1, stat)
      if (stat /= HSD_STAT_OK) buffer = ""
    end if
    if (buffer == "" .and. (.not. associated(child) .or. .not. hasInlineData(child))) then
      ctrl%tReadMDVelocities = .false.
    else
      call hsd_get_matrix(node, "Velocities", tmpVelocities, nMatRows, nMatCols, &
          & stat=stat, order="column-major")
      if (stat /= 0 .or. nMatRows /= 3) then
        call dftbp_error(child, "Invalid velocity data (expected rows of 3 reals)")
      end if
      nVelocities = nMatCols
      if (nVelocities /= nAtom) then
        call dftbp_error(node, "Incorrect number of specified velocities: " &
            & // i2c(3*nVelocities) // " supplied, " &
            & // i2c(3*nAtom) // " required.")
      end if
      if (len(modifier) > 0) then
        call convertUnitHsd(modifier, VelocityUnits, child, &
            & tmpVelocities)
      end if
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
    integer :: stat

    call hsd_get(node, "#text", temp, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required temperature value")
    if (allocated(node%attrib)) then
      modifier = node%attrib
    else
      modifier = ""
    end if
    call convertUnitHsd(modifier, energyUnits, node, temp)
    if (temp < 0.0_dp) call dftbp_error(node, "Negative temperature.")
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

    character(len=:), allocatable :: text, bufferStr
    integer :: bufferInt(1)
    real(dp) :: bufferReal(1)
    integer :: iStart, iErr, nItem, nEntries, ii
    character(len=20), allocatable :: tmpC1(:)
    logical :: success

    ! Read raw text content and parse (string, int, real) triplets
    call getFirstTextChild(node, text)

    ! First pass: count entries
    nEntries = 0
    iStart = 1
    iErr = TOKEN_OK
    do while (iErr == TOKEN_OK)
      call getNextToken(text, bufferStr, iStart, iErr)
      if (iErr /= TOKEN_OK) exit
      call getNextToken(text, bufferInt, iStart, iErr, nItem)
      if (iErr /= TOKEN_OK) then
        call dftbp_error(node, "Invalid integer value in temperature profile")
      end if
      call getNextToken(text, bufferReal, iStart, iErr, nItem)
      if (iErr /= TOKEN_OK) then
        call dftbp_error(node, "Invalid real value in temperature profile")
      end if
      nEntries = nEntries + 1
    end do

    if (nEntries < 1) then
      call dftbp_error(node, "At least one annealing step must be specified.")
    end if

    ! Second pass: read data into arrays
    allocate(tmpC1(nEntries))
    allocate(tempProfInp%tempInts(nEntries))
    allocate(tempProfInp%tempValues(nEntries))
    iStart = 1
    do ii = 1, nEntries
      call getNextToken(text, bufferStr, iStart, iErr)
      tmpC1(ii) = bufferStr
      call getNextToken(text, bufferInt, iStart, iErr, nItem)
      tempProfInp%tempInts(ii) = bufferInt(1)
      call getNextToken(text, bufferReal, iStart, iErr, nItem)
      tempProfInp%tempValues(ii) = bufferReal(1)
    end do
    allocate(tempProfInp%tempMethods(size(tmpC1)))
    do ii = 1, size(tmpC1)
      call identifyTempProfile(tempProfInp%tempMethods(ii), tmpC1(ii), success)
      if (success) then
        cycle
      end if
      call dftbp_error(node, "Invalid annealing method name '" // trim(tmpC1(ii)) // "'.")
    end do

    if (any(tempProfInp%tempInts < 0)) then
      call dftbp_error(node, "Step values must not be negative.")
    end if

    if (sum(tempProfInp%tempInts) == 0) then
      call dftbp_error(node, "Sum of steps in the profile must be greater than zero.")
    end if

    if (any(tempProfInp%tempValues < 0.0_dp)) then
      call dftbp_error(node, "Negative temperature.")
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
    integer :: stat

    call hsd_get(node, "Steps", input%steps, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'Steps'")
    call hsd_get(node, "TimeStep", input%dt, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'TimeStep'")
    call hsd_get_attrib(node, "TimeStep", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    call convertUnitHsd(modifier, timeUnits, node, input%dt)

    call hsd_get_or_set(node, "Populations", input%tPopulations, .false.)
    call hsd_get_or_set(node, "WriteFrequency", input%writeFreq, 50)
    call hsd_get_or_set(node, "Restart", input%tReadRestart, .false.)
    if (input%tReadRestart) then
      call hsd_get_or_set(node, "RestartFromAscii", input%tReadRestartAscii, .false.)
    end if
    call hsd_get_or_set(node, "WriteRestart", input%tWriteRestart, .true.)
    if (input%tWriteRestart) then
      call hsd_get_or_set(node, "WriteAsciiRestart", input%tWriteRestartAscii, .false.)
    end if
    call hsd_get_or_set(node, "RestartFrequency", input%restartFreq, max(input%Steps / 10, 1))
    call hsd_get_or_set(node, "Forces", input%tForces, .false.)
    call hsd_get_or_set(node, "WriteBondEnergy", input%tBondE, .false.)
    call hsd_get_or_set(node, "WriteBondPopulation", input%tBondP, .false.)
    call hsd_get_or_set(node, "WriteAtomicEnergies", input%tWriteAtomEnergies, .false.)
    call hsd_get_or_set(node, "Pump", input%tPump, .false.)
    call hsd_get_or_set(node, "FillingsFromFile", input%tFillingsFromFile, .false.)

    if (input%tPump) then
      call hsd_get(node, "PumpProbeFrames", input%tdPPFrames, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'PumpProbeFrames'")
      defPpRange = [0.0_dp, input%steps * input%dt]
      block
        real(dp), allocatable :: tmpPpRange(:)
        call hsd_get_or_set(node, "PumpProbeRange", tmpPpRange, defPprange)
        input%tdPpRange = tmpPpRange
      end block
      call hsd_get_attrib(node, "PumpProbeRange", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, timeUnits, node, input%tdPpRange)

      ppRangeInvalid = (input%tdPpRange(2) <= input%tdPpRange(1))&
          & .or. (input%tdPprange(1) < defPpRange(1))&
          & .or. (input%tdPpRange(2) > defPpRange(2))
      if (ppRangeInvalid) then
        call dftbp_error(node, "Wrong definition of PumpProbeRange, either incorrect order&
            & or outside of simulation time range")
      end if
    end if

    call hsd_get_or_set(node, "Probe", input%tProbe, .false.)
    if (input%tPump .and. input%tProbe) then
      call dftbp_error(node, "Pump and probe cannot be simultaneously true.")
    end if

    call hsd_get_or_set(node, "EulerFrequency", input%eulerFreq, 0)

    call hsd_get_or_set(node, "VerboseDynamics", input%tVerboseDyn, .true.)

    if ((input%eulerFreq < 50) .and. (input%eulerFreq > 0)) then
      call dftbp_error(node, "Wrong number of Euler steps, should be above 50")
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
    call hsd_get_choice(node, "Perturbation", buffer, value1, stat)
    if (stat /= HSD_STAT_OK) then
      buffer = "none"
      value1 => null()
    else
      buffer = tolower(buffer)
    end if
    select case(buffer)

    case ("kick")
      input%pertType = pertTypes%kick
      call hsd_rename_child(value1, "PolarizationDirection", "PolarisationDirection")
      call hsd_get(value1, "PolarisationDirection", buffer2, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'PolarisationDirection'")
      input%polDir = directionConversion(unquote(buffer2), value1)

      call hsd_get_or_set(value1, "SpinType", buffer2, "Singlet")
      select case(unquote(buffer2))
      case ("singlet", "Singlet")
        input%spType = tdSpinTypes%singlet
      case ("triplet", "Triplet")
        input%spType = tdSpinTypes%triplet
      case default
        call dftbp_error(value1, "Unknown spectrum spin type " // buffer2)
      end select

      defaultWrite = .false.

    case ("laser")
      input%pertType = pertTypes%laser
      call hsd_rename_child(value1, "PolarizationDirection", "PolarisationDirection")
      block
        real(dp), allocatable :: tmpArr(:)
        call hsd_get(value1, "PolarisationDirection", tmpArr, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'PolarisationDirection'")
        input%reFieldPolVec = tmpArr
      end block
      call hsd_rename_child(value1, "ImagPolarizationDirection", "ImagPolarisationDirection")
      block
        real(dp), allocatable :: tmpArr(:)
        call hsd_get_or_set(value1, "ImagPolarisationDirection", tmpArr, &
            & [0.0_dp, 0.0_dp, 0.0_dp])
        input%imFieldPolVec = tmpArr
      end block
      call hsd_get(value1, "LaserEnergy", input%omega, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'LaserEnergy'")
      call hsd_get_attrib(value1, "LaserEnergy", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, energyUnits, value1, input%omega)
      call hsd_get_or_set(value1, "Phase", input%phase, 0.0_dp)
      call hsd_get_attrib(value1, "Phase", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, angularUnits, value1, input%phase)
      call hsd_get_or_set(value1, "ExcitedAtoms", buffer, "1:-1")
      call hsd_get_table(value1, "ExcitedAtoms", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) child => value1
      call getSelectedAtomIndices(child, buffer, geom%speciesNames, geom%species,&
          & input%indExcitedAtom)

      input%nExcitedAtom = size(input%indExcitedAtom)
      if (input%nExcitedAtom == 0) then
        call error("No atoms specified for laser excitation.")
      end if

      defaultWrite = .true.

    case ("kickandlaser")
      input%pertType = pertTypes%kickAndLaser
      call hsd_get(value1, "KickPolDir", buffer2, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'KickPolDir'")
      input%polDir = directionConversion(unquote(buffer2), value1)
      call hsd_get_or_set(value1, "SpinType", input%spType, tdSpinTypes%singlet)
      block
        real(dp), allocatable :: tmpArr(:)
        call hsd_get(value1, "LaserPolDir", tmpArr, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'LaserPolDir'")
        input%reFieldPolVec = tmpArr
      end block
      block
        real(dp), allocatable :: tmpArr(:)
        call hsd_get_or_set(value1, "LaserImagPolDir", tmpArr, [0.0_dp, 0.0_dp, 0.0_dp])
        input%imFieldPolVec = tmpArr
      end block
      call hsd_get(value1, "LaserEnergy", input%omega, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'LaserEnergy'")
      call hsd_get_attrib(value1, "LaserEnergy", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, energyUnits, value1, input%omega)
      call hsd_get_or_set(value1, "Phase", input%phase, 0.0_dp)
      call hsd_get_attrib(value1, "Phase", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, angularUnits, value1, input%phase)
      call hsd_get(value1, "LaserStrength", input%tdLaserField, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'LaserStrength'")
      call hsd_get_attrib(value1, "LaserStrength", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, EFieldUnits, value1, input%tdLaserField)

      call hsd_get_or_set(value1, "ExcitedAtoms", buffer, "1:-1")
      call hsd_get_table(value1, "ExcitedAtoms", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) child => value1
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
      call dftbp_error(node, "Unknown perturbation type " // buffer)
    end select

    if (tNeedFieldStrength) then
      call hsd_get(node, "FieldStrength", input%tdfield, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'FieldStrength'")
      call hsd_get_attrib(node, "FieldStrength", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, EFieldUnits, node, input%tdfield)
    end if

    call hsd_get_or_set(node, "WriteEnergyAndCharges", input%tdWriteExtras, defaultWrite)

    !! Different envelope functions
    call hsd_get_choice(node, "EnvelopeShape", buffer, value1, stat)
    if (stat /= HSD_STAT_OK) then
      buffer = "constant"
      value1 => null()
    else
      buffer = tolower(buffer)
    end if
    select case(buffer)

    case("constant")
      input%envType = envTypes%constant

    case("gaussian")
      input%envType = envTypes%gaussian
      call hsd_get_or_set(value1, "Time0", input%time0, 0.0_dp)
      call hsd_get_attrib(value1, "Time0", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, timeUnits, value1, input%Time0)

      call hsd_get(value1, "Time1", input%time1, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'Time1'")
      call hsd_get_attrib(value1, "Time1", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, timeUnits, value1, input%Time1)

    case("sin2")
      input%envType = envTypes%sin2
      call hsd_get_or_set(value1, "Time0", input%time0, 0.0_dp)
      call hsd_get_attrib(value1, "Time0", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, timeUnits, value1, input%Time0)

      call hsd_get(value1, "Time1", input%time1, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'Time1'")
      call hsd_get_attrib(value1, "Time1", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, timeUnits, value1, input%Time1)

    case("fromfile")
      input%envType = envTypes%fromFile
      call hsd_get_or_set(value1, "Time0", input%time0, 0.0_dp)
      call hsd_get_attrib(value1, "Time0", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call convertUnitHsd(modifier, timeUnits, value1, input%Time0)

    case default
      call dftbp_error(node, "Unknown envelope shape " // buffer)
    end select

    !! Non-adiabatic molecular dynamics
    call hsd_get_or_set(node, "IonDynamics", input%tIons, .false.)
    if (input%tIons) then
      call hsd_get_or_set(node, "MovedAtoms", buffer, "1:-1")
      call hsd_get_table(node, "MovedAtoms", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) child => node
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
      call dftbp_error(node, "Wrongly specified polarisation direction " // trim(direction)&
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
    integer :: nVelocities, nMatRows, nMatCols, stat
    real(dp), allocatable :: tmpVelocities(:,:)

    call hsd_get_table(node, "Velocities", child, stat, auto_wrap=.true.)
    modifier = ""
    buffer = ""
    if (associated(child)) then
      if (allocated(child%attrib)) modifier = child%attrib
      call hsd_get_choice(child, "", buffer, value1, stat)
      if (stat /= HSD_STAT_OK) buffer = ""
    end if
    if (buffer == "" .and. (.not. associated(child) .or. .not. hasInlineData(child))) then
       input%tReadMDVelocities = .false.
    else
       call hsd_get_matrix(node, "Velocities", tmpVelocities, nMatRows, nMatCols, &
           & stat=stat, order="column-major")
       if (stat /= 0 .or. nMatRows /= 3) then
          call dftbp_error(child, "Invalid velocity data (expected rows of 3 reals)")
       end if
       nVelocities = nMatCols
       if (nVelocities /= nAtom) then
          call dftbp_error(node, "Incorrect number of specified velocities: " &
               & // i2c(3*nVelocities) // " supplied, " &
               & // i2c(3*nAtom) // " required.")
       end if
       if (len(modifier) > 0) then
          call convertUnitHsd(modifier, VelocityUnits, child, &
               & tmpVelocities)
       end if
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
    integer :: stat

    allocate(thermostatInp, tempProfileInp)
    call hsd_get_choice(node, "Thermostat", thermName, thermNode, stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required block: 'Thermostat'")
    thermName = tolower(thermName)

    select case(thermName)

    case ("berendsen")

      thermostatInp%thermostatType = thermostatTypes%berendsen
      allocate(thermostatInp%berendsen)
      associate (inp => thermostatInp%berendsen)
        call readTempOrTempProfile_(thermNode, maxRun, tempProfileInp)
        call hsd_get_table(thermNode, "CouplingStrength", child2, stat, auto_wrap=.true.)
        if (associated(child2)) then
          call hsd_get(child2, "#text", inp%coupling, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value")
          call hsd_get_table(thermNode, "Timescale", child2, stat, auto_wrap=.true.)
          if (associated(child2)) then
            call error("Only Coupling strength OR Timescale can be set for Berendsen thermostats.")
          end if
        else
          call hsd_get_table(thermNode, "Timescale", child2, stat, auto_wrap=.true.)
          if (associated(child2)) then
            call hsd_get(child2, "#text", inp%coupling, stat=stat)
            if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value")
            if (allocated(child2%attrib)) then
              modifier = child2%attrib
            else
              modifier = ""
            end if
            call convertUnitHsd(modifier, timeUnits, child2, inp%coupling)
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
        call hsd_get(thermNode, "CouplingStrength", inp%coupling, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(thermNode, "Missing required value: 'CouplingStrength'")
        call hsd_get_attrib(thermNode, "CouplingStrength", modifier, stat)
        if (stat /= HSD_STAT_OK) modifier = ""
        call convertUnitHsd(modifier, freqUnits, thermNode, inp%coupling)

        call hsd_get_or_set(thermNode, "ChainLength", inp%chainLength, 3)
        call hsd_get_or_set(thermNode, "Order", inp%expOrder, 3, child=child2)
        if (.not. any(inp%expOrder == [3, 5])) then
          call dftbp_error(child2, "Order of Nose-Hoover thermostat must be either 3 or 5")
        end if
        call hsd_get_or_set(thermNode, "IntegratorSteps", inp%nExpSteps, 1)
        call hsd_get_table(thermNode, "Restart", child2, stat, auto_wrap=.true.)
        if (associated(child2)) then
          allocate(inp%xnose(inp%chainLength))
          allocate(inp%vnose(inp%chainLength))
          allocate(inp%gnose(inp%chainLength))
          call hsd_get(child2, "x", inp%xnose, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value: 'x'")
          call hsd_get(child2, "v", inp%vnose, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value: 'v'")
          call hsd_get(child2, "g", inp%gnose, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value: 'g'")
        end if
      end associate

    case ("andersen")

      thermostatInp%thermostatType = thermostatTypes%andersen
      allocate(thermostatInp%andersen)
      associate (inp => thermostatInp%andersen)
        call readTempOrTempProfile_(thermNode, maxRun, tempProfileInp)
        call hsd_get(thermNode, "ReselectProbability", inp%rescaleProb, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(thermNode, "Missing required value: 'ReselectProbability'")
        if (inp%rescaleProb <= 0.0_dp .or. inp%rescaleProb > 1.0_dp) then
          call dftbp_error(thermNode, "ReselectProbability must be in the range (0,1]!")
        end if
        call hsd_get(thermNode, "ReselectIndividually", inp%rescaleIndiv, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(thermNode, "Missing required value: 'ReselectIndividually'")
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
      call getNodeName2(thermNode, thermName)
      call dftbp_error(node, "Invalid thermostat '" // thermName // "'")

    end select

  contains

    !> Reads the temperature or the temperature profile
    subroutine readTempOrTempProfile_(thermNode, maxRun, tempProfileInp)
      type(hsd_table), pointer, intent(in) :: thermNode
      integer, intent(inout) :: maxRun
      type(TTempProfileInput), intent(out) :: tempProfileInp

      type(hsd_table), pointer :: value, child
      character(len=:), allocatable :: buffer, modifier
      integer :: stat

      call hsd_get_table(thermNode, "Temperature", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) call dftbp_error(thermNode, "Missing required value: 'Temperature'")
      call hsd_get_attrib(thermNode, "Temperature", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      ! Determine if this is a dispatch (table child) or plain value
      value => null()
      buffer = "#text"
      block
        class(hsd_node), pointer :: tmpNode
        integer :: ic
        do ic = 1, child%num_children
          call child%get_child(ic, tmpNode)
          if (.not. associated(tmpNode)) cycle
          select type(t => tmpNode)
          type is (hsd_table)
            value => t
            call getNodeName(value, buffer)
            exit
          end select
        end do
      end block
      select case(buffer)
      case ("#text")
        call readTemperature(child, tempProfileInp)
      case ("temperatureprofile")
        call readTemperatureProfile(value, modifier, tempProfileInp)
        maxRun = sum(tempProfileInp%tempInts) - 1
      case default
        call dftbp_error(child, "Invalid method name.")
      end select

    end subroutine readTempOrTempProfile_

  end subroutine parseThermostat


end module dftbp_dftbplus_parser_driver
