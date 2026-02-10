!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Fills the derived type with the input parameters from an HSD or an XML file.
module dftbp_solvation_solvparser
  use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_support_inf, ieee_value
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_atomicrad, only : getAtomicRad
  use dftbp_common_constants, only : AA__Bohr, amu__au, Boltzmann, kg__au
  use dftbp_common_filesystem, only : findFile, getParamSearchPaths
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_unitconversion, only : energyUnits, inverseLengthUnits, lengthUnits,&
      & massDensityUnits, massUnits
  use dftbp_dftbplus_specieslist, only : readSpeciesList
  use dftbp_extlibs_lebedev, only : gridSize
  use dftbp_io_charmanip, only : tolower, unquote
  use hsd, only : hsd_rename_child, hsd_get, hsd_get_or_set, hsd_get_table, hsd_get_choice,&
      & HSD_STAT_OK, hsd_schema_t, hsd_error_t, schema_init, schema_add_field, &
      & schema_validate, schema_destroy, FIELD_OPTIONAL, FIELD_TYPE_REAL, FIELD_TYPE_TABLE
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, getNodeName
  use dftbp_io_unitconv, only : convertUnitHsd
  use hsd_data, only : hsd_table, new_table
  use dftbp_math_bisect, only : bisection
  use dftbp_solvation_born, only : fgbKernel, TGBInput
  use dftbp_solvation_cm5, only : TCM5Input
  use dftbp_solvation_cosmo, only : TCosmoInput, TDomainDecompositionInput
  use dftbp_solvation_gbsafile, only : readParamGBSA
  use dftbp_solvation_sasa, only : TSASAInput
  use dftbp_solvation_solvdata, only : getVanDerWaalsRadiusBondi, getVanDerWaalsRadiusCosmo,&
      & getVanDerWaalsRadiusD3
  use dftbp_solvation_solventdata, only : SolventFromName, TSolventData
  use dftbp_solvation_solvinput, only : TSolvationInp
  use dftbp_type_typegeometry, only : TGeometry
  implicit none

  private
  public :: readSolvation
  public :: readSolvGB, readSolvSASA, readCM5, readSolvCosmo


  real(dp), parameter :: ambientTemperature = 298.15_dp * Boltzmann


contains


  !> Reads in solvation related settings
  subroutine readSolvation(node, geo, input)

    !> Node to parse
    type(hsd_table), pointer :: node

    !> Geometry, including atomic information
    type(TGeometry), intent(in) :: geo

    !> Solvation data on exit
    type(TSolvationInp), intent(out) :: input

    type(hsd_table), pointer :: solvModel
    character(len=:), allocatable :: buffer
    integer :: stat

    call hsd_rename_child(node, "GeneralizedBorn", "GeneralisedBorn")
    call hsd_get_choice(node, "", buffer, solvModel, stat)
    if (.not. associated(solvModel)) call dftbp_error(node, "Missing required solvation model")
    call getNodeName(solvModel, buffer)

    select case (buffer)
    case default
      call dftbp_error(node, "Invalid solvation model name.")
    case ("generalisedborn")
      allocate(input%GBInp)
      call readSolvGB(solvModel, geo, input%GBInp)
    case ("cosmo")
      allocate(input%cosmoInp)
      call readSolvCosmo(solvModel, geo, input%cosmoInp)
    case ("sasa")
      allocate(input%SASAInp)
      call readSolvSASA(solvModel, geo, input%SASAInp)
    end select
  end subroutine readSolvation


  !> Reads in generalized Born related settings
  subroutine readSolvGB(node, geo, input)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Contains the input for the solvation module on exit
    type(TGBInput), intent(out) :: input

    type(TGBInput), allocatable :: defaults
    character(len=:), allocatable :: buffer, modifier
    type(hsd_table), pointer :: child, value1, field
    logical :: tHBondCorr, tALPB
    real(dp) :: temperature, shift, alphaALPB
    character(len=:), allocatable :: searchPath(:)
    type(TSolventData) :: solvent
    real(dp), parameter :: alphaDefault = 0.571412_dp
    character(len=:), allocatable :: paramFile, paramTmp
    integer :: stat

    if (geo%tPeriodic .or. geo%tHelical) then
      call dftbp_error(node, "Generalised Born model currently not available with the selected&
          & boundary conditions")
    end if

    call hsd_get_table(node, "ParamFile", value1, stat, auto_wrap=.true.)
    if (associated(value1)) then
      allocate(defaults)
      child => value1
      call hsd_get_or_set(node, "ParamFile", buffer, "")
      paramFile = trim(unquote(buffer))
      call getParamSearchPaths(searchPath)
      call findFile(searchPath, paramFile, paramTmp)
      if (allocated(paramTmp)) call move_alloc(paramTmp, paramFile)
      write(stdOut, '(a)') "Reading GBSA parameter file '" // paramFile // "'"
      call readParamGBSA(paramFile, defaults, solvent, geo%speciesNames, node=child)
    else
      call readSolvent(node, solvent)
    end if

    call hsd_get_or_set(node, "ALPB", tALPB, .false.)
    if (tALPB) then
      call hsd_get_or_set(node, "Alpha", alphaALPB, alphaDefault)
      input%alpbet = alphaALPB / solvent%dielectricConstant
    else
      input%alpbet = 0.0_dp
    end if
    input%keps = (1.0_dp / solvent%dielectricConstant - 1.0_dp) / (1.0_dp + input%alpbet)
    input%dielectricConstant = solvent%dielectricConstant

    call hsd_get_or_set(node, "Kernel", buffer, "Still")
    call hsd_get_table(node, "Kernel", child, stat, auto_wrap=.true.)
    select case(tolower(unquote(buffer)))
    case default
      call dftbp_error(child, "Unknown interaction kernel: "//buffer)
    case("still")
      input%kernel = fgbKernel%still
    case("p16")
      input%kernel = fgbKernel%p16
    end select

    ! shift value for the free energy (usually fitted)
    if (allocated(defaults)) then
      call hsd_get_or_set(node, "FreeEnergyShift", shift, defaults%freeEnergyShift)
    else
      call hsd_get(node, "FreeEnergyShift", shift, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'FreeEnergyShift'")
    end if
    call hsd_get_table(node, "FreeEnergyShift", field, stat, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
    call convertUnitHsd(modifier, energyUnits, field, shift)

    ! temperature, influence depends on the reference state
    call hsd_get_or_set(node, "Temperature", temperature, ambientTemperature)
    call hsd_get_table(node, "Temperature", field, stat, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
    call convertUnitHsd(modifier, energyUnits, field, temperature)

    ! reference state for free energy calculation
    call readReferenceState(node, solvent, temperature, shift, input%freeEnergyShift)

    if (allocated(defaults)) then
      call hsd_get_or_set(node, "BornScale", input%bornScale, defaults%bornScale)
      call hsd_get_or_set(node, "BornOffset", input%bornOffset, defaults%bornOffset)
    else
      call hsd_get(node, "BornScale", input%bornScale, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'BornScale'")
      call hsd_get(node, "BornOffset", input%bornOffset, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'BornOffset'")
    end if
    call hsd_get_table(node, "BornOffset", field, stat, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
    call convertUnitHsd(modifier, lengthUnits, field, input%bornOffset)
    block
      real(dp), allocatable :: tmpObc(:)
      call hsd_get_or_set(node, "OBCCorrection", tmpObc, [1.00_dp, 0.80_dp, 4.85_dp])
      input%obc = tmpObc
    end block

    call hsd_get_table(node, "CM5", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      allocate(input%cm5Input)
      call readCM5(child, input%cm5Input, geo)
    end if

    call readVanDerWaalsRad(node, geo, input%vdwRad)

    allocate(input%descreening(geo%nSpecies))
    if (allocated(defaults)) then
      call hsd_get_table(node, "Descreening", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) then
        block
          type(hsd_table) :: tmpTbl, defChild
          call new_table(tmpTbl, name="descreening")
          call new_table(defChild, name="defaults")
          call tmpTbl%add_child(defChild)
          call node%add_child(tmpTbl)
        end block
        call hsd_get_table(node, "Descreening", child, stat, auto_wrap=.true.)
      end if
    else
      call hsd_get_table(node, "Descreening", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) call dftbp_error(node,&
          & "Missing required block: 'Descreening'")
    end if
    call hsd_get_choice(child, "", buffer, value1, stat)
    call getNodeName(value1, buffer)
    select case(buffer)
    case default
      call dftbp_error(child, "Unknown method '"//buffer//&
          & "' to generate descreening parameters")
    case("defaults")
      if (.not.allocated(defaults)) then
        call dftbp_error(child, "No defaults available for descreening parameters")
      end if
      call readSpeciesList(value1, geo%speciesNames, input%descreening,&
          & default=defaults%descreening)
    case("unity")
      input%descreening(:) = 1.0_dp
    case("values")
      call readSpeciesList(value1, geo%speciesNames, input%descreening)
    end select

    call hsd_get_or_set(node, "Cutoff", input%rCutoff, 35.0_dp * AA__Bohr)
    call hsd_get_table(node, "Cutoff", field, stat, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
    call convertUnitHsd(modifier, lengthUnits, field, input%rCutoff)

    call hsd_get_table(node, "SASA", value1, stat, auto_wrap=.true.)
    if (associated(value1) .or. allocated(defaults)) then
      allocate(input%sasaInput)
      if (.not.associated(value1)) then
        block
          type(hsd_table) :: tmpTbl
          call new_table(tmpTbl, name="sasa")
          call node%add_child(tmpTbl)
        end block
        call hsd_get_table(node, "SASA", value1, stat)
      end if
      if (allocated(defaults)) then
        call readSolvSASA(value1, geo, input%sasaInput, defaults%sasaInput%probeRad,&
            & defaults%sasaInput%surfaceTension)
      else
        call readSolvSASA(value1, geo, input%sasaInput)
      end if

      if (allocated(defaults)) then
        call hsd_get_or_set(node, "HBondCorr", tHBondCorr, allocated(defaults%hBondPar))
      else
        call hsd_get(node, "HBondCorr", tHBondCorr, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'HBondCorr'")
      end if

      if (tHBondCorr) then
        allocate(input%hBondPar(geo%nSpecies))
        if (allocated(defaults)) then
          call hsd_get_table(node, "HBondStrength", child, stat, auto_wrap=.true.)
          if (.not. associated(child)) then
            block
              type(hsd_table) :: tmpTbl, defChild
              call new_table(tmpTbl, name="hbondstrength")
              call new_table(defChild, name="defaults")
              call tmpTbl%add_child(defChild)
              call node%add_child(tmpTbl)
            end block
            call hsd_get_table(node, "HBondStrength", child, stat, auto_wrap=.true.)
          end if
        else
          call hsd_get_table(node, "HBondStrength", child, stat, auto_wrap=.true.)
          if (.not. associated(child)) call dftbp_error(node,&
              & "Missing required block: 'HBondStrength'")
        end if
        call hsd_get_choice(child, "", buffer, value1, stat)
        call getNodeName(value1, buffer)
        select case(buffer)
        case default
          call dftbp_error(child, "Unknown method '"//buffer//&
              & "' to generate H-bond parameters")
        case("defaults")
          if (allocated(defaults)) then
            if (.not.allocated(defaults%hBondPar)) then
              call dftbp_error(child, "No defaults available for hydrogen bond strengths")
            end if
          else
            call dftbp_error(child, "No defaults available for hydrogen bond strengths")
          end if
          call readSpeciesList(value1, geo%speciesNames, input%hBondPar, default=defaults%hBondPar)
        case("values")
          call readSpeciesList(value1, geo%speciesNames, input%hBondPar)
        end select
      end if
    end if

  end subroutine readSolvGB


  !> Reads in conductor like screening model settings
  subroutine readSolvCosmo(node, geo, input)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Contains the input for the solvation module on exit
    type(TCosmoInput), intent(out) :: input

    character(len=:), allocatable :: buffer, modifier
    type(hsd_table), pointer :: child, value1, field
    real(dp) :: temperature, shift
    real(dp), allocatable :: radScale(:), radScaleSpecies(:), radScaleDefault(:)
    type(TSolventData) :: solvent
    integer :: stat

    if (geo%tPeriodic .or. geo%tHelical) then
      call dftbp_error(node, "COSMO solvation currently not available with the selected boundary&
          & conditions")
    end if

    call readSolvent(node, solvent)
    input%dielectricConst = solvent%dielectricConstant
    input%keps = 0.5_dp * (1.0_dp - 1.0_dp/solvent%dielectricConstant)

    ! shift value for the free energy (usually zero)
    call hsd_get_or_set(node, "FreeEnergyShift", shift, 0.0_dp)
    call hsd_get_table(node, "FreeEnergyShift", field, stat, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
    call convertUnitHsd(modifier, energyUnits, field, shift)

    ! temperature, influence depends on the reference state
    call hsd_get_or_set(node, "Temperature", temperature, ambientTemperature)
    call hsd_get_table(node, "Temperature", field, stat, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
    call convertUnitHsd(modifier, energyUnits, field, temperature)

    call readReferenceState(node, solvent, temperature, shift, input%freeEnergyShift)

    call readVanDerWaalsRad(node, geo, input%vdwRad)

    call hsd_get_table(node, "RadiiScaling", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      allocate(radScaleSpecies(geo%nSpecies))
      allocate(radScaleDefault(geo%nSpecies), source=1.0_dp)
      call readSpeciesList(child, geo%speciesNames, radScaleSpecies, default=radScaleDefault)
      deallocate(radScaleDefault)
      input%vdwRad(:) = input%vdwRad * radScaleSpecies
      deallocate(radScaleSpecies)
    end if

    call readAngularGrid(node, input%gridSize)

    call hsd_get_table(node, "Solver", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: tmpTbl, defChild
        call new_table(tmpTbl, name="solver")
        call new_table(defChild, name="domaindecomposition")
        call tmpTbl%add_child(defChild)
        call node%add_child(tmpTbl)
      end block
      call hsd_get_table(node, "Solver", child, stat, auto_wrap=.true.)
    end if
    call hsd_get_choice(child, "", buffer, value1, stat)
    call getNodeName(value1, buffer)
    select case(buffer)
    case default
      call dftbp_error(child, "Unknown method '"//buffer//"' to solve COSMO equation")
    case("domaindecomposition")
      call readDomainDecomposition(value1, input%ddInput)
    end select

    call hsd_get_table(node, "SASA", value1, stat, auto_wrap=.true.)
    if (associated(value1)) then
      allocate(input%sasaInput)
      call readSolvSASA(value1, geo, input%sasaInput)
    end if

  end subroutine readSolvCosmo


  !> Read domain settings for the COSMO solver
  subroutine readDomainDecomposition(node, input)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Input for the domain decomposition algorithm
    type(TDomainDecompositionInput), intent(out) :: input

    type(hsd_table), pointer :: child
    integer :: stat

    call hsd_get(node, "MaxMoment", input%lmax, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'MaxMoment'")
    call hsd_rename_child(node, "Regularization", "Regularisation")
    call hsd_get_or_set(node, "Regularisation", input%eta, 0.2_dp)
    call hsd_get(node, "Accuracy", input%conv, stat=stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'Accuracy'")

  end subroutine readDomainDecomposition


  !> Read input data for non-polar surface area solvation model.
  subroutine readSolvSASA(node, geo, input, probeRadDefault, surfaceTensionDefault)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Contains the input for the solvation module on exit
    type(TSASAInput), intent(out) :: input

    !> Default value for the probe radius
    real(dp), intent(in), optional :: probeRadDefault

    !> Default values for the surface tension
    real(dp), intent(in), optional :: surfaceTensionDefault(:)

    character(len=:), allocatable :: buffer, modifier
    type(hsd_table), pointer :: child, value1, field
    integer :: stat

    if (geo%tPeriodic .or. geo%tHelical) then
      call dftbp_error(node, "SASA model currently not available with the selected boundary&
          & conditions")
    end if

    if (present(probeRadDefault)) then
      call hsd_get_or_set(node, "ProbeRadius", input%probeRad, probeRadDefault)
    else
      call hsd_get(node, "ProbeRadius", input%probeRad, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'ProbeRadius'")
    end if
    call hsd_get_table(node, "ProbeRadius", field, stat, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
    call convertUnitHsd(modifier, lengthUnits, field, input%probeRad)

    call hsd_get_or_set(node, "Smoothing", input%smoothingPar, 0.3_dp*AA__Bohr)
    call hsd_get_table(node, "Smoothing", field, stat, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
    call convertUnitHsd(modifier, lengthUnits, field, input%smoothingPar)

    call hsd_get_or_set(node, "Tolerance", input%tolerance, 1.0e-6_dp)

    call readAngularGrid(node, input%gridSize, 230)

    call readVanDerWaalsRad(node, geo, input%vdwRad)

    allocate(input%surfaceTension(geo%nSpecies))
    if (present(surfaceTensionDefault)) then
      call hsd_get_table(node, "SurfaceTension", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) then
        block
          type(hsd_table) :: tmpTbl, defChild
          call new_table(tmpTbl, name="surfacetension")
          call new_table(defChild, name="defaults")
          call tmpTbl%add_child(defChild)
          call node%add_child(tmpTbl)
        end block
        call hsd_get_table(node, "SurfaceTension", child, stat, auto_wrap=.true.)
      end if
    else
      call hsd_get_table(node, "SurfaceTension", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) call dftbp_error(node,&
          & "Missing required block: 'SurfaceTension'")
    end if
    call hsd_get_choice(child, "", buffer, value1, stat)
    call getNodeName(value1, buffer)
    select case(buffer)
    case default
      call dftbp_error(child, "Unknown method '"//buffer//"' to generate surface tension")
    case("defaults")
      if (.not.present(surfaceTensionDefault)) then
        call dftbp_error(child, "No defaults available for surface tension values")
      end if
      call readSpeciesList(value1, geo%speciesNames, input%surfaceTension,&
          & default=surfaceTensionDefault)
    case("values")
      call readSpeciesList(value1, geo%speciesNames, input%surfaceTension)
    end select

    call hsd_get_or_set(node, "Offset", input%sOffset, 2.0_dp * AA__Bohr)
    call hsd_get_table(node, "Offset", field, stat, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
    call convertUnitHsd(modifier, lengthUnits, field, input%sOffset)

  end subroutine readSolvSASA


  !> Read settings for charge model 5.
  subroutine readCM5(node, input, geo)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Contains the input for the CM5 module on exit
    type(TCM5Input), intent(out) :: input

    type(hsd_table), pointer :: value1, dummy, child, field
    character(len=:), allocatable :: buffer, modifier
    real(dp), allocatable :: atomicRadDefault(:)
    integer :: stat

    call hsd_get_or_set(node, "Alpha", input%alpha, 2.474_dp/AA__Bohr)
    call hsd_get_table(node, "Alpha", field, stat, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
    call convertUnitHsd(modifier, inverseLengthUnits, field, input%alpha)

    allocate(input%atomicRad(geo%nSpecies))
    call hsd_get_table(node, "Radii", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: tmpTbl, defChild
        call new_table(tmpTbl, name="radii")
        call new_table(defChild, name="atomicradii")
        call tmpTbl%add_child(defChild)
        call node%add_child(tmpTbl)
      end block
      call hsd_get_table(node, "Radii", child, stat, auto_wrap=.true.)
    end if
    call hsd_get_choice(child, "", buffer, value1, stat)
    call getNodeName(value1, buffer)
    select case(buffer)
    case default
      call dftbp_error(child, "Unknown method '"//buffer//"' to generate radii")
    case("atomicradii")
      allocate(atomicRadDefault(geo%nSpecies))
      atomicRadDefault(:) = getAtomicRad(geo%speciesNames)
      call readSpeciesList(value1, geo%speciesNames, input%atomicRad, default=atomicRadDefault,&
          & units=lengthUnits)
      deallocate(atomicRadDefault)
    case("values")
      call readSpeciesList(value1, geo%speciesNames, input%atomicRad, units=lengthUnits)
    end select
    if (any(input%atomicRad <= 0.0_dp)) then
      call dftbp_error(value1, "Atomic radii must be positive for all species")
    end if

    call hsd_get_or_set(node, "Cutoff", input%rCutoff, 30.0_dp)
    call hsd_get_table(node, "Cutoff", field, stat, auto_wrap=.true.)
    if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
    call convertUnitHsd(modifier, lengthUnits, field, input%rCutoff)

    ! -- Schema validation (warnings only) --
    block
      type(hsd_schema_t) :: schema
      type(hsd_error_t), allocatable :: schemaErrors(:)
      integer :: iErr

      call schema_init(schema, name="CM5")
      call schema_add_field(schema, "Alpha", FIELD_OPTIONAL, FIELD_TYPE_REAL)
      call schema_add_field(schema, "Radii", FIELD_OPTIONAL, FIELD_TYPE_TABLE)
      call schema_add_field(schema, "Cutoff", FIELD_OPTIONAL, FIELD_TYPE_REAL)
      call schema_validate(schema, node, schemaErrors)
      if (size(schemaErrors) > 0) then
        do iErr = 1, size(schemaErrors)
          call dftbp_warning(node, "[schema] " // schemaErrors(iErr)%message)
        end do
      end if
      call schema_destroy(schema)
    end block

  end subroutine readCM5


  !> Read the input data for the solvent model
  subroutine readSolvent(node, solvent)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Data associated with the solvent
    type(TSolventData), intent(out) :: solvent

    character(len=:), allocatable :: buffer, modifier
    type(hsd_table), pointer :: child, value1, field
    logical :: found
    integer :: stat

    call hsd_get_table(node, "Solvent", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) call dftbp_error(node, "Missing required block: 'Solvent'")
    call hsd_get_choice(child, "", buffer, value1, stat)
    call getNodeName(value1, buffer)
    select case(buffer)
    case default
      call dftbp_error(child, "Invalid solvent method '" // buffer // "'")
    case('fromname')
      call hsd_get(value1, "#text", buffer, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value")
      call SolventFromName(solvent, unquote(buffer), found)
      if (.not. found) then
        call dftbp_error(value1, "Invalid solvent " // buffer)
      end if
    case('fromconstants')
      call hsd_get(value1, "Epsilon", buffer, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'Epsilon'")
      if (unquote(buffer) == "Inf") then
         if (ieee_support_inf(solvent%dielectricConstant)) then
            solvent%dielectricConstant = ieee_value(solvent%dielectricConstant, ieee_positive_inf)
         else
            solvent%dielectricConstant = huge(solvent%dielectricConstant)
         end if
      else
         call hsd_get(value1, "Epsilon", solvent%dielectricConstant, stat=stat)
         if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'Epsilon'")
      end if
      call hsd_get(value1, "MolecularMass", solvent%molecularMass, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'MolecularMass'")
      call hsd_get_table(value1, "MolecularMass", field, stat, auto_wrap=.true.)
      if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
      call convertUnitHsd(modifier, massUnits, field, solvent%molecularMass)
      call hsd_get(value1, "Density", solvent%density, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(value1, "Missing required value: 'Density'")
      call hsd_get_table(value1, "Density", field, stat, auto_wrap=.true.)
      if (allocated(field%attrib)) then; modifier = field%attrib; else; modifier = ""; end if
      call convertUnitHsd(modifier, massDensityUnits, field, solvent%density)
    end select

  end subroutine readSolvent


  !> Reference state for free energy calculation
  subroutine readReferenceState(node, solvent, temperature, shift, freeEnergyShift)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Data associated with the solvent
    type(TSolventData), intent(in) :: solvent

    !> Temperature for calculation
    real(dp), intent(in) :: temperature

    !> Shift to free energy
    real(dp), intent(in) :: shift

    !> Free energy shift includings state specific terms
    real(dp), intent(out) :: freeEnergyShift

    character(len=:), allocatable :: state
    type(hsd_table), pointer :: child
    real(dp), parameter :: referenceDensity = kg__au/(1.0e10_dp*AA__Bohr)**3
    real(dp), parameter :: referenceMolecularMass = amu__au
    real(dp), parameter :: idealGasMolVolume = 24.79_dp
    integer :: stat

    call hsd_get_or_set(node, "State", state, "gsolv")
    call hsd_get_table(node, "State", child, stat, auto_wrap=.true.)
    select case(tolower(unquote(state)))
    case default
      call dftbp_error(child, "Unknown reference state: "//state)
    case("gsolv") ! just the bare shift
      freeEnergyShift = shift
    case("reference") ! gsolv=reference option in cosmotherm
      ! RT * ln(ideal gas mol volume) + ln(rho/M)
      freeEnergyShift = shift + temperature&
          & * (log(idealGasMolVolume * temperature / ambientTemperature)&
          & + log(solvent%density/referenceDensity * referenceMolecularMass/solvent%molecularMass))
    case("mol1bar")
      ! RT * ln(ideal gas mol volume)
      freeEnergyShift = shift + temperature&
          & * log(idealGasMolVolume * temperature / ambientTemperature)
    end select
  end subroutine readReferenceState


  !> Read the atomic radii from a choice of various sources
  subroutine readVanDerWaalsRad(node, geo, vdwRad)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Van-der-Waals Radii
    real(dp), allocatable, intent(out) :: vdwRad(:)

    character(len=:), allocatable :: buffer
    type(hsd_table), pointer :: child, value1, dummy
    real(dp), allocatable :: vdwRadDefault(:)
    integer :: stat

    allocate(vdwRad(geo%nSpecies))
    call hsd_get_table(node, "Radii", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      block
        type(hsd_table) :: tmpTbl, defChild
        call new_table(tmpTbl, name="radii")
        call new_table(defChild, name="vanderwaalsradiid3")
        call tmpTbl%add_child(defChild)
        call node%add_child(tmpTbl)
      end block
      call hsd_get_table(node, "Radii", child, stat, auto_wrap=.true.)
    end if
    call hsd_get_choice(child, "", buffer, value1, stat)
    call getNodeName(value1, buffer)
    select case(buffer)
    case default
      call dftbp_error(child, "Unknown method '"//buffer//"' to generate radii")
    case("vanderwaalsradiid3")
      allocate(vdwRadDefault(geo%nSpecies))
      vdwRadDefault(:) = getVanDerWaalsRadiusD3(geo%speciesNames)
      call readSpeciesList(value1, geo%speciesNames, vdwRad, default=vdwRadDefault,&
          & units=lengthUnits)
      deallocate(vdwRadDefault)
    case("vanderwaalsradiicosmo")
      allocate(vdwRadDefault(geo%nSpecies))
      vdwRadDefault(:) = getVanDerWaalsRadiusCosmo(geo%speciesNames)
      call readSpeciesList(value1, geo%speciesNames, vdwRad, default=vdwRadDefault,&
          & units=lengthUnits)
      deallocate(vdwRadDefault)
    case("vanderwaalsradiibondi")
      allocate(vdwRadDefault(geo%nSpecies))
      vdwRadDefault(:) = getVanDerWaalsRadiusBondi(geo%speciesNames)
      call readSpeciesList(value1, geo%speciesNames, vdwRad, default=vdwRadDefault,&
          & units=lengthUnits)
      deallocate(vdwRadDefault)
    case("values")
      call readSpeciesList(value1, geo%speciesNames, vdwRad, units=lengthUnits)
    end select

  end subroutine readVanDerWaalsRad


  !> Reads settings for angular integration grid
  subroutine readAngularGrid(node, angGrid, default)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Grid identifier
    integer, intent(out) :: angGrid

    !> Default grid size
    integer, intent(in), optional :: default

    type(hsd_table), pointer :: child
    character(lc) :: errorStr
    integer :: gridPoints
    integer :: stat

    if (present(default)) then
      call hsd_get_or_set(node, "AngularGrid", gridPoints, default)
    else
      call hsd_get(node, "AngularGrid", gridPoints, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required value: 'AngularGrid'")
    end if
    call hsd_get_table(node, "AngularGrid", child, stat, auto_wrap=.true.)
    angGrid = 0
    call bisection(angGrid, gridSize, gridPoints)
    if (angGrid == 0) then
      call dftbp_error(child, "Illegal number of grid points for numerical integration")
    end if
    if (gridSize(angGrid) /= gridPoints) then
      write(errorStr, '(a, *(1x, i0, 1x, a))')&
          & "No angular integration grid with", gridPoints, "points available, using",&
          &  gridSize(angGrid), "points instead"
      call dftbp_warning(child, trim(errorStr))
    end if

  end subroutine readAngularGrid


end module dftbp_solvation_solvparser
