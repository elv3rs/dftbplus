!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Hamiltonian-related parsing subroutines extracted from the main parser module.
module dftbp_dftbplus_parser_hamiltonian
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_filesystem, only : findFile, getParamSearchPaths, joinPathsPrettyErr
  use dftbp_common_globalenv, only : abortProgram, stdout
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_status, only : TStatus
  use dftbp_common_unitconversion, only : energyUnits
  use dftbp_dftb_dftbplusu, only : plusUFunctionals
  use dftbp_dftb_elecconstraints, only : readElecConstraintInput
  use dftbp_dftb_halogenx, only : halogenXSpecies1, halogenXSpecies2
  use dftbp_dftb_hybridxc, only : THybridXcSKTag
  use dftbp_dftb_slakoeqgrid, only : skEqGridNew, skEqGridOld
  use dftbp_dftbplus_forcetypes, only : forceTypes
  use dftbp_dftbplus_inputdata, only : TControl, TSlater
  use dftbp_elecsolvers_elecsolvers, only : electronicSolverTypes
  use dftbp_extlibs_poisson, only : TPoissonInfo
  use dftbp_extlibs_tblite, only : tbliteMethod
  use dftbp_io_charmanip, only : newline, tolower, unquote
  use dftbp_io_hsdutils, only : hsd_child_list, &
      & getChild, getChildren, getChildValue, &
      & setChild, setChildValue, &
      & getLength, getItem1, destroyNodeList, removeChild
  use hsd, only : hsd_get, hsd_get_matrix, hsd_get_or_set
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning,&
      & textNodeName, getNodeName, getNodeHSDName, getNodeName2, hasInlineData
  use dftbp_io_unitconv, only : convertUnitHsd
  use hsd_data, only : hsd_table
  use dftbp_io_message, only : error, warning
  use dftbp_md_thermostats, only : thermostatTypes
  use dftbp_reks_reks, only : reksTypes
  use dftbp_solvation_solvparser, only : readSolvation
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_linkedlist, only : append, destruct, init, TListCharLc
  use dftbp_type_typegeometry, only : TGeometry
#:if WITH_TRANSPORT
  use dftbp_transport_negfvars, only : TNEGFGreenDensInfo
#:endif
  use dftbp_transport_negfvars, only : TTransPar
  use dftbp_dftbplus_parser_dispersion, only : readDispersion
  use dftbp_dftbplus_parser_kpoints, only : readKPoints
  use dftbp_dftbplus_parser_hybrid, only : parseHybridBlock, parseChimes
  use dftbp_dftbplus_parser_skfiles, only : readSKFiles
  use dftbp_dftbplus_parser_spin, only : readSpinPolarisation, readSpinOrbit, &
      & readMaxAngularMomentum, setupOrbitals, TAngShellBlocks
  use dftbp_dftbplus_parser_electrostatics, only : readElectrostatics, readMdftb
  use dftbp_dftbplus_parser_filling, only : readElectronicFilling
  use dftbp_dftbplus_parser_sccoptions, only : readSccOptions, readForceOptions, SKTruncations,&
      & readHCorrection, readDifferentiation
  use dftbp_dftbplus_parser_customisation, only : readCustomisedHubbards, readCustomReferenceOcc
  use dftbp_dftbplus_parser_solver, only : readSolver
  use dftbp_dftbplus_parser_external, only : readExternal
  implicit none

  private
  public :: readHamiltonian, readDFTBHam, readXTBHam, readLaterHamiltonian


contains

  !> Reads Hamiltonian
#:if WITH_TRANSPORT
  subroutine readHamiltonian(node, ctrl, geo, slako, tp, greendens, poisson, errStatus)
#:else
  subroutine readHamiltonian(node, ctrl, geo, slako, poisson, errStatus)
#:endif

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    !> Slater-Koster structure to be filled
    type(TSlater), intent(inout) :: slako

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens
  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    character(len=:), allocatable :: buffer
    type(hsd_table), pointer :: child

    call getNodeName(node, buffer)
    select case (buffer)
    case ("dftb")
  #:if WITH_TRANSPORT
      call readDFTBHam(node, ctrl, geo, slako, tp, greendens, poisson, errStatus)
  #:else
      call readDFTBHam(node, ctrl, geo, slako, poisson, errStatus)
  #:endif
      @:PROPAGATE_ERROR(errStatus)
    case ("xtb")
  #:if WITH_TRANSPORT
      call readXTBHam(node, ctrl, geo, tp, greendens, poisson, errStatus)
  #:else
      call readXTBHam(node, ctrl, geo, poisson, errStatus)
  #:endif
      @:PROPAGATE_ERROR(errStatus)
    case default
      call dftbp_error(node, "Invalid Hamiltonian")
    end select

  #:if WITH_API
    call getChild(node, "ASI", child, requested=.false.)
    if (associated(child)) then
      ctrl%isASICallbackEnabled = .true.
    else
      ctrl%isASICallbackEnabled = .false.
    end if
  #:endif

  end subroutine readHamiltonian


  !> Reads DFTB-Hamiltonian
#:if WITH_TRANSPORT
  subroutine readDFTBHam(node, ctrl, geo, slako, tp, greendens, poisson, errStatus)
#:else
  subroutine readDFTBHam(node, ctrl, geo, slako, poisson, errStatus)
#:endif

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Slater-Koster structure to be filled
    type(TSlater), intent(inout) :: slako

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(hsd_table), pointer :: value1, child, child2, child3
    type(hsd_child_list), pointer :: children
    character(len=:), allocatable :: buffer, buffer2, modifier
    type(TListCharLc), allocatable :: skFiles(:,:)
    integer, allocatable :: shellsTmp(:)
    character(len=:), allocatable :: strArr(:)
    type(TAngShellBlocks), allocatable :: angShells(:)
    logical, allocatable :: repPoly(:,:)
    integer :: iSp1, iSp2, ii
    character(lc) :: prefix, suffix, separator, elem1, elem2, strTmp
    character(lc) :: errorStr
    logical :: tLower
    real(dp) :: rTmp
    integer, allocatable :: iTmpN(:)
    integer :: nShell, skInterMeth
    real(dp) :: rSKCutOff
    character(len=:), allocatable :: searchPath(:)
    character(len=:), allocatable :: strOut, strJoin
    logical :: isHalogenXCorr

    !> For hybrid functional calculations
    type(THybridXcSKTag) :: hybridXcSK

    ctrl%hamiltonian = hamiltonianTypes%dftb

    call readMaxAngularMomentum(node, geo, angShells)

    ! Orbitals and angular momenta for the given shells (once the SK files contain the full
    ! information about the basis, this will be moved to the SK reading routine).
    allocate(slako%orb)
    call setupOrbitals(slako%orb, geo, angShells)

    ! Slater-Koster files
    call getParamSearchPaths(searchPath)
    strJoin = joinPathsPrettyErr(searchPath)
    allocate(skFiles(geo%nSpecies, geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
      do iSp2 = 1, geo%nSpecies
        call init(skFiles(iSp2, iSp1))
      end do
    end do
    call getChildValue(node, "SlaterKosterFiles", value1, child=child)
    call getNodeName(value1, buffer)
    select case(buffer)
    case ("type2filenames")
      call hsd_get_or_set(value1, "Prefix", buffer2, "")
      prefix = unquote(buffer2)
      call hsd_get_or_set(value1, "Suffix", buffer2, "")
      suffix = unquote(buffer2)
      call hsd_get_or_set(value1, "Separator", buffer2, "")
      separator = unquote(buffer2)
      call hsd_get_or_set(value1, "LowerCaseTypeName", tLower, .false.)
      do iSp1 = 1, geo%nSpecies
        if (tLower) then
          elem1 = tolower(geo%speciesNames(iSp1))
        else
          elem1 = geo%speciesNames(iSp1)
        end if
        do iSp2 = 1, geo%nSpecies
          if (tLower) then
            elem2 = tolower(geo%speciesNames(iSp2))
          else
            elem2 = geo%speciesNames(iSp2)
          end if
          strTmp = trim(prefix) // trim(elem1) // trim(separator) // trim(elem2) // trim(suffix)
          call findFile(searchPath, strTmp, strOut)
          if (.not. allocated(strOut)) then
            call dftbp_error(value1, "SK file with generated name '" // trim(strTmp)&
                & // "' not found." // newline // "   (search path(s): " // strJoin // ").")
          end if
          strTmp = strOut
          call append(skFiles(iSp2, iSp1), strTmp)
        end do
      end do
    case default
      if (associated(value1)) value1%processed = .false.
      call hsd_get_or_set(child, "Prefix", buffer2, "")
      prefix = unquote(buffer2)
      call getChild(child, "Suffix", child2, requested=.false.)
      if (associated(child2)) then
        call dftbp_error(child2, "Keyword requires SlaterKosterFiles = Type2Filenames {")
      end if
      call getChild(child, "Separator", child2, requested=.false.)
      if (associated(child2)) then
        call dftbp_error(child2, "Keyword requires SlaterKosterFiles = Type2Filenames {")
      end if
      call getChild(child, "LowerCaseTypeName", child2, requested=.false.)
      if (associated(child2)) then
        call dftbp_error(child2, "Keyword requires SlaterKosterFiles = Type2Filenames {")
      end if
      do iSp1 = 1, geo%nSpecies
        do iSp2 = 1, geo%nSpecies
          strTmp = trim(geo%speciesNames(iSp1)) // "-" // trim(geo%speciesNames(iSp2))
          call hsd_get(child, trim(strTmp), strArr)
          call getChild(child, trim(strTmp), child2)
          if (size(strArr) /= angShells(iSp1)%nBlocks * angShells(iSp2)%nBlocks) then
            write(errorStr, "(A,I0,A,I0)") "Incorrect number of Slater-Koster files for " //&
                & trim(strTmp) // ", expected ", angShells(iSp1)%nBlocks * angShells(iSp2)%nBlocks,&
                & " but received ", size(strArr)
            call dftbp_error(child2, errorStr)
          end if
          do ii = 1, size(strArr)
            strTmp = trim(prefix) // trim(strArr(ii))
            call findFile(searchPath, strTmp, strOut)
            if (.not. allocated(strOut)) then
              call dftbp_error(child2, "SK file '" // trim(strTmp) // "' not found." // newline&
                  & // "   (search path(s): " // strJoin // ").")
            end if
            strTmp = strOut
            call append(skFiles(iSp2, iSp1), strTmp)
          end do
        end do
      end do
    end select

    ! Which repulsive is defined by polynomial? (Default: None)
    allocate(repPoly(geo%nSpecies, geo%nSpecies))
    call getChildValue(node, "PolynomialRepulsive", value1, "", child=child, list=.true.,&
        & allowEmptyValue=.true., dummyValue=.true.)
    call getNodeName2(value1, buffer)
    if (buffer == "" .and. hasInlineData(child)) buffer = textNodeName
    select case (buffer)
    case ("")
      repPoly(:,:) = .false.
    case("setforall")
      call getChildValue(value1, "", repPoly(1,1))
      repPoly(:,:) = repPoly(1,1)
    case default
      do iSp1 = 1, geo%nSpecies
        do iSp2 = 1, geo%nSpecies
          strTmp = trim(geo%speciesNames(iSp1)) // "-" &
              &// trim(geo%speciesNames(iSp2))
          call hsd_get_or_set(child, trim(strTmp), repPoly(iSp2, iSp1), .false.)
        end do
      end do
      if (.not. all(repPoly .eqv. transpose(repPoly))) then
        call dftbp_error(value1, "Asymmetric definition (both A-B and B-A must&
            & be defined for using polynomial repulsive)")
      end if
    end select

    call parseChimes(node, ctrl%chimesRepInput)

    ! SCC
    call hsd_get_or_set(node, "SCC", ctrl%tSCC, .false.)

    call parseHybridBlock(node, ctrl%hybridXcInp, ctrl, geo, skFiles)

    if (allocated(ctrl%hybridXcInp)) then
      if (.not.ctrl%tSCC) then
        call dftbp_error(node, "Hybrid calculations require SCC = Yes")
      end if
    end if

    if (ctrl%tSCC) then
      call hsd_get_or_set(node, "ShellResolvedSCC", ctrl%tShellResolved, .false.)
    else
      ctrl%tShellResolved = .false.
    end if

    call hsd_get_or_set(node, "OldSKInterpolation", ctrl%oldSKInter, .false.)
    if (ctrl%oldSKInter) then
      skInterMeth = skEqGridOld
    else
      skInterMeth = skEqGridNew
    end if

    if (.not. allocated(ctrl%hybridXcInp)) then
      call getChild(node, "TruncateSKRange", child, requested=.false.)
      if (associated(child)) then
        call warning("Artificially truncating the SK table, this is normally a bad idea!")
        call SKTruncations(child, rSKCutOff, skInterMeth)
        call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tShellResolved,&
            & skInterMeth, repPoly, rSKCutOff)
      else
        rSKCutOff = 0.0_dp
        call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tShellResolved,&
            & skInterMeth, repPoly)
      end if
    else
      call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tShellResolved,&
          & skInterMeth, repPoly, hybridXcSK=hybridXcSK)
      ctrl%hybridXcInp%omega = hybridXcSK%omega
      ctrl%hybridXcInp%camAlpha = hybridXcSK%camAlpha
      ctrl%hybridXcInp%camBeta = hybridXcSK%camBeta
    end if

    do iSp1 = 1, geo%nSpecies
      do iSp2 = 1, geo%nSpecies
        call destruct(skFiles(iSp2, iSp1))
      end do
    end do
    deallocate(angShells)
    deallocate(skFiles)
    deallocate(repPoly)

    ! SCC parameters
    ifSCC: if (ctrl%tSCC) then

      ! get charge mixing options
      call readSccOptions(node, ctrl, geo)

      ! DFTB hydrogen bond corrections
      call readHCorrection(node, geo, ctrl)

      !> TI-DFTB varibles for Delta DFTB
      call getChild(node, "NonAufbau", child, requested=.false.)
      if (associated(child)) then
        ctrl%isNonAufbau = .true.
        call hsd_get_or_set(child, "SpinPurify", ctrl%isSpinPurify, .true.)
        call hsd_get_or_set(child, "GroundGuess", ctrl%isGroundGuess, .false.)
        ctrl%tSpin = .true.
        ctrl%t2Component = .false.
        ctrl%nrSpinPol = 0.0_dp
        ctrl%tSpinSharedEf = .false.
      else
        ctrl%isNonAufbau = .false.
      end if

    end if ifSCC

    ! Customize the reference atomic charges for virtual doping
    call readCustomReferenceOcc(node, slako%orb, slako%skOcc, geo, &
        & ctrl%customOccAtoms, ctrl%customOccFillings)

    ! Spin calculation
    if (ctrl%reksInp%reksAlg == reksTypes%noReks  .and. .not.ctrl%isNonAufbau) then
    #:if WITH_TRANSPORT
      call readSpinPolarisation(node, ctrl, geo, tp)
    #:else
      call readSpinPolarisation(node, ctrl, geo)
    #:endif
    end if

    ! temporararily removed until debugged
    !if (.not. ctrl%tscc) then
    !  !! In a non-SCC calculation it is possible to upload charge shifts
    !  !! This is useful if the calculation can jump directly to the Analysis block
    !  call getChildValue(node, "ReadShifts", ctrl%tReadShifts, .false.)
    !end if
    ctrl%tReadShifts = .false.

    ! External fields and potentials
    call readExternal(node, ctrl, geo)

    ! Non-self-consistent spin-orbit coupling
    call readSpinOrbit(node, ctrl, geo, slako%orb)

    ! Electronic solver
  #:if WITH_TRANSPORT
    call readSolver(node, ctrl, geo, tp, greendens, poisson)

    if (tp%taskUpload) then
      ! Initialise variable, but unused
      ctrl%nrChrg =  0.0_dp
    else
      ! Charge
      call hsd_get_or_set(node, "Charge", ctrl%nrChrg, 0.0_dp)
    end if
  #:else
    call readSolver(node, ctrl, geo, poisson)

    ! Charge
    call hsd_get_or_set(node, "Charge", ctrl%nrChrg, 0.0_dp)
  #:endif

    ! K-Points
    call readKPoints(node, ctrl, geo, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    if (ctrl%tscc) then

      call getChild(node, "OrbitalPotential", child, requested=.false.)
      if (associated(child)) then
        allocate(ctrl%dftbUInp)
        call hsd_get_or_set(child, "Functional", buffer, "fll")
        select case(tolower(buffer))
        case ("fll")
          ctrl%dftbUInp%iFunctional = plusUFunctionals%fll
        case ("psic")
          ctrl%dftbUInp%iFunctional = plusUFunctionals%pSic
        case default
          call dftbp_error(child,"Unknown orbital functional :"// buffer)
        end select

        allocate(ctrl%dftbUInp%nUJ(geo%nSpecies))
        ctrl%dftbUInp%nUJ(:) = 0

        ! First pass: count blocks per species and read data into temporary arrays
        do iSp1 = 1, geo%nSpecies
          call getChildren(child, trim(geo%speciesNames(iSp1)), children)
          ctrl%dftbUInp%nUJ(iSp1) = getLength(children)
          call destroyNodeList(children)
        end do

        allocate(ctrl%dftbUInp%UJ(maxval(ctrl%dftbUInp%nUJ), geo%nSpecies))
        ctrl%dftbUInp%UJ(:,:) = 0.0_dp
        allocate(ctrl%dftbUInp%niUJ(maxval(ctrl%dftbUInp%nUJ), geo%nSpecies))
        ctrl%dftbUInp%niUJ(:,:) = 0

        ! Second pass: read shells and U-J values
        do iSp1 = 1, geo%nSpecies
          call getChildren(child, trim(geo%speciesNames(iSp1)), children)
          do ii = 1, ctrl%dftbUInp%nUJ(iSp1)
            call getItem1(children, ii, child2)

            call hsd_get(child2, "Shells", shellsTmp)
            ctrl%dftbUInp%niUJ(ii, iSp1) = size(shellsTmp)
            deallocate(shellsTmp)

            call getChildValue(child2, "uj", rTmp, 0.0_dp, modifier=modifier, &
                & child=child3)
            call convertUnitHsd(modifier, energyUnits, child3, rTmp)
            if (rTmp < 0.0_dp) then
              write(errorStr,"(F12.8)")rTmp
              call dftbp_error(child2,"Negative value of U-J:"//errorStr)
            end if
            if (rTmp <= 1.0E-10_dp) then
              write(errorStr,"(F12.8)")rTmp
              call dftbp_error(child2,"Invalid value of U-J, too small: " &
                  & //errorStr)
            end if
            ctrl%dftbUInp%UJ(ii, iSp1) = rTmp
          end do
          call destroyNodeList(children)
        end do

        allocate(ctrl%dftbUInp%iUJ(maxval(ctrl%dftbUInp%niUJ),&
            & maxval(ctrl%dftbUInp%nUJ), geo%nSpecies))
        ctrl%dftbUInp%iUJ(:,:,:) = 0

        ! Third pass: read shell indices into final array
        do iSp1 = 1, geo%nSpecies
          call getChildren(child, trim(geo%speciesNames(iSp1)), children)
          do ii = 1, ctrl%dftbUInp%nUJ(iSp1)
            call getItem1(children, ii, child2)
            call hsd_get(child2, "Shells", shellsTmp)
            ctrl%dftbUInp%iUJ(1:size(shellsTmp), ii, iSp1) = shellsTmp(:)
            deallocate(shellsTmp)
          end do
          call destroyNodeList(children)
        end do

        ! check input values
        allocate(iTmpN(slako%orb%mShell))
        do iSp1 = 1, geo%nSpecies
          iTmpN = 0
          ! loop over number of blocks for that species
          do ii = 1, ctrl%dftbUInp%nUJ(iSp1)
            iTmpN(ctrl%dftbUInp%iUJ(1:ctrl%dftbUInp%niUJ(ii,iSp1),ii,iSp1)) = &
                & iTmpN(ctrl%dftbUInp%iUJ(1:ctrl%dftbUInp%niUJ(ii,iSp1),ii,iSp1)) + 1
          end do
          if (any(iTmpN(:)>1)) then
            write(stdout, *)'Multiple copies of shells present in OrbitalPotential!'
            write(stdout, "(A,A3,A,I2)") &
                & 'The count for the occurrence of shells of species ', &
                & trim(geo%speciesNames(iSp1)),' are:'
            write(stdout, *)iTmpN(1:slako%orb%nShell(iSp1))
            call abortProgram()
          end if
        end do
        deallocate(iTmpN)

      end if

      ! On-site
      call getChildValue(node, "OnSiteCorrection", value1, "", child=child, allowEmptyValue=.true.,&
          & dummyValue=.true.)
      if (associated(value1)) then
        allocate(ctrl%onSiteElements(slako%orb%mShell, slako%orb%mShell, 2, geo%nSpecies))
        do iSp1 = 1, geo%nSpecies
          call getChildValue(child, trim(geo%speciesNames(iSp1))//"uu",&
              & ctrl%onSiteElements(:slako%orb%nShell(iSp1), :slako%orb%nShell(iSp1), 1, iSp1))
          call getChildValue(child, trim(geo%speciesNames(iSp1))//"ud",&
              & ctrl%onSiteElements(:slako%orb%nShell(iSp1), :slako%orb%nShell(iSp1), 2, iSp1))
        end do
      end if

    end if

    ! Dispersion
    call getChildValue(node, "Dispersion", value1, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%dispInp)
      call readDispersion(child, geo, ctrl%dispInp, ctrl%nrChrg, ctrl%tSCC)
    end if

    ! Solvation
    call getChildValue(node, "Solvation", value1, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%solvInp)
      call readSolvation(child, geo, ctrl%solvInp)
      call hsd_get_or_set(value1, "RescaleSolvatedFields", ctrl%isSolvatedFieldRescaled, .true.)
    end if

    ! Electronic constraints
    call getChildValue(node, "ElectronicConstraints", value1, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true., list=.true.)
    if (associated(value1)) then
      allocate(ctrl%elecConstraintInp)
      call readElecConstraintInput(child, geo, ctrl%tSpin, ctrl%t2Component, ctrl%elecConstraintInp)
      if (.not. allocated(ctrl%elecConstraintInp%mullikenConstrs)) then
        call dftbp_warning(child, "No electronic constraint specified")
        deallocate(ctrl%elecConstraintInp)
      end if
    end if

    if (ctrl%tLatOpt .and. .not. geo%tPeriodic) then
      call error("Lattice optimisation only applies for periodic structures.")
    end if

    if (ctrl%tSCC) then
    #:if WITH_TRANSPORT
      call readElectrostatics(node, ctrl, geo, tp, poisson)
    #:else
      call readElectrostatics(node, ctrl, geo, poisson)
    #:endif
    end if

    ! Multipole expansion
    ctrl%isMdftb = .false.
    call readMdftb(node, ctrl, geo)

    ! Third order stuff
    ctrl%t3rd = .false.
    ctrl%t3rdFull = .false.
    if (ctrl%tSCC) then
      call hsd_get_or_set(node, "ThirdOrder", ctrl%t3rd, .false.)
      call hsd_get_or_set(node, "ThirdOrderFull", ctrl%t3rdFull, .false.)
      if (ctrl%t3rd .and. ctrl%t3rdFull) then
        call dftbp_error(node, "You must choose either ThirdOrder or&
            & ThirdOrderFull")
      end if
      if (ctrl%t3rd .and. ctrl%tShellResolved) then
        call error("Only full third-order DFTB is compatible with orbital&
            & resolved SCC")
      end if
      if (ctrl%t3rd .or. ctrl%t3rdFull) then
        call getChild(node, 'HubbardDerivs', child, requested=.true.)
        allocate(ctrl%HubDerivs(slako%orb%mShell, geo%nSpecies))
        ctrl%hubDerivs(:,:) = 0.0_dp
        do iSp1 = 1, geo%nSpecies
          nShell = slako%orb%nShell(iSp1)
          if (ctrl%tShellResolved) then
            call getChildValue(child, geo%speciesNames(iSp1),&
                & ctrl%hubDerivs(1:nShell, iSp1))
          else
            call getChildValue(child, geo%speciesNames(iSp1),&
                & ctrl%hubDerivs(1, iSp1))
            ctrl%hubDerivs(2:nShell, iSp1) = ctrl%hubDerivs(1, iSp1)
          end if
        end do
        if (ctrl%t3rd) then
          allocate(ctrl%thirdOrderOn(geo%nAtom, 2))
          ctrl%thirdOrderOn(:,1) = 0.0_dp
          ctrl%thirdOrderOn(:,2) = ctrl%hubDerivs(1, geo%species)
        end if

        ! Halogen correction to the DFTB3 model
        isHalogenXCorr =&
            & any([(any(halogenXSpecies1(ii) == geo%speciesNames), ii=1, size(halogenXSpecies1))])&
            & .and. &
            & any([(any(halogenXSpecies2(ii) == geo%speciesNames), ii=1, size(halogenXSpecies2))])

        if (isHalogenXCorr) then
          call hsd_get_or_set(node, "HalogenXCorr", ctrl%tHalogenX, .false.)
        end if

      end if
    end if

    call readDifferentiation(node, ctrl)

    if (ctrl%tSCC) then
      ! Force type
      call readForceOptions(node, ctrl)
    else
      ctrl%forceType = forceTypes%orig
    end if

    call readCustomisedHubbards(node, geo, slako%orb, ctrl%tShellResolved, ctrl%hubbU)

  end subroutine readDFTBHam


  !> Reads xTB-Hamiltonian
#:if WITH_TRANSPORT
  subroutine readXTBHam(node, ctrl, geo, tp, greendens, poisson, errStatus)
#:else
  subroutine readXTBHam(node, ctrl, geo, poisson, errStatus)
#:endif

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer
    character(len=:), allocatable :: searchPath(:)
    integer :: method
    character(len=:), allocatable :: paramFile, paramTmp
    type(TOrbitals) :: orb

    ctrl%hamiltonian = hamiltonianTypes%xtb

    allocate(ctrl%tbliteInp)
    call ctrl%tbliteInp%setupGeometry(geo%nAtom, geo%species, geo%coords, geo%speciesNames,&
        & geo%latVecs)

    call getChild(node, "Method", child, requested=.false.)
    if (associated(child)) then
      call getChildValue(child, "", buffer)
      select case(unquote(buffer))
      case default
        call dftbp_error(child, "Unknown method "//buffer//" for xTB Hamiltonian")
      case("GFN1-xTB")
        method = tbliteMethod%gfn1xtb
      case("GFN2-xTB")
        method = tbliteMethod%gfn2xtb
      case("IPEA1-xTB")
        method = tbliteMethod%ipea1xtb
      end select
      call ctrl%tbliteInp%setupCalculator(method)
      ctrl%tbliteInp%info%name = trim(unquote(buffer))
    else
      call getChildValue(node, "ParameterFile", value1, "", child=child, allowEmptyValue=.true.,&
          & dummyValue=.true.)
      if (associated(value1)) then
        call getChildValue(child, "", buffer)
        paramFile = trim(unquote(buffer))
        call getParamSearchPaths(searchPath)
        call findFile(searchPath, paramFile, paramTmp)
        if (allocated(paramTmp)) call move_alloc(paramTmp, paramFile)
        write(stdOut, '(a)') "Using parameter file '"//paramFile//"' for xTB Hamiltonian"
        call ctrl%tbliteInp%setupCalculator(paramFile)
      else
        call dftbp_error(node, "Either a Method or ParameterFile must be specified for xTB")
      end if
    end if

    call hsd_get_or_set(node, "ShellResolvedSCC", ctrl%tShellResolved, .true.)

    ! SCC parameters
    call hsd_get_or_set(node, "SCC", ctrl%tSCC, .true.)
    ifSCC: if (ctrl%tSCC) then

      ! get charge mixing options etc.
      call readSccOptions(node, ctrl, geo)

      !> TI-DFTB varibles for Delta DFTB
      call getChild(node, "NonAufbau", child, requested=.false.)
      if (associated(child)) then
        ctrl%isNonAufbau = .true.
        call hsd_get_or_set(child, "SpinPurify", ctrl%isSpinPurify, .true.)
        call hsd_get_or_set(child, "GroundGuess", ctrl%isGroundGuess, .false.)
        ctrl%tSpin = .true.
        ctrl%t2Component = .false.
        ctrl%nrSpinPol = 0.0_dp
        ctrl%tSpinSharedEf = .false.
      else
        ctrl%isNonAufbau = .false.
      end if

    end if ifSCC

    ! Spin calculation
    if (ctrl%reksInp%reksAlg == reksTypes%noReks .and. .not.ctrl%isNonAufbau .and. ctrl%tSCC) then
    #:if WITH_TRANSPORT
      call readSpinPolarisation(node, ctrl, geo, tp)
    #:else
      call readSpinPolarisation(node, ctrl, geo)
    #:endif
    end if

    ! temporararily removed until debugged
    !if (.not. ctrl%tscc) then
    !  !! In a non-SCC calculation it is possible to upload charge shifts
    !  !! This is useful if the calculation can jump directly to the Analysis block
    !  call getChildValue(node, "ReadShifts", ctrl%tReadShifts, .false.)
    !end if
    ctrl%tReadShifts = .false.

    ! External fields and potentials
    call readExternal(node, ctrl, geo)

    ! Non-self-consistent spin-orbit coupling
    call ctrl%tbliteInp%setupOrbitals(geo%species, orb)
    call readSpinOrbit(node, ctrl, geo, orb)

    ! Electronic solver
  #:if WITH_TRANSPORT
    call readSolver(node, ctrl, geo, tp, greendens, poisson)

    if (tp%taskUpload) then
      ! Initialise, but unused
      ctrl%nrChrg =  0.0_dp
    else
      ! Charge
      call hsd_get_or_set(node, "Charge", ctrl%nrChrg, 0.0_dp)
    end if
  #:else
    call readSolver(node, ctrl, geo, poisson)

    ! Charge
    call hsd_get_or_set(node, "Charge", ctrl%nrChrg, 0.0_dp)
  #:endif

    ! K-Points
    call readKPoints(node, ctrl, geo, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Dispersion
    call getChildValue(node, "Dispersion", value1, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%dispInp)
      call readDispersion(child, geo, ctrl%dispInp, ctrl%nrChrg, ctrl%tSCC)
    end if

    ! Solvation
    call getChildValue(node, "Solvation", value1, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%solvInp)
      call readSolvation(child, geo, ctrl%solvInp)
      call hsd_get_or_set(value1, "RescaleSolvatedFields", ctrl%isSolvatedFieldRescaled, .true.)
    end if

    if (ctrl%tLatOpt .and. .not. geo%tPeriodic) then
      call error("Lattice optimisation only applies for periodic structures.")
    end if

  #:if WITH_TRANSPORT
    call readElectrostatics(node, ctrl, geo, tp, poisson)
  #:else
    call readElectrostatics(node, ctrl, geo, poisson)
  #:endif

    ! Third order stuff
    ctrl%t3rd = .true.
    ctrl%t3rdFull = .false.

    call readDifferentiation(node, ctrl)

    if (ctrl%tSCC) then
      ! Force type
      call readForceOptions(node, ctrl)
    else
      ctrl%forceType = forceTypes%orig
    end if

    ! Electronic constraints
    call getChildValue(node, "ElectronicConstraints", value1, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true., list=.true.)
    if (associated(value1)) then
      allocate(ctrl%elecConstraintInp)
      call readElecConstraintInput(child, geo, ctrl%tSpin, ctrl%t2Component, ctrl%elecConstraintInp)
    end if

  end subroutine readXTBHam


  !> Read in hamiltonian settings that are influenced by those read from REKS{}, electronDynamics{}
  subroutine readLaterHamiltonian(hamNode, ctrl, driverNode, geo)

    !> Hamiltonian node to parse
    type(hsd_table), pointer :: hamNode

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> Geometry driver node to parse
    type(hsd_table), pointer :: driverNode

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    type(hsd_table), pointer :: value1, value2, child, child2
    character(len=:), allocatable :: buffer, buffer2
    real(dp), allocatable :: dynMixMatrix(:,:)
    integer :: dynMixNRows, dynMixNCols

    if (ctrl%reksInp%reksAlg == reksTypes%noReks) then

      if (ctrl%tSCC) then

        call getChildValue(hamNode, "Mixer", value1, "Broyden", child=child)
        call getNodeName(value1, buffer)
        select case(buffer)

        case ("broyden")

          allocate(ctrl%mixerInp%broydenMixerInp)
          associate (inp => ctrl%mixerInp%broydenMixerInp)
            call hsd_get_or_set(value1, "MixingParameter", inp%mixParam, 0.2_dp)
            call hsd_get_or_set(value1, "InverseJacobiWeight", inp%omega0, 0.01_dp)
            call hsd_get_or_set(value1, "MinimalWeight", inp%minWeight, 1.0_dp)
            call hsd_get_or_set(value1, "MaximalWeight", inp%maxWeight, 1.0e5_dp)
            call hsd_get_or_set(value1, "WeightFactor", inp%weightFac, 1.0e-2_dp)
          end associate

        case ("anderson")

          allocate(ctrl%mixerInp%andersonMixerInp)
          associate (inp => ctrl%mixerInp%andersonMixerInp)
            call hsd_get_or_set(value1, "MixingParameter", inp%mixParam, 0.05_dp)
            call hsd_get_or_set(value1, "Generations", inp%iGenerations, 4)
            call hsd_get_or_set(value1, "InitMixingParameter", inp%initMixParam, 0.01_dp)
            call getChildValue(value1, "DynMixingParameters", value2, "", child=child,&
                & allowEmptyValue=.true.)
            call getNodeName2(value2, buffer2)
            if (buffer2 == "" .and. .not. hasInlineData(child)) then
              inp%nConvMixParam = 0
            else
              call hsd_get_matrix(child, "", dynMixMatrix, dynMixNRows, dynMixNCols)
              if (dynMixNCols < 1) then
                call dftbp_error(child, "At least one dynamic mixing parameter must be defined.")
              end if
              inp%nConvMixParam = dynMixNCols
              allocate(inp%convMixParam(2, inp%nConvMixParam))
              inp%convMixParam(:,:) = dynMixMatrix(:,:)
              deallocate(dynMixMatrix)
            end if
            call hsd_get_or_set(value1, "DiagonalRescaling", inp%omega0, 1.0e-2_dp)
          end associate

        case ("simple")

          allocate(ctrl%mixerInp%simpleMixerInp)
          associate (inp => ctrl%mixerInp%simpleMixerInp)
            call hsd_get_or_set(value1, "MixingParameter", inp%mixParam, 0.05_dp)
          end associate

        case ("diis")

          allocate(ctrl%mixerInp%diisMixerInp)
          associate (inp => ctrl%mixerInp%diisMixerInp)
            call hsd_get_or_set(value1, "InitMixingParameter", inp%initMixParam, 0.2_dp)
            call hsd_get_or_set(value1, "Generations", inp%iGenerations, 6)
            call hsd_get_or_set(value1, "UseFromStart", inp%tFromStart, .true.)
          end associate

        case default

          call getNodeHSDName(value1, buffer)
          call dftbp_error(child, "Invalid mixer '" // buffer // "'")

        end select

      end if

      if (ctrl%tMD) then
        if (ctrl%thermostatInp%thermostatType /= thermostatTypes%dummy) then
          call getChildValue(driverNode, "Thermostat", child, child=child2)
          if (ctrl%reksInp%reksAlg == reksTypes%noReks) then
            call hsd_get_or_set(child, "AdaptFillingTemp", ctrl%tSetFillingTemp, .false.)
          end if
        end if
      end if

    end if

    hamNeedsT: if (ctrl%reksInp%reksAlg == reksTypes%noReks) then

      if (allocated(ctrl%elecDynInp)) then
        if (ctrl%elecDynInp%tReadRestart .and. .not.ctrl%elecDynInp%tPopulations) then
          exit hamNeedsT
        end if
      end if

      if (ctrl%solver%isolver /= electronicSolverTypes%GF) then
        call readElectronicFilling(hamNode, ctrl, geo)
      end if

    end if hamNeedsT

  end subroutine readLaterHamiltonian


end module dftbp_dftbplus_parser_hamiltonian
