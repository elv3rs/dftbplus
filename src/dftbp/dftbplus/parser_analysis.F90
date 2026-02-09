#:include 'common.fypp'
#:include 'error.fypp'

!> Parser routines for Analysis and related input blocks.
module dftbp_dftbplus_parser_analysis
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_unitconversion, only : energyUnits, freqUnits, lengthUnits
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_dftbplus_parser_kpoints, only : maxSelfConsIterations
  use dftbp_elecsolvers_elecsolvers, only : electronicSolverTypes, providesEigenvalues
  use dftbp_io_charmanip, only : unquote
  use hsd, only : hsd_rename_child, hsd_get_or_set, hsd_get, hsd_get_matrix, hsd_get_attrib,&
      & hsd_table_ptr, hsd_get_child_tables, hsd_get_table, hsd_set, HSD_STAT_OK
  use dftbp_io_hsdutils, only : dftbp_error, getSelectedAtomIndices, getNodeName, getNodeName2,&
      & hasInlineData
  use dftbp_io_message, only : error
  use dftbp_io_unitconv, only : convertUnitHsd
  use hsd_data, only : hsd_table
  use dftbp_math_simplealgebra, only : determinant33
  use dftbp_solvation_solvparser, only : readCM5
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_linkedlist, only : append, init
  use dftbp_type_typegeometry, only : TGeometry
#:if WITH_TRANSPORT
  use dftbp_dftbplus_parser_transport, only : readTunAndDos
  use dftbp_transport_negfvars, only : TNEGFTunDos, TTransPar
#:endif

  implicit none

  private
  public :: readAnalysis, readLaterAnalysis

contains



  !> Reads the analysis block
#:if WITH_TRANSPORT
  subroutine readAnalysis(node, ctrl, geo, orb, transpar, tundos)
#:else
  subroutine readAnalysis(node, ctrl, geo)
#:endif

    !> Node to parse
    type(hsd_table), pointer :: node

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Orbital
    type(TOrbitals), intent(in), allocatable :: orb

    !> Transport parameters
    type(TTransPar), intent(inout) :: transpar

    !> Tunneling and Dos parameters
    type(TNEGFTunDos), intent(inout) :: tundos
  #:endif

    type(hsd_table), pointer :: val, child, child2, child3
    type(hsd_table_ptr), allocatable :: children(:)
    integer, allocatable :: pTmpI1(:)
    character(len=:), allocatable :: buffer, modifier
    integer :: nReg, iReg
    character(lc) :: strTmp
    integer :: stat
    logical :: tPipekDense
    logical :: tWriteBandDatDefault, tHaveEigenDecomposition, tHaveDensityMatrix
    logical :: isEtaNeeded

    tHaveEigenDecomposition = providesEigenvalues(ctrl%solver%isolver)
    tHaveDensityMatrix = ctrl%solver%isolver /= electronicSolverTypes%OnlyTransport

    if (tHaveEigenDecomposition) then

      call hsd_get_table(node, "ProjectStates", child, stat, auto_wrap=.true.)
      if (associated(child)) then
        call hsd_get_child_tables(child, "Region", children)
        nReg = size(children)
      else
        nReg = 0
      end if
      ctrl%tProjEigenvecs = (nReg > 0)
      if (ctrl%tProjEigenvecs) then
        allocate(ctrl%tShellResInRegion(nReg))
        allocate(ctrl%tOrbResInRegion(nReg))
        allocate(ctrl%RegionLabel(nReg))
        call init(ctrl%iAtInRegion)
        do iReg = 1, nReg
          child2 => children(iReg)%ptr
          call hsd_get(child2, "Atoms", buffer, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value: 'Atoms'")
          call hsd_get_table(child2, "Atoms", child3, stat)
          if (.not. associated(child3)) child3 => child2
          call getSelectedAtomIndices(child3, buffer, geo%speciesNames, geo%species, pTmpI1)
          call append(ctrl%iAtInRegion, pTmpI1)
          call hsd_get_or_set(child2, "ShellResolved", ctrl%tShellResInRegion(iReg), .false.,&
              & child=child3)
          if (ctrl%tShellResInRegion(iReg)) then
            if (.not. all(geo%species(pTmpI1) == geo%species(pTmpI1(1)))) then
              call dftbp_error(child3, "Shell resolved PDOS only allowed for &
                  &regions where all atoms belong to the same species")
            end if
          end if
          call hsd_get_or_set(child2, "OrbitalResolved", &
              & ctrl%tOrbResInRegion(iReg), .false., child=child3)
          if (ctrl%tOrbResInRegion(iReg)) then
            if (.not. all(geo%species(pTmpI1) == geo%species(pTmpI1(1)))) then
              call dftbp_error(child3, "Orbital resolved PDOS only allowed for &
                  &regions where all atoms belong to the same species")
            end if
          end if
          deallocate(pTmpI1)
          write(strTmp, "('region',I0)") iReg
          call hsd_get_or_set(child2, "Label", buffer, trim(strTmp))
          ctrl%RegionLabel(iReg) = unquote(buffer)
        end do
      end if

      call hsd_rename_child(node, "Localize", "Localise")
      call hsd_get_table(node, "Localise", val, stat, auto_wrap=.true.)
      if (associated(val)) then
        ctrl%tLocalise = .true.
        call hsd_get_table(val, "PipekMezey", child2, stat, auto_wrap=.true.)
        if (associated(child2)) then
          allocate(ctrl%pipekMezeyInp)
          associate(inp => ctrl%pipekMezeyInp)
            call hsd_get_or_set(child2, "MaxIterations", inp%maxIter, 100)
            tPipekDense = .true.
            if (.not. geo%tPeriodic) then
              call hsd_get_or_set(child2, "Dense", tPipekDense, .false.)
              if (.not. tPipekDense) then
                call hsd_get(child2, "SparseTolerances", inp%sparseTols, stat=stat)
                if (stat == 0) then
                  if (size(inp%sparseTols) < 1) then
                    call dftbp_error(child2, "Missing values of tolerances.")
                  end if
                else
                  inp%sparseTols = [0.1_dp, 0.01_dp, 1.0E-6_dp, 1.0E-12_dp]
                  call hsd_set(child2, "SparseTolerances", inp%sparseTols)
                end if
              end if
            end if
            if (tPipekDense) then
              call hsd_get_or_set(child2, "Tolerance", inp%tolerance, 1.0E-4_dp)
            end if
          end associate
        else
          call dftbp_error(val, "No localisation method chosen")
        end if
      end if

      call hsd_get_or_set(node, "WriteEigenvectors", ctrl%tPrintEigVecs, .false.)

    #:if WITH_SOCKETS
      tWriteBandDatDefault = .not. allocated(ctrl%socketInput)
    #:else
      tWriteBandDatDefault = .true.
    #:endif

      call hsd_get_or_set(node, "WriteBandOut", ctrl%tWriteBandDat, tWriteBandDatDefault)

      ! electric field polarisability of system
      call hsd_rename_child(node, "Polarizability", "Polarisability")
      call hsd_get_table(node, "Polarisability", child, stat, auto_wrap=.true.)
      if (associated(child)) then
        if (.not.allocated(ctrl%perturbInp)) allocate(ctrl%perturbInp)
        ctrl%perturbInp%isEPerturb = .true.
        call freqRanges(child, ctrl%perturbInp%dynEFreq)
      end if

      ! Perturbation with respect to on-site potentials (related to Fukui charges)
      call hsd_get_table(node, "ResponseKernel", child, stat, auto_wrap=.true.)
      if (associated(child)) then
        if (.not.allocated(ctrl%perturbInp)) allocate(ctrl%perturbInp)
        ctrl%perturbInp%isRespKernelPert = .true.
        if (ctrl%tSCC) then
          call hsd_get_or_set(child, "RPA", ctrl%perturbInp%isRespKernelRPA, .false.)
        else
          ctrl%perturbInp%isRespKernelRPA = .true.
        end if
        call freqRanges(child, ctrl%perturbInp%dynKernelFreq)
      end if

      ! Perturbation with respect to atom positions
      call hsd_get_table(node, "CoordDerivatives", child, stat, auto_wrap=.true.)
      if (associated(child)) then
        if (.not.allocated(ctrl%perturbInp)) allocate(ctrl%perturbInp)
      #:if WITH_MPI
        call dftbp_error(node, "CoordDerivatives not currently available for MPI enabled code")
        ctrl%perturbInp%isAtomCoordPerturb = .false.
      #:else
        ctrl%perturbInp%isAtomCoordPerturb = .true.
      #:endif
      end if

      if (allocated(ctrl%perturbInp)) then
        call hsd_get_or_set(node, "PerturbDegenTol", ctrl%perturbInp%tolDegenDFTBPT, 1.0E-9_dp,&
            & child=child)
        call hsd_get_attrib(node, "PerturbDegenTol", modifier, stat)
        if (stat /= HSD_STAT_OK) modifier = ""
        call convertUnitHsd(modifier, energyUnits, child, ctrl%perturbInp%tolDegenDFTBPT)
        if (ctrl%perturbInp%tolDegenDFTBPT < epsilon(0.0_dp)) then
          call dftbp_error(child, "Perturbation degeneracy tolerance must be above machine&
              & epsilon")
        end if
        isEtaNeeded = .false.
        if (allocated(ctrl%perturbInp%dynEFreq)) then
          if (any(ctrl%perturbInp%dynEFreq /= 0.0_dp)) then
            isEtaNeeded = .true.
          end if
        end if
        if (allocated(ctrl%perturbInp%dynKernelFreq)) then
          if (any(ctrl%perturbInp%dynKernelFreq /= 0.0_dp)) then
            isEtaNeeded = .true.
          end if
        end if
        if (isEtaNeeded) then
          allocate(ctrl%perturbInp%etaFreq)
          call hsd_get_or_set(node, "PerturbEta", ctrl%perturbInp%etaFreq, 1.0E-8_dp, child=child)
          if (ctrl%perturbInp%etaFreq < epsilon(0.0_dp)) then
            call dftbp_error(child, "Imaginary constant for finite frequency perturbation too&
                & small")
          end if
        end if

      end if

      if (allocated(ctrl%perturbInp)) then
        call maxSelfConsIterations(node, ctrl, "MaxPerturbIter", ctrl%perturbInp%maxPerturbIter)
        if (ctrl%tScc) then
          call hsd_get_or_set(node, "PerturbSccTol", ctrl%perturbInp%perturbSccTol, 1.0e-5_dp)
          ! self consistency required, or not, to proceed with perturbation
          call hsd_get_or_set(node, "ConvergedPerturb", ctrl%perturbInp%isPerturbConvRequired,&
              & .true.)
        end if
      end if

    end if

    if (tHaveDensityMatrix) then

      ! Is this compatible with Poisson solver use?
      call readElectrostaticPotential(node, geo, ctrl)

      call hsd_get_or_set(node, "MullikenAnalysis", ctrl%tPrintMulliken, .true.)
      if (ctrl%tPrintMulliken) then
        call hsd_get_or_set(node, "WriteNetCharges", ctrl%tPrintNetAtomCharges, .false.)
        if (ctrl%tPrintNetAtomCharges) then
          ctrl%tNetAtomCharges = .true.
        end if
        call hsd_get_table(node, "CM5", child, stat, auto_wrap=.true.)
        if (associated(child)) then
          allocate(ctrl%cm5Input)
          call readCM5(child, ctrl%cm5Input, geo)
        end if
      end if
      call hsd_get_or_set(node, "AtomResolvedEnergies", ctrl%tAtomicEnergy, .false.)

      if (allocated(ctrl%solvInp)) then
        call hsd_get_or_set(node, "writeCosmoFile", ctrl%tWriteCosmoFile, &
            & allocated(ctrl%solvInp%cosmoInp), child=child)
        if (ctrl%tWriteCosmoFile .and. .not.allocated(ctrl%solvInp%cosmoInp)) then
          call dftbp_error(child, "Cosmo file can only be written for Cosmo calculations")
        end if
      end if

      call hsd_get_or_set(node, "PrintForces", ctrl%tPrintForces, .false.)

    else

      ctrl%tPrintMulliken = .false.
      ctrl%tAtomicEnergy = .false.
      ctrl%tPrintForces = .false.

    end if


  #:if WITH_TRANSPORT
    call hsd_get_table(node, "TunnelingAndDOS", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      if (.not.transpar%defined) then
        call error("Block TunnelingAndDos requires Transport block.")
      end if
      if (.not.transpar%taskUpload) then
        call error("Block TunnelingAndDos not compatible with task=contactHamiltonian")
      end if
      if (.not. allocated(orb)) then
        call error("Orbital information from SK-files missing (xTB Hamiltonian not compatible&
            & with transport yet)")
      end if
      call readTunAndDos(child, orb, geo, tundos, transpar, ctrl%tempElec)
    else
      if (ctrl%solver%isolver == electronicSolverTypes%OnlyTransport) then
        call dftbp_error(node, "The TransportOnly solver requires a TunnelingAndDos block to be&
            & present.")
      end if
    endif
  #:endif

  end subroutine readAnalysis


  !> Frequency ranges for response calculations
  subroutine freqRanges(node, frequencies)

    !> Node to parse
    type(hsd_table), pointer :: node

    !> Frequencies, 0 being static
    real(dp), allocatable, intent(inout) :: frequencies(:)

    type(hsd_table), pointer :: child
    character(len=:), allocatable :: modifier
    real(dp), allocatable :: tmpFreqs(:)
    integer :: nFreq, iFreq, jFreq
    real(dp) :: tmp3R(3)
    logical :: isStatic
    integer :: stat

    call hsd_get_or_set(node, "Static", isStatic, .true.)
    if (isStatic) then
      call growFreqArray(frequencies, 1)
      ! should already be zero, but just in case:
      frequencies(:) = 0.0_dp
    end if

    call hsd_get_table(node, "Frequencies", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      call hsd_get_attrib(node, "Frequencies", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call hsd_get(node, "Frequencies", tmpFreqs)
      nFreq = size(tmpFreqs)
      if (nFreq > 0) then
        if (allocated(frequencies)) then
          iFreq = size(frequencies)
        else
          iFreq = 0
        end if
        call growFreqArray(frequencies, nFreq)
        frequencies(iFreq+1:iFreq+nFreq) = tmpFreqs
        call convertUnitHsd(modifier, freqUnits, child, frequencies(iFreq+1:iFreq+nFreq))
      end if
      if (any(frequencies < 0.0_dp)) then
        call dftbp_error(child, "Negative driving frequency requested")
      end if
    end if

    call hsd_get_table(node, "FrequencyRange", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      call hsd_get_attrib(node, "FrequencyRange", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      call hsd_get(node, "FrequencyRange", tmpFreqs)
      if (size(tmpFreqs) == 3) then
        tmp3R = tmpFreqs
        call convertUnitHsd(modifier, freqUnits, child, tmp3R)
        if (any(tmp3R(:2) < 0.0_dp)) then
          call dftbp_error(child, "Negative values in dynamic frequency range.")
        end if
        if (abs(tmp3R(3)) <= epsilon(0.0_dp)) then
          call dftbp_error(child, "Increase step size in dynamic frequency range.")
        end if
        ! how many frequencies in the specified range?
        nFreq = max(int((tmp3R(2)-tmp3R(1))/tmp3R(3))+1,0)
        if (allocated(frequencies)) then
          iFreq = size(frequencies)
        else
          iFreq = 0
        end if
        call growFreqArray(frequencies, nFreq)
        do jFreq = 1, nFreq
          frequencies(iFreq+jFreq) = tmp3R(1) + (jFreq-1) * tmp3R(3)
        end do
      else
        call dftbp_error(child,"Malformed frequency range.")
      end if
      if (any(frequencies < 0.0_dp)) then
        call dftbp_error(child, "Negative driving frequency requested")
      end if
    end if

  end subroutine freqRanges


  !> Resize array, retaining values at start
  subroutine growFreqArray(freq, nFreq)

    !> Array to expand
    real(dp), allocatable, intent(inout) :: freq(:)

    !> Number of extra elements
    integer, intent(in) :: nFreq

    real(dp), allocatable :: tmpFreq(:)
    integer :: nElem

    if (allocated(freq)) then
      nElem = size(freq)
      call move_alloc(freq, tmpFreq)
    else
      nElem =0
    end if
    allocate(freq(nElem + nFreq))
    if (nElem > 0) then
      freq(:nElem) = tmpFreq
    end if
    freq(nElem+1:) = 0.0_dp

  end subroutine growFreqArray


  !> Read in settings that are influenced by those read from Options{} but belong in Analysis{}
  subroutine readLaterAnalysis(node, ctrl)

    !> Node to parse
    type(hsd_table), pointer :: node

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl


    logical :: tPrintEigVecs

    tPrintEigVecs = ctrl%tPrintEigVecs
    if (allocated(ctrl%lrespini)) tPrintEigvecs = tPrintEigvecs .or. ctrl%lrespini%tPrintEigVecs
    if (tPrintEigVecs) then
      call hsd_get_or_set(node, "EigenvectorsAsText", ctrl%tPrintEigVecsTxt, .false.)
    end if

  end subroutine readLaterAnalysis



  !> Reads the settings for electrostatic potential plotting
  subroutine readElectrostaticPotential(node, geo, ctrl)

    !> Node containing optional electrostatic settings
    type(hsd_table), pointer, intent(in) :: node

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Control structure
    type(TControl), intent(inout) :: ctrl

    type(hsd_table), pointer :: child, child2, child3
    character(len=:), allocatable :: buffer, modifier
    integer :: nMatRows, nMatCols
    integer :: stat

    call hsd_get_table(node, "ElectrostaticPotential", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      return
    end if

    if (.not. ctrl%tSCC) then
      call error("Electrostatic potentials only available in an SCC calculation")
    end if
    allocate(ctrl%elStatPotentialsInp)
    call hsd_get_or_set(child, "OutputFile", buffer, "ESP.dat")
    ctrl%elStatPotentialsInp%espOutFile = unquote(buffer)
    ctrl%elStatPotentialsInp%tAppendEsp = .false.
    if (ctrl%isGeoOpt .or. ctrl%tMD) then
      call hsd_get_or_set(child, "AppendFile", ctrl%elStatPotentialsInp%tAppendEsp, .false.)
    end if
    ! discrete points
    call hsd_get_table(child, "Points", child3, stat, auto_wrap=.true.)
    call hsd_get_attrib(child, "Points", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    if (associated(child3)) then
      call hsd_get_matrix(child, "Points", ctrl%elStatPotentialsInp%espGrid, nMatRows, nMatCols,&
          & order="column-major", stat=stat)
      if (stat == HSD_STAT_OK .and. nMatRows > 0 .and. nMatCols > 0) then
        if (geo%tPeriodic .and. (modifier == "F" .or. modifier == "f")) then
          ctrl%elStatPotentialsInp%espGrid = matmul(geo%latVecs, ctrl%elStatPotentialsInp%espGrid)
        else
          call convertUnitHsd(modifier, lengthUnits, child3,&
              & ctrl%elStatPotentialsInp%espGrid)
        end if
      end if
    end if

    ! grid specification for points instead
    call hsd_get_table(child, "Grid", child2, stat, auto_wrap=.true.)
    if (associated(child2)) then
      call hsd_get_attrib(child, "Grid", modifier, stat)
      if (stat /= HSD_STAT_OK) modifier = ""
      if (allocated(ctrl%elStatPotentialsInp%espGrid)) then
        call error("Both grid and point specification not both currently possible")
      end if
      if (geo%tPeriodic) then
        call readGrid(ctrl%elStatPotentialsInp%espGrid, child2, modifier,&
            & latVecs=geo%latVecs, nPoints=ctrl%elStatPotentialsInp%gridDimensioning,&
            & origin=ctrl%elStatPotentialsInp%origin,&
            & axes=ctrl%elStatPotentialsInp%axes)
      else
        call readGrid(ctrl%elStatPotentialsInp%espGrid, child2, modifier,&
            & nPoints=ctrl%elStatPotentialsInp%gridDimensioning,&
            & origin=ctrl%elStatPotentialsInp%origin,&
            & axes=ctrl%elStatPotentialsInp%axes)
      end if
    end if
    if (.not.allocated(ctrl%elStatPotentialsInp%espGrid)) then
      call dftbp_error(child,"Either a grid or set of points must be specified")
    end if
    call hsd_get_or_set(child, "Softening", ctrl%elStatPotentialsInp%softenESP, 1.0E-6_dp,&
        & child=child2)
    call hsd_get_attrib(child, "Softening", modifier, stat)
    if (stat /= HSD_STAT_OK) modifier = ""
    call convertUnitHsd(modifier, lengthUnits, child2, ctrl%elStatPotentialsInp%softenEsp)

  end subroutine readElectrostaticPotential


  !> Read in a grid specification
  subroutine readGrid(points, node, modifier, latVecs, nPoints, origin, axes)

    !> Points in the grid
    real(dp), allocatable, intent(out) :: points(:,:)

    !> input data to parse
    type(hsd_table), pointer, intent(in) :: node

    !> unit modifier for the grid
    character(len=*), intent(in) :: modifier

    !> geometry of the system
    real(dp), intent(in), optional :: latVecs(:,:)

    !> Number of grid points in each direction, if required
    integer, intent(out), optional :: nPoints(3)

    !> origin of grid if required
    real(dp), intent(out), optional :: origin(3)

    !> axes of the grid if required
    real(dp), intent(out), optional :: axes(3,3)

    type(hsd_table), pointer :: child
    real(dp) :: r3Tmp(3), r3Tmpb(3)
    integer :: i3Tmp(3), iPt, ii, jj, kk
    logical :: tPeriodic
    real(dp) :: axes_(3,3), r33Tmp(3,3)
    integer :: stat

    tPeriodic = present(latvecs)

    if (.not.tPeriodic .and. (modifier == "F" .or. modifier == "f")) then
      call dftbp_error(node, "Fractional grid specification only available for periodic&
          & geometries")
    end if

    block
      real(dp), allocatable :: tmpR(:)
      call hsd_get(node, "Spacing", tmpR, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required array: 'Spacing'")
      r3Tmp(:min(size(tmpR),3)) = tmpR(:min(size(tmpR),3))
    end block
    block
      real(dp), allocatable :: tmpR(:)
      call hsd_get(node, "Origin", tmpR, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required array: 'Origin'")
      r3Tmpb(:min(size(tmpR),3)) = tmpR(:min(size(tmpR),3))
    end block
    block
      integer, allocatable :: tmpI(:)
      call hsd_get(node, "GridPoints", tmpI, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required array: 'GridPoints'")
      i3Tmp(:min(size(tmpI),3)) = tmpI(:min(size(tmpI),3))
    end block
    call hsd_get_table(node, "GridPoints", child, stat)
    if (.not. associated(child)) child => node
    if (any(i3Tmp < 1)) then
      call dftbp_error(child,"Grid must be at least 1x1x1")
    end if
    if (any(abs(r3Tmp) < epsilon(1.0_dp) .and. i3Tmp > 1)) then
      call dftbp_error(child,"Grid spacings must be non-zero")
    end if
    allocate(points(3,product(i3Tmp)))
    if (present(nPoints)) then
      nPoints = i3Tmp
    end if

    !  length not fraction modifier
    if (.not.(tPeriodic .and. (modifier == "F" .or. modifier == "f"))) then
      call convertUnitHsd(modifier, lengthUnits, child, r3Tmp)
      call convertUnitHsd(modifier, lengthUnits, child, r3Tmpb)
    end if

    points = 0.0_dp
    iPt = 0
    do ii = 0, i3Tmp(1)-1
      do jj = 0, i3Tmp(2)-1
        do kk = 0, i3Tmp(3)-1
          iPt = iPt + 1
          points(1,iPt) = ii * r3Tmp(1) + r3Tmpb(1)
          points(2,iPt) = jj * r3Tmp(2) + r3Tmpb(2)
          points(3,iPt) = kk * r3Tmp(3) + r3Tmpb(3)
        end do
      end do
    end do

    ! transformation matrix on directions, could use a 4x4 homogeneous coordinate transform instead
    if (.not.(modifier == "F" .or. modifier == "f") .or. .not.tPeriodic) then
      r33Tmp = reshape([1,0,0,0,1,0,0,0,1],[3,3])
      block
        real(dp), allocatable :: tmpMat(:,:)
        integer :: nRows, nCols, nr, nc
        call hsd_get_matrix(node, "Directions", tmpMat, nRows, nCols, stat=stat, &
            & order="column-major")
        if (stat /= HSD_STAT_OK) then
          axes_ = r33Tmp
          call hsd_set(node, "Directions", reshape(r33Tmp, [size(r33Tmp)]))
        else
          nr = min(size(tmpMat, 1), size(axes_, 1))
          nc = min(size(tmpMat, 2), size(axes_, 2))
          axes_(:nr, :nc) = tmpMat(:nr, :nc)
        end if
      end block
      call hsd_get_table(node, "Directions", child, stat)
      if (.not. associated(child)) child => node
      if (abs(determinant33(axes_)) < epsilon(1.0_dp)) then
        call dftbp_error(child, "Dependent axis directions")
      end if
      do ii = 1, 3
        axes_(:,ii) = axes_(:,ii) / sqrt(sum(axes_(:,ii)**2))
      end do
      points = matmul(axes_,points)
      if (present(axes)) then
        axes = axes_*spread(r3Tmp,2,3)
      end if
    end if

    if (present(origin)) then
      origin = r3Tmpb
    end if

    ! Fractional specification of points
    if (tPeriodic .and. (modifier == "F" .or. modifier == "f")) then
      points = matmul(latVecs,points)
      if (present(origin)) then
        origin = matmul(latVecs,origin)
      end if
      if (present(axes)) then
        axes = latVecs * spread(r3Tmp,2,3)
      end if
    end if

  end subroutine readGrid



end module dftbp_dftbplus_parser_analysis
