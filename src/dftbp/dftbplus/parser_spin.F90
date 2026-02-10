!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Spin-related parsing routines extracted from the main parser module.
module dftbp_dftbplus_parser_spin
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_constants, only : maxL, shellNames
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_unitconversion, only : energyUnits
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_io_charmanip, only : i2c, tolower, unquote
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, getSelectedAtomIndices
  use dftbp_io_message, only : error
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_reks_reks, only : reksTypes
  use dftbp_type_commontypes, only : TOrbitals

  use dftbp_type_typegeometry, only : TGeometry
  use hsd, only : hsd_get, hsd_rename_child, hsd_get_or_set, hsd_table_ptr, hsd_get_child_tables,&
      & hsd_get_table, hsd_get_choice, hsd_get_attrib, HSD_STAT_OK, hsd_schema_t, hsd_error_t, &
      & schema_init, schema_add_field, schema_validate, schema_destroy, FIELD_OPTIONAL, &
      & FIELD_REQUIRED, FIELD_TYPE_INTEGER, FIELD_TYPE_REAL, FIELD_TYPE_LOGICAL, FIELD_TYPE_TABLE, &
      & hsd_get_name, hsd_has_value_children
  use hsd_data, only : hsd_table
#:if WITH_TRANSPORT
  use dftbp_transport_negfvars, only : TTransPar
#:endif

  implicit none

  private
  public :: readSpinPolarisation, readSpinOrbit, readMaxAngularMomentum
  public :: setupOrbitals, getInitialSpins, getInitialCharges, readSpinConstants

  !> Storage for angular shell blocks (replaces TListIntR1)
  type, public :: TAngShellBlocks
    integer :: nBlocks = 0
    integer :: capacity = 0
    integer, allocatable :: blockSize(:)
    integer, allocatable :: angMom(:,:)  ! (maxL+1, capacity)
  end type TAngShellBlocks

  public :: initAngShellBlocks, appendAngShellBlock, getAngShellBlock

contains


  !> Reads in settings for spin orbit enabled calculations
  subroutine readSpinOrbit(node, ctrl, geo, orb)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Information about the orbitals of the species/atoms in the system
    class(TOrbitals), intent(in) :: orb

    type(hsd_table), pointer :: child, child2
    character(len=:), allocatable :: modifier
    integer :: iSp
    integer :: stat

    call hsd_get_table(node, "SpinOrbit", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) then
      ctrl%tSpinOrbit = .false.
      allocate(ctrl%xi(0,0))
    else
      if (ctrl%tSpin .and. .not. ctrl%t2Component) then
        call error("Spin-orbit coupling incompatible with collinear spin.")
      end if

      ctrl%tSpinOrbit = .true.
      ctrl%t2Component = .true.

      call hsd_get_or_set(child, "Dual", ctrl%tDualSpinOrbit, .true.)

      allocate(ctrl%xi(orb%mShell,geo%nSpecies), source = 0.0_dp)
      do iSp = 1, geo%nSpecies
        block
          real(dp), allocatable :: tmpArr(:)
          call hsd_get(child, geo%speciesNames(iSp), tmpArr, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing required value: '" &
              & // trim(geo%speciesNames(iSp)) // "'")
          ctrl%xi(:orb%nShell(iSp),iSp) = tmpArr(:orb%nShell(iSp))
        end block
        call hsd_get_attrib(child, geo%speciesNames(iSp), modifier, stat)
        if (stat /= HSD_STAT_OK) modifier = ""
        call hsd_get_table(child, geo%speciesNames(iSp), child2, stat, auto_wrap=.true.)
        if (.not. associated(child2)) child2 => child
        call convertUnitHsd(modifier, energyUnits, child2,&
            & ctrl%xi(:orb%nShell(iSp),iSp))
      end do

      ! -- Schema validation (warnings only) --
      block
        type(hsd_schema_t) :: schema
        type(hsd_error_t), allocatable :: schemaErrors(:)
        integer :: iErr

        call schema_init(schema, name="SpinOrbit")
        call schema_add_field(schema, "Dual", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_validate(schema, child, schemaErrors)
        if (size(schemaErrors) > 0) then
          do iErr = 1, size(schemaErrors)
            call dftbp_warning(child, "[schema] " // schemaErrors(iErr)%message)
          end do
        end if
        call schema_destroy(schema)
      end block
    end if

  end subroutine readSpinOrbit


  !> Read in maximal angular momenta or selected shells
  subroutine readMaxAngularMomentum(node, geo, angShells)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> List containing the angular momenta of the shells
    type(TAngShellBlocks), allocatable, intent(out) :: angShells(:)

    type(hsd_table), pointer :: value1, child, child2
    character(len=:), allocatable :: buffer
    integer :: iSp1, ii, jj, kk
    character(lc) :: strTmp
    integer :: nShell
    character(1) :: tmpCh
    character(len=:), allocatable :: strArr(:)
    integer :: stat
    logical :: tShellIncl(4), tFound
    integer :: angShell(maxL+1), angShellOrdered(maxL+1)

    do ii = 1, maxL+1
      angShellOrdered(ii) = ii - 1
    end do
    call hsd_get_table(node, "MaxAngularMomentum", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) call dftbp_error(node, "Missing required block: 'MaxAngularMomentum'")
    allocate(angShells(geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
      call initAngShellBlocks(angShells(iSp1))
      call hsd_get_table(child, geo%speciesNames(iSp1), child2, stat, auto_wrap=.true.)
      if (.not. associated(child2)) call dftbp_error(child, "Missing '" &
          & // trim(geo%speciesNames(iSp1)) // "'")
      call hsd_get_choice(child2, "", buffer, value1, stat)
      if (.not. associated(value1)) buffer = "#text"
      select case(buffer)
      case("selectedshells")
        call hsd_get(value1, "#text", strArr, stat=stat)
        if (stat /= 0) call dftbp_error(value1, "Error reading shell selections")
        do ii = 1, size(strArr)
          strTmp = tolower(unquote(trim(strArr(ii))))
          if (len_trim(strTmp) > 4 .or. len_trim(strTmp) < 1) then
            call dftbp_error(value1, "Invalid shell selection '" &
                &// trim(strTmp) &
                &// "'. Nr. of selected shells must be between 1 and 4.")
          end if
          tShellIncl(:) = .false.
          nShell = len_trim(strTmp)
          do jj = 1, nShell
            tmpCh = strTmp(jj:jj)
            tFound = .false.
            do kk = 1, size(shellNames)
              if (tmpCh == trim(shellNames(kk))) then
                if (tShellIncl(kk)) then
                  call dftbp_error(value1, "Double selection of the same shell&
                      & '" // tmpCh // "' in shell selection block '" &
                      &// trim(strTmp) // "'")
                end if
                tShellIncl(kk) = .true.
                angShell(jj) = kk - 1
                tFound = .true.
                exit
              end if
            end do
            if (.not. tFound) then
              call dftbp_error(value1, "Invalid shell name '" // tmpCh // "'")
            end if
          end do
          call appendAngShellBlock(angShells(iSp1), angShell(1:nShell))
        end do

      case("#text")
        call hsd_get(child2, "#text", buffer, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value")
        strTmp = unquote(buffer)
        do jj = 1, size(shellNames)
          if (trim(strTmp) == trim(shellNames(jj))) then
            call appendAngShellBlock(angShells(iSp1), angShellOrdered(:jj))
          end if
        end do
        if (angShells(iSp1)%nBlocks < 1) then
          call dftbp_error(child2, "Invalid orbital name '" // &
              &trim(strTmp) // "'")
        end if

      case default
        call hsd_get_name(value1, buffer)
        call dftbp_error(child2, "Invalid shell specification method '" //&
            & buffer // "'")
      end select
    end do

  end subroutine readMaxAngularMomentum


  !> Setup information about the orbitals of the species/atoms from angShell lists
  subroutine setupOrbitals(orb, geo, angShells)

    !> Information about the orbitals of the species/atoms in the system
    class(TOrbitals), intent(out) :: orb

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> List containing the angular momenta of the shells
    type(TAngShellBlocks), intent(inout) :: angShells(:)

    integer :: nShell, iSp1, iSh1, ii, jj, ind
    integer :: angShell(maxL+1)

    allocate(orb%nShell(geo%nSpecies))
    allocate(orb%nOrbSpecies(geo%nSpecies))
    allocate(orb%nOrbAtom(geo%nAtom))
    orb%mOrb = 0
    orb%mShell = 0
    do iSp1 = 1, geo%nSpecies
      orb%nShell(iSp1) = 0
      orb%nOrbSpecies(iSp1) = 0
      do ii = 1, angShells(iSp1)%nBlocks
        call getAngShellBlock(angShells(iSp1), ii, angShell, nShell)
        orb%nShell(iSp1) = orb%nShell(iSp1) + nShell
        do jj = 1, nShell
          orb%nOrbSpecies(iSp1) = orb%nOrbSpecies(iSp1) &
              &+ 2 * angShell(jj) + 1
        end do
      end do
    end do
    orb%mShell = maxval(orb%nShell)
    orb%mOrb = maxval(orb%nOrbSpecies)
    orb%nOrbAtom(:) = orb%nOrbSpecies(geo%species(:))
    orb%nOrb = sum(orb%nOrbAtom)

    allocate(orb%angShell(orb%mShell, geo%nSpecies))
    allocate(orb%iShellOrb(orb%mOrb, geo%nSpecies))
    allocate(orb%posShell(orb%mShell+1, geo%nSpecies))
    orb%angShell(:,:) = 0
    do iSp1 = 1, geo%nSpecies
      ind = 1
      iSh1 = 1
      do ii = 1, angShells(iSp1)%nBlocks
        call getAngShellBlock(angShells(iSp1), ii, angShell, nShell)
        do jj = 1, nShell
          orb%posShell(iSh1, iSp1) = ind
          orb%angShell(iSh1, iSp1) = angShell(jj)
          orb%iShellOrb(ind:ind+2*angShell(jj), iSp1) = iSh1
          ind = ind + 2 * angShell(jj) + 1
          iSh1 = iSh1 + 1
        end do
        orb%posShell(iSh1, iSp1) = ind
      end do
    end do

  end subroutine setupOrbitals


  !> Spin calculation
#:if WITH_TRANSPORT
  subroutine readSpinPolarisation(node, ctrl, geo, tp)
#:else
  subroutine readSpinPolarisation(node, ctrl, geo)
#:endif

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp
  #:endif

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer
    integer :: stat

    call hsd_rename_child(node, "SpinPolarization", "SpinPolarisation")
    call hsd_get_table(node, "SpinPolarisation", child, stat, auto_wrap=.true.)
    value1 => null()
    buffer = ""
    if (associated(child)) then
      call hsd_get_choice(child, "", buffer, value1, stat)
      if (len_trim(buffer) == 0 .and. hsd_has_value_children(child)) buffer = "#text"
    end if
    select case(tolower(buffer))
    case ("")
      ctrl%tSpin = .false.
      ctrl%t2Component = .false.
      ctrl%nrSpinPol = 0.0_dp

    case ("colinear", "collinear")
      ctrl%tSpin = .true.
      ctrl%t2Component = .false.
      call hsd_get_or_set(value1, 'UnpairedElectrons', ctrl%nrSpinPol, 0.0_dp)
      call hsd_get_or_set(value1, 'RelaxTotalSpin', ctrl%tSpinSharedEf, .false.)
      if (.not. ctrl%tReadChrg) then
        call getInitialSpins(value1, geo, 1, ctrl%initialSpins)
      end if

    case ("noncolinear", "noncollinear")
      ctrl%tSpin = .true.
      ctrl%t2Component = .true.
      if (.not. ctrl%tReadChrg) then
        call getInitialSpins(value1, geo, 3, ctrl%initialSpins)
      end if

    case default
      call dftbp_error(child, "Invalid spin polarisation type '" //&
          & buffer // "'")
    end select

    ! -- Schema validation for spin polarisation choice table (warnings only) --
    if (associated(value1)) then
      block
        type(hsd_schema_t) :: schema
        type(hsd_error_t), allocatable :: schemaErrors(:)
        integer :: iErr

        call schema_init(schema, name="SpinPolarisation")
        call schema_add_field(schema, "UnpairedElectrons", FIELD_OPTIONAL, FIELD_TYPE_REAL)
        call schema_add_field(schema, "RelaxTotalSpin", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_add_field(schema, "InitialSpins", FIELD_OPTIONAL, FIELD_TYPE_TABLE)
        call schema_validate(schema, value1, schemaErrors)
        if (size(schemaErrors) > 0) then
          do iErr = 1, size(schemaErrors)
            call dftbp_warning(value1, "[schema] " // schemaErrors(iErr)%message)
          end do
        end if
        call schema_destroy(schema)
      end block
    end if

  end subroutine readSpinPolarisation


  !> Reads initial charges
  subroutine getInitialCharges(node, geo, initCharges)

    !> relevant node in input tree
    type(hsd_table), pointer :: node

    !> geometry, including atomic type information
    type(TGeometry), intent(in) :: geo

    !> initial atomic charges
    real(dp), allocatable :: initCharges(:)

    type(hsd_table), pointer :: child, child2, child3, val
    type(hsd_table_ptr), allocatable :: children(:)
    integer, allocatable :: pTmpI1(:)
    character(len=:), allocatable :: buffer
    real(dp) :: rTmp
    integer :: ii, jj, iAt
    integer :: stat

    call hsd_get_table(node, "InitialCharges", child, stat, auto_wrap=.true.)

    ! Read either all atom charges, or individual atom specifications
    child2 => null()
    if (associated(child)) then
      call hsd_get_table(child, "AllAtomCharges", child2, stat, auto_wrap=.true.)
    end if
    if (associated(child2)) then
      call hsd_get(child2, "#text", initCharges, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Error reading charges")
    else
      if (associated(child)) then
        call hsd_get_child_tables(child, "AtomCharge", children)
      else
        allocate(children(0))
      end if
      if (size(children) > 0) then
        allocate(initCharges(geo%nAtom))
        initCharges = 0.0_dp
      end if
      do ii = 1, size(children)
        child2 => children(ii)%ptr
        call hsd_get(child2, "Atoms", buffer, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value: 'Atoms'")
        call hsd_get_table(child2, "Atoms", child3, stat, auto_wrap=.true.)
        if (.not. associated(child3)) child3 => child2
        call getSelectedAtomIndices(child3, buffer, geo%speciesNames, geo%species, pTmpI1)
        call hsd_get(child2, "ChargePerAtom", rTmp, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value: 'ChargePerAtom'")
        do jj = 1, size(pTmpI1)
          iAt = pTmpI1(jj)
          if (initCharges(iAt) /= 0.0_dp) then
            call dftbp_warning(child3, "Previous setting for the charge &
                &of atom" // i2c(iAt) // " overwritten")
          end if
          initCharges(iAt) = rTmp
        end do
        deallocate(pTmpI1)
      end do
    end if

  end subroutine getInitialCharges


  !> Reads initial spins
  subroutine getInitialSpins(node, geo, nSpin, initSpins)

    !> relevant node in input data
    type(hsd_table), pointer :: node

    !> geometry, including atomic information
    type(TGeometry), intent(in) :: geo

    !> number of spin channels
    integer, intent(in) :: nSpin

    !> initial spins on return
    real(dp), allocatable :: initSpins(:,:)

    type(hsd_table), pointer :: child, child2, child3, val
    type(hsd_table_ptr), allocatable :: children(:)
    integer, allocatable :: pTmpI1(:)
    character(len=:), allocatable :: buffer
    real(dp), allocatable :: rTmp(:)
    integer :: ii, jj, iAt
    integer :: stat

    @:ASSERT(nSpin == 1 .or. nSpin == 3)

    call hsd_get_table(node, "InitialSpins", child, stat, auto_wrap=.true.)

    ! Read either all atom spins, or individual spin specifications
    child2 => null()
    if (associated(child)) then
      call hsd_get_table(child, "AllAtomSpins", child2, stat, auto_wrap=.true.)
    end if
    if (associated(child2)) then
      block
        real(dp), allocatable :: flat(:)
        call hsd_get(child2, "#text", flat, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Error reading initial spins")
        allocate(initSpins(nSpin, geo%nAtom))
        initSpins = reshape(flat, [nSpin, geo%nAtom])
      end block
    else
      if (associated(child)) then
        call hsd_get_child_tables(child, "AtomSpin", children)
      else
        allocate(children(0))
      end if
      if (size(children) > 0) then
        allocate(initSpins(nSpin, geo%nAtom))
        initSpins = 0.0_dp
      end if
      allocate(rTmp(nSpin))
      do ii = 1, size(children)
        child2 => children(ii)%ptr
        call hsd_get(child2, "Atoms", buffer, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value: 'Atoms'")
        call hsd_get_table(child2, "Atoms", child3, stat, auto_wrap=.true.)
        if (.not. associated(child3)) child3 => child2
        call getSelectedAtomIndices(child3, buffer, geo%speciesNames, geo%species, pTmpI1)
        call hsd_get(child2, "SpinPerAtom", rTmp, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required value: 'SpinPerAtom'")
        do jj = 1, size(pTmpI1)
          iAt = pTmpI1(jj)
          if (any(initSpins(:,iAt) /= 0.0_dp)) then
            call dftbp_warning(child3, "Previous setting for the spin of atom" // i2c(iAt) //&
                & " overwritten")
          end if
          initSpins(:,iAt) = rTmp
        end do
        deallocate(pTmpI1)
      end do
      deallocate(rTmp)
    end if

  end subroutine getInitialSpins


  !> Reads W values if required by settings in the Hamiltonian or the excited state
  subroutine readSpinConstants(hamNode, geo, orb, ctrl)

    !> node for Hamiltonian data
    type(hsd_table), pointer :: hamNode

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> control structure
    type(TControl), intent(inout) :: ctrl

    type(hsd_table), pointer :: child
    logical :: tLRNeedsSpinConstants, tShellResolvedW
    integer :: iSp1, nConstants
    real(dp), allocatable :: realArr(:)
    integer :: stat
    character(lc) :: strTmp

    tLRNeedsSpinConstants = .false.

    if (allocated(ctrl%lrespini)) then
      select case (ctrl%lrespini%sym)
      case ("T", "B", " ")
        tLRNeedsSpinConstants = .true.
      case ("S")
        tLRNeedsSpinConstants = .false.
      case default
      end select
    end if

    if (tLRNeedsSpinConstants .or. ctrl%tSpin .or. &
        & ctrl%reksInp%reksAlg /= reksTypes%noReks) then
      allocate(ctrl%spinW(orb%mShell, orb%mShell, geo%nSpecies))
      ctrl%spinW(:,:,:) = 0.0_dp

      call hsd_get_table(hamNode, "SpinConstants", child, stat, auto_wrap=.true.)
      if (.not. associated(child)) call dftbp_error(hamNode, "Missing required block: 'SpinConstants'")
      if (ctrl%hamiltonian == hamiltonianTypes%xtb) then
        call hsd_get_or_set(child, "ShellResolvedSpin", tShellResolvedW, .true.)
      else
        if (.not.ctrl%tShellResolved) then
          call hsd_get_or_set(child, "ShellResolvedSpin", tShellResolvedW, .false.)
        else
          tShellResolvedW = .true.
        end if
      end if

      do iSp1 = 1, geo%nSpecies
        call hsd_get(child, geo%speciesNames(iSp1), realArr, stat=stat)
        if (stat /= 0) call dftbp_error(child, "Error reading spin constants for " &
            &// trim(geo%speciesNames(iSp1)))
        nConstants = size(realArr)
        if (tShellResolvedW) then
          if (nConstants == orb%nShell(iSp1)**2) then
            ctrl%spinW(:orb%nShell(iSp1), :orb%nShell(iSp1), iSp1) =&
                & reshape(realArr(:orb%nShell(iSp1)**2), [orb%nShell(iSp1), orb%nShell(iSp1)])
          else
            write(strTmp, "(A,I0,A,I0,A,A,A)")'Expecting a ', orb%nShell(iSp1), ' x ',&
                & orb%nShell(iSp1), ' spin constant matrix for "', trim(geo%speciesNames(iSp1)),&
                & '", as ShellResolvedSpin enabled.'
            call dftbp_error(child, trim(strTmp))
          end if
        else
          if (nConstants == 1) then
            ! only one value for all atom spin constants
            ctrl%spinW(:orb%nShell(iSp1), :orb%nShell(iSp1), iSp1) = realArr(1)
          else
            write(strTmp, "(A,A,A)")'Expecting a single spin constant for "',&
                & trim(geo%speciesNames(iSp1)),'", as ShellResolvedSpin not enabled.'
            call dftbp_error(child, trim(strTmp))
          end if
        end if
      end do

      ! -- Schema validation (warnings only) --
      block
        type(hsd_schema_t) :: schema
        type(hsd_error_t), allocatable :: schemaErrors(:)
        integer :: iErr

        call schema_init(schema, name="SpinConstants")
        call schema_add_field(schema, "ShellResolvedSpin", FIELD_OPTIONAL, FIELD_TYPE_LOGICAL)
        call schema_validate(schema, child, schemaErrors)
        if (size(schemaErrors) > 0) then
          do iErr = 1, size(schemaErrors)
            call dftbp_warning(child, "[schema] " // schemaErrors(iErr)%message)
          end do
        end if
        call schema_destroy(schema)
      end block
    end if

  end subroutine readSpinConstants


  !> Initialise an angular shell block list
  subroutine initAngShellBlocks(self)

    !> Angular shell blocks to initialise
    type(TAngShellBlocks), intent(out) :: self

    self%nBlocks = 0
    self%capacity = 4
    allocate(self%blockSize(4))
    allocate(self%angMom(maxL + 1, 4))
    self%blockSize(:) = 0
    self%angMom(:,:) = 0

  end subroutine initAngShellBlocks


  !> Append a block of angular momentum values
  subroutine appendAngShellBlock(self, shells)

    !> Angular shell blocks to append to
    type(TAngShellBlocks), intent(inout) :: self

    !> Angular momentum values for this block
    integer, intent(in) :: shells(:)

    integer, allocatable :: tmpSize(:)
    integer, allocatable :: tmpMom(:,:)
    integer :: newCap

    if (self%nBlocks >= self%capacity) then
      newCap = self%capacity * 2
      allocate(tmpSize(newCap))
      allocate(tmpMom(maxL + 1, newCap))
      tmpSize(:self%capacity) = self%blockSize(:self%capacity)
      tmpMom(:, :self%capacity) = self%angMom(:, :self%capacity)
      tmpSize(self%capacity + 1:) = 0
      tmpMom(:, self%capacity + 1:) = 0
      call move_alloc(tmpSize, self%blockSize)
      call move_alloc(tmpMom, self%angMom)
      self%capacity = newCap
    end if
    self%nBlocks = self%nBlocks + 1
    self%blockSize(self%nBlocks) = size(shells)
    self%angMom(:size(shells), self%nBlocks) = shells

  end subroutine appendAngShellBlock


  !> Retrieve a block of angular momentum values (replaces intoArray)
  subroutine getAngShellBlock(self, idx, shells, nShells)

    !> Angular shell blocks to query
    type(TAngShellBlocks), intent(in) :: self

    !> Index of the block (1-based)
    integer, intent(in) :: idx

    !> Angular momentum values (output, must be at least maxL+1 in size)
    integer, intent(out) :: shells(:)

    !> Number of shells in this block
    integer, intent(out) :: nShells

    nShells = self%blockSize(idx)
    shells(:nShells) = self%angMom(:nShells, idx)

  end subroutine getAngShellBlock


end module dftbp_dftbplus_parser_spin
