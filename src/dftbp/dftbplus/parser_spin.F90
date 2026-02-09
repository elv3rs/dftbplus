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
  use dftbp_io_hsdutils, only : hsd_child_list, &
      & getChild, getChildren, getChildValue, &
      & getLength, getItem1, destroyNodeList
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, getSelectedAtomIndices, &
      & textNodeName, getNodeName, getNodeHSDName, getNodeName2, hasInlineData
  use dftbp_io_message, only : error
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_reks_reks, only : reksTypes
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_linkedlist, only : append, asArray, destruct, get, init, intoArray, len, &
      & TListIntR1, TListReal, TListString
  use dftbp_type_typegeometry, only : TGeometry
  use hsd, only : hsd_rename_child
  use hsd_data, only : hsd_table
#:if WITH_TRANSPORT
  use dftbp_transport_negfvars, only : TTransPar
#:endif

  implicit none

  private
  public :: readSpinPolarisation, readSpinOrbit, readMaxAngularMomentum
  public :: setupOrbitals, getInitialSpins, getInitialCharges, readSpinConstants

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

    call getChild(node, "SpinOrbit", child, requested=.false.)
    if (.not. associated(child)) then
      ctrl%tSpinOrbit = .false.
      allocate(ctrl%xi(0,0))
    else
      if (ctrl%tSpin .and. .not. ctrl%t2Component) then
        call error("Spin-orbit coupling incompatible with collinear spin.")
      end if

      ctrl%tSpinOrbit = .true.
      ctrl%t2Component = .true.

      call getChildValue(child, "Dual", ctrl%tDualSpinOrbit, .true.)

      allocate(ctrl%xi(orb%mShell,geo%nSpecies), source = 0.0_dp)
      do iSp = 1, geo%nSpecies
        call getChildValue(child, geo%speciesNames(iSp), &
            & ctrl%xi(:orb%nShell(iSp),iSp), modifier=modifier, child=child2 )
        call convertUnitHsd(modifier, energyUnits, child2,&
            & ctrl%xi(:orb%nShell(iSp),iSp))
      end do
    end if

  end subroutine readSpinOrbit


  !> Read in maximal angular momenta or selected shells
  subroutine readMaxAngularMomentum(node, geo, angShells)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> List containing the angular momenta of the shells
    type(TListIntR1), allocatable, intent(out) :: angShells(:)

    type(hsd_table), pointer :: value1, child, child2
    character(len=:), allocatable :: buffer
    integer :: iSp1, ii, jj, kk
    character(lc) :: strTmp
    integer :: nShell
    character(1) :: tmpCh
    type(TListString) :: lStr
    logical :: tShellIncl(4), tFound
    integer :: angShell(maxL+1), angShellOrdered(maxL+1)

    do ii = 1, maxL+1
      angShellOrdered(ii) = ii - 1
    end do
    call getChild(node, "MaxAngularMomentum", child)
    allocate(angShells(geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
      call init(angShells(iSp1))
      call getChildValue(child, geo%speciesNames(iSp1), value1, child=child2)
      call getNodeName(value1, buffer)
      select case(buffer)
      case("selectedshells")
        call init(lStr)
        call getChildValue(value1, "", lStr)
        do ii = 1, len(lStr)
          call get(lStr, strTmp, ii)
          strTmp = tolower(unquote(trim(strTmp)))
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
          call append(angShells(iSp1), angShell(1:nShell))
        end do
        call destruct(lStr)

      case(textNodeName)
        call getChildValue(child2, "", buffer)
        strTmp = unquote(buffer)
        do jj = 1, size(shellNames)
          if (trim(strTmp) == trim(shellNames(jj))) then
            call append(angShells(iSp1), angShellOrdered(:jj))
          end if
        end do
        if (len(angShells(iSp1)) < 1) then
          call dftbp_error(child2, "Invalid orbital name '" // &
              &trim(strTmp) // "'")
        end if

      case default
        call getNodeHSDName(value1, buffer)
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

    !> List containing the angular momenta of the shells,
    !> must be inout, since intoArray requires inout arguments
    type(TListIntR1), intent(inout) :: angShells(:)

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
      do ii = 1, len(angShells(iSp1))
        call intoArray(angShells(iSp1), angShell, nShell, ii)
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
      do ii = 1, len(angShells(iSp1))
        call intoArray(angShells(iSp1), angShell, nShell, ii)
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

    call hsd_rename_child(node, "SpinPolarization", "SpinPolarisation")
    call getChildValue(node, "SpinPolarisation", value1, "", child=child, allowEmptyValue=.true.)
    call getNodeName2(value1, buffer)
    if (buffer == "" .and. hasInlineData(child)) buffer = textNodeName
    select case(buffer)
    case ("")
      ctrl%tSpin = .false.
      ctrl%t2Component = .false.
      ctrl%nrSpinPol = 0.0_dp

    case ("colinear", "collinear")
      ctrl%tSpin = .true.
      ctrl%t2Component = .false.
      call getChildValue(value1, 'UnpairedElectrons', ctrl%nrSpinPol, 0.0_dp)
      call getChildValue(value1, 'RelaxTotalSpin', ctrl%tSpinSharedEf, .false.)
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
      call getNodeHSDName(value1, buffer)
      call dftbp_error(child, "Invalid spin polarisation type '" //&
          & buffer // "'")
    end select

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
    type(hsd_child_list), pointer :: children
    integer, allocatable :: pTmpI1(:)
    character(len=:), allocatable :: buffer
    real(dp) :: rTmp
    integer :: ii, jj, iAt

    call getChildValue(node, "InitialCharges", val, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! Read either all atom charges, or individual atom specifications
    call getChild(child, "AllAtomCharges", child2, requested=.false.)
    if (associated(child2)) then
      allocate(initCharges(geo%nAtom))
      call getChildValue(child2, "", initCharges)
    else
      call getChildren(child, "AtomCharge", children)
      if (getLength(children) > 0) then
        allocate(initCharges(geo%nAtom))
        initCharges = 0.0_dp
      end if
      do ii = 1, getLength(children)
        call getItem1(children, ii, child2)
        call getChildValue(child2, "Atoms", buffer, child=child3, multiple=.true.)
        call getSelectedAtomIndices(child3, buffer, geo%speciesNames, geo%species, pTmpI1)
        call getChildValue(child2, "ChargePerAtom", rTmp)
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
      call destroyNodeList(children)
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
    type(hsd_child_list), pointer :: children
    integer, allocatable :: pTmpI1(:)
    character(len=:), allocatable :: buffer
    real(dp), allocatable :: rTmp(:)
    integer :: ii, jj, iAt

    @:ASSERT(nSpin == 1 .or. nSpin == 3)

    call getChildValue(node, "InitialSpins", val, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! Read either all atom spins, or individual spin specifications
    call getChild(child, "AllAtomSpins", child2, requested=.false.)
    if (associated(child2)) then
      allocate(initSpins(nSpin, geo%nAtom))
      call getChildValue(child2, "", initSpins)
    else
      call getChildren(child, "AtomSpin", children)
      if (getLength(children) > 0) then
        allocate(initSpins(nSpin, geo%nAtom))
        initSpins = 0.0_dp
      end if
      allocate(rTmp(nSpin))
      do ii = 1, getLength(children)
        call getItem1(children, ii, child2)
        call getChildValue(child2, "Atoms", buffer, child=child3, multiple=.true.)
        call getSelectedAtomIndices(child3, buffer, geo%speciesNames, geo%species, pTmpI1)
        call getChildValue(child2, "SpinPerAtom", rTmp)
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
      call destroyNodeList(children)
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
    type(TListReal) :: realBuffer
    character(lc) :: strTmp
    real(dp) :: rWork(maxval(orb%nShell)**2)

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

      call getChild(hamNode, "SpinConstants", child)
      if (ctrl%hamiltonian == hamiltonianTypes%xtb) then
        call getChildValue(child, "ShellResolvedSpin", tShellResolvedW, .true.)
      else
        if (.not.ctrl%tShellResolved) then
          call getChildValue(child, "ShellResolvedSpin", tShellResolvedW, .false.)
        else
          tShellResolvedW = .true.
        end if
      end if

      do iSp1 = 1, geo%nSpecies
        call init(realBuffer)
        call getChildValue(child, geo%speciesNames(iSp1), realBuffer)
        nConstants = len(realBuffer)
        if (tShellResolvedW) then
          if (nConstants == orb%nShell(iSp1)**2) then
            call asArray(realBuffer, rWork(:orb%nShell(iSp1)**2))
            ctrl%spinW(:orb%nShell(iSp1), :orb%nShell(iSp1), iSp1) =&
                & reshape(rWork(:orb%nShell(iSp1)**2), [orb%nShell(iSp1), orb%nShell(iSp1)])
          else
            write(strTmp, "(A,I0,A,I0,A,A,A)")'Expecting a ', orb%nShell(iSp1), ' x ',&
                & orb%nShell(iSp1), ' spin constant matrix for "', trim(geo%speciesNames(iSp1)),&
                & '", as ShellResolvedSpin enabled.'
            call dftbp_error(child, trim(strTmp))
          end if
        else
          if (nConstants == 1) then
            call asArray(realBuffer, rWork(:1))
            ! only one value for all atom spin constants
            ctrl%spinW(:orb%nShell(iSp1), :orb%nShell(iSp1), iSp1) = rWork(1)
          else
            write(strTmp, "(A,A,A)")'Expecting a single spin constant for "',&
                & trim(geo%speciesNames(iSp1)),'", as ShellResolvedSpin not enabled.'
            call dftbp_error(child, trim(strTmp))
          end if
        end if
        call destruct(realBuffer)
      end do
    end if

  end subroutine readSpinConstants


end module dftbp_dftbplus_parser_spin
