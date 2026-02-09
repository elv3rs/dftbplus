!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Reads in dispersion related settings from the HSD input.
module dftbp_dftbplus_parser_dispersion
  use dftbp_common_accuracy, only : dp, mc
  use dftbp_common_constants, only : symbolToNumber
  use dftbp_common_status, only : TStatus
  use dftbp_common_unitconversion, only : chargeUnits, energyUnits, lengthUnits, volumeUnits
  use dftbp_dftb_coordnumber, only : cnType, getD3Radius, getElectronegativity, TCNInput
  use dftbp_dftb_dftd4param, only : getEeqChi, getEeqGam, getEeqKcn, getEeqRad
  use dftbp_dftb_dispersions, only : getUffValues, TDispDftD4Inp, TDispersionInp, TDispSlaKirkInp,&
      & TDispUffInp, TSimpleDftD3Input
  use dftbp_dftb_encharges, only : TEeqInput
  use dftbp_dftb_periodic, only : getCellTranslations, TNeighbourList, TNeighbourlist_init,&
      & updateNeighbourList
  use dftbp_dftbplus_specieslist, only : readSpeciesList
  use dftbp_extlibs_sdftd3, only : dampingFunction, TSDFTD3Input
  use dftbp_io_charmanip, only : tolower, unquote
  use hsd, only : hsd_rename_child
  use hsd_data, only : hsd_table
  use dftbp_io_hsdutils, only : getChild, getChildValue
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, textNodeName, getNodeName,&
      & getNodeHSDName, splitModifier
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_io_message, only : error, warning
  use dftbp_type_typegeometry, only : TGeometry
#:if WITH_MBD
  use dftbp_dftb_dispmbd, only : TDispMbdInp
#:endif
  implicit none

  private
  public :: readDispersion

contains


  !> Reads in dispersion related settings
  subroutine readDispersion(node, geo, input, nrChrg, tSCC)

    !> Node to parse
    type(hsd_table), pointer :: node

    !> geometry, including atomic information
    type(TGeometry), intent(in) :: geo

    !> dispersion data on exit
    type(TDispersionInp), intent(out) :: input

    !> net charge
    real(dp), intent(in) :: nrChrg

    !> SCC calculation?
    logical, intent(in) :: tScc

    type(hsd_table), pointer :: dispModel
    character(len=:), allocatable :: buffer

    call getChildValue(node, "", dispModel)
    call getNodeName(dispModel, buffer)
    select case (buffer)
    case ("slaterkirkwood")
      allocate(input%slakirk)
      call readDispSlaKirk(dispModel, geo, input%slakirk)
    case ("lennardjones")
      allocate(input%uff)
      call readDispVdWUFF(dispModel, geo, input%uff)
    case ("dftd3")
      allocate(input%dftd3)
      call readDFTD3(dispModel, geo, input%dftd3)
    case ("simpledftd3")
      allocate(input%sdftd3)
      call readSimpleDFTD3(dispModel, geo, input%sdftd3)
    case ("dftd4")
      allocate(input%dftd4)
      call readDispDFTD4(dispModel, geo, input%dftd4, nrChrg)
    case ("ts")
  #:if WITH_MBD
      allocate(input%mbd)
      call readDispTs(dispModel, input%mbd)
  #:else
      call dftbp_error(node, "Program must be compiled with the mbd library for TS-dispersion")
  #:endif
    case ("mbd")
  #:if WITH_MBD
      allocate(input%mbd)
      call readDispMbd(dispModel, input%mbd)
  #:else
      call dftbp_error(node, "Program must be compiled with the mbd library for MBD-dispersion")
  #:endif
    case default
      call dftbp_error(node, "Invalid dispersion model name.")
    end select

  end subroutine readDispersion


  !> Reads in the dispersion input data for the Slater-Kirkwood dispersion model
  subroutine readDispSlaKirk(node, geo, input)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Contains the input for the dispersion module on exit
    type(TDispSlaKirkInp), intent(out) :: input

    type(hsd_table), pointer :: value1, value2, child, child2, child3
    character(len=:), allocatable :: buffer, modifier, modifier2
    character(len=mc) :: modifiers(3)
    real(dp), allocatable :: tmpR2(:,:), tmp2R2(:,:), rCutoffs(:)
    real(dp) :: mCutoff, rTmp
    integer :: iAt1, iAt2f, iSp1, iSp2, iNeigh
    integer, allocatable :: nNeighs(:)
    real(dp), allocatable :: cellVec(:,:), rCellVec(:,:)
    real(dp), allocatable :: coords(:,:)
    integer, allocatable :: img2CentCell(:), iCellVec(:)
    integer :: nAllAtom
    type(TNeighbourList) :: neighs
    type(TStatus) :: errStatus

    allocate(tmpR2(3, geo%nAtom))
    allocate(input%polar(geo%nAtom))
    allocate(input%rWaals(geo%nAtom))
    allocate(input%charges(geo%nAtom))
    call getChildValue(node, "PolarRadiusCharge", value1, child=child, modifier=modifier)
    call getNodeName(value1, buffer)
    select case (buffer)
    case (textNodeName)
      call getChildValue(child, "", tmpR2, modifier=modifier)
      if (len(modifier) > 0) then
        call splitModifier(modifier, child, modifiers)
        call convertUnitHsd(modifiers(1), volumeUnits, child, tmpR2(1,:),&
            &.false.)
        call convertUnitHsd(modifiers(2), lengthUnits, child, tmpR2(2,:),&
            &.false.)
        call convertUnitHsd(modifiers(3), chargeUnits, child, tmpR2(3,:),&
            &.false.)
      end if

    case ("hybriddependentpol")
      if (len(modifier) > 0) then
        call dftbp_error(child, "PolarRadiusCharge is not allowed to carry &
            &a modifier, if the HybridDependentPol method is used.")
      end if
      allocate(rCutoffs(geo%nSpecies))
      allocate(tmp2R2(13, geo%nSpecies))
      do iSp1 = 1, geo%nSpecies
        call getChildValue(value1, geo%speciesNames(iSp1), value2, &
            &child=child2, dummyValue=.true.)
        call getChildValue(child2, "CovalentRadius", rCutoffs(iSp1), &
            &modifier=modifier2, child=child3)
        call convertUnitHsd(modifier2, lengthUnits, child3, &
            &rCutoffs(iSp1))
        call hsd_rename_child(child2, "HybridPolarizations", "HybridPolarisations")
        call getChildValue(child2, "HybridPolarisations", tmp2R2(:, iSp1), &
            &modifier=modifier2, child=child3)
        if (len(modifier2) > 0) then
          call splitModifier(modifier2, child, modifiers)
          call convertUnitHsd(modifiers(1), volumeUnits, child, &
              &tmp2R2(1:6, iSp1), .false.)
          call convertUnitHsd(modifiers(2), lengthUnits, child, &
              &tmp2R2(7:12, iSp1), .false.)
          call convertUnitHsd(modifiers(3), chargeUnits, child, &
              &tmp2R2(13, iSp1), .false.)
        end if
      end do
      mCutoff = 2.0_dp * maxval(rCutoffs)
      if (geo%tPeriodic) then
        call getCellTranslations(cellVec, rCellVec, geo%latVecs, geo%recVecs2p, mCutoff)
      else
        allocate(cellVec(3, 1))
        allocate(rCellVec(3, 1))
        cellVec(:, 1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
        rCellVec(:, 1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
      end if
      call TNeighbourlist_init(neighs, geo%nAtom, 10)
      if (geo%tPeriodic) then
        ! Make some guess for the nr. of all interacting atoms
        nAllAtom = int((real(geo%nAtom, dp)**(1.0_dp/3.0_dp) + 3.0_dp)**3)
      else
        nAllAtom = geo%nAtom
      end if
      allocate(coords(3, nAllAtom))
      allocate(img2CentCell(nAllAtom))
      allocate(iCellVec(nAllAtom))
      call updateNeighbourList(coords, img2CentCell, iCellVec, neighs, nAllAtom, geo%coords,&
          & mCutoff, rCellVec, errStatus)
      if (errStatus%hasError()) then
        call error(errStatus%message)
      end if
      allocate(nNeighs(geo%nAtom))
      nNeighs(:) = 0
      do iAt1 = 1, geo%nAtom
        iSp1 = geo%species(iAt1)
        do iNeigh = 1, neighs%nNeighbour(iAt1)
          iAt2f = img2CentCell(neighs%iNeighbour(iNeigh, iAt1))
          iSp2 = geo%species(iAt2f)
          rTmp = rCutoffs(iSp1) + rCutoffs(iSp2)
          if (neighs%neighDist2(iNeigh, iAt1) <= rTmp**2) then
            nNeighs(iAt1) = nNeighs(iAt1) + 1
            nNeighs(iAt2f) = nNeighs(iAt2f) + 1
          end if
        end do
      end do
      do iAt1 = 1, geo%nAtom
        iSp1 = geo%species(iAt1)
        if (nNeighs(iAt1) <= 4 ) then
          tmpR2(1, iAt1) = tmp2R2(1+nNeighs(iAt1), iSp1)
          tmpR2(2, iAt1) = tmp2R2(7+nNeighs(iAt1), iSp1)
        else
          tmpR2(1, iAt1) = tmp2R2(6, iSp1)
          tmpR2(2, iAt1) = tmp2R2(12, iSp1)
        end if
        tmpR2(3, iAt1) = tmp2R2(13, iSp1)
      end do

    case default
      call dftbp_error(value1, "Invalid method for PolarRadiusCharge.")
    end select

    input%polar(:) = tmpR2(1,:)
    input%rWaals(:) = tmpR2(2,:)
    input%charges(:) = tmpR2(3,:)

  end subroutine readDispSlaKirk


  !> Reads in initialization data for the UFF dispersion model
  subroutine readDispVdWUFF(node, geo, input)

    !> Node to process
    type(hsd_table), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit
    type(TDispUffInp), intent(out) :: input

    character(len=:), allocatable :: buffer
    type(hsd_table), pointer :: child, value1, child2
    integer :: iSp
    logical :: found

    call getChildValue(node, "Parameters", value1, child=child)
    allocate(input%distances(geo%nSpecies))
    allocate(input%energies(geo%nSpecies))
    call getNodeName(value1, buffer)
    select case(buffer)
    case("uffparameters")
      do iSp = 1, geo%nSpecies
        call getUffValues(geo%speciesNames(iSp), input%distances(iSp), &
            &input%energies(iSp), found)
        if (.not. found) then
          call dftbp_error(value1, "UFF parameters for species '" // geo&
              &%speciesNames(iSp) // "' not found.")
        end if
      end do
    case default
      if (associated(value1)) value1%processed = .false.
      do iSp = 1, geo%nSpecies
        call getChild(child, geo%speciesNames(iSp), child2)
        call getChildValue(child2, "Distance", input%distances(iSp), &
            &modifier=buffer)
        call convertUnitHsd(buffer, lengthUnits, child, &
            &input%distances(iSp))
        call getChildValue(child2, "Energy", input%energies(iSp), &
            &modifier=buffer)
        call convertUnitHsd(buffer, energyUnits, child, &
            &input%energies(iSp))
      end do
    end select

  end subroutine readDispVdWUFF


  !> Reads in initialization data for the DFTD3 dispersion module
  subroutine readDFTD3(node, geo, input)

    !> Node to process.
    type(hsd_table), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit.
    type(TSDFTD3Input), intent(out) :: input

    integer :: iSp
    integer, allocatable :: izpDefault(:)
    type(hsd_table), pointer :: child, childval
    character(len=:), allocatable :: buffer
    integer, parameter :: d3MaxNum = 94
    logical :: unknownSpecies, threebody

    call getChildValue(node, "Damping", childval, child=child)
    call getNodeName(childval, buffer)
    select case (buffer)
    case ("beckejohnson")
      input%dampingFunction = dampingFunction%rational
      call getChildValue(childval, "a1", input%a1)
      call getChildValue(childval, "a2", input%a2)
    case ("zerodamping")
      input%dampingFunction = dampingFunction%zero
      call getChildValue(childval, "sr6", input%sr6)
      call getChildValue(childval, "alpha6", input%alpha6, default=14.0_dp)
    case ("modifiedzerodamping")
      input%dampingFunction = dampingFunction%mzero
      call getChildValue(childval, "sr6", input%sr6)
      call getChildValue(childval, "beta", input%beta)
      call getChildValue(childval, "alpha6", input%alpha6, default=14.0_dp)
    case default
      call getNodeHSDName(childval, buffer)
      call dftbp_error(child, "Invalid damping method '" // buffer // "'")
    end select
    call getChildValue(node, "s6", input%s6)
    call getChildValue(node, "s8", input%s8)
    call getChildValue(node, "cutoff", input%cutoff, default=sqrt(9000.0_dp), &
        & modifier=buffer, child=child)
    call convertUnitHsd(buffer, lengthUnits, child, input%cutoff)
    call getChildValue(node, "cutoffcn", input%cutoffCN, default=40.0_dp, &
        & modifier=buffer, child=child)
    call convertUnitHsd(buffer, lengthUnits, child, input%cutoffCN)
    call getChildValue(node, "threebody", threebody, default=.false.)
    input%s9 = merge(1.0_dp, 0.0_dp, threebody)
    ! D3H5 - additional H-H repulsion
    call getChildValue(node, "hhrepulsion", input%hhrepulsion, default=.false.)

    ! Initialize default atomic numbers
    allocate(izpDefault(size(geo%speciesNames)))
    do iSp = 1, size(geo%speciesNames)
      izpDefault(iSp) = symbolToNumber(geo%speciesNames(iSp))
    end do

    ! See if we find user specified overwrites for atomic numbers
    call getChild(node, "AtomicNumbers", child, requested=.false.)
    if (associated(child)) then
      allocate(input%izp(size(geo%speciesNames)))
      call readSpeciesList(child, geo%speciesNames, input%izp, default=izpDefault)
      deallocate(izpDefault)
    else
      call move_alloc(izpDefault, input%izp)
    end if

    unknownSpecies = .false.
    do iSp = 1, size(geo%speciesNames)
      if (input%izp(iSp) <= 0 .or. input%izp(iSp) > d3MaxNum) then
        unknownSpecies = .true.
        call warning("Species '"//trim(geo%speciesNames(iSp))// &
          & "' is not supported by DFT-D3")
      end if
    end do
    if (unknownSpecies) then
      call dftbp_error(node, "DFT-D3 does not support all species present")
    end if

  end subroutine readDFTD3


  !> Reads in initialization data for the simple D3 dispersion model.
  subroutine readSimpleDFTD3(node, geo, input)

    !> Node to process.
    type(hsd_table), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit.
    type(TSimpleDftD3Input), intent(out) :: input

    type(hsd_table), pointer :: child
    character(len=:), allocatable :: buffer

    call getChildValue(node, "s6", input%s6, default=1.0_dp)
    call getChildValue(node, "s8", input%s8)
    call getChildValue(node, "s10", input%s10, default=0.0_dp)
    call getChildValue(node, "a1", input%a1)
    call getChildValue(node, "a2", input%a2)
    call getChildValue(node, "alpha", input%alpha, default=14.0_dp)
    call getChildValue(node, "weightingFactor", input%weightingFactor, default=4.0_dp)
    call getChildValue(node, "cutoffInter", input%cutoffInter, default=64.0_dp, modifier=buffer,&
        & child=child)
    call convertUnitHsd(buffer, lengthUnits, child, input%cutoffInter)

    call readCoordinationNumber(node, input%cnInput, geo, "exp", 0.0_dp)

  end subroutine readSimpleDFTD3


  !> Reads in initialization data for the D4 dispersion model.
  !>
  !> The D4 dispersion model is usually constructed in a failsafe way, so
  !> it only requires to know the damping parameters s8, a1 and a2.
  !> Here we additionally require a s9, since the non-addititive contributions
  !> tend to be expensive especially in the tight-binding context, s9 = 0.0_dp
  !> will disable the calculation.
  subroutine readDispDFTD4(node, geo, input, nrChrg)

    !> Node to process.
    type(hsd_table), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit.
    type(TDispDftD4Inp), intent(out) :: input

    !> Net charge of the system.
    real(dp), intent(in) :: nrChrg

    integer :: iSp
    integer, allocatable :: izpDefault(:)
    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer
    real(dp), allocatable :: d4Chi(:), d4Gam(:), d4Kcn(:), d4Rad(:)
    integer, parameter :: d4MaxNum = 86
    logical :: unknownSpecies

    call getChildValue(node, "s6", input%s6, default=1.0_dp)
    call getChildValue(node, "s8", input%s8)
    call getChildValue(node, "s9", input%s9)
    call getChildValue(node, "s10", input%s10, default=0.0_dp)
    call getChildValue(node, "a1", input%a1)
    call getChildValue(node, "a2", input%a2)
    call getChildValue(node, "alpha", input%alpha, default=16.0_dp)
    call getChildValue(node, "WeightingFactor", input%weightingFactor, default=6.0_dp)
    call getChildValue(node, "ChargeSteepness", input%chargeSteepness, default=2.0_dp)
    call getChildValue(node, "ChargeScale", input%chargeScale, default=3.0_dp)
    call getChildValue(node, "CutoffInter", input%cutoffInter, default=64.0_dp, modifier=buffer,&
        & child=child)
    call convertUnitHsd(buffer, lengthUnits, child, input%cutoffInter)
    call getChildValue(node, "CutoffThree", input%cutoffThree, default=40.0_dp, modifier=buffer,&
        & child=child)
    call convertUnitHsd(buffer, lengthUnits, child, input%cutoffThree)

    call getChildValue(node, "ChargeModel", value1, "EEQ", child=child)
    call getNodeName(value1, buffer)
    select case(buffer)
    case default
      call dftbp_error(value1, "Unknown method '"//buffer//"' for ChargeModel")
    case ("selfconsistent")
      input%selfConsistent = .true.
    case ("eeq")
      allocate(input%eeqInput)
      allocate(d4Chi(geo%nSpecies))
      d4Chi(:) = getEeqChi(geo%speciesNames)
      allocate(d4Gam(geo%nSpecies))
      d4Gam(:) = getEeqGam(geo%speciesNames)
      allocate(d4Kcn(geo%nSpecies))
      d4Kcn(:) = getEeqKcn(geo%speciesNames)
      allocate(d4Rad(geo%nSpecies))
      d4Rad(:) = getEeqRad(geo%speciesNames)
      call readEeqModel(value1, input%eeqInput, geo, nrChrg, d4Chi, d4Gam, d4Kcn, d4Rad)
    end select

    ! Initialize default atomic numbers
    allocate(izpDefault(size(geo%speciesNames)))
    do iSp = 1, size(geo%speciesNames)
      izpDefault(iSp) = symbolToNumber(geo%speciesNames(iSp))
    end do

    ! See if we find user specified overwrites for atomic numbers
    call getChild(node, "AtomicNumbers", child, requested=.false.)
    if (associated(child)) then
      allocate(input%izp(size(geo%speciesNames)))
      call readSpeciesList(child, geo%speciesNames, input%izp, default=izpDefault)
      deallocate(izpDefault)
    else
      call move_alloc(izpDefault, input%izp)
    end if

    call readCoordinationNumber(node, input%cnInput, geo, "Cov", 0.0_dp)

    unknownSpecies = .false.
    do iSp = 1, size(geo%speciesNames)
      if (input%izp(iSp) <= 0 .or. input%izp(iSp) > d4MaxNum) then
        unknownSpecies = .true.
        call warning("Species '"//trim(geo%speciesNames(iSp))// &
          & "' is not supported by DFT-D4")
      end if
    end do
    if (unknownSpecies) then
      call dftbp_error(node, "DFT-D4 does not support all species present")
    end if

  end subroutine readDispDFTD4


  !> Read settings regarding the EEQ charge model
  subroutine readEeqModel(node, input, geo, nrChrg, kChiDefault, kGamDefault, &
      & kKcnDefault, kRadDefault)

    !> Node to process.
    type(hsd_table), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit.
    type(TEeqInput), intent(out) :: input

    !> Net charge of the system.
    real(dp), intent(in) :: nrChrg

    !> Electronegativities default values
    real(dp), intent(in) :: kChiDefault(:)

    !> Chemical hardnesses default values
    real(dp), intent(in) :: kGamDefault(:)

    !> CN scaling default values
    real(dp), intent(in) :: kKcnDefault(:)

    !> Charge widths default values
    real(dp), intent(in) :: kRadDefault(:)

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer

    input%nrChrg = nrChrg

    allocate(input%chi(geo%nSpecies))
    allocate(input%gam(geo%nSpecies))
    allocate(input%kcn(geo%nSpecies))
    allocate(input%rad(geo%nSpecies))

    call getChildValue(node, "Chi", value1, "Defaults", child=child)
    call getNodeName(value1, buffer)
    select case(buffer)
    case default
      call dftbp_error(child, "Unknown method '"//buffer//"' for chi")
    case ("defaults")
      call readSpeciesList(value1, geo%speciesNames, input%chi, default=kChiDefault)
    case ("values")
      call readSpeciesList(value1, geo%speciesNames, input%chi)
    end select

    call getChildValue(node, "Gam", value1, "Defaults", child=child)
    call getNodeName(value1, buffer)
    select case(buffer)
    case default
      call dftbp_error(child, "Unknown method '"//buffer//"' for gam")
    case ("defaults")
      call readSpeciesList(value1, geo%speciesNames, input%gam, default=kGamDefault)
    case ("values")
      call readSpeciesList(value1, geo%speciesNames, input%gam)
    end select

    call getChildValue(node, "Kcn", value1, "Defaults", child=child)
    call getNodeName(value1, buffer)
    select case(buffer)
    case default
      call dftbp_error(child, "Unknown method '"//buffer//"' for kcn")
    case ("defaults")
      call readSpeciesList(value1, geo%speciesNames, input%kcn, default=kKcnDefault)
    case ("values")
      call readSpeciesList(value1, geo%speciesNames, input%kcn)
    end select

    call getChildValue(node, "Rad", value1, "Defaults", child=child)
    call getNodeName(value1, buffer)
    select case(buffer)
    case default
      call dftbp_error(child, "Unknown method '"//buffer//"' for rad")
    case ("defaults")
      call readSpeciesList(value1, geo%speciesNames, input%rad, default=kRadDefault)
    case ("values")
      call readSpeciesList(value1, geo%speciesNames, input%rad)
    end select

    call getChildValue(node, "Cutoff", input%cutoff, default=40.0_dp, modifier=buffer,&
        & child=child)
    call convertUnitHsd(buffer, lengthUnits, child, input%cutoff)

    call getChildValue(node, "EwaldParameter", input%parEwald, 0.0_dp)
    call getChildValue(node, "EwaldTolerance", input%tolEwald, 1.0e-9_dp)

    call readCoordinationNumber(node, input%cnInput, geo, "Erf", 8.0_dp)

  end subroutine readEeqModel


  !> Read in coordination number settings
  subroutine readCoordinationNumber(node, input, geo, cnDefault, cutDefault)

    !> Node to get the information from
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TCNInput), intent(inout) :: input

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Default value for the coordination number type
    character(len=*), intent(in) :: cnDefault

    !> Default value for the maximum CN used for cutting (0 turns it off)
    real(dp), intent(in) :: cutDefault

    type(hsd_table), pointer :: value1, value2, child, child2, field
    character(len=:), allocatable :: buffer, modifier
    real(dp), allocatable :: kENDefault(:), kRadDefault(:)

    call getChildValue(node, "CoordinationNumber", value1, cnDefault, child=child)
    call getNodeName(value1, buffer)

    select case(buffer)
    case default
      call dftbp_error(child, "Invalid coordination number type specified")
    case("exp")
      input%cnType = cnType%exp
    case("erf")
      input%cnType = cnType%erf
    case("cov")
      input%cnType = cnType%cov
    case("gfn")
      input%cnType = cnType%gfn
    end select

    call getChildValue(value1, "CutCN", input%maxCN, cutDefault, &
        & child=child2)

    call getChildValue(value1, "Cutoff", input%rCutoff, 40.0_dp, &
        & modifier=modifier, child=field)
    call convertUnitHsd(modifier, lengthUnits, field, input%rCutoff)

    allocate(input%en(geo%nSpecies))
    if (input%cnType == cnType%cov) then
      call getChildValue(value1, "Electronegativities", value2, "PaulingEN", child=child2)
      call getNodeName(value2, buffer)
      select case(buffer)
      case default
        call dftbp_error(child2, "Unknown method '" // buffer //&
            & "' to generate electronegativities")
      case("paulingen")
        allocate(kENDefault(geo%nSpecies))
        kENDefault(:) = getElectronegativity(geo%speciesNames)
        call readSpeciesList(value2, geo%speciesNames, input%en, default=kENDefault)
        deallocate(kENDefault)
      case("values")
        call readSpeciesList(value2, geo%speciesNames, input%en)
      end select
      if (any(input%en <= 0.0_dp)) then
        call dftbp_error(value1, "Electronegativities are not defined for all species")
      end if
    else
      ! array is not used, but should still be populated with dummies
      input%en(:) = 0.0_dp
    end if

    allocate(input%covRad(geo%nSpecies))
    call getChildValue(value1, "Radii", value2, "CovalentRadiiD3", child=child2)
    call getNodeName(value2, buffer)
    select case(buffer)
    case default
      call dftbp_error(child2, "Unknown method '"//buffer//"' to generate radii")
    case("covalentradiid3")
      allocate(kRadDefault(geo%nSpecies))
      kRadDefault(:) = getD3Radius(geo%speciesNames)
      call readSpeciesList(value2, geo%speciesNames, input%covRad, default=kRadDefault)
      deallocate(kRadDefault)
    case("values")
      call readSpeciesList(value2, geo%speciesNames, input%covRad)
    end select

    if (any(input%covRad <= 0.0_dp)) then
      call dftbp_error(value1, "Covalent radii are not defined for all species")
    end if

  end subroutine readCoordinationNumber


#:if WITH_MBD

  !> Reads in settings for Tkatchenko-Scheffler dispersion
  subroutine readDispTs(node, input)

    !> data to parse
    type(hsd_table), pointer, intent(in) :: node

    !> control data coming back
    type(TDispMbdInp), intent(out) :: input

    character(len=:), allocatable :: buffer
    type(hsd_table), pointer :: child

    input%method = 'ts'
    call getChild(node, "EnergyAccuracy", child, requested=.false.)
    if (associated(child)) then
      call dftbp_warning(child, "The energy accuracy setting will be ignored as it is not&
          & supported/need by libMBD any more")
    end if
    call getChild(node, "ForceAccuracy", child, requested=.false.)
    if (associated(child)) then
      call dftbp_warning(child, "The force accuracy setting will be ignored as it is not&
          & supported/need by libMBD any more")
    end if
    call getChildValue(node, "Damping", input%ts_d, default=(input%ts_d))
    call getChildValue(node, "RangeSeparation", input%ts_sr, default=(input%ts_sr))
    call getChildValue(node, "ReferenceSet", buffer, 'ts', child=child)
    input%vdw_params_kind = tolower(unquote(buffer))
    call checkManyBodyDispRefName(input%vdw_params_kind, child)
    call getChildValue(node, "LogLevel", input%log_level, default=(input%log_level))

  end subroutine readDispTs


  !> Reads in many-body dispersion settings
  subroutine readDispMbd(node, input)

    !> data to parse
    type(hsd_table), pointer, intent(in) :: node

    !> control data coming back
    type(TDispMbdInp), intent(out) :: input

    character(len=:), allocatable :: buffer
    type(hsd_table), pointer :: child

    input%method = 'mbd-rsscs'
    call getChildValue(node, "Beta", input%mbd_beta, input%mbd_beta)
    call getChildValue(node, "NOmegaGrid", input%n_omega_grid, default=(input%n_omega_grid))
    call getChildValue(node, "KGrid", input%k_grid)
    call getChildValue(node, "KGridShift", input%k_grid_shift, default=(input%k_grid_shift))
    call getChildValue(node, "ReferenceSet", buffer, 'ts', child=child)
    input%vdw_params_kind = tolower(unquote(buffer))
    call checkManyBodyDispRefName(input%vdw_params_kind, child)
    call getChildValue(node, "LogLevel", input%log_level, default=(input%log_level))

  end subroutine readDispMbd


  !> Check the dispersion label matches allowed cases
  subroutine checkManyBodyDispRefName(name, node)

    !> Label
    character(*), intent(in) :: name

    !> data tree for error usage
    type(hsd_table), pointer, intent(in) :: node

    if (name /= 'ts' .and. name /= 'tssurf') then
      call dftbp_error(node, 'Invalid reference set name for TS/MBD-dispersion')
    end if

  end subroutine checkManyBodyDispRefName

#:endif


end module dftbp_dftbplus_parser_dispersion
