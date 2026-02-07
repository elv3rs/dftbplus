!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Reads Slater-Koster files and performs consistency checks.
module dftbp_dftbplus_parser_skfiles
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_constants, only : maxL
  use dftbp_common_globalenv, only : stdOut, stdout
  use dftbp_dftb_hybridxc, only : THybridXcSKTag
  use dftbp_dftb_repulsive_polyrep, only : TPolyRep, TPolyRepInp
  use dftbp_dftb_repulsive_splinerep, only : TSplineRep, TSplineRepInp
  use dftbp_dftb_slakocont, only : addTable, init
  use dftbp_dftb_slakoeqgrid, only : init, skEqGridNew, skEqGridOld, TSlakoEqGrid
  use dftbp_dftbplus_inputdata, only : TSlater
  use dftbp_io_message, only : error
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_linkedlist, only : get, init, intoArray, len, TListCharLc, TListIntR1
  use dftbp_type_oldskdata, only : readFromFile, TOldSKData
  implicit none

  private
  public :: readSKFiles

contains


  !> Reads Slater-Koster files.
  !> Should be replaced with a more sophisticated routine, once the new SK-format has been
  !> established.
  subroutine readSKFiles(skFiles, nSpecies, slako, orb, angShells, orbRes, skInterMeth, repPoly,&
      & truncationCutOff, hybridXcSK)

    !> List of SK file names to read in for every interaction
    type(TListCharLc), intent(inout) :: skFiles(:,:)

    !> Nr. of species in the system
    integer, intent(in) :: nSpecies

    !> Data type for slako information
    type(TSlater), intent(inout) :: slako

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> For every species, a list of rank one arrays. Each array contains the angular momenta to pick
    !> from the appropriate SK-files.
    type(TListIntR1), intent(inout) :: angShells(:)

    !> Are the Hubbard Us different for each l-shell?
    logical, intent(in) :: orbRes

    !> Method of the sk interpolation
    integer, intent(in) :: skInterMeth

    !> is this a polynomial or spline repulsive?
    logical, intent(in) :: repPoly(:,:)

    !> Distances to artificially truncate tables of SK integrals
    real(dp), intent(in), optional :: truncationCutOff

    !> If calculation uses a hybrid functional, try to read extra data from end of SK file and
    !! confirm it is a suitable parameterisation
    type(THybridXcSKTag), intent(inout), optional :: hybridXcSK

    integer :: iSp1, iSp2, nSK1, nSK2, iSK1, iSK2, ind, nInteract, iSh1
    integer :: angShell(maxL+1), nShell
    logical :: readRep, readAtomic
    character(lc) :: fileName
    real(dp), allocatable, target :: skHam(:,:), skOver(:,:)
    real(dp) :: dist
    type(TOldSKData), allocatable :: skData12(:,:), skData21(:,:)
    type(TSlakoEqGrid), allocatable :: pSlakoEqGrid1, pSlakoEqGrid2
    type(TSplineRepInp) :: repSplineIn1, repSplineIn2
    type(TPolyRepInp) :: repPolyIn1, repPolyIn2

    ! if artificially cutting the SK tables
    integer :: nEntries

    @:ASSERT(size(skFiles, dim=1) == size(skFiles, dim=2))
    @:ASSERT((size(skFiles, dim=1) > 0) .and. (size(skFiles, dim=1) == nSpecies))
    @:ASSERT(all(shape(repPoly) == shape(skFiles)))

    allocate(slako%skSelf(orb%mShell, nSpecies))
    allocate(slako%skHubbU(orb%mShell, nSpecies))
    allocate(slako%skOcc(orb%mShell, nSpecies))
    allocate(slako%mass(nSpecies))
    slako%skSelf(:,:) = 0.0_dp
    slako%skHubbU(:,:) = 0.0_dp
    slako%skOcc(:,:) = 0.0_dp

    allocate(slako%skHamCont)
    call init(slako%skHamCont, nSpecies)
    allocate(slako%skOverCont)
    call init(slako%skOverCont, nSpecies)
    allocate(slako%pairRepulsives(nSpecies, nSpecies))

    write(stdout, "(A)") "Reading SK-files:"
    lpSp1: do iSp1 = 1, nSpecies
      nSK1 = len(angShells(iSp1))
      lpSp2: do iSp2 = iSp1, nSpecies
        nSK2 = len(angShells(iSp2))
        allocate(skData12(nSK2, nSK1))
        allocate(skData21(nSK1, nSK2))
        ind = 1
        do iSK1 = 1, nSK1
          do iSK2 = 1, nSK2
            readRep = (iSK1 == 1 .and. iSK2 == 1)
            readAtomic = (iSp1 == iSp2 .and. iSK1 == iSK2)
            call get(skFiles(iSp2, iSp1), fileName, ind)
            write(stdOut, "(a)") trim(fileName)
            if (.not. present(hybridXcSK)) then
              if (readRep .and. repPoly(iSp2, iSp1)) then
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, polyRepInp=repPolyIn1)
              elseif (readRep) then
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, iSp1, iSp2,&
                    & splineRepInp=repSplineIn1)
              else
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic)
              end if
            else
              if (readRep .and. repPoly(iSp2, iSp1)) then
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, polyRepInp=repPolyIn1,&
                    & hybridXcSK=hybridXcSK)
              elseif (readRep) then
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, iSp1, iSp2,&
                    & splineRepInp=repSplineIn1, hybridXcSK=hybridXcSK)
              else
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, hybridXcSK=hybridXcSK)
              end if
            end if
            ind = ind + 1
          end do
        end do
        if (iSp1 == iSp2) then
          skData21 = skData12
          if (repPoly(iSp1, iSp2)) then
            repPolyIn2 = repPolyIn1
          else
            repSplineIn2 = repSplineIn1
          end if
          ind = 1
          do iSK1 = 1, nSK1
            call intoArray(angShells(iSp1), angShell, nShell, iSK1)
            do iSh1 = 1, nShell
              slako%skSelf(ind, iSp1) = &
                  &skData12(iSK1,iSK1)%skSelf(angShell(iSh1)+1)
              slako%skOcc(ind:ind, iSp1) = &
                  &skData12(iSK1,iSK1)%skOcc(angShell(iSh1)+1)
              slako%skHubbU(ind, iSp1) = &
                  &skData12(iSK1,iSK1)%skHubbU(angShell(iSh1)+1)
              ind = ind + 1
            end do
          end do
          if (.not. orbRes) then
            slako%skHubbU(2:,iSp1) = slako%skHubbU(1,iSp1)
          end if
          slako%mass(iSp1) = skData12(1,1)%mass
        else
          ind = 1
          do iSK2 = 1, nSK2
            do iSK1 = 1, nSK1
              readRep = (iSK1 == 1 .and. iSK2 == 1)
              call get(sKFiles(iSp1, iSp2), fileName, ind)
              if (readRep .and. repPoly(iSp1, iSp2)) then
                call readFromFile(skData21(iSK1,iSK2), fileName, readAtomic, &
                    &polyRepInp=repPolyIn2)
              elseif (readRep) then
                call readFromFile(skData21(iSK1,iSK2), fileName, readAtomic, &
                    &iSp2, iSp1, splineRepInp=repSplineIn2)
              else
                call readFromFile(skData21(iSK1,iSK2), fileName, readAtomic)
              end if
              ind = ind + 1
            end do
          end do
        end if

        ! Check for SK and repulsive consistentcy
        call checkSKCompElec(skData12, skData21, iSp1, iSp2)
        if (repPoly(iSp1, iSp2)) then
          call checkSKCompRepPoly(repPolyIn1, repPolyIn2, iSp1, iSp2)
        else
          call checkSKCompRepSpline(repSplineIn1, repSplineIn2, iSp1, iSp2)
        end if

        ! Create full H/S table for all interactions of iSp1-iSp2
        nInteract = getNSKIntegrals(iSp1, iSp2, orb)
        allocate(skHam(size(skData12(1,1)%skHam, dim=1), nInteract))
        allocate(skOver(size(skData12(1,1)%skOver, dim=1), nInteract))
        call getFullTable(skHam, skOver, skData12, skData21, angShells(iSp1),&
            & angShells(iSp2))

        ! Add H/S tables to the containers for iSp1-iSp2
        dist = skData12(1,1)%dist
        if (present(truncationCutOff)) then
          nEntries = floor(truncationCutOff / dist)
          nEntries = min(nEntries, size(skData12(1,1)%skHam, dim=1))
        else
          nEntries = size(skData12(1,1)%skHam, dim=1)
        end if
        allocate(pSlakoEqGrid1, pSlakoEqGrid2)
        call init(pSlakoEqGrid1, dist, skHam(:nEntries,:), skInterMeth)
        call init(pSlakoEqGrid2, dist, skOver(:nEntries,:), skInterMeth)
        call addTable(slako%skHamCont, pSlakoEqGrid1, iSp1, iSp2)
        call addTable(slako%skOverCont, pSlakoEqGrid2, iSp1, iSp2)
        deallocate(skHam)
        deallocate(skOver)
        if (iSp1 /= iSp2) then
          ! Heteronuclear interactions: the same for the reverse interaction
          allocate(skHam(size(skData12(1,1)%skHam, dim=1), nInteract))
          allocate(skOver(size(skData12(1,1)%skOver, dim=1), nInteract))
          call getFullTable(skHam, skOver, skData21, skData12, angShells(iSp2),&
              & angShells(iSp1))
          allocate(pSlakoEqGrid1, pSlakoEqGrid2)
          call init(pSlakoEqGrid1, dist, skHam(:nEntries,:), skInterMeth)
          call init(pSlakoEqGrid2, dist, skOver(:nEntries,:), skInterMeth)
          call addTable(slako%skHamCont, pSlakoEqGrid1, iSp2, iSp1)
          call addTable(slako%skOverCont, pSlakoEqGrid2, iSp2, iSp1)
          deallocate(skHam)
          deallocate(skOver)
        end if
        deallocate(skData12)
        deallocate(skData21)

        ! Create repulsive container

        ! Add repulsives to the containers.
        if (repPoly(iSp2, iSp1)) then
          slako%pairRepulsives(iSp2, iSp1)%item = TPolyRep(repPolyIn1)
        else
          slako%pairRepulsives(iSp2, iSp1)%item = TSplineRep(repSplineIn1)
          deallocate(repSplineIn1%xStart)
          deallocate(repSplineIn1%spCoeffs)
        end if
        if (iSp1 /= iSp2) then
          if (repPoly(iSp1, iSp2)) then
            slako%pairRepulsives(iSp1, iSp2)%item = TPolyRep(repPolyIn2)
          else
            slako%pairRepulsives(iSp1, iSp2)%item = TSplineRep(repSplineIn2)
            deallocate(repSplineIn2%xStart)
            deallocate(repSplineIn2%spCoeffs)
          end if
        end if
      end do lpSp2
    end do lpSp1
    write(stdout, "(A)") "Done."

  end subroutine readSKFiles


  !> Checks if the provided set of SK-tables for a the interactions A-B and B-A are consistent
  subroutine checkSKCompElec(skData12, skData21, sp1, sp2)

    !> Slater-Koster integral set for the interaction A-B
    type(TOldSKData), intent(in), target :: skData12(:,:)

    !> Slater-Koster integral set for the interaction B-A
    type(TOldSKData), intent(in), target :: skData21(:,:)

    !> Species number for A (for error messages)
    integer, intent(in) :: sp1

    !> Species number for B (for error messages)
    integer, intent(in) :: sp2

    integer :: iSK1, iSK2, nSK1, nSK2
    integer :: nGrid
    real(dp) :: dist
    type(TOldSKData), pointer :: pSK12, pSK21
    character(lc) :: errorStr

    nSK1 = size(skData12, dim=2)
    nSK2 = size(skData12, dim=1)

    @:ASSERT(size(skData21, dim=1) == nSK1)
    @:ASSERT(size(skData21, dim=2) == nSK2)

    nGrid = skData12(1,1)%nGrid
    dist = skData12(1,1)%dist

    ! All SK files should have the same grid separation and table length
    nGrid = skData12(1,1)%nGrid
    dist = skData12(1,1)%dist
    do iSK1 = 1, nSK1
      do iSK2 = 1, nSK2
        pSK12 => skData12(iSK2, iSK1)
        pSK21 => skData21(iSK1, iSK2)

        if (pSK12%dist /= dist .or. pSK21%dist /= dist) then
          write (errorStr, "(A,I2,A,I2)") "Incompatible SK grid separations &
              &for species ", sp1, ", ", sp2
          call error(errorStr)
        end if
        if (pSK12%nGrid /= nGrid .or. pSK21%nGrid /= nGrid) then
          write (errorStr, "(A,I2,A,I2)") "Incompatible SK grid lengths for &
              &species pair ", sp1, ", ", sp2
          call error(errorStr)
        end if
      end do
    end do

  end subroutine checkSKCompElec


  !> Checks if the provided repulsive splines for A-B and B-A are compatible
  subroutine checkSKCompRepSpline(repIn1, repIn2, sp1, sp2)

    !> Repulsive spline for interaction A-B
    type(TSplineRepInp), intent(in) :: repIn1

    !> Repulsive spline for interaction B-A
    type(TSplineRepInp), intent(in) :: repIn2

    !> Number of species A (for error messages only)
    integer, intent(in) :: sp1

    !> Number of species B (for error messages only)
    integer, intent(in) :: sp2


    !> Tolerance for the agreement in the repulsive data
    real(dp), parameter :: tolRep = 1.0e-8_dp


    !> string for error return
    character(lc) :: errorStr

    ! Repulsives for A-B and B-A should be the same
    if (size(repIn1%xStart) /= size(repIn2%xStart)) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible nr. of repulsive &
          &intervals for species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if
    if (maxval(abs(repIn1%xStart - repIn2%xStart)) > tolRep) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive spline &
          &intervals for species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if
    if (maxval(abs(repIn1%spCoeffs - repIn2%spCoeffs)) > tolRep &
        &.or. maxval(abs(repIn1%spLastCoeffs - repIn2%spLastCoeffs)) &
        &> tolRep) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive spline &
          &coefficients for species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if
    if (maxval(abs(repIn1%expCoeffs - repIn2%expCoeffs)) > tolRep) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive spline &
          &exp. coefficients for species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if
    if (abs(repIn1%cutoff - repIn2%cutoff) > tolRep) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive spline &
          &cutoffs for species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if

  end subroutine checkSKCompRepSpline


  !> Checks if repulsive polynomials for A-B and B-A are compatible
  subroutine checkSKCompRepPoly(repIn1, repIn2, sp1, sp2)

    !> Repulsive polynomial for interaction A-B
    type(TPolyRepInp), intent(in) :: repIn1

    !> Repulsive polynomial for interaction B-A
    type(TPolyRepInp), intent(in) :: repIn2

    !> Number of species A (for error messages only)
    integer, intent(in) :: sp1

    !> Number of species B (for error messages only)
    integer, intent(in) :: sp2


    !> for error string return
    character(lc) :: errorStr

    if (any(repIn1%polyCoeffs /= repIn2%polyCoeffs)) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive polynomial &
          &coefficients  for the species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if
    if (repIn1%cutoff /= repIn2%cutoff) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive cutoffs  &
          &for the species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if

  end subroutine checkSKCompRepPoly


  !> Returns the nr. of Slater-Koster integrals necessary to describe the interactions between two
  !> species
  pure function getNSKIntegrals(sp1, sp2, orb) result(nInteract)

    !> Index of the first species
    integer, intent(in) :: sp1

    !> Index of the second species
    integer, intent(in) :: sp2

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Nr. of Slater-Koster interactions
    integer :: nInteract

    integer :: iSh1, iSh2

    nInteract = 0
    do iSh1 = 1, orb%nShell(sp1)
      do iSh2 = 1, orb%nShell(sp2)
        nInteract = nInteract + min(orb%angShell(iSh2, sp2), orb%angShell(iSh1, sp1)) + 1
      end do
    end do

  end function getNSKIntegrals


  !> Creates from the columns of the Slater-Koster files for A-B and B-A a full table for A-B,
  !> containing all integrals.
  subroutine getFullTable(skHam, skOver, skData12, skData21, angShells1, angShells2)

    !> Resulting table of H integrals
    real(dp), intent(out) :: skHam(:,:)

    !> Resulting table of S integrals
    real(dp), intent(out) :: skOver(:,:)

    !> Contains all SK files describing interactions for A-B
    type(TOldSKData), intent(in), target :: skData12(:,:)

    !> Contains all SK files describing interactions for B-A
    type(TOldSKData), intent(in), target :: skData21(:,:)

    !> Angular momenta to pick from the SK-files for species A
    type(TListIntR1), intent(inout) :: angShells1

    !> Angular momenta to pick from the SK-files for species B
    type(TListIntR1), intent(inout) :: angShells2

    integer :: ind, iSK1, iSK2, iSh1, iSh2, nSh1, nSh2, l1, l2, lMin, lMax, mm
    integer :: angShell1(maxL+1), angShell2(maxL+1)
    real(dp), pointer :: pHam(:,:), pOver(:,:)


    !> Maps (mm, l1, l2 ) onto an element in the SK table.
    !> l2 >= l1 (l1 = 0, 1, ...; l2 = 0, 1, ...), m <= l1.
    integer, parameter :: skMap(0:maxL, 0:maxL, 0:maxL) &
        &= reshape((/&
        &20, 0,  0,  0,  19,  0,  0,  0,  18,  0,  0,  0,  17,  0,  0,  0,&
        & 0, 0,  0,  0,  15, 16,  0,  0,  13, 14,  0,  0,  11, 12,  0,  0,&
        & 0, 0,  0,  0,   0,  0,  0,  0,   8,  9, 10,  0,   5,  6,  7,  0,&
        & 0, 0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   1,  2,  3,  4/),&
        &(/maxL + 1, maxL + 1, maxL + 1/))

    ind = 1
    do iSK1 = 1, len(angShells1)
      call intoArray(angShells1, angShell1, nSh1, iSK1)
      do iSh1 = 1, nSh1
        l1 = angShell1(iSh1)
        do iSK2 = 1, len(angShells2)
          call intoArray(angShells2, angShell2, nSh2, iSK2)
          do iSh2 = 1, nSh2
            l2 = angShell2(iSh2)
            if (l1 <= l2) then
              pHam => skData12(iSK2,iSK1)%skHam
              pOver => skData12(iSK2,iSK1)%skOver
              lMin = l1
              lMax = l2
            else
              pHam => skData21(iSK1,iSK2)%skHam
              pOver => skData21(iSK1,iSK2)%skOver
              lMin = l2
              lMax = l1
            end if
            do mm = 0, lMin
              ! Safety check, if array size are appropriate
              @:ASSERT(all(shape(skHam) >= (/ size(pHam, dim=1), ind /)))
              @:ASSERT(all(shape(skOver) >= (/ size(pOver, dim=1), ind /)))
              @:ASSERT(size(pHam, dim=1) == size(pOver, dim=1))
              skHam(:,ind) = pHam(:,skMap(mm,lMax,lMin))
              skOver(:,ind) = pOver(:,skMap(mm,lMax,lMin))
              ind = ind + 1
            end do
          end do
        end do
      end do
    end do

  end subroutine getFullTable


end module dftbp_dftbplus_parser_skfiles
