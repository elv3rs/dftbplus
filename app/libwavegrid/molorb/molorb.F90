!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to calculate the value of one or more molecular orbitals composed from STOs on
!! an equidistant grid.
module libwavegrid_molorb
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : imag
  use dftbp_dftb_boundarycond, only : TBoundaryConds
  use dftbp_dftb_periodic, only : getCellTranslations
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_typegeometry, only : TGeometry
  use libwavegrid_molorb_parallel, only : evaluateParallel
  use libwavegrid_molorb_types, only : TSpeciesBasis, TSystemParams, TPeriodicParams, TBasisParams
  use libwavegrid_slater, only : TSlaterOrbital

  implicit none

  private
  !> Data type containing information for molecular orbital calculator.
  type TMolecularOrbital
    !> System geometry and composition
    type(TSystemParams) :: system
    !> Periodic boundary conditions data
    type(TPeriodicParams) :: periodic
    !> Basis set data in SoA format
    type(TBasisParams) :: basis
    !> If it is initialised
    logical :: tInitialised = .false.
  contains
    private
    procedure, pass(this) :: initSystem
    procedure, pass(this) :: initBasis
    procedure, pass(this) :: initPeriodic
    procedure, pass(this) :: initCoords
  end type TMolecularOrbital


  !> Returns the value of one or more molecular orbitals on a grid
  interface getValue
    module procedure TMolecularOrbital_getValue_real
    module procedure TMolecularOrbital_getValue_cmpl
  end interface

  public :: TSpeciesBasis
  public :: TMolecularOrbital, TMolecularOrbital_init, getValue

contains


  !> Initialises MolecularOrbital instance.
  subroutine TMolecularOrbital_init(this, geometry, boundaryCond, basisInput)
    type(TMolecularOrbital), intent(out) :: this
    type(TGeometry), intent(in) :: geometry
    type(TBoundaryConds), intent(in) :: boundaryCond
    type(TSpeciesBasis), intent(in) :: basisInput(:)
    real(dp) :: maxCutoff

    @:ASSERT(.not. this%tInitialised)

    call this%initSystem(geometry, basisInput)
    call this%initBasis(basisInput)
    maxCutoff = sqrt(maxval(this%basis%cutoffsSq))
    call this%initPeriodic(geometry, maxCutoff)
    call this%initCoords(geometry, boundaryCond)

    this%tInitialised = .true.
  end subroutine TMolecularOrbital_init


  !> Initialises system parameters (non-periodic)
  subroutine initSystem(this, geometry, basis)
    class(TMolecularOrbital), intent(inout) :: this
    type(TGeometry), intent(in) :: geometry
    type(TSpeciesBasis), intent(in) :: basis(:)

    integer :: nOrbTotal, nOrbSpecies, iSpec, iAtom, ind, iOrb, angMom, nStosTotal

    @:ASSERT(geometry%nSpecies == size(basis))

    this%system%nAtom = geometry%nAtom
    this%system%nSpecies = geometry%nSpecies
    allocate(this%system%species(this%system%nAtom))
    this%system%species(:) = geometry%species

    ! Get total number of STOs
    nStosTotal = 0
    do iSpec = 1, this%system%nSpecies
      nStosTotal = nStosTotal + basis(iSpec)%nOrb
    end do
    allocate(this%system%iStos(this%system%nSpecies + 1))
    this%basis%nStos = nStosTotal

    ! Create STO index map: iStos(i) points to the first STO of species i
    ind = 1
    do iSpec = 1, this%system%nSpecies
      this%system%iStos(iSpec) = ind
      ind = ind + basis(iSpec)%nOrb
    end do
    this%system%iStos(iSpec) = ind

    ! Count total number of orbitals (including m-dependence)
    nOrbTotal = 0
    do iAtom = 1, this%system%nAtom
      iSpec = this%system%species(iAtom)
      do iOrb = 1, basis(iSpec)%nOrb
        angMom = basis(iSpec)%stos(iOrb)%angMom
        nOrbTotal = nOrbTotal + 1 + 2 * angMom
      end do
    end do
    this%system%nOrb = nOrbTotal
  end subroutine initSystem


  !> Converts basis from AoS to SoA format for performance.
  subroutine initBasis(this, basisInput)
    class(TMolecularOrbital), intent(inout) :: this
    type(TSpeciesBasis), intent(in) :: basisInput(:)
    integer :: iSpec, iOrb, ind
    type(TSlaterOrbital), allocatable :: stos_flat(:)

    ! Flatten the basis array for easier iteration
    allocate(stos_flat(this%basis%nStos))
    ind = 1
    do iSpec = 1, this%system%nSpecies
      stos_flat(ind:ind+basisInput(iSpec)%nOrb-1) = basisInput(iSpec)%stos(:)
      ind = ind + basisInput(iSpec)%nOrb
    end do

    ! Allocate SoA arrays
    allocate(this%basis%angMoms(this%basis%nStos))
    allocate(this%basis%sto_nPows(this%basis%nStos))
    allocate(this%basis%sto_nAlphas(this%basis%nStos))
    allocate(this%basis%cutoffsSq(this%basis%nStos))

    ! Populate SoA arrays
    do iOrb = 1, this%basis%nStos
      this%basis%angMoms(iOrb) = stos_flat(iOrb)%angMom
      this%basis%cutoffsSq(iOrb) = stos_flat(iOrb)%cutoff ** 2
      this%basis%sto_nPows(iOrb) = stos_flat(iOrb)%nPow
      this%basis%sto_nAlphas(iOrb) = stos_flat(iOrb)%nAlpha
    end do
    this%basis%maxNPows = maxval(this%basis%sto_nPows)
    this%basis%maxNAlphas = maxval(this%basis%sto_nAlphas)

    ! Allocate and populate coefficient/alpha matrices
    allocate(this%basis%sto_coeffs(this%basis%maxNPows, this%basis%maxNAlphas, this%basis%nStos))
    allocate(this%basis%sto_alphas(this%basis%maxNAlphas, this%basis%nStos))
    do iOrb = 1, this%basis%nStos
      this%basis%sto_coeffs(1:stos_flat(iOrb)%nPow, 1:stos_flat(iOrb)%nAlpha, iOrb) = stos_flat(iOrb)%aa
      this%basis%sto_alphas(1:stos_flat(iOrb)%nAlpha, iOrb) = stos_flat(iOrb)%alpha
    end do
  end subroutine initBasis


  !> Initializes periodic parameters
  subroutine initPeriodic(this, geometry, maxCutoff)
    class(TMolecularOrbital), intent(inout) :: this
    type(TGeometry), intent(in) :: geometry
    real(dp), intent(in) :: maxCutoff
    real(dp), allocatable :: rCellVec(:,:)
  
    this%periodic%isPeriodic = geometry%tPeriodic
    if (this%periodic%isPeriodic) then
      allocate(this%periodic%latVecs(3, 3))
      allocate(this%periodic%recVecs2pi(3, 3))

      this%periodic%latVecs(:,:) = geometry%latVecs
      call invert33(this%periodic%recVecs2pi, this%periodic%latVecs)
      this%periodic%recVecs2pi(:,:) = transpose(this%periodic%recVecs2pi)
      call getCellTranslations(this%periodic%fCellVec, this%periodic%rCellVec, this%periodic%latVecs, &
          & this%periodic%recVecs2pi, maxCutoff)
      this%periodic%nCell = size(this%periodic%fCellVec, dim=2)
    else
      this%periodic%nCell = 1
      allocate(this%periodic%latVecs(3, 0))
      allocate(this%periodic%recVecs2pi(3, 0))
      allocate(this%periodic%fCellVec(3, 1))
      allocate(this%periodic%rCellVec(3, 1))
      this%periodic%fCellVec(:,:) = 0.0_dp
      this%periodic%rCellVec(:,:) = 0.0_dp
    end if
  end subroutine initPeriodic


  !> Initializes coordinates including periodic images.
  ! Depends on initPeriodic having run first.
  subroutine initCoords(this, geometry, boundaryCond)
      class(TMolecularOrbital), intent(inout) :: this
      type(TGeometry), intent(in) :: geometry
      type(TBoundaryConds), intent(in) :: boundaryCond
      integer :: iCell, iAtom

      allocate(this%system%coords(3, this%system%nAtom, this%periodic%nCell))
      this%system%coords(:,:,1) = geometry%coords
      call boundaryCond%foldCoordsToCell(this%system%coords(:,:,1), this%periodic%latVecs)

      if (this%periodic%isPeriodic) then
        do iCell = 2, this%periodic%nCell
          do iAtom = 1, this%system%nAtom
            this%system%coords(:, iAtom, iCell) = this%system%coords(:, iAtom, 1) + this%periodic%rCellVec(:, iCell)
          end do
        end do
      end if
  end subroutine initCoords



  !> Returns molecular orbitals on a grid.
  subroutine TMolecularOrbital_getValue_real(this, origin, gridVecs, eigVecsReal, &
      & valueOnGrid, addDensities, preferCPU, occupationVec)

    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this
    !> Origin of the grid
    real(dp), intent(in) :: origin(:)
    !> Grid vectors
    real(dp), intent(in) :: gridVecs(:,:)
    !> Summation coefficients for the STOs
    real(dp), intent(in) :: eigVecsReal(:,:)
    !> Molecular orbitals on a grid
    real(dp), intent(out) :: valueOnGrid(:,:,:,:)
    !> Add densities instead of wave functions
    logical, intent(in), optional :: addDensities
    !> Whether to prefer CPU for calculation
    logical, intent(in), optional :: preferCPU
    !> if present, calculate total charge. Coefficients for each squared state
    real(dp), intent(in), optional :: occupationVec(:)

    integer :: kIndexes(0)
    complex(dp) :: valueCmpl(0, 0, 0, 0), eigVecsCmpl(0, 0), phases(0,0)
    logical :: isDensityCalc, doPreferCPU
    logical, parameter :: isRealInput = .true.

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(origin) == 3)
    @:ASSERT(all(shape(gridVecs) == [3, 3]))
    @:ASSERT(size(eigVecsReal, dim=1) == this%system%nOrb)
    @:ASSERT(all(shape(valueOnGrid) > [1, 1, 1, 0]))
    @:ASSERT(size(eigVecsReal, dim=2) == size(valueOnGrid, dim=4))

    isDensityCalc = .false.
    if (present(addDensities)) then
      isDensityCalc = addDensities
    end if
    doPreferCPU = .false.
    if (present(preferCPU)) then
      doPreferCPU = preferCPU
    end if

    if (present(occupationVec)) then
      call evaluateParallel(origin, gridVecs, this%system, this%periodic, kIndexes, phases, this%basis, &
        & isRealInput, isDensityCalc, doPreferCPU, eigVecsReal, eigVecsCmpl, &
        & valueOnGrid, valueCmpl, occupationVec)
    else
      call evaluateParallel(origin, gridVecs, this%system, this%periodic, kIndexes, phases, this%basis, &
        & isRealInput, isDensityCalc, doPreferCPU, eigVecsReal, eigVecsCmpl, &
        & valueOnGrid, valueCmpl)
    end if

  end subroutine TMolecularOrbital_getValue_real


  !> Returns molecular orbitals on a grid.
  subroutine TMolecularOrbital_getValue_cmpl(this, origin, gridVecs, eigVecsCmpl, kPoints,&
      & kIndexes, valueOnGrid, preferCPU)

    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this
    !> Origin of the grid
    real(dp), intent(in) :: origin(:)
    !> Grid vectors
    real(dp), intent(in) :: gridVecs(:,:)
    !> Summation coefficients for the STOs
    complex(dp), intent(in) :: eigVecsCmpl(:,:)
    !> Array of k-points
    real(dp), intent(in) :: kPoints(:,:)
    !> Index of the k-points in kPoints for every mol.orbital
    integer, intent(in) :: kIndexes(:)
    !> Molecular orbitals on grid on exit.
    complex(dp), intent(out) :: valueOnGrid(:,:,:,:)
    !> Whether to prefer CPU for calculation
    logical, intent(in), optional :: preferCPU

    real(dp) :: valueReal(0,0,0,0), eigVecsReal(0,0)
    complex(dp), allocatable :: phases(:,:)
    logical, parameter :: isRealInput = .false.
    logical, parameter :: isDensityCalc = .false.
    logical :: doPreferCPU

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(origin) == 3)
    @:ASSERT(all(shape(gridVecs) == [3, 3]))
    @:ASSERT(size(eigVecsCmpl, dim=1) == this%system%nOrb)
    @:ASSERT(all(shape(valueOnGrid) > [0, 0, 0, 0]))
    @:ASSERT(size(eigVecsCmpl, dim=2) == size(valueOnGrid, dim=4))
    @:ASSERT(size(kPoints, dim=1) == 3)
    @:ASSERT(size(kPoints, dim=2) > 0)
    @:ASSERT(size(kIndexes) == size(eigVecsCmpl, dim=2))
    @:ASSERT(maxval(kIndexes) <= size(kPoints, dim=2))
    @:ASSERT(minval(kIndexes) > 0)
    doPreferCPU = .false.
    if (present(preferCPU)) then
      doPreferCPU = preferCPU
    end if

    allocate(phases(this%periodic%nCell, size(kPoints, dim =2)))
    if (this%periodic%isPeriodic) then
      phases(:,:) = exp(imag * matmul(transpose(this%periodic%fCellVec), kPoints))
    else
      phases(1,:) = (1.0_dp, 0.0_dp)
    end if

    call evaluateParallel(origin, gridVecs, this%system, this%periodic, kIndexes, phases, this%basis, &
      & isRealInput, isDensityCalc, preferCPU, eigVecsReal, eigVecsCmpl, &
      & valueReal, valueOnGrid)

  end subroutine TMolecularOrbital_getValue_cmpl

end module libwavegrid_molorb
