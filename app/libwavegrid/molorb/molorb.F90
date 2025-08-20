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
  subroutine TMolecularOrbital_init(this, geometry, boundaryCond, basis)

    !> Molecular Orbital
    type(TMolecularOrbital), intent(out) :: this

    !> Geometrical information.
    type(TGeometry), intent(in) :: geometry

    !> Boundary condition
    type(TBoundaryConds), intent(in) :: boundaryCond

    !> Basis for each species.
    type(TSpeciesBasis), intent(in) :: basis(:)

    integer :: nOrb, ii, jj, ind, iSp, iOrb
    real(dp) :: maxCutoff
    real(dp), allocatable :: rCellVec(:,:)
    type(TSlaterOrbital), allocatable :: stos_temp(:)

    @:ASSERT(.not. this%tInitialised)
    @:ASSERT(geometry%nSpecies == size(basis))

    ! Populate system parameters
    this%system%nAtom = geometry%nAtom
    this%system%nSpecies = geometry%nSpecies
    allocate(this%system%species(this%system%nAtom))
    this%system%species(:) = geometry%species

    ! Create an sto index 
    nOrb = 0
    do ii = 1, this%system%nSpecies
      nOrb = nOrb + basis(ii)%nOrb
    end do
    allocate(this%system%iStos(this%system%nSpecies+1))
    allocate(stos_temp(nOrb))
    this%basis%nStos = nOrb

    ind = 1
    do ii = 1, this%system%nSpecies
      this%system%iStos(ii) = ind
      nOrb = basis(ii)%nOrb
      stos_temp(ind:ind+nOrb-1) = basis(ii)%stos(1:nOrb)
      ind = ind + nOrb
    end do
    this%system%iStos(ii) = ind

    ! Count all orbitals (including atom, m-dependence)
    nOrb = 0
    do ii = 1, this%system%nAtom
      iSp = this%system%species(ii)
      do jj = this%system%iStos(iSp), this%system%iStos(iSp+1)-1
        nOrb = nOrb + 1 + 2 * stos_temp(jj)%angMom
      end do
    end do
    this%system%nOrb = nOrb

    ! Convert basis from AoS to SoA format
    allocate(this%basis%angMoms(this%basis%nStos), source=0)
    allocate(this%basis%sto_nPows(this%basis%nStos), source=0)
    allocate(this%basis%sto_nAlphas(this%basis%nStos), source=0)
    allocate(this%basis%cutoffsSq(this%basis%nStos), source=0.0_dp)
    do iOrb = 1, this%basis%nStos
      this%basis%angMoms(iOrb) = stos_temp(iOrb)%angMom
      this%basis%cutoffsSq(iOrb) = stos_temp(iOrb)%cutoff ** 2
      this%basis%sto_nPows(iOrb) = stos_temp(iOrb)%nPow
      this%basis%sto_nAlphas(iOrb) = stos_temp(iOrb)%nAlpha
    end do
    this%basis%maxNPows = maxval(this%basis%sto_nPows)
    this%basis%maxNAlphas = maxval(this%basis%sto_nAlphas)
    maxCutoff = sqrt(maxval(this%basis%cutoffsSq))
    allocate(this%basis%sto_coeffs(this%basis%maxNPows, this%basis%maxNAlphas, this%basis%nStos))
    allocate(this%basis%sto_alphas(this%basis%maxNAlphas, this%basis%nStos))
    do iOrb = 1, this%basis%nStos
      this%basis%sto_coeffs(1:this%basis%sto_nPows(iOrb), 1:this%basis%sto_nAlphas(iOrb), iOrb) = stos_temp(iOrb)%aa
      this%basis%sto_alphas(1:this%basis%sto_nAlphas(iOrb), iOrb) = stos_temp(iOrb)%alpha
    end do




    ! Populate periodic parameters
    this%periodic%isPeriodic = geometry%tPeriodic
    if (this%periodic%isPeriodic) then
      allocate(this%periodic%latVecs(3, 3))
      allocate(this%periodic%recVecs2p(3, 3))
      this%periodic%latVecs(:,:) = geometry%latVecs
      call invert33(this%periodic%recVecs2p, this%periodic%latVecs)
      this%periodic%recVecs2p(:,:) = reshape(this%periodic%recVecs2p, [3, 3], order=[2, 1])
      call getCellTranslations(this%periodic%cellVec, rCellVec, this%periodic%latVecs, &
          & this%periodic%recVecs2p, maxCutoff)
      this%periodic%nCell = size(this%periodic%cellVec,dim=2)
    else
      allocate(this%periodic%latVecs(3, 0))
      allocate(this%periodic%recVecs2p(3, 0))
      allocate(this%periodic%cellVec(3, 1))
      this%periodic%cellVec(:,:) = 0.0_dp
      allocate(rCellVec(3, 1))
      rCellVec(:,:) = 0.0_dp
      this%periodic%nCell = 1
    end if

    ! Create coordinates for central cell and periodic images
    allocate(this%system%coords(3, this%system%nAtom, this%periodic%nCell))
    this%system%coords(:,:,1) = geometry%coords
    call boundaryCond%foldCoordsToCell(this%system%coords(:,:,1), this%periodic%latVecs)
    if (this%periodic%isPeriodic) then
      do ii = 2, this%periodic%nCell
        do jj = 1, this%system%nAtom
          this%system%coords(:, jj, ii) = this%system%coords(:, jj, 1) + rCellVec(:, ii)
        end do
      end do
    end if



    this%tInitialised = .true.

  end subroutine TMolecularOrbital_init


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
      phases(:,:) = exp(imag * matmul(transpose(this%periodic%cellVec), kPoints))
    else
      phases(1,:) = (1.0_dp, 0.0_dp)
    end if

    call evaluateParallel(origin, gridVecs, this%system, this%periodic, kIndexes, phases, this%basis, &
      & isRealInput, isDensityCalc, preferCPU, eigVecsReal, eigVecsCmpl, &
      & valueReal, valueOnGrid)

  end subroutine TMolecularOrbital_getValue_cmpl

end module libwavegrid_molorb
