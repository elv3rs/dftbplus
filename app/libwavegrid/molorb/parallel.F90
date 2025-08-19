!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------! 

#:include 'common.fypp'

module libwavegrid_molorb_parallel
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : imag
  use dftbp_io_message, only : error
  use libwavegrid_slater, only : TSlaterOrbital, realTessY, getVerboseRadial
  use omp_lib, only : omp_is_initial_device, omp_get_num_devices
  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_double_complex, c_bool
  implicit none  
  
  public :: evaluateParallel

contains

  !> Returns the values of several molecular orbitals on grids.
  !! Caveat: The flag tPeriodic decides if the complex or the real version is read/written for the
  !! various parameters.
  subroutine evaluateParallel(origin, gridVecs, eigVecsReal, eigVecsCmpl, nAtom, nOrb, coords,&
      & species, iStos, stos, tPeriodic, tReal, latVecs, recVecs2p, kPoints,&
      & kIndexes, nCell, cellVec, tAddDensities, valueReal, valueCmpl)

    !> Origin of the grid
    real(dp), intent(in) :: origin(:)

    !> Grid vectors
    real(dp), intent(in) :: gridVecs(:,:)

    !> Real eigenvectors, or null-array
    real(dp), intent(in) :: eigVecsReal(:,:)

    !> Complex eigenvectors, or null-array
    complex(dp), intent(in) :: eigVecsCmpl(:,:)

    !> Nr. of atoms
    integer, intent(in) :: nAtom

    !> Nr. of orbitals
    integer, intent(in) :: nOrb

    !> Coordinates of the atoms
    real(dp), intent(in) :: coords(:,:,:)

    !> Species for each atom
    integer, intent(in) :: species(:)

    !> Starting position of the STOs for each species
    integer, intent(in) :: iStos(:)

    !> Array containing the STOs
    type(TSlaterOrbital), intent(in) :: stos(:)

    !> If the system is periodic
    logical, intent(in) :: tPeriodic

    !> If the system is real
    logical, intent(in) :: tReal

    !> Lattice vectors or null-array
    real(dp), intent(in) :: latVecs(:,:)

    !> Reciprocal vectors divided by 2pi (periodic) or a null-array (molecular)
    real(dp), intent(in) :: recVecs2p(:,:)

    !> Kpoints or null-array
    real(dp), intent(in) :: kPoints(:,:)

    !> Index of the k-points for each orbital in KPoints
    integer, intent(in) :: kIndexes(:)

    !> Nr. of cells to consider
    integer, intent(in) :: nCell

    !> Translation vector of the considered cells
    real(dp), intent(in) :: cellVec(:,:)

    !> If densities should be added instead of wave funcs
    logical, intent(in) :: tAddDensities

    !> Contains the real grid on exit
    real(dp), intent(out) :: valueReal(:,:,:,:)

    !> Contains the complex grid on exit
    complex(dp), intent(out) :: valueCmpl(:,:,:,:)

    !!Complex
    !complex(dp), allocatable :: atomOrbValCmpl(:)
    complex(dp) :: phases(nCell, size(kPoints, dim=2))
    integer :: nPoints(4)



    integer :: nEigOut

    !! SOA for stos contents
    integer, allocatable :: angMoms(:)
    real(dp), allocatable :: cutoffsSq(:)

    integer, allocatable :: sto_nPows(:)
    integer, allocatable :: sto_nAlphas(:)
    real(dp), allocatable :: sto_coeffs(:,:,:)
    real(dp), allocatable :: sto_alphas(:,:)
    integer :: maxNPows, maxNAlphas, iOrb

    !! Create arrays out of sto object data to simplify OMP parallelization
    allocate(angMoms(size(stos)), source=0)
    allocate(sto_nPows(size(stos)), source=0)
    allocate(sto_nAlphas(size(stos)), source=0)
    allocate(cutoffsSq(size(stos)), source=0.0_dp)


    do iOrb = 1, size(stos)
      angMoms(iOrb) = stos(iOrb)%angMom
      cutoffsSq(iOrb) = stos(iOrb)%cutoff ** 2
      sto_nPows(iOrb) = stos(iOrb)%nPow
      sto_nAlphas(iOrb) = stos(iOrb)%nAlpha
    end do


    !! Determine max power/ coefficient array size
    maxNPows = maxval(sto_nPows)
    maxNAlphas = maxval(sto_nAlphas)

    allocate(sto_coeffs(maxNPows, maxNAlphas, size(stos)))
    allocate(sto_alphas(maxNAlphas, size(stos)))

    ! Copy data from the sto objects to the arrays
    do iOrb = 1, size(stos)
      sto_coeffs(1:sto_nPows(iOrb), 1:sto_nAlphas(iOrb), iOrb) = stos(iOrb)%aa
      sto_alphas(1:sto_nAlphas(iOrb), iOrb) = stos(iOrb)%alpha
    end do


    ! Array for the contribution of each orbital (and its periodic images)
    !allocate(atomOrbValReal(nOrb))
    nPoints = shape(valueReal)
    !valueReal(:,:,:,:) = 0.0_dp


    !print *, "Devices:", omp_get_num_devices()
    !print *, "maxNPows:", maxNPows, "maxNAlphas:", maxNAlphas
    !print *, "Max sto_nPows in data:", maxval(sto_nPows)
    !print *, "iStos", iStos
    !print *, "Stos:", size(stos), "nOrb:", nOrb
    !print *, "species", species
    !print *, "sto_alphas:", sto_alphas
    !print *, "sto_nalphas:", sto_nAlphas
    !print *, "sto_nPows:", sto_nPows

    !print *, "EV:", sum(eigVecsReal(1, :)), sum(eigVecsReal(:,1))

    ! Phase factors for the periodic image cell. Note: This will be conjugated in the scalar product
    ! below. This is fine as, in contrast to what was published, DFTB+ implicitly uses exp(-ikr) as
    ! a phase factor, as the unpack routines assemble the lower triangular matrix with exp(ikr) as
    ! factor.
    phases(:,:) = exp(imag * matmul(transpose(cellVec), kPoints))
    nEigOut = nPoints(4) 
    if (tAddDensities) then
      nEigOut = 1
    end if

    
    #: set VARIANT = 'CUDA' if WITH_CUDA else 'OMP'
    !#: set VARIANT = 'OMP'
    print *, "Running molorb using ${VARIANT}$ kernel."

    call evaluate${VARIANT}$(nPointsX=nPoints(1), nPointsY=nPoints(2), nPointsZ=nPoints(3), &
        & nEigIn=size(eigVecsReal, dim=2), nEigOut=nEigOut, nOrb=nOrb, nStos=size(stos), &
        & maxNPows=maxNPows, maxNAlphas=maxNAlphas, &
        & nAtom=nAtom, nCell=nCell, nSpecies=size(iStos), isRealInput=tReal, isPeriodic=tPeriodic, &
        & isDensityCalc=tAddDensities, origin=origin, gridVecs=gridVecs, eigVecsReal=eigVecsReal, eigVecsCmpl=eigVecsCmpl, &
        & coords=coords, species=species, iStos=iStos, &
        & latVecs=latVecs, recVecs2p=recVecs2p, kIndexes=kIndexes, phases=phases, &
        & sto_angMoms=angMoms, sto_nPows=sto_nPows, sto_nAlphas=sto_nAlphas, &
        & sto_cutoffsSq=cutoffsSq, sto_coeffs=sto_coeffs, sto_alphas=sto_alphas, &
        & valueReal=valueReal, valueCmpl=valueCmpl)

  end subroutine evaluateParallel



#:if WITH_CUDA
  subroutine evaluateCuda(nPointsX, nPointsY, nPointsZ, nEigIn, nEigOut, nOrb, nStos, maxNPows, &
      & maxNAlphas, nAtom, nCell, nSpecies, &
      & isRealInput, isPeriodic, isDensityCalc, origin, gridVecs, eigVecsReal, eigVecsCmpl, coords, species, iStos, &
      & latVecs, recVecs2p, kIndexes, phases, &
      & sto_angMoms, sto_nPows, sto_nAlphas, sto_cutoffsSq, sto_coeffs, sto_alphas, &
      & valueReal, valueCmpl)
    use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_loc, c_double, c_double_complex


    integer, intent(in) :: nPointsX ! number of grid points in x direction
    integer, intent(in) :: nPointsY ! number of grid points in y direction
    integer, intent(in) :: nPointsZ ! number of grid points in z direction
    integer, intent(in) :: nEigIn ! number of eigenvalues in input
    integer, intent(in) :: nEigOut ! number of eigenvalues in output
    integer, intent(in) :: nOrb ! total number of orbitals
    integer, intent(in) :: nStos ! num of unique stos
    integer, intent(in) :: maxNPows ! max number of powers in a sto
    integer, intent(in) :: maxNAlphas ! max number of alphas in a sto
    integer, intent(in) :: nAtom ! number of atoms
    integer, intent(in) :: nCell ! number of adjacent cells
    integer, intent(in) :: nSpecies ! number of different atom kinds


    !> System parameters
    logical, intent(in) :: isRealInput  ! if the system is real
    logical, intent(in) :: isPeriodic ! if the system is periodic
    logical, intent(in) :: isDensityCalc ! if the calculation is for density instead of wave functions
    real(dp), intent(in), target :: origin(3)
    real(dp), intent(in), target :: gridVecs(3, 3)
    real(dp), intent(in), target :: eigVecsReal(nOrb, nEigIn)
    complex(dp), intent(in), target :: eigVecsCmpl(nOrb, nEigIn)
    real(dp), intent(in), target :: coords(3, nAtom, nCell)
    integer, intent(in), target :: species(nAtom)
    integer, intent(in), target :: iStos(nSpecies + 1)

    !> Additional Periodic system parameters
    real(dp), intent(in), target :: latVecs(3, 3) ! lattice vectors
    real(dp), intent(in), target :: recVecs2p(3, 3) ! reciprocal vectors divided by 2pi, or null-array (molecular)
    integer, intent(in), target :: kIndexes(nEigIn) ! index of the k-points for each orbital in KPoints
    complex(dp), intent(in), target :: phases(nCell, nEigIn) ! phases for the periodic system

    ! STO data
    integer, intent(in), target :: sto_angMoms(nStos)
    integer, intent(in), target :: sto_nPows(nStos)
    integer, intent(in), target :: sto_nAlphas(nStos)
    real(dp), intent(in), target :: sto_cutoffsSq(nStos)
    real(dp), intent(in), target :: sto_coeffs(maxNPows, maxNAlphas, nStos)
    real(dp), intent(in), target :: sto_alphas(maxNAlphas, nStos)

    !> Contains the real grid on exit
    real(dp), intent(out), target :: valueReal(nPointsX, nPointsY, nPointsZ, nEigOut)
    !> Contains the complex grid on exit
    complex(dp), intent(out), target :: valueCmpl(nPointsX, nPointsY, nPointsZ, nEigOut)

    type, bind(c) :: TGridParams
      integer(c_int) :: nPointsX, nPointsY, nPointsZ
      type(c_ptr) :: origin, gridVecs
    end type

    type, bind(c) :: TSystemParams
      integer(c_int) :: nAtom, nCell, nSpecies, nOrb
      type(c_ptr) :: coords, species, iStos
    end type

    type, bind(c) :: TPeriodicParams
      integer(c_int) :: isPeriodic
      type(c_ptr) :: latVecs, recVecs2p, kIndexes, phases
    end type

    type, bind(c) :: TBasisParams
      integer(c_int) :: nStos, maxNPows, maxNAlphas
      type(c_ptr) :: sto_angMoms, sto_nPows, sto_nAlphas
      type(c_ptr) :: sto_cutoffsSq, sto_coeffs, sto_alphas
    end type

    type, bind(c) :: TCalculationParams
      integer(c_int) :: nEigIn, nEigOut, isRealInput, isDensityCalc, accDensity
      type(c_ptr) :: eigVecsReal, eigVecsCmpl
      type(c_ptr) :: valueReal_out, valueCmpl_out
    end type


    ! Interface to the C-function defined in kernel.cu
    interface

      ! Now, define the subroutine with the new signature
      subroutine evaluate_on_device_c(grid, system, periodic, basis, calc) bind(C, name='evaluate_on_device_c')
        import
        type(TGridParams), intent(in) :: grid
        type(TSystemParams), intent(in) :: system
        type(TPeriodicParams), intent(in) :: periodic
        type(TBasisParams), intent(in) :: basis
        type(TCalculationParams), intent(in) :: calc
      end subroutine evaluate_on_device_c

    end interface

    type(TGridParams) :: grid_p
    type(TSystemParams) :: system_p
    type(TPeriodicParams) :: periodic_p
    type(TBasisParams) :: sto_basis_p
    type(TCalculationParams) :: calc_p
    type(c_ptr) :: valueReal_out, valueCmpl_out

    ! Populate the structs
    ! Grid parameters
    grid_p%nPointsX = nPointsX
    grid_p%nPointsY = nPointsY
    grid_p%nPointsZ = nPointsZ
    grid_p%origin   = c_loc(origin)
    grid_p%gridVecs = c_loc(gridVecs)
    ! System parameters
    system_p%nAtom    = nAtom
    system_p%nCell    = nCell
    system_p%nSpecies = nSpecies
    system_p%nOrb     = nOrb
    system_p%coords   = c_loc(coords)
    system_p%species  = c_loc(species)
    system_p%iStos    = c_loc(iStos)
    ! Periodic system parameters
    periodic_p%isPeriodic = merge(1, 0, isPeriodic)
    periodic_p%latVecs   = c_loc(latVecs)
    periodic_p%recVecs2p = c_loc(recVecs2p)
    periodic_p%kIndexes  = c_loc(kIndexes)
    periodic_p%phases    = c_loc(phases)
    ! Basis parameters
    sto_basis_p%nStos       = nStos
    sto_basis_p%maxNPows    = maxNPows
    sto_basis_p%maxNAlphas  = maxNAlphas
    sto_basis_p%sto_angMoms = c_loc(sto_angMoms)
    sto_basis_p%sto_nPows   = c_loc(sto_nPows)
    sto_basis_p%sto_nAlphas = c_loc(sto_nAlphas)
    sto_basis_p%sto_cutoffsSq = c_loc(sto_cutoffsSq)
    sto_basis_p%sto_coeffs  = c_loc(sto_coeffs)
    sto_basis_p%sto_alphas  = c_loc(sto_alphas)
    ! Calculation parameters
    calc_p%nEigIn         = nEigIn
    calc_p%nEigOut        = nEigOut
    calc_p%isRealInput    = merge(1, 0, isRealInput)
    calc_p%isDensityCalc  = merge(1, 0, isDensityCalc)
    calc_p%accDensity     = merge(1, 0, isDensityCalc)
    calc_p%eigVecsReal    = c_loc(eigVecsReal)
    calc_p%eigVecsCmpl    = c_loc(eigVecsCmpl)
    calc_p%valueReal_out  = c_loc(valueReal)
    calc_p%valueCmpl_out  = c_loc(valueCmpl)

    call evaluate_on_device_c(grid_p, system_p, periodic_p, sto_basis_p, calc_p)

  end subroutine evaluateCuda
#:endif




  subroutine evaluateOMP(nPointsX, nPointsY, nPointsZ, nEigIn, nEigOut, nOrb, nStos, maxNPows, maxNAlphas, nAtom, nCell, nSpecies, &
      & isRealInput, isPeriodic, isDensityCalc, origin, gridVecs, eigVecsReal, eigVecsCmpl, coords, species, iStos, &
      & latVecs, recVecs2p, kIndexes, phases, &
      & sto_angMoms, sto_nPows, sto_nAlphas, sto_cutoffsSq, sto_coeffs, sto_alphas, &
      & valueReal, valueCmpl)

    integer, intent(in) :: nPointsX ! number of grid points in x direction
    integer, intent(in) :: nPointsY ! number of grid points in y direction
    integer, intent(in) :: nPointsZ ! number of grid points in z direction
    integer, intent(in) :: nEigIn ! number of eigenvalues in input
    integer, intent(in) :: nEigOut ! number of eigenvalues in output
    integer, intent(in) :: nOrb ! total number of orbitals
    integer, intent(in) :: nStos ! num of unique stos 
    integer, intent(in) :: maxNPows ! max number of powers in a sto
    integer, intent(in) :: maxNAlphas ! max number of alphas in a sto
    integer, intent(in) :: nAtom ! number of atoms
    integer, intent(in) :: nCell ! number of adjacent cells
    integer, intent(in) :: nSpecies ! number of different atom kinds
    

    !> System parameters
    logical, intent(in) :: isRealInput  ! if the system is real
    logical, intent(in) :: isPeriodic  ! if the system is periodic
    logical, intent(in) :: isDensityCalc ! if the calculation is for density. Implies isRealInput.
    real(dp), intent(in) :: origin(3)
    real(dp), intent(in) :: gridVecs(3, 3)
    real(dp), intent(in) :: eigVecsReal(nOrb, nEigIn)
    complex(dp), intent(in) :: eigVecsCmpl(nOrb, nEigIn)
    real(dp), intent(in) :: coords(3, nAtom, nCell)
    integer, intent(in) :: species(nAtom)
    integer, intent(in) :: iStos(nSpecies + 1)

    !> Additional Periodic system parameters
    real(dp), intent(in) :: latVecs(3, 3) ! lattice vectors
    real(dp), intent(in) :: recVecs2p(3, 3) ! reciprocal vectors divided by 2pi, or null-array (molecular)
    integer, intent(in) :: kIndexes(nEigIn) ! index of the k-points for each orbital in KPoints
    complex(dp), intent(in) :: phases(nCell, nEigIn) ! phases for the periodic system

    ! STO data
    integer, intent(in) :: sto_angMoms(nStos)
    integer, intent(in) :: sto_nPows(nStos)
    integer, intent(in)  :: sto_nAlphas(nStos)
    real(dp), intent(in) :: sto_cutoffsSq(nStos)
    real(dp), intent(in) :: sto_coeffs(maxNPows, maxNAlphas, nStos)
    real(dp), intent(in) :: sto_alphas(maxNAlphas, nStos)

    !> Contains the real grid on exit
    real(dp), intent(out) :: valueReal(nPointsX, nPointsY, nPointsZ, nEigOut)
    !> Contains the complex grid on exit
    complex(dp), intent(out) :: valueCmpl(nPointsX, nPointsY, nPointsZ, nEigOut)

    

    !! Thread private variables
    integer ::  ind, iSpecies
    real(dp) :: sto_tmp_pows(16), xyz(3), diff(3)
    real(dp) :: rSq, r, val, radialVal, sto_tmp_rexp, frac(3)

    !! Loop Variables
    integer :: i1, i2, i3, iEig, iAtom, iOrb, iM, iL, iCell, ii, jj
    
    !$omp parallel do collapse(3) &
    !$omp&    private(i1, i2, i3, iCell, iAtom, iOrb, iEig, iL, iM, ii, jj, xyz, diff, &
    !$omp&              r, val, radialVal, sto_tmp_pows, sto_tmp_rexp, ind, iSpecies, rSq) &
    !$omp&    shared(gridVecs, origin, species, nPointsX, nPointsY, nPointsZ, &
    !$omp&              sto_angMoms, sto_nPows, sto_cutoffsSq, sto_nAlphas, sto_coeffs, &
    !$omp&              sto_alphas, eigVecsReal, eigVecsCmpl, nEigIn, nEigOut, nOrb, nStos, &
    !$omp&              coords, iStos, nAtom, nCell, isPeriodic, recVecs2p, latVecs, &
    !$omp&              phases, isRealInput, isDensityCalc, valueReal, valueCmpl)
    lpI3: do i3 = 1, nPointsZ
      lpI2: do i2 = 1, nPointsY
        lpI1: do i1 = 1, nPointsX
            valueReal(i1, i2, i3, :) = 0.0_dp
            xyz(:) = origin(:) + real(i1 - 1, dp) * gridVecs(:, 1) &                                                               
                             & + real(i2 - 1, dp) * gridVecs(:, 2) &
                             & + real(i3 - 1, dp) * gridVecs(:, 3)

            ! Fold coordinates into unit cell
            if (isPeriodic) then
              frac(:) = matmul(xyz, recVecs2p)
              xyz(:) = matmul(latVecs, frac - real(floor(frac), dp))
            end if

            ! Get contribution from every atom in every cell for current point
            lpCell: do iCell = 1, nCell
              ind = 0
              lpAtom: do iAtom = 1, nAtom
                iSpecies = species(iAtom)
                diff(:) = xyz - coords(:, iAtom, iCell)
                rSq = dot_product(diff, diff)

                lpOrb: do iOrb = iStos(iSpecies), iStos(iSpecies + 1) - 1
                  iL = sto_angMoms(iOrb)
                  ! Calculate wave function only if atom is inside the cutoff
                  if (rSq > sto_cutoffsSq(iOrb)) then
                    ! Skip this orbital
                    ind = ind + 2*iL + 1
                    cycle lpOrb
                  end if
                  r = sqrt(rSq)

                call getVerboseRadial(iL, &
                                    & sto_nPows(iOrb), &
                                    & sto_nAlphas(iOrb), &
                                    & sto_coeffs(1:sto_nPows(iOrb), 1:sto_nAlphas(iOrb), iOrb), &
                                    & sto_alphas(1:sto_nAlphas(iOrb), iOrb), &
                                    & r, &
                                    & radialVal)

                  lpM : do iM = -iL, iL
                    ind = ind + 1
                    val =  radialVal * realTessY(iL, iM, diff, r)

                    if (isRealInput) then
                      do iEig = 1, nEigIn
                        valueReal(i1, i2, i3, iEig) = valueReal(i1, i2, i3, iEig) + val * eigVecsReal(ind, iEig)
                      end do
                    else ! Complex
                      do iEig = 1, nEigIn
                        valueCmpl(i1, i2, i3, iEig) = valueCmpl(i1, i2, i3, iEig) + val &
                            & * phases(iCell, kIndexes(iEig)) *  eigVecsCmpl(ind, iEig)
                      end do
                    end if
                  end do lpM
                end do lpOrb
              end do lpAtom
            end do lpCell
        end do lpI1
      end do lpI2
    end do lpI3
    !$omp end parallel do
    end subroutine evaluateOMP

end module libwavegrid_molorb_parallel
