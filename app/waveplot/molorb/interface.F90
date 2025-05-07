!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include 'common.fypp'

!> Module to interface with the C++ implementation of pointwise MO evaluation (REAL ONLY).
module waveplot_molorb_cpp_interface
  use iso_c_binding
  use dftbp_common_accuracy, only : dp
  use waveplot_slater, only : TSlaterOrbital
  implicit none

  private

  public :: evaluatePointwise_cpp_real_wrapper


  ! Interface to the C++ function (REAL ONLY version)
  interface
    subroutine evaluatePointwise_cpp_real( &
        origin, gridVecs, nGridX, nGridY, nGridZ, & ! Grid definition
        eigVecsReal, nEig, &                         ! Eigenvectors (Real)
        nAtom, nOrb, coords, species, iStos, nSpecies, & ! System Geometry & Basis
        sto_angMom, sto_cutoff, sto_nPow, sto_nAlpha, & ! STO Data (flat)
        sto_aa_flat, sto_alpha_flat, sto_aa_offsets, sto_alpha_offsets, & ! STO Data (flat)
        nTotalStos, &                                 ! STO Data (flat)
        tPeriodic, latVecs, recVecs2p, nCell, &       ! Periodicity
        tAddDensities, &                              ! Calculation Type
        valueReal &                                   ! Output Arrays (Real)
      ) bind(C, name='evaluatePointwise_cpp_real')

      import :: c_double, c_int, c_bool

      real(c_double), intent(in) :: origin(3)
      real(c_double), intent(in) :: gridVecs(3, 3)
      integer(c_int), value, intent(in) :: nGridX, nGridY, nGridZ
      real(c_double), intent(in) :: eigVecsReal(nOrb, nEig)
      integer(c_int), value, intent(in) :: nEig
      integer(c_int), value, intent(in) :: nAtom, nOrb
      real(c_double), intent(in) :: coords(3, nAtom, nCell)
      integer(c_int), intent(in) :: species(nAtom)
      integer(c_int), intent(in) :: iStos(*)
      integer(c_int), value, intent(in) :: nSpecies
      integer(c_int), intent(in) :: sto_angMom(*)
      real(c_double), intent(in) :: sto_cutoff(*)
      integer(c_int), intent(in) :: sto_nPow(*)
      integer(c_int), intent(in) :: sto_nAlpha(*)
      real(c_double), intent(in) :: sto_aa_flat(*)
      real(c_double), intent(in) :: sto_alpha_flat(*)
      integer(c_int), intent(in) :: sto_aa_offsets(*)
      integer(c_int), intent(in) :: sto_alpha_offsets(*)
      integer(c_int), value, intent(in) :: nTotalStos
      logical(c_bool), value, intent(in) :: tPeriodic
      real(c_double), intent(in) :: latVecs(3, 3)
      real(c_double), intent(in) :: recVecs2p(3, 3)
      integer(c_int), value, intent(in) :: nCell
      logical(c_bool), value, intent(in) :: tAddDensities
      real(c_double), intent(out) :: valueReal(nGridX, nGridY, nGridZ, nEig)

    end subroutine evaluatePointwise_cpp_real
  end interface


contains


  !> Fortran wrapper to prepare data and call the C++ evaluation routine (REAL ONLY).
  !! Arguments mostly match the original Fortran evaluatePointwise, but complex args removed.
  subroutine evaluatePointwise_cpp_real_wrapper(origin, gridVecs, eigVecsReal, & ! REMOVED eigVecsCmpl
      & nAtom, nOrb, coords, species, iStos, stos, tPeriodic, &
      & latVecs, recVecs2p, &
      & nCell,  tAddDensities, valueReal & 
     )

    real(dp), intent(in) :: origin(3)
    real(dp), intent(in) :: gridVecs(3, 3)
    real(dp), intent(in), target :: eigVecsReal(:,:)
    integer, intent(in) :: nAtom
    integer, intent(in) :: nOrb
    real(dp), intent(in), target :: coords(:,:,:)
    integer, intent(in), target :: species(:)
    integer, intent(in), target :: iStos(:)
    type(TSlaterOrbital), intent(in), target :: stos(:)
    logical, intent(in) :: tPeriodic
    real(dp), intent(in), target :: latVecs(:,:)
    real(dp), intent(in), target :: recVecs2p(:,:)
    integer, intent(in) :: nCell
    logical, intent(in) :: tAddDensities
    real(dp), intent(out) :: valueReal(:,:,:,:)

    ! Local variables for flattened STO data
    integer :: nTotalStos, nSpecies
    integer, allocatable :: sto_angMom(:)
    real(dp), allocatable :: sto_cutoff(:)
    integer, allocatable :: sto_nPow(:)
    integer, allocatable :: sto_nAlpha(:)
    real(dp), allocatable :: sto_aa_flat(:)
    real(dp), allocatable :: sto_alpha_flat(:)
    integer, allocatable :: sto_aa_offsets(:)
    integer, allocatable :: sto_alpha_offsets(:)

    integer :: i, nEig
    integer :: nGridX, nGridY, nGridZ
    integer :: total_aa_size, total_alpha_size
    integer :: current_aa_offset, current_alpha_offset
    integer :: istat


    ! --- Determine sizes ---
    nTotalStos = size(stos)
    nSpecies = size(iStos) - 1
    nGridX = size(valueReal, 1)
    nGridY = size(valueReal, 2)
    nGridZ = size(valueReal, 3)
    nEig = size(valueReal, 4)

    ! --- Prepare flattened STO data (Same as before) ---
    allocate(sto_angMom(nTotalStos), sto_cutoff(nTotalStos), sto_nPow(nTotalStos), &
             sto_nAlpha(nTotalStos), sto_aa_offsets(nTotalStos), sto_alpha_offsets(nTotalStos))
    total_aa_size = 0
    total_alpha_size = 0
    do i = 1, nTotalStos
      if (.not. allocated(stos(i)%aa) .or. .not. allocated(stos(i)%alpha)) then
         write(*,*) 'ERROR: STO components aa or alpha not allocated for STO index: ', i
         stop 1
      end if
      sto_angMom(i) = stos(i)%angMom
      sto_cutoff(i) = stos(i)%cutoff
      sto_nPow(i) = stos(i)%nPow
      sto_nAlpha(i) = stos(i)%nAlpha
      total_aa_size = total_aa_size + int(stos(i)%nPow) * int(stos(i)%nAlpha)
      total_alpha_size = total_alpha_size + stos(i)%nAlpha
    end do
    allocate(sto_aa_flat(total_aa_size), sto_alpha_flat(total_alpha_size))
    current_aa_offset = 0
    current_alpha_offset = 0
    do i = 1, nTotalStos
      sto_aa_offsets(i) = current_aa_offset     ! 0-based offset
      sto_alpha_offsets(i) = current_alpha_offset ! 0-based offset
      if (stos(i)%nPow > 0 .and. stos(i)%nAlpha > 0) then
         sto_aa_flat(current_aa_offset + 1 : current_aa_offset + stos(i)%nPow * stos(i)%nAlpha) = &
              reshape(stos(i)%aa, (/ stos(i)%nPow * stos(i)%nAlpha /))
         current_aa_offset = current_aa_offset + stos(i)%nPow * stos(i)%nAlpha
      end if
      if (stos(i)%nAlpha > 0) then
        sto_alpha_flat(current_alpha_offset + 1 : current_alpha_offset + stos(i)%nAlpha) = stos(i)%alpha
        current_alpha_offset = current_alpha_offset + stos(i)%nAlpha
      end if
    end do

    ! --- Call the C++ function (Real only version) ---
    call evaluatePointwise_cpp_real( &
        origin, gridVecs, nGridX, nGridY, nGridZ, &
        eigVecsReal, nEig, &
        nAtom, nOrb, coords, species, iStos, nSpecies, &
        sto_angMom, sto_cutoff, sto_nPow, sto_nAlpha, &
        sto_aa_flat, sto_alpha_flat, sto_aa_offsets, sto_alpha_offsets, &
        nTotalStos, &
        logical(tPeriodic, kind=c_bool), latVecs, recVecs2p, nCell,  &
        logical(tAddDensities, kind=c_bool), &
        valueReal &
        )

    ! --- Deallocate temporary arrays ---
    deallocate(sto_angMom, sto_cutoff, sto_nPow, sto_nAlpha, sto_aa_flat, sto_alpha_flat, &
               sto_aa_offsets, sto_alpha_offsets, stat=istat)
    if (istat /= 0) then
       write(*,*) 'ERROR: Failed to deallocate temporary STO arrays in wrapper.'
    end if

  end subroutine evaluatePointwise_cpp_real_wrapper

end module waveplot_molorb_cpp_interface
