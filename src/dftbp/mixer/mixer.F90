!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set FLAVOURS = [('real', 'real', 'Real'),('cmplx', 'complex', 'Cmplx')]

!> Provides a general mixer which contains the desired actual mixer.
module dftbp_mixer_mixer
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
#:for NAME, TYPE, LABEL in FLAVOURS
!  use dftbp_mixer_andersonmixer, only : TAndersonMixer${LABEL}$, TAndersonMixer${LABEL}$_mix,&
!      & TAndersonMixer${LABEL}$_reset
!  use dftbp_mixer_broydenmixer, only : TBroydenMixer${LABEL}$, TBroydenMixer${LABEL}$_mix,&
!    & TBroydenMixer${LABEL}$_reset
! use dftbp_mixer_diismixer, only : TDiisMixer${LABEL}$, TDiisMixer${LABEL}$_mix,&
!      & TDiisMixer${LABEL}$_reset
!  use dftbp_mixer_simplemixer, only : TSimpleMixer${LABEL}$, TSimpleMixer${LABEL}$_mix,&
!      & TSimpleMixer${LABEL}$_reset
#:endfor
  implicit none

  private
#:for NAME, TYPE, LABEL in FLAVOURS
  public :: TMixer${LABEL}$
  !public :: TMixer${LABEL}$_init,
  public::  TMixer${LABEL}$_reset, TMixer${LABEL}$_mix
#:endfor
  public :: TMixerReal_hasInverseJacobian, TMixerReal_getInverseJacobian
  public :: mixerTypes


#:for NAME, TYPE, LABEL in FLAVOURS





  type, abstract :: TMixer${LABEL}$
    contains
        !> Mixer Initialisation. 
        !-- Since init has different signatures, each mixer has its unique init
        !-- Todo: create a input type to unify the signature
        !--procedure(Initialise), deferred :: init
        procedure(IReset${LABEL}$), deferred :: reset

        !> Actual mixing routine
        procedure(IMix1D${LABEL}$), deferred :: mix1D
        
        !> Inverse Jacobian (only for Broyden mixer)
        procedure(IhasInverseJacobian${LABEL}$), deferred :: hasInverseJacobian
        procedure(IgetInverseJacobian${LABEL}$), deferred :: getInverseJacobian

        !> Interface declared below, allowing to mix 1d, 3d and 6d arrays 
        !--procedure, non-overridable :: TMixer${LABEL}$_mix
  end type TMixer${LABEL}$


  abstract interface 
    subroutine IReset${LABEL}$(this, nElem)
      import :: TMixer${LABEL}$
      class(TMixer${LABEL}$), intent(inout) :: this
      integer, intent(in) :: nElem
    end subroutine IReset${LABEL}$

    logical function IhasInverseJacobian${LABEL}$(this)
      import :: TMixer${LABEL}$
      class(TMixer${LABEL}$), intent(in) :: this
    end function IhasInverseJacobian${LABEL}$

    subroutine IMix1D${LABEL}$(this, qInpRes, qDiff)
      import :: TMixer${LABEL}$, dp
      class(TMixer${LABEL}$), intent(inout) :: this
      ${TYPE}$(dp), intent(inout) :: qInpRes(:)
      ${TYPE}$(dp), intent(in) :: qDiff(:)
    end subroutine IMix1D${LABEL}$

    subroutine IgetInverseJacobian${LABEL}$(this, invJac)
      import :: TMixer${LABEL}$, dp
      class(TMixer${LABEL}$), intent(in) :: this
      ${TYPE}$(dp), intent(out) :: invJac(:,:)
    end subroutine IgetInverseJacobian${LABEL}$
  end interface



  !> Initialises specific mixer in use
  !interface TMixer${LABEL}$_init
   ! module procedure TMixer${LABEL}$_initSimple
   ! module procedure TMixer${LABEL}$_initAnderson
   ! module procedure TMixer${LABEL}$_initBroyden
  !  module procedure TMixer${LABEL}$_initDiis
  !end interface TMixer${LABEL}$_init


  !> Does the actual mixing
  interface TMixer${LABEL}$_mix
    module procedure TMixer${LABEL}$_mix1D
    module procedure TMixer${LABEL}$_mix3D
    module procedure TMixer${LABEL}$_mix6D
  end interface TMixer${LABEL}$_mix
#:endfor


  type :: TMixerTypesEnum
    integer :: simple = 1
    integer :: anderson = 2
    integer :: broyden = 3
    integer :: diis = 4
  end type TMixerTypesEnum

  !> Contains mixer types
  type(TMixerTypesEnum), parameter :: mixerTypes = TMixerTypesEnum()


contains

#:for NAME, TYPE, LABEL in FLAVOURS
  !> Resets the mixer.
  subroutine TMixer${LABEL}$_reset(this, nElem)

    !> Mixer instance
    class(TMixer${LABEL}$), intent(inout) :: this

    !> Size of the vectors to mix
    integer, intent(in) :: nElem

    call this%reset(nElem)

  end subroutine TMixer${LABEL}$_reset


  !> Mixes two vectors.
  subroutine TMixer${LABEL}$_mix1D(this, qInpRes, qDiff)

    !> Mixer instance
    class(TMixer${LABEL}$), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    ${TYPE}$(dp), intent(inout) :: qInpRes(:)

    !> Difference between input and output vectors (measure of lack of convergence)
    ${TYPE}$(dp), intent(in) :: qDiff(:)

    call this%mix1D(qInpRes, qDiff)

  end subroutine TMixer${LABEL}$_mix1D


  !> Mixes two 3D matrices.
  subroutine TMixer${LABEL}$_mix3D(this, qInpResSqr, qDiffSqr)

    !> Mixer instance
    class(TMixer${LABEL}$), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    ${TYPE}$(dp), intent(inout), contiguous, target :: qInpResSqr(:,:,:)

    !> Difference between input and output vectors (measure of lack of convergence)
    ${TYPE}$(dp), intent(in), contiguous, target :: qDiffSqr(:,:,:)

    !! Difference between input and output vectors (1D pointer)
    ${TYPE}$(dp), pointer :: qDiff(:)

    !! Input vector on entry, result vector on exit (1D pointer)
    ${TYPE}$(dp), pointer :: qInpRes(:)

    qDiff(1:size(qDiffSqr)) => qDiffSqr
    qInpRes(1:size(qInpResSqr)) => qInpResSqr

    call TMixer${LABEL}$_mix1D(this, qInpRes, qDiff)

  end subroutine TMixer${LABEL}$_mix3D


  !> Mixes two 6D matrices.
  subroutine TMixer${LABEL}$_mix6D(this, qInpResSqr, qDiffSqr)

    !> Mixer instance
    class(TMixer${LABEL}$), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    ${TYPE}$(dp), intent(inout), contiguous, target :: qInpResSqr(:,:,:,:,:,:)

    !> Difference between input and output vectors (measure of lack of convergence)
    ${TYPE}$(dp), intent(in), contiguous, target :: qDiffSqr(:,:,:,:,:,:)

    !!Difference between input and output vectors (1D pointer)
    ${TYPE}$(dp), pointer :: qDiff(:)

    !! Input vector on entry, result vector on exit (1D pointer)
    ${TYPE}$(dp), pointer :: qInpRes(:)

    qDiff(1:size(qDiffSqr)) => qDiffSqr
    qInpRes(1:size(qInpResSqr)) => qInpResSqr

    call TMixer${LABEL}$_mix1D(this, qInpRes, qDiff)

  end subroutine TMixer${LABEL}$_mix6D


  !-- By default, no Jacobian is provided. As of right now, only the Real Broyden mixer provides an inverse Jacobian.

  !> Tells whether the mixer is able to provide the inverse Jacobian.
  logical function TMixer${LABEL}$_hasInverseJacobian(this) result(has)
    class(TMixer${LABEL}$), intent(in) :: this
    has = this%hasInverseJacobian()
  end function TMixer${LABEL}$_hasInverseJacobian

  !> Returns an inverse Jacobian if possible, halting if not.
  subroutine TMixer${LABEL}$_getInverseJacobian(this, invJac)
    class(TMixer${LABEL}$), intent(in) :: this
    ${TYPE}$(dp), intent(out) :: invJac(:,:)
    call this%getInverseJacobian(invJac)
  end subroutine TMixer${LABEL}$_getInverseJacobian
#:endfor

end module dftbp_mixer_mixer
