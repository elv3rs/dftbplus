!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Simple mixer for mixing charges
module dftbp_mixer_simplemixer
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
  use dftbp_mixer_mixer, only: TMixerReal, TMixerCmplx
  implicit none

#:set FLAVOURS = [('cmplx', 'complex', 'Cmplx'), ('real', 'real', 'Real')]


!  private
#:for NAME, TYPE, LABEL in FLAVOURS
  public :: TSimpleMixer${LABEL}$
  public :: TSimpleMixer${LABEL}$_init, TSimpleMixer${LABEL}$_reset, TSimpleMixer${LABEL}$_mix
#:endfor


#:for NAME, TYPE, LABEL in FLAVOURS
  !> Contains data for a simple mixer
  type, extends (TMixer${LABEL}$) :: TSimpleMixer${LABEL}$
    !private

    !> Mixing parameter
    real(dp) :: mixParam
    contains
      procedure :: SimpleInit => TSimpleMixer${LABEL}$_init
      procedure :: reset => TSimpleMixer${LABEL}$_reset
      procedure :: mix1D => TSimpleMixer${LABEL}$_mix
      procedure :: hasInverseJacobian => TSimpleMixer${LABEL}$_hasInverseJacobian
      procedure :: getInverseJacobian => TSimpleMixer${LABEL}$_getInverseJacobian
  end type TSimpleMixer${LABEL}$
#:endfor


contains

#:for NAME, TYPE, LABEL in FLAVOURS
  !> Creates a simple mixer.
  subroutine TSimpleMixer${LABEL}$_init(this, mixParam)

    !> Simple mixer instance on exit
    class(TSimpleMixer${LABEL}$), intent(out) :: this

    !> Mixing parameter
    real(dp), intent(in) :: mixParam

    this%mixParam = mixParam

  end subroutine TSimpleMixer${LABEL}$_init


  !> Resets the mixer.
  subroutine TSimpleMixer${LABEL}$_reset(this, nElem)

    !> Simple mixer instance
    class(TSimpleMixer${LABEL}$), intent(inout) :: this

    !> Length of the vectors to mix
    integer, intent(in) :: nElem

    @:ASSERT(nElem > 0)

    continue

  end subroutine TSimpleMixer${LABEL}$_reset


  !> Does the actual mixing.
  subroutine TSimpleMixer${LABEL}$_mix(this, qInpRes, qDiff)

    !> SimpleMixer instance
    class(TSimpleMixer${LABEL}$), intent(inout) :: this

    !> Input charge on entry, mixed charge on exit
    ${TYPE}$(dp), intent(inout) :: qInpRes(:)

    !> Charge difference
    ${TYPE}$(dp), intent(in) :: qDiff(:)

    @:ASSERT(size(qInpRes) == size(qDiff))

    qInpRes(:) = qInpRes + this%mixParam * qDiff

  end subroutine TSimpleMixer${LABEL}$_mix


  !> Currently, an inverse Jacobian is only provided by the Broyden mixer.
  logical function TSimpleMixer${LABEL}$_hasInverseJacobian(this) result(hasInverseJacobian)
        class(TSimpleMixer${LABEL}$), intent(in) :: this
        hasInverseJacobian = .false.
  end function TSimpleMixer${LABEL}$_hasInverseJacobian

        

  !> Throw an error if the inverse Jacobian is requested.
  subroutine TSimpleMixer${LABEL}$_getInverseJacobian(this, invJac)
      class(TSimpleMixer${LABEL}$), intent(in) :: this
      ${TYPE}$(dp), intent(out) :: invJac(:,:)
      !stop "SimpleMixer does not provide an inverse Jacobian. Perhaps you meant to use the Broyden mixer?"
      !call error("SimpleMixer does not provide an inverse Jacobian. Perhaps you meant to use the Broyden mixer?")
  end subroutine TSimpleMixer${LABEL}$_getInverseJacobian

#:endfor

end module dftbp_mixer_simplemixer
