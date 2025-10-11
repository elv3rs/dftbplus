
!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set pure = "" if defined('WITH_ASSERT') else "pure"
!> Holds TGaussianOrbital, a concrete implementation of TOrbital.
module dftbp_wavegrid_basis_gaussian
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
  use dftbp_wavegrid_basis_orbital, only : TOrbital
  implicit none

  private

  public :: TGaussianOrbital
  
  !> Represents a contracted Gaussian-Type Orbital (GTO):
  !> R(r) = r^l * Sum_p [ coeff_p * exp(-alpha_p * r^2) ]
  type, extends(TOrbital) :: TGaussianOrbital
    !> Summation coefficients [nPrimitives]
    real(dp), allocatable :: coeff(:) 
    !> Exponents of primitive Gaussians [nPrimitives]
    real(dp), allocatable :: alpha(:) 
  contains
    procedure :: getRadial => TGaussianOrbital_getRadial
    procedure :: init => TGaussianOrbital_init
    procedure, pass(lhs) :: assign => TGaussianOrbital_assign
  end type TGaussianOrbital



contains

  !> Store the parameters for the TGaussianOrbital.
  subroutine TGaussianOrbital_init(this, coeff, alpha, angMom, cutoff)
    class(TGaussianOrbital), intent(out) :: this
    !> Summation coefficients [nPrimitives] 
    real(dp), intent(in) :: coeff(:)

    !> Exponential coefficients [nPrimitives]
    real(dp), intent(in) :: alpha(:)

    !> Angular momentum of the orbital (l)
    integer, intent(in) :: angMom

    !> Cutoff, after which orbital is assumed to be zero
    real(dp), intent(in) :: cutoff

    @:ASSERT(size(coeff) == size(alpha))

    this%angMom = angMom
    this%cutoffSq = cutoff**2

    allocate(this%coeff, source=coeff)
    allocate(this%alpha, source=alpha)
  end subroutine TGaussianOrbital_init


  !> Calculates the value of a contracted GTO analytically.
  ${pure}$ function TGaussianOrbital_getRadial(this, r) result(val)
    class(TGaussianOrbital), intent(in) :: this
    real(dp), intent(in) :: r
    real(dp) :: val

    integer :: iPrim
    real(dp) :: weighedSum

    ! Sum and weigh all primitive Gaussians
    weighedSum = 0.0_dp
    do iPrim = 1, size(this%coeff)
      weighedSum = weighedSum + this%coeff(iPrim) * exp(-this%alpha(iPrim) * r**2)
    end do

    ! Multiply with r^l (Avoid 0**0 on pre 2008 compilers)
    if (r < epsilon(1.0_dp)) then
      val = weighedSum
      if (this%angMom > 0) val = 0.0_dp
    else
      val = r**this%angMom * weighedSum
    end if

  end function TGaussianOrbital_getRadial




  !> Assignment operator
  subroutine TGaussianOrbital_assign(lhs, rhs)
    class(TGaussianOrbital), intent(out) :: lhs
    class(TOrbital), intent(in) :: rhs
    select type (rhs)
      type is (TGaussianOrbital)
        lhs%angMom = rhs%angMom
        lhs%cutoffSq = rhs%cutoffSq

        if (allocated(lhs%alpha)) deallocate(lhs%alpha)
        if (allocated(lhs%coeff)) deallocate(lhs%coeff)

        allocate(lhs%alpha, source=rhs%alpha)
        allocate(lhs%coeff, source=rhs%coeff)
      class default
        call error("Cannot assign non-Gaussian orbital to Gaussian orbital.")
    end select

  end subroutine TGaussianOrbital_assign
end module dftbp_wavegrid_basis_gaussian
