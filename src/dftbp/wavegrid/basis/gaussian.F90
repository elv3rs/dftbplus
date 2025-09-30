
!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set pure = "" if defined('WITH_ASSERT') else "pure"
!> Holds TSlaterOrbital, a concrete implementation of TOrbital.
module dftbp_wavegrid_basis_gaussian
  use dftbp_common_accuracy, only : dp
  use dftbp_wavegrid_basis_orbital, only : TOrbital
  implicit none

  private

  public :: TSlaterOrbital
  
  !> Concrete class for a contracted Gaussian-Type Orbital (GTO).
  !> R(r) = r^l * Sum_p [ coeff_p * exp(-alpha_p * r^2) ]
  type, extends(TOrbital) :: TGaussianOrbital
    integer :: nprim = 0
    !> Exponents of primitive Gaussians
    real(dp), allocatable :: alpha(:) 
    !> Contraction coefficients
    real(dp), allocatable :: coeff(:) 
  contains
    procedure :: getRadial => TGaussianOrbital_getRadial
    procedure :: init => TGaussianOrbital_init
  end type TGaussianOrbital



contains

  !> Initializes a TGaussianOrbital with its analytical parameters.
  subroutine TGaussianOrbital_init(this, angMom, coeff, alpha)
    class(TGaussianOrbital), intent(out) :: this
    integer, intent(in) :: angMom
    real(dp), intent(in) :: coeff(:)
    real(dp), intent(in) :: alpha(:)

    this%angMom = angMom
    this%nprim = size(coeff)
    ! ASSERT(this%nprim == size(alpha))

    allocate(this%coeff(this%nprim), source=coeff)
    allocate(this%alpha(this%nprim), source=alpha)
  end subroutine TGaussianOrbital_init


  !> Calculates the value of a contracted GTO analytically.
  ${pure}$ function TGaussianOrbital_getRadial(this, r) result(val)
    class(TGaussianOrbital), intent(in) :: this
    real(dp), intent(in) :: r
    real(dp) :: val

    integer :: iPrim
    real(dp) :: rToTheL, sumOfPrims
    real(dp), parameter :: zeroTol = 1.0e-12_dp

    ! Sum of primitive Gaussians: Sum_p[coeff_p * exp(-alpha_p * r^2)]
    sumOfPrims = 0.0_dp
    do iPrim = 1, this%nprim
      sumOfPrims = sumOfPrims + this%coeff(iPrim) * exp(-this%alpha(iPrim) * r**2)
    end do

    ! The full radial part is R(r) = r^l * Sum_p[...]
    ! Handle r=0 case carefully to avoid 0**0.
    if (r < epsilon(1.0_dp)) then
      rToTheL = 1.0_dp
      if (this%angMom > 0) rToTheL = 0.0_dp
    else
      val = r**this%angMom * sumOfPrims
    end if

  end function TGaussianOrbital_getRadial




end module dftbp_wavegrid_basis_gaussian
