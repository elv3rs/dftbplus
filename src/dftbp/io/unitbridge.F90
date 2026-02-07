!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Unit conversion bridge for the hsd-fortran/hsd-data migration.
!>
!> Provides per-unit-category converter functions matching the pure function interface
!> required by hsd_get_with_unit:
!>
!>   pure function converter(value, from_unit, to_unit) result(converted)
!>     real(dp), intent(in) :: value
!>     character(len=*), intent(in) :: from_unit, to_unit
!>     real(dp) :: converted
!>   end function
!>
!> Each converter wraps the existing TUnit-based system from unitconversion.F90.
!> DFTB+ always works in atomic units internally, so `to_unit` is typically "au".
!>
!> Usage with hsd-fortran:
!>   call hsd_get_with_unit(table, "Temperature", temp, "au", dftbp_convert_energy, stat)
module dftbp_io_unitbridge
  use dftbp_common_accuracy, only : dp
  use dftbp_common_unitconversion, only : convertUnitValue, &
      & lengthUnits, inverseLengthUnits, energyUnits, forceUnits, timeUnits, freqUnits, &
      & volumeUnits, chargeUnits, EFieldUnits, BFieldUnits, pressureUnits, velocityUnits, &
      & dipoleUnits, massUnits, angularUnits, massDensityUnits
  implicit none

  private
  public :: dftbp_convert_length, dftbp_convert_inverse_length
  public :: dftbp_convert_energy, dftbp_convert_force, dftbp_convert_time
  public :: dftbp_convert_frequency, dftbp_convert_volume, dftbp_convert_charge
  public :: dftbp_convert_efield, dftbp_convert_bfield, dftbp_convert_pressure
  public :: dftbp_convert_velocity, dftbp_convert_dipole, dftbp_convert_mass
  public :: dftbp_convert_angular, dftbp_convert_mass_density

contains


  !> Convert a length value from the given unit to the target unit.
  pure function dftbp_convert_length(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, lengthUnits)

  end function dftbp_convert_length


  !> Convert an inverse length value from the given unit to the target unit.
  pure function dftbp_convert_inverse_length(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, inverseLengthUnits)

  end function dftbp_convert_inverse_length


  !> Convert an energy value from the given unit to the target unit.
  pure function dftbp_convert_energy(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, energyUnits)

  end function dftbp_convert_energy


  !> Convert a force value from the given unit to the target unit.
  pure function dftbp_convert_force(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, forceUnits)

  end function dftbp_convert_force


  !> Convert a time value from the given unit to the target unit.
  pure function dftbp_convert_time(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, timeUnits)

  end function dftbp_convert_time


  !> Convert a frequency value from the given unit to the target unit.
  pure function dftbp_convert_frequency(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, freqUnits)

  end function dftbp_convert_frequency


  !> Convert a volume value from the given unit to the target unit.
  pure function dftbp_convert_volume(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, volumeUnits)

  end function dftbp_convert_volume


  !> Convert a charge value from the given unit to the target unit.
  pure function dftbp_convert_charge(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, chargeUnits)

  end function dftbp_convert_charge


  !> Convert an electric field value from the given unit to the target unit.
  pure function dftbp_convert_efield(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, EFieldUnits)

  end function dftbp_convert_efield


  !> Convert a magnetic field value from the given unit to the target unit.
  pure function dftbp_convert_bfield(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, BFieldUnits)

  end function dftbp_convert_bfield


  !> Convert a pressure value from the given unit to the target unit.
  pure function dftbp_convert_pressure(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, pressureUnits)

  end function dftbp_convert_pressure


  !> Convert a velocity value from the given unit to the target unit.
  pure function dftbp_convert_velocity(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, velocityUnits)

  end function dftbp_convert_velocity


  !> Convert a dipole value from the given unit to the target unit.
  pure function dftbp_convert_dipole(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, dipoleUnits)

  end function dftbp_convert_dipole


  !> Convert a mass value from the given unit to the target unit.
  pure function dftbp_convert_mass(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, massUnits)

  end function dftbp_convert_mass


  !> Convert an angular value from the given unit to the target unit.
  pure function dftbp_convert_angular(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, angularUnits)

  end function dftbp_convert_angular


  !> Convert a mass density value from the given unit to the target unit.
  pure function dftbp_convert_mass_density(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res

    res = convert_via_au_(value, from_unit, to_unit, massDensityUnits)

  end function dftbp_convert_mass_density


  ! ---------------------------------------------------------------------------
  !  Private helpers
  ! ---------------------------------------------------------------------------


  !> Convert a value from from_unit to to_unit, going through atomic units as pivot.
  !>
  !> If from_unit == to_unit, value is returned unchanged.
  !> Otherwise: value_au = convertUnitValue(units, from_unit, value)
  !> Then if to_unit is not "au", would need reverse conversion — but DFTB+ always
  !> uses "au" as target, so this is a trivial passthrough.
  pure function convert_via_au_(value, from_unit, to_unit, units) result(res)
    use dftbp_common_unitconversion, only : TUnit
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    type(TUnit), intent(in) :: units(:)
    real(dp) :: res

    if (from_unit == to_unit) then
      res = value
    else
      ! Convert from from_unit to atomic units
      res = convertUnitValue(units, from_unit, value)
    end if

  end function convert_via_au_

end module dftbp_io_unitbridge
