program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_libwavegrid_simple, only : simple_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      simple_tests()&
    ])&
  )
  ! Cases to handle:
  ! Both cpu, gpu.
  ! real : default, totChrg, addDensities, lut
  ! cmplx: default, totChrg, lut


end program testapp
