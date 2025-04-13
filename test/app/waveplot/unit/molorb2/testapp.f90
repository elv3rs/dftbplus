program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_waveplot_cache, only : cache_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      cache_tests()&
    ])&
  )

end program testapp
