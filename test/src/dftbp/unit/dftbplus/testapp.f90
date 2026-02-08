program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_dftbplus_parser, only : parser_tests => tests
  use test_dftbplus_oldcompat, only : oldcompat_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      parser_tests(),&
      oldcompat_tests()&
    ])&
  )

end program testapp
