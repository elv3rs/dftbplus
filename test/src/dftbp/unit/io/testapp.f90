program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_io_hsd_smoke, only : hsd_smoke_tests => tests
  use test_io_indexselection, only : indexselection_tests => tests
  use test_io_tokenreader, only : tokenreader_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      hsd_smoke_tests(),&
      indexselection_tests(),&
      tokenreader_tests()&
    ])&
  )

end program testapp
