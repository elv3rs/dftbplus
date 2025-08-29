add_subdirectory(dftb+)
add_subdirectory(modes)
add_subdirectory(libwavegrid)
if(WITH_TRANSPORT)
  add_subdirectory(phonons)
endif()

