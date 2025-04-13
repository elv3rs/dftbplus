add_subdirectory(dftb+)
add_subdirectory(modes)
add_subdirectory(waveplot)
if(WITH_TRANSPORT)
  add_subdirectory(phonons)
endif()

