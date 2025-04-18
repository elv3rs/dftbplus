set(builddir ${CMAKE_CURRENT_BINARY_DIR})
set(srcdir ${CMAKE_CURRENT_SOURCE_DIR})
set(projectdir ${PROJECT_SOURCE_DIR})

execute_process(
  COMMAND ${projectdir}/utils/test/testlist_to_fypp
  INPUT_FILE ${srcdir}/tests
  OUTPUT_FILE ${builddir}/_modes_tests.fypp)

# Make sure, the line-marker option is not set in fypp
set(fypp_flags ${FYPP_CONFIG_FLAGS})
list(APPEND fypp_flags -I${projectdir}/src/dftbp/include)

execute_process(
  COMMAND ${FYPP} ${fypp_flags} -DMPI_PROCS=${TEST_MPI_PROCS} -DOMP_THREADS=${TEST_OMP_THREADS}
  INPUT_FILE ${builddir}/_modes_tests.fypp
  OUTPUT_FILE ${builddir}/_modes_tests)

file(STRINGS ${builddir}/_modes_tests modes_tests_raw)
foreach(line IN LISTS modes_tests_raw)
  string(STRIP "${line}" testname)
  if(NOT "${testname}" STREQUAL "")
    list(APPEND tests ${testname})
  endif()
endforeach()

foreach(test IN LISTS tests)
  add_test(
    NAME app/modes/${test}
    COMMAND ${srcdir}/bin/autotest2
        -r ${srcdir} -w ${builddir} -d ${srcdir}/bin/tagdiff
        -P "${MODES_RUNNER}" -p "$<TARGET_FILE:modes>"
        -s P,R,C,S ${test})
  set_tests_properties(
    app/modes/${test}
    PROPERTIES
    ENVIRONMENT "DFTBPLUS_PARAM_DIR=${PROJECT_SOURCE_DIR}/external")
endforeach()
