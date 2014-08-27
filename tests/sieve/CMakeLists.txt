# CMake 2.8.4 and newer has a WORKING_DIRECTORY property, but we may want
# to stay compatible with older versions

# Build makefb and las before running check
add_dependencies(check makefb las las_tracek las_kbucket las_mbucket)

# Get paths to the makefb and las binaries and store in variables
# makefb_path and las_path, resp.
get_property(makefb_path TARGET makefb PROPERTY LOCATION)
get_property(las_path TARGET las PROPERTY LOCATION)
get_property(las_tracek_path TARGET las_tracek PROPERTY LOCATION)
get_property(las_kbucket_path TARGET las_kbucket PROPERTY LOCATION)
get_property(las_mbucket_path TARGET las_mbucket PROPERTY LOCATION)

set(F9_roots "${CMAKE_CURRENT_BINARY_DIR}/F9.roots")
add_test(F9_makefbtest "${makefb_path}" -poly "${CADO_NFS_SOURCE_DIR}/params/F9.poly" -alim 2300000 -maxbits 10 -out "${F9_roots}")

add_test(F9_sievetest bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh\" \"${F9_roots}\" \"${las_path}\" \"${CADO_NFS_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}/F9.checksums\"")
# Before running F9_sievetest, remove the file with the captured checksums
# We have to make this a "test" so that ctest runs it, even though it does
# not actually test anything
add_test(F9_sievetest_cleanup "${CMAKE_COMMAND}" -E remove -f "${CMAKE_CURRENT_BINARY_DIR}/F9.checksums")
SET_TESTS_PROPERTIES(F9_sievetest PROPERTIES DEPENDS "F9_sievetest_cleanup;F9_makefbtest")

# We make these depend on sievetest, for checking the checksum file.
# Since sievetest depends on makefbtest, we don't have to specify it
# again.
add_test(F9_sievetest_unsieve1 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh\" \"${F9_roots}\" \"${las_path}\" \"${CADO_NFS_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}/F9.checksums\" -unsievethresh 1")
SET_TESTS_PROPERTIES(F9_sievetest_unsieve1 PROPERTIES DEPENDS F9_sievetest)
add_test(F9_kbucket_test bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh\" \"${F9_roots}\" \"${las_kbucket_path}\" \"${CADO_NFS_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}/F9.checksums\"")
SET_TESTS_PROPERTIES(F9_kbucket_test PROPERTIES DEPENDS F9_sievetest)
add_test(F9_mbucket_test bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh\" \"${F9_roots}\" \"${las_mbucket_path}\" \"${CADO_NFS_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}/F9.checksums\"")
SET_TESTS_PROPERTIES(F9_mbucket_test PROPERTIES DEPENDS F9_sievetest)

# The TRACE_K test does only a single special-q, so we don't compare output
# relations or checksums
add_test(F9_tracektest bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.tracektest.sh\" \"${F9_roots}\" \"${las_tracek_path}\" \"${CADO_NFS_SOURCE_DIR}\" -traceab -8517,584707")
SET_TESTS_PROPERTIES(F9_tracektest PROPERTIES DEPENDS F9_makefbtest)

# Try a too-large lambda value
add_test(F9_rlambdatest bash -c "rlambda=5 \"${CMAKE_CURRENT_SOURCE_DIR}/F9.tracektest.sh\" \"${F9_roots}\" \"${las_path}\" \"${CADO_NFS_SOURCE_DIR}\" -regex \"rlambda>.* does not make sense\"")
SET_TESTS_PROPERTIES(F9_rlambdatest PROPERTIES DEPENDS F9_makefbtest)
add_test(F9_alambdatest bash -c "alambda=5 \"${CMAKE_CURRENT_SOURCE_DIR}/F9.tracektest.sh\" \"${F9_roots}\" \"${las_path}\" \"${CADO_NFS_SOURCE_DIR}\" -regex \"alambda>.* does not make sense\"")
SET_TESTS_PROPERTIES(F9_alambdatest PROPERTIES DEPENDS F9_makefbtest)

# Try with default TRIALDIV_MAXLEN, whatever it may be
add_executable(test_trialdiv test_trialdiv.c "${CADO_NFS_SOURCE_DIR}/sieve/trialdiv.c")
target_link_libraries(test_trialdiv utils tests)
run_test(test_trialdiv test_trialdiv)

macro (add_trialdiv_test MAXLEN)
  add_executable(test_trialdiv${MAXLEN} test_trialdiv.c)
  set_property(TARGET test_trialdiv${MAXLEN} PROPERTY COMPILE_DEFINITIONS TRIALDIV_MAXLEN=${MAXLEN})
  target_link_libraries(test_trialdiv${MAXLEN} trialdiv utils tests)
  run_test(test_trialdiv${MAXLEN} test_trialdiv${MAXLEN})
endmacro(add_trialdiv_test)

macro (add_init_norms_bucket_region)
  add_executable(test_init_norms_bucket_region test_init_norms_bucket_region.c)
  set_property(TARGET test_init_norms_bucket_region PROPERTY COMPILE_DEFINITIONS)
  target_link_libraries(test_init_norms_bucket_region las-norms utils tests)
  run_test(test_init_norms_bucket_region test_init_norms_bucket_region)
endmacro(add_init_norms_bucket_region)

# Try specific TRIALDIV_MAXLEN from 1 to 6
add_trialdiv_test(1)
add_trialdiv_test(2)
add_trialdiv_test(3)
add_trialdiv_test(4)
add_trialdiv_test(5)
add_trialdiv_test(6)

add_init_norms_bucket_region()

add_executable(test_fb_make_step test_fb_make_step.c "${CADO_NFS_SOURCE_DIR}/sieve/fb.c")
target_link_libraries(test_fb_make_step utils tests)

add_subdirectory (ecm)
