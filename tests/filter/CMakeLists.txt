cado_define_test(SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/testsm.sh
    PROVIDE_TEMPORARY_WDIR
    TARGET_DEPENDENCIES sm)

cado_define_test(SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_filter_galois.sh
    PROVIDE_TEMPORARY_WDIR
    TARGET_DEPENDENCIES filter_galois)

cado_define_test(SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_dup1.sh
    PROVIDE_TEMPORARY_WDIR
    TARGET_DEPENDENCIES dup1)

set(sage_opt)
set(deps)
if(SAGE)
    set(sage_opt _with_sagemath)
    set(deps explain_indexed_relation)
endif()

cado_define_test(test_dup2${sage_opt}
   NO_DEFAULT_RUN
   SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_dup2.sh
   PROVIDE_TEMPORARY_WDIR
   TARGET_DEPENDENCIES freerel debug_renumber dup2 ${deps}
   )

if (SAGE)
    add_dependencies(all_sagemath_test_dependencies test_dup2${sage_opt}_dependencies)
endif()

if(DEFINED ENV{CHECKS_EXPENSIVE})
    # rels found by:
    # ./build/sieve/makefb   -poly ./tests/misc/test_renumber.data/c105_30003.poly   -lim 1600000 -side 1 -out /tmp/r1.gz
    # ./build/sieve/las  -I 14    -poly ./tests/misc/test_renumber.data/c105_30003.poly   -sqside 1 -q0 300000 -q1 300020 -lim0 1000000 -lim1 1600000 -lpb0 22 -lpb1 23 -mfb0 50 -mfb1 52 -fb1 /tmp/r1.gz -v
cado_divert_test(test_dup2${sage_opt} c105_30003
    APPEND_ARGUMENTS
    -poly ${CADO_NFS_SOURCE_DIR}/tests/misc/test_renumber.data/c105_30003.poly
    -rels ${CMAKE_CURRENT_SOURCE_DIR}/c105_30003.rels.gz
    -lpbs 22,23
    TIMEOUT 90
    ${AVOID_CONCURRENT_UNDER_COVERAGE}
    )
endif()

cado_divert_test(test_dup2${sage_opt} c105_30004_with_sides
    APPEND_ARGUMENTS
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/p40_30004.poly
    -rels ${CMAKE_CURRENT_SOURCE_DIR}/p40_30004_with_sides.rels.gz
    -lpbs 14,14
    -dl
    -S daa017e17f50607820e7896dcc3a05370364b3b4
    -R 2959973f1
    TIMEOUT 90
    )

# dup2 with -large-ab should give the same results as without -large-ab when
# rels have small a,b values
cado_divert_test(test_dup2${sage_opt} c105_30004_largeab
    APPEND_ARGUMENTS
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/p40_30004.poly
    -rels ${CMAKE_CURRENT_SOURCE_DIR}/p40_30004.rels.gz
    -lpbs 14,14
    -dl
    -S daa017e17f50607820e7896dcc3a05370364b3b4
    -R 2959973f1
    -large-ab
    TIMEOUT 90
    )

cado_define_test(test_sm_with_small_ell
    PROGRAM
    $<TARGET_FILE:sm_append>
    -ell 79
    -in   ${CMAKE_CURRENT_SOURCE_DIR}/c30_sm_small_ell.purged
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/c30_sm_small_ell.poly
    TARGET_DEPENDENCIES sm_append
)

cado_define_test(test_purge_torture
    PROGRAM
    $<TARGET_FILE:purge> -keep 0 -col-min-index 8
    -required_excess 0.0 ${CMAKE_CURRENT_SOURCE_DIR}/purge_torturetest1.rels
    TARGET_DEPENDENCIES purge
)

if (NOT DEFINED ENV{NO_PYTHON_CHECK})
cado_define_test(test_purge_heavy_torture
    PROVIDE_TEMPORARY_WDIR
    PROGRAM
    ${CMAKE_CURRENT_SOURCE_DIR}/heavy_torture_test_purge.sh
    TARGET_DEPENDENCIES purge
)
endif()

add_executable(create_fake_relation_file create_fake_relation_file.cpp)
target_link_libraries(create_fake_relation_file utils)

# under valgrind, this takes easily 200 seconds to run only 16 tests.
# Let's limit the number of trials in this case. Anyway it doesn't make
# much sense to test multiple times like this under valgrind.
set(nrepeats REPEATS 32)
if (VALGRIND)
    set(nrepeats REPEATS 2)
endif()
cado_define_test(dup1_stress_30138
    PROVIDE_TEMPORARY_WDIR_ARG -out
    ${nrepeats}
    PROGRAM
    $<TARGET_FILE:dup1> --mkdir -prefix dup1.0 -n 3
    -filelist ${CMAKE_CURRENT_SOURCE_DIR}/dup1_stress_30138.data/dup1_stress_30138.files
    -basepath ${CMAKE_CURRENT_SOURCE_DIR}/dup1_stress_30138.data
)
