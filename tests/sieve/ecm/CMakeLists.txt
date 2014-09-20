# Build testbench before running check
add_dependencies(check testbench)

# Get path to the testbench binary and store in variable testbench_path
# see comment in tests/sieve/CMakeLists.txt
# get_property(testbench_path TARGET testbench PROPERTY LOCATION)
set(testbench_path ${CADO_NFS_BINARY_DIR}/sieve/ecm/testbench)

add_test(test_factor_pm1 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/test_testbench.sh\" \"${testbench_path}\" \"${CMAKE_CURRENT_SOURCE_DIR}/test_factor_pm1_2_100_1000.inp\" -q -vf -vcf -pm1 100 1000")
add_test(test_factor_pp1_27 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/test_testbench.sh\" \"${testbench_path}\" \"${CMAKE_CURRENT_SOURCE_DIR}/test_factor_pp1_27_100_1000.inp\" -q -vf -vcf -pp1_27 100 1000")
add_test(test_factor_pp1_65 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/test_testbench.sh\" \"${testbench_path}\" \"${CMAKE_CURRENT_SOURCE_DIR}/test_factor_pp1_65_100_1000.inp\" -q -vf -vcf -pp1_65 100 1000")
add_test(test_factor_ecm_10 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/test_testbench.sh\" \"${testbench_path}\" \"${CMAKE_CURRENT_SOURCE_DIR}/test_factor_ecm_10_100_1000.inp\" -q -vf -vcf -ecm 100 1000 10")
add_test(test_factor_ecm_11 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/test_testbench.sh\" \"${testbench_path}\" \"${CMAKE_CURRENT_SOURCE_DIR}/test_factor_ecm_11_100_1000.inp\" -q -vf -vcf -ecm 100 1000 11")
add_test(test_factor_ecmm12_2 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/test_testbench.sh\" \"${testbench_path}\" \"${CMAKE_CURRENT_SOURCE_DIR}/test_factor_ecmm12_2_100_1000.inp\" -q -vf -vcf -ecmm12 100 1000 2")
add_test(test_factor_ecmm12_3 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/test_testbench.sh\" \"${testbench_path}\" \"${CMAKE_CURRENT_SOURCE_DIR}/test_factor_ecmm12_4_100_1000.inp\" -q -vf -vcf -ecmm12 100 1000 4")
add_test(test_factor_ecmm16_1 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/test_testbench.sh\" \"${testbench_path}\" \"${CMAKE_CURRENT_SOURCE_DIR}/test_factor_ecmm16_1_100_1000.inp\" -q -vf -vcf -ecmm16 100 1000 1")

macro(add_order_test METHOD SIGMA PMIN PMAX)
  add_test("test_order_ecm${METHOD}_${SIGMA}_${PMIN}_${PMAX}" bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/test_testbench.sh\" \"${testbench_path}\" \"${CMAKE_CURRENT_SOURCE_DIR}/test_order_ecm${METHOD}_${SIGMA}_${PMIN}_${PMAX}.inp\" -q -vf \"-po${METHOD}\" ${SIGMA}")
endmacro(add_order_test)


add_order_test("" 10 1000 2000)
add_order_test("" 10 1000000 1001000)
# Only one test with large primes as the sieve of Eratosthenes take a while
# to initalised for large starting prime
add_order_test("" 10 1000000000 1000001000)
add_order_test("" 11 1000 2000)
add_order_test("" 11 1000000 1001000)
add_order_test("" 12 1000 2000)
add_order_test("" 12 1000000 1001000)
add_order_test(m12 4 1000 2000)
add_order_test(m12 5 1000 2000)
add_order_test(m12 4 1000000 1001000)
add_order_test(m12 5 1000000 1001000)
add_order_test(m12 6 1000 2000)
add_order_test(m12 6 1000000 1001000)

# Point counting with torsion 16 curves is currently broken
# add_order_test(m16 1 1000 2000)
# add_order_test(m16 1 1000000 1001000)
