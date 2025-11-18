include_directories(${CADO_NFS_SOURCE_DIR}/linalg)
include_directories(${CADO_NFS_SOURCE_DIR}/linalg/bwc)
include_directories(${CADO_NFS_SOURCE_DIR}/utils)
include_directories(${CADO_NFS_SOURCE_DIR}/linalg/bblas)

add_custom_target(all_gdb_printer_tests)

add_compile_options(-g -O0)
if(CMAKE_C_COMPILER_ID MATCHES "^Clang")
add_compile_options(-fno-limit-debug-info)
endif()

add_executable(test_gdb_printers test_gdb_printers.cpp)
target_link_libraries(test_gdb_printers utils)
add_dependencies(all_gdb_printer_tests test_gdb_printers)
if(NOT ENABLE_SHARED)
    configure_file(${CADO_NFS_SOURCE_DIR}/utils/hook_cado_nfs_gdb.py.in
        test_gdb_printers-gdb.py ESCAPE_QUOTES @ONLY)
endif()
cado_define_test(
    test_gdb_printers
    PROGRAM
    ${GDB_EXECUTABLE}
    -ex "add-auto-load-safe-path ${CADO_NFS_BINARY_DIR}"
    -ex "file ${CMAKE_CURRENT_BINARY_DIR}/test_gdb_printers"
    -x ${CMAKE_CURRENT_SOURCE_DIR}/test_gdb_printers.gdb
    --batch
    TARGET_DEPENDENCIES
    test_gdb_printers
    # 79cddcf25b16b8faaa4bc0e3ad26c2fb64c6d5d1 is without the poly
    # printers
    EXPECT_SHA1 79cddcf25b16b8faaa4bc0e3ad26c2fb64c6d5d1,bf57e3cb5e93d755beb3bf148114e55570422931
    SHA1_ON_REGEXP_LINES "^\\$"
    ENVIRONMENT
    CADO_NFS_GDB_SHOW_MATRIX_COMPACT=1
)

# primes are defined in tests/linalg/bwc/CMakeLists.txt

set(BWC_GFP_ARITHMETIC_BACKENDS_FOR_TESTS)
foreach(layer ${BWC_GFP_ARITHMETIC_BACKENDS})
    parse_arith_layer(${layer} dummy)
    set(layer ${canonical_arith_layer})

    string(REGEX MATCH "^(p1|p3|pz)$" t "${layer}")
    if(t)
        list(APPEND BWC_GFP_ARITHMETIC_BACKENDS_FOR_TESTS ${t})
    else()
        message(STATUS "Note: tests for ${layer} are skipped, because we restrict the routine tests")
    endif()
endforeach()

set(gdb_sha1_p1 a0fa9faf88f15e2429ea2c3f4eedf5dce721a606)
set(gdb_sha1_p3 a0fa9faf88f15e2429ea2c3f4eedf5dce721a606)
set(gdb_sha1_pz a0fa9faf88f15e2429ea2c3f4eedf5dce721a606)
set(gdb_sha1_b64 ef3c0bfb4bbdcb7eed1bdc52f58a1ade65a87b81)

foreach(layer ${BWC_GFP_ARITHMETIC_BACKENDS_FOR_TESTS} b64)
    parse_arith_layer(${layer} our_compile_flags)

    string(REGEX MATCH "p([0-9]+)" t "${layer}")

    set(extra)
    if(t)
        set(width ${CMAKE_MATCH_1})
        set(p "${prime_${width}}")
        set(extra "--prime ${p}")
    elseif(layer STREQUAL "pz")
        set(t TRUE)
        set(width "z")
        set(p "${prime_${width}}")
        set(extra "--prime ${p}")
    elseif(layer STREQUAL "b64")
        link_directories (${CADO_NFS_BINARY_DIR}/gf2x/fft/.libs)
        link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)
        include_directories(${CADO_NFS_BINARY_DIR}/gf2x/fft)
        include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft)
        list(APPEND our_compile_flags -DLINGEN_BINARY)
    else()
        message(FATAL_ERROR "bad arithmetic backend ${layer}")
    endif()

    add_executable(test_gdb_printers_bwc_${layer}
         test_gdb_printers_bwc.cpp)
    string_join(our_compile_flags_string " " ${our_compile_flags})
    set_target_properties(test_gdb_printers_bwc_${layer} PROPERTIES
        COMPILE_FLAGS ${our_compile_flags_string})

    target_link_libraries(test_gdb_printers_bwc_${layer}
        lingen_${layer}_support)
    # There should be a way to inherit this. Anyway, since
    # lingen_${layer}_support is mpi, this binary has to be mpi too!
    mark_targets_as_mpi(test_gdb_printers_bwc_${layer})

    add_dependencies(all_gdb_printer_tests test_gdb_printers_bwc_${layer})

    set(s ${gdb_sha1_${layer}})
    if(NOT s)
        set(s please_fill_sha1_value)
    endif()

    cado_define_test(
        test_gdb_printers_bwc_${layer}
        PROGRAM
        ${GDB_EXECUTABLE}
        -ex "add-auto-load-safe-path ${CADO_NFS_BINARY_DIR}"
        -ex "file ${CMAKE_CURRENT_BINARY_DIR}/test_gdb_printers_bwc_${layer}"
        -ex "set args ${extra}"
        -x ${CMAKE_CURRENT_SOURCE_DIR}/test_gdb_printers_bwc.gdb
        --batch
        TARGET_DEPENDENCIES
        test_gdb_printers_bwc_${layer}
        SHA1_ON_REGEXP_LINES "^\\$"
        EXPECT_SHA1 ${s}
        ENVIRONMENT
        CADO_NFS_GDB_SHOW_MATRIX_COMPACT=1
    )

    if(NOT ENABLE_SHARED)
        configure_file(${CADO_NFS_SOURCE_DIR}/utils/hook_cado_nfs_gdb.py.in
            test_gdb_printers_bwc_${layer}-gdb.py ESCAPE_QUOTES @ONLY)
    endif()
endforeach()
