if(HAVE_AVX512F AND (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang" OR
    CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_CXX_COMPILER_ID
    MATCHES "Clang") OR (
    CMAKE_CXX_COMPILER_ID MATCHES "^Clang" AND
    NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "13"
    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15"))

    message(STATUS "Testing if bug 30057 (llvm bug 53842) is present")
    execute_process(COMMAND
        ${CMAKE_CXX_COMPILER} -O3 -mllvm -mattr=+avx512f
        INPUT_FILE ${PROJECT_SOURCE_DIR}/config/bug-llvm53842.bc
        OUTPUT_QUIET
        ERROR_VARIABLE llc_err)
    # try_compile(bug_llvm53842_compiles ${PROJECT_BINARY_DIR}/config ${PROJECT_SOURCE_DIR}/config/bug-llvm53842.bc)
    if(NOT llc_err)
        message(STATUS "Testing if llvm bug 30057 (llvm bug 53842) is present -- No")
    else()
        message(STATUS "Testing if llvm bug 30057 (llvm bug 53842) is present -- Yes, disabling avx512f. Use a newer compiler to fix. See https://gitlab.inria.fr/cado-nfs/cado-nfs/issues/30057")
        string(REPLACE "-mavx512f" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
        string(REPLACE "-mavx512f" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        set (HAVE_AVX512F 0)
    endif()
endif()
