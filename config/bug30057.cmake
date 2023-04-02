if(HAVE_AVX512F AND (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang" OR
    CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_CXX_COMPILER_ID
    MATCHES "Clang") OR (
    CMAKE_CXX_COMPILER_ID MATCHES "^Clang" AND
    NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "13"
    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15"))

    message(STATUS "Testing if bug 30057 (llvm bug 53842) is present")
    execute_process(COMMAND
        ${CMAKE_CXX_COMPILER} --print-prog-name llc
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE LLC_PATH
        ERROR_VARIABLE llc_path_err)
    if(llc_path_err OR NOT LLC_PATH)
        message(STATUS ${llc_path_err})
        message(FATAL_ERROR "${CMAKE_CXX_COMPILER} won't tell us where llc is!")
    endif()
    # message(STATUS "${LLC_PATH} -mattr=+avx512f < ${PROJECT_SOURCE_DIR}/config/bug30057.bc")
    execute_process(COMMAND
        ${LLC_PATH} -mattr=+avx512f
        INPUT_FILE ${PROJECT_SOURCE_DIR}/config/bug30057.bc
        OUTPUT_QUIET
        ERROR_VARIABLE llc_err)
    if(NOT llc_err)
        message(STATUS "Testing if llvm bug 30057 (llvm bug 53842) is present -- No")
    else()
        message(STATUS "Testing if llvm bug 30057 (llvm bug 53842) is present -- Yes, disabling avx512. Use a newer compiler to fix. See https://gitlab.inria.fr/cado-nfs/cado-nfs/issues/30057")
        string(REPLACE "-mavx512f" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
        string(REPLACE "-mavx512f" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        set (HAVE_AVX512F 0)
    endif()
endif()
