
# avx512dq
message(STATUS "Testing whether avx512dq code can be used")
if (HAVE_AVX512F)
    try_run(avx512dq_runs avx512dq_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/avx512dq.c)
    if(avx512dq_compiles)
        if (avx512dq_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether avx512dq code can be used -- No")
            set (HAVE_AVX512DQ 0)
        else()
            message(STATUS "Testing whether avx512dq code can be used -- Yes")
            set (HAVE_AVX512DQ 1)
        endif()
    else()
        try_run(avx512dq_runs avx512dq_compiles
            ${PROJECT_BINARY_DIR}/config
            ${PROJECT_SOURCE_DIR}/config/avx512dq.c
            COMPILE_DEFINITIONS -mavx512dq)
        if(avx512dq_compiles)
            if (avx512dq_runs MATCHES FAILED_TO_RUN)
                message(STATUS "Testing whether avx512dq code can be used -- No")
                set (HAVE_AVX512DQ 0)
            else()
                message(STATUS "Testing whether avx512dq code can be used -- Yes, with -mavx512dq")
                set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx512dq")
                set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512dq")
                set (HAVE_AVX512DQ 1)
            endif()
        else()
            message(STATUS "Testing whether avx512dq code can be used -- No")
            set (HAVE_AVX512DQ 0)
        endif()
    endif()
else()
message(STATUS "Testing whether avx512dq code can be used -- skipped")
endif()
