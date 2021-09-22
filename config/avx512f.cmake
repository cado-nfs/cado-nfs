
# avx512f
message(STATUS "Testing whether avx512f code can be used")
if (HAVE_AVX)
    try_run(avx512f_runs avx512f_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/avx512f.c)
    if(avx512f_compiles)
        if (avx512f_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether avx512f code can be used -- No")
            set (HAVE_AVX512F 0)
        else()
            message(STATUS "Testing whether avx512f code can be used -- Yes")
            set (HAVE_AVX512F 1)
        endif()
    else()
        try_run(avx512f_runs avx512f_compiles
            ${PROJECT_BINARY_DIR}/config
            ${PROJECT_SOURCE_DIR}/config/avx512f.c
            COMPILE_DEFINITIONS -mavx512f)
        if(avx512f_compiles)
            if (avx512f_runs MATCHES FAILED_TO_RUN)
                message(STATUS "Testing whether avx512f code can be used -- No")
                set (HAVE_AVX512F 0)
            else()
                message(STATUS "Testing whether avx512f code can be used
                -- Yes, with -mavx512f")
                set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx512f")
                set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f")
                set (HAVE_AVX512F 1)
            endif()
        else()
            message(STATUS "Testing whether avx512f code can be used -- No")
            set (HAVE_AVX512F 0)
        endif()
    endif()
else()
message(STATUS "Testing whether avx512f code can be used -- skipped")
endif()
