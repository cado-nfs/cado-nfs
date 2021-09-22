
# avx512vl
message(STATUS "Testing whether avx512vl code can be used")
if (HAVE_AVX512F)
    try_run(avx512vl_runs avx512vl_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/avx512vl.c)
    if(avx512vl_compiles)
        if (avx512vl_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether avx512vl code can be used -- No")
            set (HAVE_AVX512VL 0)
        else()
            message(STATUS "Testing whether avx512vl code can be used -- Yes")
            set (HAVE_AVX512VL 1)
        endif()
    else()
        try_run(avx512vl_runs avx512vl_compiles
            ${PROJECT_BINARY_DIR}/config
            ${PROJECT_SOURCE_DIR}/config/avx512vl.c
            COMPILE_DEFINITIONS -mavx512vl)
        if(avx512vl_compiles)
            if (avx512vl_runs MATCHES FAILED_TO_RUN)
                message(STATUS "Testing whether avx512vl code can be used -- No")
                set (HAVE_AVX512VL 0)
            else()
                message(STATUS "Testing whether avx512vl code can be used
                -- Yes, with -mavx512vl")
                set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx512vl")
                set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512vl")
                set (HAVE_AVX512VL 1)
            endif()
        else()
            message(STATUS "Testing whether avx512vl code can be used -- No")
            set (HAVE_AVX512VL 0)
        endif()
    endif()
else()
message(STATUS "Testing whether avx512vl code can be used -- skipped")
endif()
