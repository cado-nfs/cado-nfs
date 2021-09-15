
# SSE-2
message(STATUS "Testing whether sse-2 code can be used")
try_run(sse2_runs sse2_compiles
            ${PROJECT_BINARY_DIR}/config
            ${PROJECT_SOURCE_DIR}/config/sse2.c)
if(sse2_compiles)
    if (sse2_runs MATCHES FAILED_TO_RUN)
        message(STATUS "Testing whether sse-2 code can be used -- No")
        set (HAVE_SSE2 0)
    else()
        message(STATUS "Testing whether sse-2 code can be used -- Yes")
        set (HAVE_SSE2 1)
    endif()
elseif(CMAKE_C_FLAGS MATCHES "-march")
   message(STATUS "Testing whether sse code can be used -- No (not testing -msse because -march is already present)")
else()
    try_run(sse2_runs sse2_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/sse2.c
        COMPILE_DEFINITIONS -msse2)
    if(sse2_compiles)
        if (sse2_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether sse-2 code can be used -- No")
	    set (HAVE_SSE2 0)
        else()
            message(STATUS "Testing whether sse-2 code can be used -- Yes, with -msse2")
            set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
            set (HAVE_SSE2 1)
        endif()
    else()
        message(STATUS "Testing whether sse-2 code can be used -- No (cannot compile)")
        set (HAVE_SSE2 0)
    endif()
endif()
