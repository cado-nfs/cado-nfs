
# ARM NEON
message(STATUS "Testing whether ARM NEON code can be used")
try_run(neon_runs neon_compiles
            ${PROJECT_BINARY_DIR}/config
            ${PROJECT_SOURCE_DIR}/config/neon.c)
if(neon_compiles)
    if (neon_runs MATCHES FAILED_TO_RUN)
        message(STATUS "Testing whether ARM NEON code can be used -- No")
        set (HAVE_ARM_NEON 0)
    else()
        message(STATUS "Testing whether ARM NEON code can be used -- Yes")
        set (HAVE_ARM_NEON 1)
    endif()
else()
    try_run(neon_runs neon_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/neon.c
        COMPILE_DEFINITIONS -mfpu=neon)
    if(neon_compiles)
        if (neon_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether ARM NEON code can be used -- No")
	    set (HAVE_ARM_NEON 0)
        else()
            message("${neon_runs}")
            message(STATUS "Testing whether ARM NEON code can be used -- Yes, with -mfpu=neon")
            set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfpu=neon")
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpu=neon")
            set (HAVE_ARM_NEON 1)
        endif()
    else()
        message(STATUS "Testing whether ARM NEON code can be used -- No (cannot compile)")
        set (HAVE_ARM_NEON 0)
    endif()
endif()
