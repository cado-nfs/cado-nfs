
# MMX
message(STATUS "Testing whether mmx code can be used")
try_run(mmx_runs mmx_compiles
            ${PROJECT_BINARY_DIR}/config
            ${PROJECT_SOURCE_DIR}/config/mmx.c)
if(mmx_compiles)
    if (mmx_runs MATCHES FAILED_TO_RUN)
        message(STATUS "Testing whether mmx code can be used -- No")
        set (HAVE_MMX 0)
    else()
        message(STATUS "Testing whether mmx code can be used -- Yes")
        set (HAVE_MMX 1)
    endif()
else()
    try_run(mmx_runs mmx_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/mmx.c
        COMPILE_DEFINITIONS -mmmx)
    if(mmx_compiles)
        if (mmx_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether mmx code can be used -- No")
	    set (HAVE_MMX 0)
        else()
            message("${mmx_runs}")
            message(STATUS "Testing whether mmx code can be used -- Yes, with -mmmx")
            set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mmmx")
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmmx")
            set (HAVE_MMX 1)
        endif()
    else()
        message(STATUS "Testing whether mmx code can be used -- No (cannot compile)")
        set (HAVE_MMX 0)
    endif()
endif()
