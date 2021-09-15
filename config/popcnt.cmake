
# popcnt
message(STATUS "Testing whether popcnt code can be used")
if (HAVE_SSSE3)
    try_run(popcnt_runs popcnt_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/popcnt.c
        )
    if(popcnt_compiles)
        if (popcnt_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether popcnt code can be used -- No")
            set (HAVE_POPCNT 0)
        else()
            message(STATUS "Testing whether popcnt code can be used -- Yes")
            set (HAVE_POPCNT 1)
        endif()
    elseif(CMAKE_C_FLAGS MATCHES "-march")
       message(STATUS "Testing whether popcnt code can be used -- No (not testing -mpopcnt because -march is already present)")
    else()
        try_run(popcnt_runs popcnt_compiles
            ${PROJECT_BINARY_DIR}/config
            ${PROJECT_SOURCE_DIR}/config/popcnt.c
            COMPILE_DEFINITIONS -mpopcnt
            )
        if(popcnt_compiles)
            if (popcnt_runs MATCHES FAILED_TO_RUN)
                message(STATUS "Testing whether popcnt code can be used -- No")
                set (HAVE_POPCNT 0)
            else()
                message(STATUS "Testing whether popcnt code can be used -- Yes, with -mpopcnt")
                set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mpopcnt")
                set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
                set (HAVE_POPCNT 1)
            endif()
        else()
            message(STATUS "Testing whether popcnt code can be used -- No")
            set (HAVE_POPCNT 0)
        endif()
    endif()
else()
message(STATUS "Testing whether popcnt code can be used -- skipped")
endif()
