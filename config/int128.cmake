# __int128 types
message(STATUS "Testing whether __int128 is supported")
try_run(int128_runs int128_compiles
            ${PROJECT_BINARY_DIR}/config
            ${PROJECT_SOURCE_DIR}/config/int128.c)
if(int128_compiles AND (NOT (int128_runs MATCHES FAILED_TO_RUN)))
    message(STATUS "Testing whether __int128 is supported -- Yes")
    set (HAVE_INT128 1)
else()
    message(STATUS "Testing whether __int128 is supported -- No")
    set (HAVE_INT128 0)
endif()
