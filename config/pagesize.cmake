
message(STATUS "Checking the OS page size")

try_run(pagesize_runs 
    pagesize_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/pagesize.c
    RUN_OUTPUT_VARIABLE OS_PAGE_SIZE
    )

if(pagesize_runs EQUAL 0)
    string(STRIP "${OS_PAGE_SIZE}" OS_PAGE_SIZE)
    message(STATUS "Checking the OS page size -- ${OS_PAGE_SIZE}")
else()
    message(FATAL_ERROR "Checking the OS page size -- Cannot run check (fatal)")
endif()
