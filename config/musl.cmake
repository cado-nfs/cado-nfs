message(STATUS "Testing if the C library is musl libc")
execute_process(
    COMMAND ${PROJECT_SOURCE_DIR}/config/musl.sh
    RESULT_VARIABLE MUSL_TEST_EXITCODE
    OUTPUT_VARIABLE MUSL_TEST_OUTPUT
)
set(HAVE_MUSL)
if(MUSL_TEST_EXITCODE EQUAL 0)
    set(HAVE_MUSL 1)
    STRING(STRIP ${MUSL_TEST_OUTPUT} MUSL_VERSION)
    message(STATUS "Testing if the C library is musl libc -- yes, version ${MUSL_VERSION}")
else()
    message(STATUS "Testing if the C library is musl libc -- no")
endif()

