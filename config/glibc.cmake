message(STATUS "Testing if the C library is GNU libc")
try_compile(HAVE_GLIBC
    ${PROJECT_BINARY_DIR}/config
    ${PROJECT_SOURCE_DIR}/config/glibc.c)
if(HAVE_GLIBC)
    message(STATUS "Testing if the C library is GNU libc -- Success")
else()
    message(STATUS "Testing if the C library is GNU libc -- Failed")
endif()

