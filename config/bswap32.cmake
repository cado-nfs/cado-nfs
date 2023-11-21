
if(HAVE_SYS_ENDIAN_H)
    message(STATUS "Testing if we have bswap32 in sys/endian.h")
    try_compile(HAVE_BSWAP32_IN_SYS_ENDIAN_H
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/bswap32.c)
    if(HAVE_BSWAP32_IN_SYS_ENDIAN_H)
        message(STATUS "Testing if we have bswap32 in sys/endian.h -- Success")
    else()
        message(STATUS "Testing if we have bswap32 in sys/endian.h -- Failed")
    endif()
endif()

