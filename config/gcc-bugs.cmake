# test a bug appeared in Ubuntu/Linaro 4.6.3-1ubuntu5
# Since it is related to the treatment of amd64 asm constraints, we may
# skip it in other cases (or we get a spurious "Error (cannot compile)"
# message).
message(STATUS "Testing known compiler bugs")
set(GCC_BUGS_LIST)

if(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
    # Test for Ubuntu bug
    try_run(gcc-ubuntu-bug_runs gcc-ubuntu-bug_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/gcc-ubuntu-bug.c)

    if(gcc-ubuntu-bug_compiles)
        if (gcc-ubuntu-bug_runs EQUAL 0)
            set(VOLATILE_IF_GCC_UBUNTU_BUG 0)
        else()
            list(APPEND GCC_BUGS_LIST "Ubuntu/Linaro 4.6.3-1ubuntu5")
            set(VOLATILE_IF_GCC_UBUNTU_BUG 1)
        endif()
    else()
        message(STATUS "Testing known compiler bugs -- Error (cannot compile)")
        set(VOLATILE_IF_GCC_UBUNTU_BUG 0)
    endif()
else()
    set(VOLATILE_IF_GCC_UBUNTU_BUG 0)
endif()


if(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
    # Test for gcc bug 58805
    try_run(gcc-bug-58805_runs gcc-bug-58805_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/gcc-bug-58805.c)

    if(gcc-bug-58805_compiles)
        if (gcc-bug-58805_runs EQUAL 0)
            set(VOLATILE_IF_GCC_58805_BUG 0)
        else()
            list(APPEND GCC_BUGS_LIST "gcc bug 58805")
            set(VOLATILE_IF_GCC_58805_BUG 1)
        endif()
    else()
        message(STATUS "Testing known compiler bugs -- Error (cannot compile)")
        set(VOLATILE_IF_GCC_58805_BUG 0)
    endif()
else()
    set(VOLATILE_IF_GCC_58805_BUG 0)
endif()


# test for libstdc++ bug 114153
try_run(libstdc++-bug-114153_runs libstdc++-bug-114153_compiles
    ${PROJECT_BINARY_DIR}/config
    ${PROJECT_SOURCE_DIR}/config/libstdc++-bug-114153.cpp)
if(libstdc++-bug-114153_compiles)
    if (libstdc++-bug-114153_runs EQUAL 0)
        set(HAVE_LIBSTDCXX_BUG_114153 0)
    else()
        list(APPEND GCC_BUGS_LIST "libstdc++ bug 114153")
        set(HAVE_LIBSTDCXX_BUG_114153 1)
    endif()
else()
    message(STATUS "Testing known compiler bugs -- Error (cannot compile)")
    set(HAVE_LIBSTDCXX_BUG_114153 0)
endif()


if(NOT GCC_BUGS_LIST)
    set(GCC_BUGS_LIST_STRING "None found")
else()
    list(JOIN GCC_BUGS_LIST ", " GCC_BUGS_LIST_STRING)
endif()

message(STATUS "Testing known compiler bugs -- ${GCC_BUGS_LIST_STRING}")

