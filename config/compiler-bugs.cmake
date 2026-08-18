# test a bug appeared in Ubuntu/Linaro 4.6.3-1ubuntu5
# Since it is related to the treatment of amd64 asm constraints, we may
# skip it in other cases (or we get a spurious "Error (cannot compile)"
# message).
message(STATUS "Testing known compiler bugs")
set(COMPILER_BUGS_LIST)

if(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
    # Test for Ubuntu bug
    try_run(gcc-ubuntu-bug_runs gcc-ubuntu-bug_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/gcc-ubuntu-bug.c)

    if(gcc-ubuntu-bug_compiles)
        if (gcc-ubuntu-bug_runs EQUAL 0)
            set(VOLATILE_IF_GCC_UBUNTU_BUG 0)
        else()
            list(APPEND COMPILER_BUGS_LIST "Ubuntu/Linaro 4.6.3-1ubuntu5")
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
            list(APPEND COMPILER_BUGS_LIST "gcc bug 58805")
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
        list(APPEND COMPILER_BUGS_LIST "libstdc++ bug 114153")
        set(HAVE_LIBSTDCXX_BUG_114153 1)
    endif()
else()
    message(STATUS "Testing known compiler bugs -- Error (cannot compile)")
    set(HAVE_LIBSTDCXX_BUG_114153 0)
endif()

if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    # Test for llvm bug 201065
    try_run(llvm-bug-201065_runs llvm-bug-201065_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/llvm-bug-201065.c)

    if(llvm-bug-201065_compiles)
        if (llvm-bug-201065_runs EQUAL 0)
        else()
            list(APPEND COMPILER_BUGS_LIST "llvm-bug-201065")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mllvm -unroll-add-parallel-reductions=false")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mllvm -unroll-add-parallel-reductions=false")
        endif()
    else()
        message(STATUS "Testing known compiler bugs -- Error (cannot compile)")
    endif()
endif()



if(NOT COMPILER_BUGS_LIST)
    set(COMPILER_BUGS_LIST_STRING "None found")
else()
    list(JOIN COMPILER_BUGS_LIST ", " COMPILER_BUGS_LIST_STRING)
endif()

message(STATUS "Testing known compiler bugs -- ${COMPILER_BUGS_LIST_STRING}")

