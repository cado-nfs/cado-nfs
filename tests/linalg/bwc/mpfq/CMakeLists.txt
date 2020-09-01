
include_directories(${PROJECT_SOURCE_DIR}/linalg/bwc/mpfq)

# because the source files we compile here are generated, we need to tell
# cmake that the current source dir does still serve a purpose...
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

foreach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS})

    add_custom_command(OUTPUT
        ${CMAKE_CURRENT_BINARY_DIR}/mpfq_test_${gfp_layer}.cpp
        COMMAND sed
            -e '1 s/^.*toplevel include.*/\#include \"cado.h\"/'
            -e 's/TAG/${gfp_layer}/g'
            -e 's,PATH_TO_SRC_TEST_DIRECTORY,${CMAKE_CURRENT_SOURCE_DIR}/,'
            -e 's/^\#warning.*//'
            -e 's/^\#define _POSIX_C_SOURCE.*//'
            -e 's/^\#define _GNU_SOURCE.*//'
            -e 's/^\#define _BSD_SOURCE.*//'
            -e 's/^\#define _DEFAULT_SOURCE.*//'
            < mpfq_test.cpp.meta
            > ${CMAKE_CURRENT_BINARY_DIR}/mpfq_test_${gfp_layer}.cpp
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DEPENDS mpfq_test.cpp.meta
        )

    set(toggle_variable_size_flag 0)
    string(REGEX MATCH "p_([0-9]+)" t "${gfp_layer}")
    if(t)
        set(width ${CMAKE_MATCH_1})
    else()
        if(gfp_layer STREQUAL "pz")
            set(t TRUE)
            set(width "z")
            set(toggle_variable_size_flag 1)
        endif()
    endif()
    set(p "${prime_${width}}")

    cado_define_test(${CMAKE_CURRENT_BINARY_DIR}/mpfq_test_${gfp_layer}.cpp
                ${PROJECT_SOURCE_DIR}/linalg/bwc/mpfq/mpfq_${gfp_layer}.c
                LIBRARIES ${gmp_libname}
                ARGUMENTS -s 1 -p ${p})
    if(toggle_variable_size_flag)
        set_target_properties(mpfq_test_${gfp_layer} PROPERTIES COMPILE_FLAGS "-DVARIABLE_SIZE_PRIME")
    endif()

endforeach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS})

# mpfq_test_fixmp is not seeded as of now, because it uses random2, and
# there is no mpn_rrandomb as of gmp-6.0.0
cado_define_test(mpfq_test_fixmp.c LIBRARIES ${gmp_libname} ARGUMENTS 40)
