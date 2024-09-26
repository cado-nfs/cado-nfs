# check_type_equality does _two_ things.
#
# 1 - check if the given type is exactly equal to one of the provided
# types, and define the corresponding macro FOO_IS_EXACTLY_BAR. It is
# not required that any of these succeeds. This is used to avoid
# exposing templates that would be ambiguous.
#
# 2 - check if the given type is compatible (size and
# signedness-wise) to one of the provided types, and define the
# corresponding macro FOO_IS_COMPATIBLE_WITH_BAR.  It *is* required
# that one of these succeeds. This is chiefly used to define MPI type
# aliases.


# include(CheckCXXSourceCompiles)
include(${CADO_NFS_SOURCE_DIR}/config/cado_check_cxx_source_compiles.cmake)

macro(testcode_type_exact_eq type1 type2)
    set(test_code "
#include <type_traits>
#include <stdint.h>
#include <stdio.h>
#include <gmp.h>
int main()
{
    static_assert(std::is_same<${type1}, ${type2}>::value, \"not this type\");
    return 0;
}
"
)
endmacro()

macro(testcode_type_compatible_eq type1 type2)
    set(test_code "
#include <type_traits>
#include <stdint.h>
#include <stdio.h>
#include <gmp.h>
int main()
{
    static_assert(std::is_same<${type1}, ${type2}>::value ||
        (sizeof(${type1}) == sizeof(${type2})
        && std::is_signed<${type1}>::value == std::is_signed<${type2}>::value)
        , \"not this type\");
    return 0;
}
"
)
endmacro()

set(CADO_C_NAME_UINT64_T "uint64_t")
set(CADO_C_NAME_UINT32_T "uint32_t")
set(CADO_C_NAME_INT64_T "int64_t")
set(CADO_C_NAME_INT32_T "int32_t")
set(CADO_C_NAME_UNSIGNED "unsigned int")
set(CADO_C_NAME_INT "int")
set(CADO_C_NAME_UNSIGNED_LONG "unsigned long")
set(CADO_C_NAME_UNSIGNED_LONG_LONG "unsigned long long")
set(CADO_C_NAME_LONG "long")
set(CADO_C_NAME_LONG_LONG "long long")
set(CADO_C_NAME_SIZE_T "size_t")
set(CADO_C_NAME_SSIZE_T "ssize_t")
set(CADO_C_NAME_MP_LIMB_T "mp_limb_t")
set(CADO_C_NAME_MP_SIZE_T "mp_size_t")
set(CADO_C_NAME_MPZ_INTERNAL_SIZE_T "decltype(__mpz_struct::_mp_size)")

set(CMAKE_REQUIRED_LINK_OPTIONS)
set(CMAKE_REQUIRED_LIBDIRS ${GMP_LIBDIR})
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_INCLUDES ${GMP_INCDIR})
set(CMAKE_REQUIRED_LIBRARIES ${gmp_libname})

macro(check_type_equality basetype)
    set(found_compatible)
    set(found_exact)
    set(must_list)
    foreach(t ${ARGN})
        if (${t} STREQUAL ${basetype})
            # do not rant about stuff like "long is exactly long"
            set(found_compatible 1)
        else()
            testcode_type_exact_eq("${CADO_C_NAME_${basetype}}" "${CADO_C_NAME_${t}}")
            CADO_CHECK_CXX_SOURCE_COMPILES("${test_code}" ${basetype}_IS_EXACTLY_${t})
            if(${basetype}_IS_EXACTLY_${t})
                message(STATUS "${CADO_C_NAME_${basetype}} == ${CADO_C_NAME_${t}} (exact match)")
                if(NOT found_compatible)
                    set(${basetype}_IS_COMPATIBLE_WITH_${t})
                    set(found_compatible 1)
                    set(CADO_MPI_${basetype} MPI_${t})
                endif()
            else()
                testcode_type_compatible_eq("${CADO_C_NAME_${basetype}}" "${CADO_C_NAME_${t}}")
                CADO_CHECK_CXX_SOURCE_COMPILES("${test_code}" ${basetype}_IS_COMPATIBLE_WITH_${t})
                if(${basetype}_IS_COMPATIBLE_WITH_${t})
                    message(STATUS "${CADO_C_NAME_${basetype}} == ${CADO_C_NAME_${t}} (compatibility match)")
                    if(NOT found_compatible)
                        # keep only the first match.
                        set(found_compatible 1)
                        set(CADO_MPI_${basetype} MPI_${t})
                    endif()
                endif()
            endif()
        endif()
        list(APPEND must_list "${CADO_C_NAME_${t}}")
    endforeach()
    if(NOT found_compatible)
        string_join(must " or " ${must_list})
        message(FATAL_ERROR "${CADO_C_NAME_${basetype}} must be either unsigned long or unsigned long long")
    endif()
endmacro()

check_type_equality(UINT64_T  UNSIGNED_LONG UNSIGNED_LONG_LONG)
check_type_equality(UNSIGNED_LONG_LONG  UNSIGNED_LONG UNSIGNED_LONG_LONG)
check_type_equality(UNSIGNED_LONG  UNSIGNED UNSIGNED_LONG)
check_type_equality(SIZE_T UNSIGNED_LONG UNSIGNED_LONG_LONG)
check_type_equality(MP_LIMB_T UNSIGNED_LONG UNSIGNED_LONG_LONG)
check_type_equality(UINT32_T UNSIGNED_LONG UNSIGNED)
check_type_equality(INT64_T LONG LONG_LONG)
check_type_equality(LONG_LONG  LONG LONG_LONG)
check_type_equality(LONG INT LONG)
check_type_equality(SSIZE_T LONG LONG_LONG)
check_type_equality(INT32_T LONG INT)
check_type_equality(MP_SIZE_T LONG_LONG LONG INT)
check_type_equality(MPZ_INTERNAL_SIZE_T LONG_LONG LONG INT)
