
set(cado_precompiled_header_set
    "${CADO_NFS_BINARY_DIR}/cado_config.h"
    "${CADO_NFS_SOURCE_DIR}/cado.h"
    "$<$<COMPILE_LANGUAGE:CXX>:<algorithm$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<array$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<atomic$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<barrier$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cassert$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cctype$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cerrno$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cfenv$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cfloat$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<charconv$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cinttypes$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<climits$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<clocale$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cmath$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<compare$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<complex$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<concepts$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<condition_variable$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<csignal$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cstdarg$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cstddef$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cstdint$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cstdio$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cstdlib$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<cstring$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<ctime$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<deque$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<exception$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<filesystem$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<forward_list$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<fstream$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<functional$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<future$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<gmp.h$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<initializer_list$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<iomanip$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<ios$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<iosfwd$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<iostream$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<istream$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<iterator$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<limits$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<list$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<locale$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<map$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<math.h$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<memory$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<mutex$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<new$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<numbers$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<numeric$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<ostream$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<queue$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<ranges$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<sstream$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<streambuf$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<string$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<system_error$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<thread$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<tuple$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<type_traits$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<typeinfo$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<unordered_map$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<utility$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<valarray$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<variant$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<vector$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<version$<ANGLE-R>>"
)

if (ENABLE_SHARED)
    set(is_pic_choices OFF ON)
else()
    set(is_pic_choices OFF)
endif()

if (OPENMP_FOUND)
    set(is_openmp_choices OFF ON)
else()
    set(is_openmp_choices OFF)
endif()

if (WITH_MPI)
    set(is_mpi_choices OFF ON)
else()
    set(is_mpi_choices OFF)
endif()

# Some build configs are known to cause compiler errors. We're not that
# interested in chasing down those bugs. Better discard the PCH setting
# in such cases.
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
        AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0)
    message(STATUS "Disabling precompiled headers for this build configuration")
    set(WITH_PCH 0)
else()
    set(WITH_PCH 1)
endif()

if (WITH_PCH)
set(PCH_DUMMY_C "${CMAKE_CURRENT_BINARY_DIR}/cado_pch_dummy.c")
file(WRITE "${PCH_DUMMY_C}" "// Dummy file for C PCH generation\n")

set(PCH_DUMMY_CXX "${CMAKE_CURRENT_BINARY_DIR}/cado_pch_dummy.cpp")
file(WRITE "${PCH_DUMMY_CXX}" "// Dummy file for CXX PCH generation\n")

foreach(is_pic IN ITEMS ${is_pic_choices})
    foreach(is_omp IN ITEMS ${is_openmp_choices})
        foreach(is_mpi IN ITEMS ${is_mpi_choices})
            set(pch_name "cado_pch")
            if(is_pic)
                string(APPEND pch_name "_pic")
            endif()
            if(is_omp)
                string(APPEND pch_name "_omp")
            endif()
            if(is_mpi)
                string(APPEND pch_name "_mpi")
            endif()

            add_library(${pch_name} OBJECT "${PCH_DUMMY_CXX}" "${PCH_DUMMY_C}")
            set_target_properties(${pch_name} PROPERTIES POSITION_INDEPENDENT_CODE ${is_pic})
            target_precompile_headers(${pch_name} PRIVATE ${cado_precompiled_header_set})

            # Apply the flags to the PCH itself so the AST matches!
            if(is_omp)
                mark_targets_as_openmp(${pch_name})
            endif()
            if(is_mpi)
                mark_targets_as_mpi(${pch_name})
            endif()
        endforeach()
    endforeach()
endforeach()

function(cado_apply_pch target)
    if(NOT TARGET ${target})
        return()
    endif()

    if("${target}" MATCHES "^cado_pch")
        return()
    endif()

    get_target_property(target_type ${target} TYPE)
    set(is_pic OFF)
    if(target_type STREQUAL "SHARED_LIBRARY")
        set(is_pic ON)
    elseif(target_type STREQUAL "OBJECT_LIBRARY")
        get_target_property(pic_prop ${target} POSITION_INDEPENDENT_CODE)
        if(pic_prop)
            set(is_pic ON)
        endif()
    endif()

    get_target_property(is_omp ${target} CADO_USES_OPENMP)
    if(NOT is_omp)
        set(is_omp OFF)
    endif()

    get_target_property(is_mpi ${target} CADO_USES_MPI)
    if(NOT is_mpi)
        set(is_mpi OFF)
    endif()

    set(pch_name "cado_pch")
    if(is_pic)
        string(APPEND pch_name "_pic")
    endif()
    if(is_omp)
        string(APPEND pch_name "_omp")
    endif()
    if(is_mpi)
        string(APPEND pch_name "_mpi")
    endif()

    target_precompile_headers(${target} REUSE_FROM ${pch_name})
endfunction()
else()
    function(cado_apply_pch target)
    endfunction()
endif()
