
set(_compdefs -DULONG_BITS=${ULONG_BITS})
if(HAVE_GMP)
    set(_my_incdir ${GMP_INCDIR})
    set(_my_libdir ${GMP_LIBDIR})
    set(_compdefs ${_compdefs} -DHAVE_GMP=${HAVE_GMP})
elseif(HAVE_MPIR)
    set(_my_incdir ${MPIR_INCDIR})
    set(_my_libdir ${MPIR_LIBDIR})
    set(_compdefs ${_compdefs} -DHAVE_MPIR=${HAVE_MPIR})
endif()
message(STATUS "Testing for gmp_random predictability")
try_run(gmp_randstate_runs gmp_randstate_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/gmp_randstate.c
        CMAKE_FLAGS
        -DINCLUDE_DIRECTORIES=${_my_incdir}
        -DLINK_DIRECTORIES=${_my_libdir}
        COMPILE_DEFINITIONS ${_compdefs}
        LINK_LIBRARIES ${gmp_libname}
        COMPILE_OUTPUT_VARIABLE compilevar
        RUN_OUTPUT_VARIABLE outvar
        )


if(NOT gmp_randstate_runs EQUAL 0)
    message(FATAL_ERROR "Cannot run a simple gmp program\nCompile output (${gmp_randstate_compiles}): ${compilevar}\nRun output (${gmp_randstate_runs}): ${outvar}")
endif()

string(STRIP "${outvar}" outvar)
if("${outvar}" STREQUAL "d41c91186caf806b_45558c7335696741_71096848fde90ec7_7b34411325e1217a")
    message(STATUS "Testing for gmp_random predictability -- ok")
    set(HAVE_KNOWN_GMP_RANDOM_BEHAVIOUR 1)
else()
    message(STATUS "Testing for gmp_random predictability -- failed")
    message(STATUS "gmp_randstate test code returns [[${outvar}]]")
    message(STATUS "Warning: ${gmp_libname}'s internal random behaves differently from what we expect. This is not a big issue, but it is definitely going to ruin tests that rely on it giving a deterministic and known output once the seed is set")
    set(HAVE_KNOWN_GMP_RANDOM_BEHAVIOUR 0)
endif()
