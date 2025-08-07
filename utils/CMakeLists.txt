
set(fmt_sources)
if(NOT HAVE_FMT)
    set(fmt_sources
        embedded/fmt/format.cc
        embedded/fmt/os.cc
        )
endif()

set (sources 
        cado_poly.c timing.cpp getprime.c
        relation.cpp
        relation-tools.c
        arith/mod_ul.c arith/modredc_ul.c arith/modredc_15ul.c arith/modredc_2ul2.c arith/modul_poly.c
        arith/mod_mpz.c
        gmp_aux.c
        arithxx/mod64.cpp
        arithxx/modredc64.cpp
        arithxx/modredc96.cpp
        arithxx/modredc126.cpp
        arithxx/mod_mpz_new.cpp
        gmp_aux2.cpp
        double_poly.cpp rootfinder.cpp params.cpp params-extra.cpp gcd.c
        gzip.cpp misc.c misc2.cpp utils_cxx.cpp
        mpz_mat.cpp mpz_poly.cpp crc.c purgedfile.c bit_vector.c
        barrier.c usp.c fix-endianness.c memusage.cpp
        cachesize_cpuid.c cachesize_guess.c roots_mod.c
        renumber.cpp
        renumber_proxy.cpp
        ringbuf.c sm_utils.cpp galois_action.cpp memalloc.cpp
        special-q.cpp
	mpz_vector.c memory.c stats.c lll.c
        verbose.c 
        polyroots.c 
        tdict.cpp
        version_info.c
	modified_files.c threadpool.cpp mpz_poly_bivariate.cpp
        gmpxx.cpp
        rho.c gpf.cpp
        sha1.c
        fd_streambuf.cpp
        cado_pipe_streambuf.cpp
        numbertheory/number_field.cpp
        numbertheory/number_field_element.cpp
        numbertheory/number_field_fractional_ideal.cpp
        numbertheory/number_field_order.cpp
        numbertheory/number_field_order_element.cpp
        numbertheory/number_field_prime_ideal.cpp
        numbertheory/numbertheory_internals.cpp
        numbertheory/all_valuations_above_p.cpp
        mpz_mat_accessors.cpp
        badideals.cpp
        json.cpp
        indexed_relation.cpp
        cado-sighandlers.cpp
        mmap_allocator.cpp
        random_distributions.cpp
        imaginary_quadratic_class_groups.cpp
        ${fmt_sources}
)
if (NOT MINGW)
    list(APPEND sources cado_popen.c)
endif()
if(HAVE_HWLOC)
    list(APPEND sources hwloc-aux.c)
endif()
if(HAVE_MPFR)
    list(APPEND sources cxx_mpfr.cpp)
endif()

set (sources_io filter_io.cpp)

############################################################
# Get a list of modified files in the working directory,
# with their SHA1 sums, and write as a string constant to
# modified_files.c

# How we do this mechanism of generating the modified files is very quirky,
# and has varied with cmake versions. See http://public.kitware.com/Bug/view.php?id=11332
#
# The following seems to work also with pre-2.8.4 versions

add_custom_target(list_modified_files
    COMMAND "${CADO_NFS_SOURCE_DIR}/scripts/list_modified_files.sh" "${CMAKE_CURRENT_BINARY_DIR}/modified_files.c"
    WORKING_DIRECTORY "${CADO_NFS_SOURCE_DIR}"
    COMMENT "Generating list of modified files in working tree"
)
configure_file(version_info.c.in version_info.c ESCAPE_QUOTES @ONLY)
set_source_files_properties("${CMAKE_CURRENT_BINARY_DIR}/version_info.c" PROPERTIES GENERATED TRUE)

# the utils library is used by some shared objects (prominently the
# matmul shared libs in bwc), and therefore must be shared.
add_library (utils ${CADO_SHARED} ${sources})
set_source_files_properties("${CMAKE_CURRENT_BINARY_DIR}/modified_files.c" PROPERTIES GENERATED TRUE)
add_dependencies(utils list_modified_files)
mark_targets_as_openmp(utils)
target_link_libraries(utils pthread ${gmp_libname} m)
if(HAVE_HWLOC)
target_link_libraries(utils hwloc)
endif()
if(HAVE_FMT)
target_link_libraries(utils fmt)
endif()
if(HAVE_MPC)
target_link_libraries(utils mpc)
endif()
if(HAVE_MPFR)
target_link_libraries(utils mpfr)
endif()

add_library (utils_with_io ${CADO_SHARED} ${sources_io})
target_link_libraries(utils_with_io utils)

if(CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
    message(STATUS "Disabling antebuffer on ${CMAKE_SYSTEM_NAME}, see #21410")
    # https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues/21410
else()
    add_executable(antebuffer antebuffer.c timing.cpp memusage.cpp
        ${fmt_sources})
    target_link_libraries(antebuffer pthread)
    if(HAVE_FMT)
        target_link_libraries(antebuffer fmt)
    endif()
    install(TARGETS antebuffer RUNTIME DESTINATION ${LIBSUFFIX}/utils)
    add_dependencies(utils_with_io antebuffer)
endif()

if(ENABLE_SHARED)
    install(TARGETS utils utils_with_io DESTINATION ${LIBSUFFIX}/utils)
endif()

add_executable(numbertheory_tool numbertheory_tool.cpp)
target_link_libraries (numbertheory_tool utils)

install(TARGETS numbertheory_tool RUNTIME DESTINATION ${LIBSUFFIX}/utils)

if(ENABLE_SHARED AND HAVE_GDB_PYTHON_PRETTY_PRINT)
configure_file(hook_cado_nfs_gdb.py.in libutils.so-gdb.py ESCAPE_QUOTES @ONLY)
endif()

include(CMakeLists-nodist.txt OPTIONAL)

