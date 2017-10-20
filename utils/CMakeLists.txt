
set (sources 
        cado_poly.c timing.cpp getprime.c mod_ul.c relation.cpp
        relation-tools.c
        modredc_ul.c modredc_15ul.c modredc_2ul2.c modul_poly.c gmp_aux.c
        double_poly.cpp rootfinder.c params.c gcd.c
        gzip.c misc.c
        mpz_mat.cpp mpz_poly.cpp crc.c purgedfile.c bit_vector.c
        barrier.c usp.c fix-endianness.c memusage.c
        cachesize_cpuid.c cachesize_guess.c roots_mod.c mod_mpz.c renumber.c
        ringbuf.c sm_utils.c galois_utils.c memalloc.c 
	mpz_vector.c memory.c stats.c lll.c
        verbose.c 
        polyroots.c 
        tdict.cpp
        version_info.c
	modified_files.c threadpool.cpp mpz_poly_bivariate.c
        gmpxx.cpp
        rho.c
        minimum_spanning_tree.cpp
)
if (NOT MINGW)
    set (sources ${sources} cado_popen.c)
endif()

# gmpxx.cpp now included unconditionally above.
#if (NOT HAVE_GMPXX AND NOT HAVE_MPIRXX)
#    set (sources ${sources} gmpxx.cpp)
#endif()

set (sources_io
        filter_io.cpp
)

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
target_link_libraries(utils pthread ${gmp_libname} m)
if (HAVE_GMPXX OR HAVE_MPIRXX)
    target_link_libraries(utils "${gmp_libname}xx")
endif()

add_library (utils_with_io ${CADO_SHARED} ${sources_io})
target_link_libraries(utils_with_io utils)

if(CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
    message(STATUS "Disabling antebuffer on ${CMAKE_SYSTEM_NAME}, see #21410")
    # http://cado-nfs.gforge.inria.fr/bug.php?21410
else()
    add_executable(antebuffer antebuffer.c timing.cpp memusage.c)
    target_link_libraries(antebuffer pthread)
    install(TARGETS antebuffer RUNTIME DESTINATION ${LIBSUFFIX}/utils)
    add_dependencies(utils_with_io antebuffer)
endif()

if(ENABLE_SHARED)
    install(TARGETS utils utils_with_io DESTINATION ${LIBSUFFIX}/utils)
endif()


include(CMakeLists-nodist.txt OPTIONAL)

