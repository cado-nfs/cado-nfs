
set (sources 
        cado_poly.c timing.cpp getprime.c mod_ul.c relation.cpp
        relation-tools.c
        modredc_ul.c modredc_15ul.c modredc_2ul2.c modul_poly.c gmp_aux.c
        double_poly.c rootfinder.c params.c gcd.c
        mpz_array.c gzip.c misc.c
        mpz_poly.c crc.c purgedfile.c bit_vector.c
        barrier.c usp.c fix-endianness.c memusage.c
        cachesize_cpuid.c cachesize_guess.c roots_mod.c mod_mpz.c renumber.c
        ringbuf.c sm_utils.c memalloc.c mpz_vector.c memory.c stats.c lll.c
        verbose.c 
        version_info.c
	modified_files.c threadpool.cpp
)
if (NOT MINGW)
    set (sources ${sources} cado_popen.c)
endif()

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

add_library (utils_with_io ${CADO_SHARED} ${sources_io})
target_link_libraries(utils_with_io utils)
add_dependencies(utils_with_io antebuffer)

if(ENABLE_SHARED)
install(TARGETS utils utils_with_io DESTINATION ${LIBSUFFIX}/utils)
endif()

# To avoid the warning: implicit declaration of function nanosleep. It's POSIX, non c99.
set (CMAKE_OLD_C_FLAGS "${CMAKE_C_FLAGS}")
set (CMAKE_C_FLAGS "prout -D _POSIX_C_SOURCE=200809L ${CMAKE_C_FLAGS}")
add_executable(antebuffer antebuffer.c timing.cpp memusage.c)
set (CMAKE_C_FLAGS "${CMAKE_OLD_C_FLAGS}")

target_link_libraries(antebuffer pthread)

install(TARGETS antebuffer RUNTIME DESTINATION ${LIBSUFFIX}/utils)

include(CMakeLists-nodist.txt OPTIONAL)

