
include_directories(. ..)
# include_directories(mpfq)

# Use with care. Normally, binaries not needing MPI for real should not
# emit MPI symbols since they don't do anything depending on this config
# variable. However this mechanism is somewhat fragile.
if (${WITH_MPI})
add_definitions(-DWITH_MPI)
endif()

### bwc-specific configuration.

# unset(BUILD_DYNAMICALLY_LINKABLE_BWC)
# Dynamically linkable bwc needs shared libs.
set(BUILD_DYNAMICALLY_LINKABLE_BWC ${ENABLE_SHARED})
# having both shared and static is possible but messy. One needs two
# add_library calls, and then I'm not sure whether
# target_link_libraries picks the static or shared ones...

set(BWC_GF2_ARITHMETIC_BACKENDS m128 u64k1 u64k2 u64k3 u64k4
    CACHE STRING
    "which low-level arithmetic backends to use for GF(2) linear algebra")
set(BWC_GF2_MATMUL_BACKENDS bucket basic sliced
    CACHE STRING
    "which low-level matrix multiplication backends to use for GF(2) linear algebra")
set(BWC_GFP_ARITHMETIC_BACKENDS p_1 p_2 p_3 p_4 p_5 p_6 p_7 p_8 p_9 p_10 p_11 p_12 p_13 p_14 p_15 pz
    CACHE STRING
    "which GF(p) backends to compile for DLOG mode (pz excluded)")
if(HAVE_CXX11)
    set(bmb basicp zone)
else()
    set(bmb basicp)
endif()
set(BWC_GFP_MATMUL_BACKENDS ${bmb}
    CACHE STRING
    "which low-level matrix multiplication backends to use for GF(2) linear algebra")

# In the fflas branch, this is used to plug in fflas.
set(BWC_EXTRA_BACKENDS "" CACHE STRING "list of pairs (arithmetic,matmul) backends for linear algebra")

# In case there's a user override, just make mention of the situation.
if(NOT BWC_GFP_ARITHMETIC_BACKENDS OR NOT HAVE_GCC_STYLE_AMD64_INLINE_ASM)
    # The code in linalg/bwc/mpfq/ is 64-bit only, and relies on GNU
    # CC idioms. We don't want to expose it on 32-bit platforms.
    set(BWC_GFP_ARITHMETIC_BACKENDS CACHE STRING
        "which GF(p) backends to compile for DLOG mode" FORCE)
    message(STATUS "Warning: DL mode forcibly disabled for linear algebra (amd64-only)")
endif()

### end bwc-specific configuration.


# These are just shorthands. Some targets also like to define compile and
# link flags for their own, so it's really messy, because I don't know
# how to have cmake *append* to the target's link flags.

# Note that whenever a target (at least very certainly an executable)
# links an mpi compiled library, it must be marked as mpi as well.

macro(mark_mpi_target target)
    set_target_properties(${target} PROPERTIES
        COMPILE_FLAGS "--mpi"
        LINK_FLAGS "--mpi")
endmacro()
macro(mpi_add_library target)
    add_library(${target} ${ARGN})
    mark_mpi_target(${target})
endmacro()
macro(mpi_add_executable target)
    add_executable(${target} ${ARGN})
    mark_mpi_target(${target})
endmacro()

### bwc_base

# Note that because the balancing code contains balancing_get_matrix_u32,
# which calls the parallelizing_info layer, it's difficult to make them
# two distinct libraries. So far it's the only stopper though.
set(balancing_sources
    balancing.c
    balancing_workhorse.c
    balancing_file_source.c
    balancing_mpi_source.c)

if (HAVE_CURL)
set(balancing_sources ${balancing_sources} balancing_curl_source.c)
endif()

set(bwc_base_sources
    parallelizing_info.c
    intersections.c
    xvectors.c
    bw-common.c
    async.c
    rolling.c
    matmul.c
    matmul_top.c
    random_matrix.c
    mf.c 
    mf_bal.c
    ../rowset_heap.cpp
    ${balancing_sources}
)

if(HAVE_CXX11 AND HAVE_HWLOC)
set(bwc_base_sources
    ${bwc_base_sources}
    cpubinding.cpp
    )
endif()

mpi_add_library(bwc_base ${CADO_SHARED} ${bwc_base_sources})

target_link_libraries(bwc_base ${gmp_libname})
if (HAVE_CURL)
target_link_libraries(bwc_base curl)
endif()

if(HAVE_CXX11 AND HAVE_HWLOC)
target_link_libraries(bwc_base hwloc)
endif()


if(ENABLE_SHARED)
    install(TARGETS bwc_base DESTINATION ${LIBSUFFIX}/linalg/bwc)
endif()

### end bwc_base

mpi_add_executable(mf_scan ${EXCLUDE_FROM_RELEASE_BUILD} mf_scan.c)
target_link_libraries(mf_scan bwc_base)
if (HAVE_OPENMP)
mpi_add_executable(mf_scan2 ${EXCLUDE_FROM_RELEASE_BUILD} mf_scan2.cpp)
target_link_libraries(mf_scan2 bwc_base)
set_target_properties(mf_scan2 PROPERTIES
    COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
    LINK_FLAGS "${OpenMP_CXX_FLAGS}")
endif()


include_directories(${CADO_NFS_SOURCE_DIR}/linalg)

mpi_add_executable(mf_bal mf_bal_main.c)
target_link_libraries(mf_bal bwc_base m)
install(TARGETS mf_bal RUNTIME DESTINATION ${LIBSUFFIX}/linalg/bwc)


configure_file(bwc.pl ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl ESCAPE_QUOTES @ONLY)

add_executable(random_matrix ${EXCLUDE_FROM_RELEASE_BUILD} random_matrix.c)
set_property(TARGET random_matrix PROPERTY COMPILE_DEFINITIONS WANT_MAIN)
target_link_libraries(random_matrix utils m)
set_target_properties(random_matrix PROPERTIES
        COMPILE_FLAGS "--mpi"
        LINK_FLAGS "--mpi")

### bwc_mpfq

set(bwc_mpfq_sources mpfq/mpfq_vbase.c)
set(bwc_mpfq_flags)
foreach(v ${BWC_GF2_ARITHMETIC_BACKENDS})
    if (v STREQUAL "m128" AND NOT HAVE_SSE41)
        message(STATUS "GF(2) arithmetic backend m128 not compiled in, sse-4.1 support missing")
    else()
        set(bwc_mpfq_sources ${bwc_mpfq_sources} mpfq/mpfq_${v}.c mpfq/mpfq_${v}.cpp mpfq/mpfq_${v}_t.c)
        set(COMPILE_MPFQ_BINARY_FIELD_${v} 1 CACHE INTERNAL "${v} is compiled in")
        set(bwc_mpfq_flags "${bwc_mpfq_flags} -DCOMPILE_MPFQ_BINARY_FIELD_${v}")
        message(STATUS "GF(2) arithmetic backend enabled for factoring linear algebra: ${v}")
    endif()
endforeach(v)

foreach(v ${BWC_GFP_ARITHMETIC_BACKENDS})
    set(bwc_mpfq_sources ${bwc_mpfq_sources} mpfq/mpfq_${v}.c mpfq/mpfq_${v}.cpp mpfq/mpfq_${v}_t.c)
    set(COMPILE_MPFQ_PRIME_FIELD_${v} 1 CACHE INTERNAL "${v} is compiled in")
    set(bwc_mpfq_flags "${bwc_mpfq_flags} -DCOMPILE_MPFQ_PRIME_FIELD_${v}")
    message(STATUS "GF(p) arithmetic backend enabled for DL linear algebra: ${v}")
endforeach(v)

add_library(bwc_mpfq ${CADO_SHARED} ${bwc_mpfq_sources})

if(ENABLE_SHARED)
    install(TARGETS bwc_mpfq DESTINATION ${LIBSUFFIX}/linalg/bwc)
endif()

target_link_libraries(bwc_mpfq ${gmp_libname} utils)

set_target_properties(bwc_mpfq PROPERTIES COMPILE_FLAGS "--mpi ${bwc_mpfq_flags}")
set_target_properties(bwc_mpfq PROPERTIES LINK_FLAGS "--mpi")

### end bwc_mpfq


### link dependencies

IF(UNIX)
    FIND_PROGRAM(CMAKE_UNAME uname /bin /usr/bin /usr/local/bin )
    IF(CMAKE_UNAME)
        EXEC_PROGRAM(uname ARGS -s OUTPUT_VARIABLE CMAKE_SYSTEM_NAME)
    ENDIF()
ENDIF()
IF(CMAKE_SYSTEM_NAME MATCHES BSD)
    target_link_libraries(bwc_base utils bwc_mpfq ${pthread_libs})
ELSE()
    if(BUILD_DYNAMICALLY_LINKABLE_BWC)
        target_link_libraries(bwc_base utils bwc_mpfq ${pthread_libs} dl)
    else()
        target_link_libraries(bwc_base utils bwc_mpfq ${pthread_libs})
    endif()
ENDIF()

# matmul_common is a dependency of all the impl-specific libs. When doing
# dynamic linking, there's no real need to make matmul_common *also*
# shared, as embarking the code within the impl-specific .so is fine.
# Unfortunately, we do want matmul_common to be compiled with -fPIC,
# then. Which is easiest to trigger by making it a shared object.
mpi_add_library(matmul_common ${CADO_SHARED} matmul-common.c)
if(ENABLE_SHARED)
    install(TARGETS matmul_common DESTINATION ${LIBSUFFIX}/linalg/bwc)
endif()


set(matmul_basic_sources    matmul-basic.c)
set(matmul_basicp_sources   matmul-basicp.cpp)
set(matmul_zone_sources   matmul-zone.cpp)
set(matmul_sliced_sources   matmul-sliced.cpp)
set(matmul_bucket_sources   matmul-bucket.cpp)
set(matmul_threaded_sources matmul-threaded.c worker-threads.c)

macro(CONFIGURE_MATMUL_LIB v i)
    if(v STREQUAL "pz")
        # message(STATUS "Information: pz matmul backends are now disabled")
    else()
    set(COOKED_BWC_BACKENDS ${COOKED_BWC_BACKENDS} "DO(${v}, ${i})")
    string(COMPARE EQUAL "${v}" "u64k1" doing_u64)
    string(COMPARE EQUAL "${i}" "bucket" doing_bucket)
    set(sources matmul_facade.c ${matmul_${i}_sources})
    if(doing_u64 AND doing_bucket)
        if(HAVE_GAS_SYNTAX_ASSEMBLY_SOURCES)
            set(sources ${sources}
                    matmul-sub-small1.S
                    matmul-sub-small2.S
                    matmul-sub-large-fbi.S
                    matmul-sub-large-fbd.S
                    matmul-sub-vsc-combine.S
                    matmul-sub-vsc-dispatch.S
                    )
            set_property(SOURCE matmul-sub-small1.S PROPERTY LANGUAGE C)
            set_property(SOURCE matmul-sub-small2.S PROPERTY LANGUAGE C)
            set_property(SOURCE matmul-sub-large-fbi.S PROPERTY LANGUAGE C)
            set_property(SOURCE matmul-sub-large-fbd.S PROPERTY LANGUAGE C)
            set_property(SOURCE matmul-sub-vsc-dispatch.S PROPERTY LANGUAGE C)
            set_property(SOURCE matmul-sub-vsc-combine.S PROPERTY LANGUAGE C)
        endif()
    endif()

        add_library(matmul_${v}_${i} ${CADO_SHARED} ${sources})
        target_link_libraries(matmul_${v}_${i} matmul_common bwc_mpfq)
        set_target_properties(matmul_${v}_${i} PROPERTIES
            COMPILE_FLAGS "--mpi -DSELECT_MPFQ_LAYER_${v} -DMM_IMPL=${i} -DMM_MPFQ_LAYER=${v}"
                LINK_FLAGS "--mpi")
        set(all_matmul_libs ${all_matmul_libs} "matmul_${v}_${i}")
        if(ENABLE_SHARED)
            install(TARGETS matmul_${v}_${i} DESTINATION ${LIBSUFFIX}/linalg/bwc)
        endif()
    endif()
endmacro()

string_join(FOO ", " ${BWC_GF2_MATMUL_BACKENDS})
message(STATUS "GF(2) SpMV backends for linear algebra: ${FOO}")
if (BWC_GFP_ARITHMETIC_BACKENDS)
string_join(FOO ", " ${BWC_GFP_MATMUL_BACKENDS})
message(STATUS "GF(p) SpMV backends for linear algebra: ${FOO}")
endif()



foreach(v ${BWC_GF2_ARITHMETIC_BACKENDS})
    if (v STREQUAL "m128" AND NOT HAVE_SSE41)
    else()
        foreach(w ${BWC_GF2_MATMUL_BACKENDS})
            CONFIGURE_MATMUL_LIB(${v} ${w})
        endforeach(w)
    endif()
endforeach(v)
foreach(v ${BWC_GFP_ARITHMETIC_BACKENDS})
    foreach(w ${BWC_GFP_MATMUL_BACKENDS})
        CONFIGURE_MATMUL_LIB(${v} ${w})
    endforeach(w)
endforeach(v)

foreach(vw ${BWC_EXTRA_BACKENDS})
    string_split("${vw}" "+" FOO)
    message(STATUS "output variable: ${FOO}")
    list(GET FOO 0 v)
    list(GET FOO 1 w)
    message(STATUS "Extra bwc backend: arithmetic=${v}, matmul=${w}")
    CONFIGURE_MATMUL_LIB(${v} ${w})
endforeach(vw ${BWC_EXTRA_BACKENDS})

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/bwc_config.h)
    message(FATAL_ERROR "Please remove the file ${CMAKE_CURRENT_SOURCE_DIR}/bwc_config.h before proceeding")
endif()
configure_file(bwc_config_h.in bwc_config.h ESCAPE_QUOTES @ONLY)


if(BUILD_DYNAMICALLY_LINKABLE_BWC)
set(MATMUL_LIBS_PREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
set(MATMUL_LIBS_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
endif()

configure_file(matmul-libnames.h.in matmul-libnames.h ESCAPE_QUOTES @ONLY)
include_directories (${CMAKE_CURRENT_BINARY_DIR})

if(NOT BUILD_DYNAMICALLY_LINKABLE_BWC)
target_link_libraries(bwc_base ${all_matmul_libs})
endif()

### matmul_mf
mpi_add_library(matmul_mf ${CADO_SHARED} matmul-mf.c)
target_link_libraries(matmul_mf bwc_base)
if(ENABLE_SHARED)
    install(TARGETS matmul_mf DESTINATION ${LIBSUFFIX}/linalg/bwc)
endif()
### end matmul_mf

macro(bwc_program_base v)
    mpi_add_executable(${v} ${ARGN})
    set_target_properties(${v} PROPERTIES
            LINKER_LANGUAGE CXX)
        # -Wl,-rpath,${CMAKE_CURRENT_BINARY_DIR}")
    target_link_libraries(${v} bwc_base)
endmacro()

# The macro bwc_program_extra is mostly used by the binaries specific to
# the git tree.
macro(bwc_program_extra v)
    bwc_program_base(${v} EXCLUDE_FROM_ALL ${ARGN})
endmacro()

macro(bwc_program v)
    bwc_program_base(${v} ${ARGN})
    install(TARGETS ${v} RUNTIME DESTINATION ${LIBSUFFIX}/linalg/bwc)
endmacro()

bwc_program(acollect acollect.c)
bwc_program(prep prep.cpp)
target_link_libraries(prep bitlinalg)
bwc_program(krylov krylov.cpp xdotprod.c)
bwc_program(bwccheck bwccheck.cpp xdotprod.c)
bwc_program(dispatch dispatch.cpp xdotprod.c)
bwc_program(mksol mksol.cpp xdotprod.c)
bwc_program(gather gather.cpp)
bwc_program(secure secure.cpp)
bwc_program(blocklanczos blocklanczos.c)
target_link_libraries(blocklanczos bitlinalg)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl DESTINATION
    ${LIBSUFFIX}/linalg/bwc)

# Also lingen.

### This is the old binary lingen code. It is meant to go away.
link_directories (${CADO_NFS_BINARY_DIR}/gf2x/fft/.libs)
link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)
include_directories(${CADO_NFS_BINARY_DIR}/gf2x/fft)
include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft)
add_executable(lingen lingen_binary.cpp tree_stats.cpp logline.cpp
    lingen_round_operand_size.cpp
    lingen_qcode_binary.cpp
    # matpoly isn't really used by the legacy lingen, but the basecase code in
    # qcode_binary is shared between the two implementations now, and some of
    # the interface there relies on matpoly being available.
    # the two other files are just pulled by this.
    lingen_matpoly_binary.cpp
    lingen_memory_pool.cpp
    lingen_bmstatus.cpp
    lingen_expected_pi_length.cpp
    )
# This **MUST** match the #define in gf2x/Makefile.am
set_source_files_properties(lingen_qcode_binary.cpp COMPILE_FLAGS
    "-DCANTOR_BASE_FIELD_SIZE=128")
set_source_files_properties(lingen_binary.cpp COMPILE_FLAGS
    "-DCANTOR_BASE_FIELD_SIZE=128")
# for debug only
# set_source_files_properties(lingen_binary.cpp COMPILE_FLAGS
# "-DCANTOR_BASE_FIELD_SIZE=128 -DLINGEN_BINARY_TRACE_MALLOCS=20")
install(TARGETS lingen RUNTIME DESTINATION ${LIBSUFFIX}/linalg/bwc)
target_link_libraries(lingen bwc_base utils gf2x-fft gf2x)
add_dependencies(lingen gf2x-build)
# It's ugly. FindOpenMP.cmake doesn't know how to differentiate between
# required compile-time flags and required link-time flags...
set_target_properties(lingen PROPERTIES
    COMPILE_FLAGS
    "-DSELECT_MPFQ_LAYER_u64k1 --mpi ${OpenMP_CXX_FLAGS}"
    LINK_FLAGS "--mpi ${OpenMP_CXX_FLAGS}")




add_executable(cleanup cleanup.c)
target_link_libraries(cleanup bitlinalg utils m)

install(TARGETS cleanup RUNTIME DESTINATION ${LIBSUFFIX}/linalg/bwc)

if(BWC_GFP_ARITHMETIC_BACKENDS)
    # This is a convenience library which is used **ONLY** by one binary,
    # so there's absolutely no point in making it shared.
    file(GLOB flint_fft_files flint-fft/*.c)
    add_library(flint-fft ${flint_fft_files})
    set_target_properties(flint-fft PROPERTIES
        COMPILE_FLAGS "--mpi ${OpenMP_C_FLAGS}"
        LINK_FLAGS "--mpi ${OpenMP_C_FLAGS}")
    add_library(flint-fft-debug EXCLUDE_FROM_ALL ${flint_fft_files})
    set_target_properties(flint-fft-debug PROPERTIES
        COMPILE_FLAGS "--mpi ${OpenMP_C_FLAGS} -DDEBUG_FFT"
        LINK_FLAGS "--mpi ${OpenMP_C_FLAGS}")
endif()

# These files are not dependent on the mpfq layer, therefore it's ok tu
# put all of them in a common convenience library.
add_library(lingen_common_support 
    lingen_memory_pool.cpp
    logline.cpp
    check_mpi_problems.c
    tree_stats.cpp
    lingen_round_operand_size.cpp
    bw-common.c
    lingen_platform.cpp
    lingen_hints.cpp
    )
target_link_libraries(lingen_common_support ${gmp_libname} utils)
set_target_properties(lingen_common_support PROPERTIES
    COMPILE_FLAGS "--mpi ${OpenMP_CXX_FLAGS}"
    LINK_FLAGS "--mpi ${OpenMP_CXX_FLAGS}"
    )

foreach(v u64k1 ${BWC_GFP_ARITHMETIC_BACKENDS})
    set(lingen_${v}_sources
        lingen.cpp
        lingen_tuning.cpp
        lingen_tuning_cache.cpp
        lingen_expected_pi_length.cpp
        lingen_bmstatus.cpp
        )

    set(lingen_${v}_support_sources
            lingen_matpoly_ft.cpp
            )
    if (v STREQUAL "u64k1")
        set(lingen_${v}_sources
            ${lingen_${v}_sources}
            lingen_qcode_binary.cpp
            )
        set(lingen_${v}_support_sources
            ${lingen_${v}_support_sources}
            lingen_matpoly_binary.cpp
            )
    else()
        set(lingen_${v}_sources
            ${lingen_${v}_sources}
            lingen_tune_cutoffs.cpp
            lingen_qcode_prime.cpp
            )
        set(lingen_${v}_support_sources
            ${lingen_${v}_support_sources}
            mpfq/mpfq_${v}.c
            mpfq/mpfq_${v}.cpp
            # I don't have the courage to ditch lingen_polymat entirely
            # right now, because it's a bit tied inside
            # lingen_tune_cutoffs.cpp (which I seldom use, if ever). But
            # really, it should go. On the other hand, some of the
            # functionality in lingen_tune_cutoffs.cpp should probably be
            # made more of a first-class citizen
            lingen_polymat.cpp
            # lingen_bigpolymat.c # deleted 20150826
            lingen_matpoly.cpp
            lingen_bigmatpoly.cpp
            lingen_bigmatpoly_ft.cpp
            )
    endif()

    # make lingen_pz a default target, but not the others.
    if(v STREQUAL "pz")
        add_library(lingen_${v}_support ${lingen_${v}_support_sources})
        add_executable(lingen_${v} ${lingen_${v}_sources})
        install(TARGETS lingen_${v} RUNTIME DESTINATION ${LIBSUFFIX}/linalg/bwc)
    else()
        add_library(lingen_${v}_support EXCLUDE_FROM_ALL ${lingen_${v}_support_sources})
        add_executable(lingen_${v} EXCLUDE_FROM_ALL ${lingen_${v}_sources})
    endif()

    if (v STREQUAL "u64k1")
        add_dependencies(lingen_${v}_support gf2x-build)
        target_link_libraries(lingen_${v}_support gf2x-fft gf2x)
    else()
        target_link_libraries(lingen_${v}_support flint-fft)
    endif()

    target_link_libraries(lingen_${v}_support lingen_common_support)
    target_link_libraries(lingen_${v} lingen_${v}_support)

    # we also need link flags, and config/openmp.cmake won't provide them
    # separately.

    # set_source_files_properties(lingen.cpp PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
    set_target_properties(lingen_${v}_support PROPERTIES
        COMPILE_FLAGS "--mpi ${OpenMP_CXX_FLAGS} -DSELECT_MPFQ_LAYER_${v}"
        LINK_FLAGS "--mpi ${OpenMP_CXX_FLAGS}"
        )
    set_target_properties(lingen_${v} PROPERTIES
        COMPILE_FLAGS "--mpi ${OpenMP_CXX_FLAGS} -DSELECT_MPFQ_LAYER_${v}"
        LINK_FLAGS "--mpi ${OpenMP_CXX_FLAGS}"
        )
endforeach(v)


add_custom_target(bwc_matmul)
add_dependencies(bwc_matmul ${all_matmul_libs})

add_custom_target(bwc_full)
add_dependencies(bwc_full mf_bal random_matrix acollect lingen
    cleanup prep krylov dispatch mksol gather secure bwc_base
    bwc_matmul bwccheck)

# Those are used only for tests.
bwc_program_extra(build_matcache build_matcache.c)
target_link_libraries(build_matcache matmul_mf)
bwc_program_extra(bench_matcache bench_matcache.c worker-threads.c)
target_link_libraries(bench_matcache matmul_mf)


include(CMakeLists-nodist.txt OPTIONAL)
