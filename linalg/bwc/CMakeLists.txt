
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

# This is a default setting, but in reality we don't have much latitude
# with it. We *have* to compile all the u64k* which are in the set of
# sources, because of cross dependencies which are presently not
# correctly guarded.
set(BWC_GF2_ARITHMETIC_BACKENDS u64k1 u64k2 u64k3 u64k4
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
if(NOT BWC_GFP_ARITHMETIC_BACKENDS)
    message(STATUS "DL mode disabled for linear algebra")
else()
    if(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
        string_join(FOO ", " ${BWC_GFP_ARITHMETIC_BACKENDS})
        message(STATUS "GF(p) implementations enabled for DL linear algebra: ${FOO}")
        string_join(FOO ", " ${BWC_GFP_MATMUL_BACKENDS})
        message(STATUS "GF(p) SpMV backends for DL linear algebra: ${FOO}")
    else()
        # The code in linalg/bwc/mpfq/ is 64-bit only, and relies on GNU
        # CC idioms. As it is not used in production yet, we don't want
        # to expose it on 32-bit platforms.
        set(BWC_GFP_ARITHMETIC_BACKENDS CACHE STRING
                "which GF(p) backends to compile for DLOG mode" FORCE)
        message(STATUS "Warning: DL mode forcibly disabled for linear algebra (amd64-only)")
    endif()
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

# FIXME: does this have to be the same as BWC_GF2_ARITHMETIC_BACKENDS ?
# Is there any reason for having two variables ??
set(bwc_abases u64k1 u64k2 u64k3 u64k4)

set(bwc_mpfq_sources mpfq/mpfq_vbase.c)
foreach(v ${bwc_abases})
    set(bwc_mpfq_sources ${bwc_mpfq_sources} mpfq/mpfq_${v}.c mpfq/mpfq_${v}_t.c )
endforeach(v)

foreach(v ${BWC_GFP_ARITHMETIC_BACKENDS})
    set(bwc_mpfq_sources ${bwc_mpfq_sources} mpfq/mpfq_${v}.c mpfq/mpfq_${v}_t.c )
endforeach(v)

add_library(bwc_mpfq ${CADO_SHARED} ${bwc_mpfq_sources})

set(bwc_mpfq_flags)
foreach(v ${BWC_GFP_ARITHMETIC_BACKENDS})
    set(bwc_mpfq_flags "${bwc_mpfq_flags} -DCOMPILE_MPFQ_PRIME_FIELD_${v}")
endforeach(v)

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

string_join(FOO ", " ${BWC_GF2_ARITHMETIC_BACKENDS})
message(STATUS "GF(2) implementations enabled for linear algebra: ${FOO}")
string_join(FOO ", " ${BWC_GF2_MATMUL_BACKENDS})
message(STATUS "GF(2) SpMV backends for linear algebra: ${FOO}")



foreach(v ${BWC_GF2_ARITHMETIC_BACKENDS})
    foreach(w ${BWC_GF2_MATMUL_BACKENDS})
        CONFIGURE_MATMUL_LIB(${v} ${w})
    endforeach(w)
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

link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)

include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft)
include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft/mpfq)
add_executable(lingen lingen_binary.cpp tree_stats.cpp logline.c
    lingen_qcode.c lingen_qcode_do.cpp)
# This **MUST** match the #define in gf2x/Makefile.am
set_source_files_properties(lingen_qcode_do.cpp COMPILE_FLAGS
    "-DCANTOR_BASE_FIELD_SIZE=128")
set_source_files_properties(lingen_binary.cpp COMPILE_FLAGS
    "-DCANTOR_BASE_FIELD_SIZE=128")
# for debug only
# set_source_files_properties(lingen_binary.cpp COMPILE_FLAGS
# "-DCANTOR_BASE_FIELD_SIZE=128 -DLINGEN_BINARY_TRACE_MALLOCS=20")
install(TARGETS lingen RUNTIME DESTINATION ${LIBSUFFIX}/linalg/bwc)
target_link_libraries(lingen bwc_base utils gf2x)
add_dependencies(lingen gf2x-build)

# It's ugly. FindOpenMP.cmake doesn't know how to differentiate between
# required compile-time flags and required link-time flags...
set_target_properties(lingen PROPERTIES
    COMPILE_FLAGS "--mpi ${OpenMP_CXX_FLAGS}"
    LINK_FLAGS "--mpi ${OpenMP_CXX_FLAGS}")

add_executable(cleanup cleanup.c)
target_link_libraries(cleanup bitlinalg utils m)

install(TARGETS cleanup RUNTIME DESTINATION ${LIBSUFFIX}/linalg/bwc)

if(WITH_MPIR)
    # This is a convenience library which is used **ONLY** by one binary,
    # so there's absolutely no point in making it shared.
    file(GLOB flint_fft_files flint-fft/*.c)
    add_library(flint-fft ${flint_fft_files})
    add_executable(test-flint test-flint.c)
    target_link_libraries(test-flint flint-fft utils ${gmp_libname})
endif()

foreach(v ${BWC_GFP_ARITHMETIC_BACKENDS})
    set(plingen_${v}_sources
        bw-common.c
        plingen.cpp
        plingen-tuning.cpp
        tree_stats.cpp
        logline.c
        lingen-matpoly.c
        lingen-bigmatpoly.c
        lingen-polymat.c
        # lingen-bigpolymat.c # deleted 20150826
        mpfq/mpfq_${v}.c
        )
    if(WITH_MPIR)
        set(plingen_${v}_sources
            ${plingen_${v}_sources}
            lingen-matpoly-ft.c
            lingen-bigmatpoly-ft.c)
    endif()

    # make plingen_pz a default target, but not the others.
    if(v STREQUAL "pz")
        add_executable(plingen_${v} ${plingen_${v}_sources})
        install(TARGETS plingen_${v} RUNTIME DESTINATION ${LIBSUFFIX}/linalg/bwc)
    else()
        add_executable(plingen_${v} EXCLUDE_FROM_ALL ${plingen_${v}_sources})
    endif()

    if(WITH_MPIR)
        target_link_libraries(plingen_${v} flint-fft ${gmp_libname} utils)
    else()
        target_link_libraries(plingen_${v} ${gmp_libname} utils)
    endif()
    set_target_properties(plingen_${v} PROPERTIES
        COMPILE_FLAGS "--mpi -DSELECT_MPFQ_LAYER_${v}"
        LINK_FLAGS "--mpi"
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
