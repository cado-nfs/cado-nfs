
include_directories(. .. ../bblas)

# Use with care. Normally, binaries not needing MPI for real should not
# emit MPI symbols since they don't do anything depending on this config
# variable. However this mechanism is somewhat fragile.
if (${WITH_MPI})
add_definitions(-DWITH_MPI)
endif()

### bwc-specific configuration.

# Dynamically linkable bwc needs shared libs.
if(ENABLE_SHARED)
    set(BUILD_DYNAMICALLY_LINKABLE_BWC 1)
else()
    unset(BUILD_DYNAMICALLY_LINKABLE_BWC)
endif()

# having both shared and static is possible but messy. One needs two
# add_library calls, and then I'm not sure whether
# target_link_libraries picks the static or shared ones...

set(BWC_GF2_ARITHMETIC_BACKENDS b64 b128 b256 bz
    CACHE STRING
    "which low-level arithmetic backends to use for GF(2) linear algebra")
set(BWC_GF2_MATMUL_BACKENDS bucket basic sliced
    CACHE STRING
    "which low-level matrix multiplication backends to use for GF(2) linear algebra")
set(BWC_GFP_ARITHMETIC_BACKENDS p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 pz
    CACHE STRING
    "which GF(p) backends to compile for DLOG mode")
set(bmb basicp zone)
set(BWC_GFP_MATMUL_BACKENDS ${bmb}
    CACHE STRING
    "which low-level matrix multiplication backends to use for GF(2) linear algebra")

# In the fflas branch, this is used to plug in fflas.
set(BWC_EXTRA_BACKENDS "" CACHE STRING "list of pairs (arithmetic,matmul) backends for linear algebra")

### end bwc-specific configuration.

# These are just shorthands. Some targets also like to define compile and
# link flags for their own, so it's really messy, because I don't know
# how to have cmake *append* to the target's link flags.

# Note that whenever a target (at least very certainly an executable)
# links an mpi compiled library, it must be marked as mpi as well.

macro(mpi_add_library target)
    add_library(${target} ${ARGN})
    mark_targets_as_mpi(${target})
endmacro()
macro(mpi_add_executable target)
    add_executable(${target} ${ARGN})
    mark_targets_as_mpi(${target})
endmacro()

### bwc_base

set(bwc_parallelization_sources
    parallelizing_info.cpp
    )
if(HAVE_HWLOC)
    list(APPEND bwc_parallelization_sources
    cpubinding.cpp
    )
endif()

# The ${CADO_SHARED} flag is ignored here, really. As of current cmake
# (3.26) it seems that we can't give the SHARED flag to an OBJECT
# library. This results in object files in the object library which are,
# by default, compiled without -fPIC. Eventually, this is likely to yield
# issues like:
# /usr/bin/ld: CMakeFiles/bwc_parallelization.dir/parallelizing_info.cpp.o: relocation R_X86_64_PC32 against symbol `_Z17pi_go_helper_funcP14pi_go_helper_s' can not be used when making a shared object; recompile with -fPIC
#
# This is the exact same problem as https://stackoverflow.com/questions/50600708/combining-cmake-object-libraries-with-shared-libraries
# More info here as well: https://gitlab.kitware.com/cmake/cmake/-/issues/17356
#
# Therefore we get away with the exact same solution, although it feels a
# bit cumbersome.
mpi_add_library(bwc_parallelization ${CADO_SHARED} OBJECT ${bwc_parallelization_sources})
if(ENABLE_SHARED)
    set_property(TARGET bwc_parallelization PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

set(bwc_base_sources
    intersections.c
    xvectors.cpp
    bw-common.c
    async.cpp
    rolling.cpp
    matmul.cpp
    matmul_top.cpp
    matmul_top_vec.cpp
    matmul_top_comm.cpp
    random_matrix.cpp
    mmt_vector_pair.cpp
    mf_bal.cpp
    ../rowset_heap.cpp
    balancing.cpp
    balancing_workhorse.cpp
    arith-cross.cpp
    matrix_u32.cpp
    bwc_filenames.cpp
)

mpi_add_library(bwc_base ${CADO_SHARED} ${bwc_base_sources})

target_link_libraries(bwc_base bwc_parallelization)

target_link_libraries(bwc_base ${gmp_libname})
if (HAVE_CURL)
target_link_libraries(bwc_base curl)
endif()

target_link_libraries(bwc_base bblas)

if(HAVE_HWLOC)
target_link_libraries(bwc_base hwloc)
endif()


if(ENABLE_SHARED)
    install(TARGETS bwc_base DESTINATION ${LIBSUFFIX}/linalg/bwc)
endif()

### end bwc_base

mpi_add_executable(mf_scan2 ${EXCLUDE_FROM_RELEASE_BUILD} mf_scan2.cpp)
target_link_libraries(mf_scan2 bwc_base)


include_directories(${CADO_NFS_SOURCE_DIR}/linalg)

mpi_add_executable(mf_bal mf_bal_main.cpp)
target_link_libraries(mf_bal bwc_base m)
install(TARGETS mf_bal RUNTIME DESTINATION ${LIBSUFFIX}/linalg/bwc)


configure_file(bwc.pl ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl ESCAPE_QUOTES @ONLY)

add_executable(random_matrix ${EXCLUDE_FROM_RELEASE_BUILD}
    random_matrix.cpp)
set_property(TARGET random_matrix PROPERTY COMPILE_DEFINITIONS WANT_MAIN)
target_link_libraries(random_matrix utils)
mark_targets_as_mpi(random_matrix)

if(BUILD_DYNAMICALLY_LINKABLE_BWC)
    set(SOLIB_PREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
    set(SOLIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
    configure_file(solib-naming.h.in solib-naming.h ESCAPE_QUOTES @ONLY)
endif()

include_directories (${CMAKE_CURRENT_BINARY_DIR})


macro(parse_arith_layer v outlist)
    set(${outlist})
    list(APPEND ${outlist} -DARITH_LAYER_${v})
    list(APPEND ${outlist} -DARITH_LAYER=${v})
    string(REGEX REPLACE "^([bp])([0-9]*|z)$" "\\1;\\2" v_items ${v})
    list(GET v_items 0 vbp) # binary or prime
    if(vbp STREQUAL "b")
        list(GET v_items 1 vgroup)
        set(vwidth 0)
        if(vgroup STREQUAL "z")
            set(vgroup 0)
        endif()
        list(APPEND ${outlist} -DARITH_MOD2)
        list(APPEND ${outlist} -DARITH_SIMD_GROUPSIZE=${vgroup})
        set(cooked_macro "DO_${vbp}(${vgroup})")
        set(canonical_arith_layer "${v}")
    elseif(vbp STREQUAL "p")
        set(vgroup 1)
        list(GET v_items 1 vwidth)
        if(vwidth STREQUAL "z")
            set(vwidth 0)
        endif()
        list(APPEND ${outlist} -DARITH_MODP)
        list(APPEND ${outlist} -DARITH_SIMD_GROUPSIZE=1)
        list(APPEND ${outlist} -DARITH_PRIME_WIDTH=${vwidth})
        set(cooked_macro "DO_${vbp}(${vwidth})")
        set(canonical_arith_layer "${v}")
    elseif(${v} MATCHES "^u64k[0-9]*$")
        string(REGEX REPLACE "^u64k([0-9]*)$" "\\1" k ${v})
        math(EXPR vgroup "64*${k}")
        set(vwidth 0)
        set(vbp b)
        list(APPEND ${outlist} -DARITH_MOD2)
        list(APPEND ${outlist} -DARITH_SIMD_GROUPSIZE=${vgroup})
        set(cooked_macro "DO_${vbp}(${vgroup})")
        if(NOT warned_obsolete_${v})
            set(warned_obsolete_${v} 1)
            message(WARNING "${v} is an old backend name. Change your local.sh to use ${vbp}${vgroup} instead")
        endif()
        set(canonical_arith_layer "${vbp}${vgroup}")
    elseif(${v} STREQUAL "m128")
        set(vgroup 128)
        set(vbp b)
        list(APPEND ${outlist} -DARITH_MOD2)
        list(APPEND ${outlist} -DARITH_SIMD_GROUPSIZE=${vgroup})
        set(cooked_macro "DO_${vbp}(${vgroup})")
        if(NOT warned_obsolete_${v})
            set(warned_obsolete_${v} 1)
            message(WARNING "${v} is an old backend name. Change your local.sh to use ${vbp}${vgroup} instead")
        endif()
        set(canonical_arith_layer "${vbp}${vgroup}")
    elseif(${v} MATCHES "^p_[0-9]*$")
        string(REGEX REPLACE "^p_([0-9]*)$" "\\1" vwidth ${v})
        set(vbp p)
        list(APPEND ${outlist} -DARITH_MODP)
        list(APPEND ${outlist} -DARITH_SIMD_GROUPSIZE=1)
        list(APPEND ${outlist} -DARITH_PRIME_WIDTH=${vwidth})
        set(cooked_macro "DO_${vbp}(${vwidth})")
        if(NOT warned_obsolete_${v})
            set(warned_obsolete_${v} 1)
            message(WARNING "${v} is an old backend name. Change your local.sh to use ${vbp}${vwidth} instead")
        endif()
        set(canonical_arith_layer "${vbp}${vwidth}")
    elseif(${v} STREQUAL "pz")
        set(vwidth 0)
        set(vbp p)
        list(APPEND ${outlist} -DARITH_MODP)
        list(APPEND ${outlist} -DARITH_SIMD_GROUPSIZE=1)
        list(APPEND ${outlist} -DARITH_PRIME_WIDTH=${vwidth})
        set(cooked_macro "DO_${vbp}(${vwidth})")
        if(NOT warned_obsolete_${v})
            set(warned_obsolete_${v} 1)
            message(WARNING "${v} is an old backend name. Change your local.sh to use ${vbp}z instead")
        endif()
        set(canonical_arith_layer "${vbp}z")
    else()
        # Try to understand the old formats.
        message(FATAL_ERROR "unsupported arithmetic backend ${v}")
    endif()
endmacro()


add_library(bwc_arith ${CADO_SHARED} arith-generic.cpp)
target_link_libraries(bwc_arith utils)
mark_targets_as_mpi(bwc_arith)



set(internal_arithmetic_backends_gfp)
foreach(v pz ${BWC_GFP_ARITHMETIC_BACKENDS})
    set(dummy)
    parse_arith_layer(${v} dummy)
    list(APPEND internal_arithmetic_backends_gfp ${canonical_arith_layer})
endforeach()
list(REMOVE_DUPLICATES internal_arithmetic_backends_gfp)

set(internal_arithmetic_backends_gf2)
foreach(v b64 bz ${BWC_GF2_ARITHMETIC_BACKENDS})
    set(dummy)
    parse_arith_layer(${v} dummy)
    list(APPEND internal_arithmetic_backends_gf2 ${canonical_arith_layer})
endforeach()
list(REMOVE_DUPLICATES internal_arithmetic_backends_gf2)

set(internal_arithmetic_backends_all ${internal_arithmetic_backends_gf2} ${internal_arithmetic_backends_gfp})


set(all_arithmetic_gfp_libs)
set(all_arithmetic_gf2_libs)
set(COOKED_ARITHMETIC_BACKENDS)

if(BUILD_DYNAMICALLY_LINKABLE_BWC)
    foreach(v ${internal_arithmetic_backends_all})
        parse_arith_layer(${v} our_compile_flags)

        # To be honest, we don't really need the compile flags here. All
        # code is pulled in from arith-generic.cpp, which does not
        # include arith-hard.hpp

        add_library(arithmetic_${v} ${CADO_SHARED} arith-generic.cpp)
        target_link_libraries(arithmetic_${v} utils)
        mark_targets_as_mpi(arithmetic_${v})
        set_property(TARGET arithmetic_${v}
            APPEND PROPERTY COMPILE_OPTIONS ${our_compile_flags})
        install(TARGETS arithmetic_${v} DESTINATION ${LIBSUFFIX}/linalg/bwc)
    endforeach()

    foreach(v ${internal_arithmetic_backends_gf2})
        list(APPEND all_arithmetic_gf2_libs arithmetic_${v})
    endforeach()

    foreach(v ${internal_arithmetic_backends_gfp})
        list(APPEND all_arithmetic_gfp_libs arithmetic_${v})
    endforeach()

    # This is mandatory if we want bwc_arith to find our small .so
    # files for arithmetic.
    # Note that we do not need to add utils, it's added
    # automatically.
    set_target_properties(bwc_arith PROPERTIES BUILD_RPATH
        "${CMAKE_CURRENT_BINARY_DIR}")

    target_link_libraries(bwc_arith dl)

    install(TARGETS bwc_arith DESTINATION ${LIBSUFFIX}/linalg/bwc)
else()
    foreach(v ${internal_arithmetic_backends_all})
        parse_arith_layer(${v} dummy)
        set(COOKED_ARITHMETIC_BACKENDS ${COOKED_ARITHMETIC_BACKENDS} "${cooked_macro};")
    endforeach()
    # DO NOT USE our_compile_flags here, because those only make sense
    # when a single implementation is compiled.
endif()

### link dependencies

# we will probably need to link bwc_base with bwc_mpfq
target_link_libraries(bwc_base utils ${pthread_libs} bwc_arith)

# matmul_common is a dependency of all the impl-specific libs. When doing
# dynamic linking, there's no real need to make matmul_common *also*
# shared, as embarking the code within the impl-specific .so is fine.
# Unfortunately, we do want matmul_common to be compiled with -fPIC,
# then. Which is easiest to trigger by making it a shared object.
mpi_add_library(matmul_common ${CADO_SHARED} matmul-common.cpp)
target_link_libraries(matmul_common utils)
if(ENABLE_SHARED)
    install(TARGETS matmul_common DESTINATION ${LIBSUFFIX}/linalg/bwc)
endif()


set(matmul_basic_sources    matmul-basic.cpp)
set(matmul_basicp_sources   matmul-basicp.cpp)
set(matmul_zone_sources   matmul-zone.cpp)
set(matmul_sliced_sources   matmul-sliced.cpp)
set(matmul_bucket_sources   matmul-bucket.cpp)



macro(CONFIGURE_MATMUL_LIB v i)
        parse_arith_layer(${v} our_compile_flags)
        list(APPEND our_compile_flags -DMM_IMPL=${i})

        set(COOKED_BWC_BACKENDS ${COOKED_BWC_BACKENDS} "DO(${v}, ${i})")
        string(COMPARE EQUAL "${v}" "b64" doing_b64)
        string(COMPARE EQUAL "${i}" "bucket" doing_bucket)
        set(sources matmul_facade.cpp ${matmul_${i}_sources})
        if(doing_b64 AND doing_bucket)
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
        target_link_libraries(matmul_${v}_${i} matmul_common)#  bwc_mpfq)
        mark_targets_as_mpi(matmul_${v}_${i})
        set_property(TARGET matmul_${v}_${i}
            APPEND PROPERTY COMPILE_OPTIONS ${our_compile_flags})
        if(ENABLE_SHARED)
            install(TARGETS matmul_${v}_${i} DESTINATION ${LIBSUFFIX}/linalg/bwc)
        endif()
endmacro()

string_join(FOO ", " ${internal_arithmetic_backends_gf2})
message(STATUS "GF(2) SpMV backends for linear algebra: ${FOO}")
string_join(FOO ", " ${internal_arithmetic_backends_gfp})
message(STATUS "GF(p) SpMV backends for linear algebra: ${FOO}")



foreach(v ${internal_arithmetic_backends_gf2})
    foreach(w ${BWC_GF2_MATMUL_BACKENDS})
        CONFIGURE_MATMUL_LIB(${v} ${w})
        list(APPEND all_matmul_gf2_libs "matmul_${v}_${w}")
    endforeach(w)
endforeach(v)
foreach(v ${internal_arithmetic_backends_gfp})
    foreach(w ${BWC_GFP_MATMUL_BACKENDS})
        CONFIGURE_MATMUL_LIB(${v} ${w})
        list(APPEND all_matmul_gfp_libs "matmul_${v}_${w}")
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


# If we use static libraries, then unfortunately we _have_ to build all
# possible implementation for just any single run of bwc. Define
# ENABLE_SHARED=1 in local.sh to avoid this pitfall.
if(NOT BUILD_DYNAMICALLY_LINKABLE_BWC)
target_link_libraries(bwc_base ${all_matmul_gf2_libs} ${all_matmul_gfp_libs})
endif()

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
bwc_program(prep prep.cpp xdotprod.cpp)
target_link_libraries(prep bitlinalg)
target_link_libraries(prep bblas)
bwc_program(krylov krylov.cpp xdotprod.cpp)
bwc_program(bwccheck bwccheck.cpp xdotprod.cpp bwc_filenames.cpp)
bwc_program(dispatch dispatch.cpp xdotprod.cpp)
bwc_program(mksol mksol.cpp xdotprod.cpp)
bwc_program(det_from_lingen det_from_lingen.cpp)
bwc_program(gather gather.cpp)
bwc_program(secure secure.cpp)
bwc_program(blocklanczos blocklanczos.cpp blocklanczos_extraction.cpp)
target_link_libraries(blocklanczos bitlinalg)
target_link_libraries(blocklanczos bblas)

# bwc_program_base is good for programs that we don't want to see when we
# do "make install"
bwc_program_base(bench_cpu_bwc bench_cpu_bwc.cpp)
bwc_program_base(bench_comm_bwc bench_comm_bwc.cpp)
set_target_properties(bench_cpu_bwc PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(bench_comm_bwc PROPERTIES EXCLUDE_FROM_ALL TRUE)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl DESTINATION
    ${LIBSUFFIX}/linalg/bwc)

# Also lingen.

### This is the old binary lingen code. It is meant to go away.
if(NOT HAVE_GF2X)
    include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft)
    include_directories(${CADO_NFS_BINARY_DIR}/gf2x/fft)
    link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)
    link_directories (${CADO_NFS_BINARY_DIR}/gf2x/fft/.libs)
else()
    include_directories(${GF2X_INCDIR})
    link_directories(${GF2X_LIBDIR})
endif()
# This **MUST** match the #define in gf2x/Makefile.am
set_property(SOURCE lingen_qcode_binary.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS GF2X_CANTOR_BASE_FIELD_SIZE=128)

add_executable(cleanup cleanup.cpp)
target_link_libraries(cleanup bitlinalg utils)

install(TARGETS cleanup RUNTIME DESTINATION ${LIBSUFFIX}/linalg/bwc)

# This is a convenience library which is used **ONLY** by one binary,
# so there's absolutely no point in making it shared.
file(GLOB flint_fft_files flint-fft/*.c)
add_library(flint-fft OBJECT ${flint_fft_files})
if(ENABLE_SHARED)
    set_property(TARGET flint-fft PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
add_library(flint-fft-debug OBJECT EXCLUDE_FROM_ALL ${flint_fft_files})
mark_targets_as_mpi(flint-fft flint-fft-debug)
mark_targets_as_openmp(flint-fft flint-fft-debug)
set_property(TARGET flint-fft-debug
    APPEND PROPERTY COMPILE_DEFINITIONS DEBUG_FFT)

# These files are not dependent on the mpfq layer, therefore it's ok to
# put all of them in a common convenience library.
add_library(lingen_common_support 
    OBJECT
    lingen_memory_pool.cpp
    logline.cpp
    check_mpi_problems.c
    tree_stats.cpp
    lingen_round_operand_size.cpp
    bw-common.c
    lingen_platform.cpp
    lingen_hints.cpp
    lingen_call_companion.cpp
    )
if(ENABLE_SHARED)
    set_property(TARGET lingen_common_support PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
# not-for-OBJECT target_link_libraries(lingen_common_support ${gmp_libname} utils)
mark_targets_as_openmp(lingen_common_support)
mark_targets_as_mpi(lingen_common_support)

foreach(v ${internal_arithmetic_backends_gfp} b64)
    parse_arith_layer(${v} our_compile_flags)

    set(lingen_${v}_support_sources
        lingen_io_matpoly.cpp
        lingen_checkpoints.cpp
        lingen_average_matsize.cpp
        lingen_bmstatus.cpp
        lingen_bigmatpoly.cpp
        lingen_matpoly_ft.cpp
        )

    set(lingen_${v}_core_sources
        lingen.cpp
        lingen_expected_pi_length.cpp
        lingen_bigmatpoly_ft.cpp
        lingen_tuning.cpp
        lingen_tuning_cache.cpp
        )

    set(lingen_${v}_main_sources
        lingen_main.cpp
        lingen_io_wrappers.cpp
        )

    if (v STREQUAL "b64")
        list(APPEND lingen_${v}_core_sources
            lingen_qcode_binary.cpp
            )
        list(APPEND lingen_${v}_support_sources
            lingen_matpoly_binary.cpp
            )
        list(APPEND our_compile_flags -DLINGEN_BINARY)
    else()
        list(APPEND lingen_${v}_core_sources
            lingen_tune_cutoffs.cpp
            lingen_qcode_prime.cpp
            )

        list(APPEND lingen_${v}_support_sources
            ## mpfq/mpfq_${v}.c
            ## mpfq/mpfq_${v}.cpp
            # I don't have the courage to ditch lingen_polymat entirely
            # right now, because it's a bit tied inside
            # lingen_tune_cutoffs.cpp (which I seldom use, if ever). But
            # really, it should go. On the other hand, some of the
            # functionality in lingen_tune_cutoffs.cpp should probably be
            # made more of a first-class citizen
            lingen_polymat.cpp
            # lingen_bigpolymat.c # deleted 20150826
            lingen_matpoly.cpp
            )
    endif()

    add_executable(lingen_verify_checkpoints_${v}
        lingen_verify_checkpoints.cpp
        )
    mark_targets_as_mpi(lingen_verify_checkpoints_${v})
    mark_targets_as_openmp(lingen_verify_checkpoints_${v})
    set_property(TARGET lingen_verify_checkpoints_${v}
        APPEND PROPERTY COMPILE_OPTIONS ${our_compile_flags})
    target_link_libraries(lingen_verify_checkpoints_${v} utils)

    # lingen_bz could make sense for testing, but we certainly do not
    # want to  use it.
    if(v STREQUAL "pz" OR v STREQUAL "b64" OR v STREQUAL "bz")
        add_library(lingen_${v}_support ${CADO_SHARED} ${lingen_${v}_support_sources})
        add_library(lingen_${v}_core OBJECT ${lingen_${v}_core_sources})
        add_executable(lingen_${v} ${lingen_${v}_main_sources})
        install(TARGETS lingen_${v} RUNTIME DESTINATION ${LIBSUFFIX}/linalg/bwc)
    else()
        add_library(lingen_${v}_support ${CADO_SHARED} EXCLUDE_FROM_ALL ${lingen_${v}_support_sources})
        add_library(lingen_${v}_core OBJECT EXCLUDE_FROM_ALL ${lingen_${v}_core_sources})
        add_executable(lingen_${v} EXCLUDE_FROM_ALL ${lingen_${v}_main_sources})
    endif()
    install(TARGETS lingen_${v}_support DESTINATION ${LIBSUFFIX}/linalg/bwc)

    target_link_libraries(lingen_${v}_support lingen_common_support)
    target_link_libraries(lingen_${v}_core lingen_${v}_support)

    target_link_libraries(lingen_${v}_support 
        lingen_common_support utils)

    if (v STREQUAL "b64")
        if(NOT HAVE_GF2X)
            add_dependencies(lingen_${v}_support gf2x-build)
        endif()
        target_link_libraries(lingen_${v}_support gf2x bblas)
    else()
        # In fact the dependency could be put on lingen_${v} instead ;
        # however time_matpoly_parallel relies on _support to provide the
        # full thing.
        target_link_libraries(lingen_${v}_support flint-fft)
    endif()

    target_link_libraries(lingen_${v}
        lingen_${v}_core
        lingen_${v}_support
    )
    target_link_libraries(lingen_verify_checkpoints_${v}
        lingen_${v}_support
    )

    # we also need link flags, and config/openmp.cmake won't provide them
    # separately.

    mark_targets_as_mpi(lingen_${v}_support lingen_${v}_core lingen_${v})
    mark_targets_as_openmp(lingen_${v}_support lingen_${v}_core lingen_${v})
    set_property(TARGET lingen_${v}_support lingen_${v}_core lingen_${v}
        APPEND PROPERTY COMPILE_OPTIONS ${our_compile_flags})
endforeach(v)

add_custom_target(bwc_matmul_gf2)
add_custom_target(bwc_matmul_gfp)
if(BUILD_DYNAMICALLY_LINKABLE_BWC)
add_dependencies(bwc_matmul_gf2 ${all_arithmetic_gf2_libs})
add_dependencies(bwc_matmul_gf2 ${all_matmul_gf2_libs})
add_dependencies(bwc_matmul_gfp ${all_arithmetic_gfp_libs})
add_dependencies(bwc_matmul_gfp ${all_matmul_gfp_libs})
endif()

add_custom_target(bwc_full_common)
add_dependencies(bwc_full_common mf_bal random_matrix acollect
    prep krylov dispatch mksol gather secure bwc_base
    bwccheck)

add_custom_target(bwc_full_gf2)
add_dependencies(bwc_full_gf2 bwc_full_common bwc_matmul_gf2 lingen_b64 cleanup)
add_custom_target(bwc_full_gfp)
add_dependencies(bwc_full_gfp bwc_full_common bwc_matmul_gfp)
if(BWC_GFP_ARITHMETIC_BACKENDS)
    add_dependencies(bwc_full_gfp lingen_pz)
    if(pz IN_LIST BWC_GFP_ARITHMETIC_BACKENDS)
        add_dependencies(bwc_full_gfp lingen_pz)
    endif()
endif()

add_custom_target(bwc_full)
add_dependencies(bwc_full bwc_full_gf2 bwc_full_gfp)

# Those are used only for tests.
bwc_program_extra(build_matcache build_matcache.cpp)
target_link_libraries(build_matcache bwc_base)
bwc_program_extra(bench_matcache bench_matcache.cpp worker-threads.c)
target_link_libraries(bench_matcache bwc_base)

add_executable(bench_polmatmul EXCLUDE_FROM_ALL bench_polmatmul.cpp)
set_property(SOURCE bench_polmatmul.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS GF2X_CANTOR_BASE_FIELD_SIZE=128)
mark_targets_as_openmp(bench_polmatmul)

if(NOT HAVE_GF2X)
add_dependencies(bench_polmatmul gf2x-build)
endif()
target_link_libraries(bench_polmatmul utils gf2x)


include(CMakeLists-nodist.txt OPTIONAL)
