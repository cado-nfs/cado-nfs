
include_directories(..)

# Use with care. Normally, binaries not needed MPI for real should not
# emit MPI symbols since they don't do anything depending on this config
# variable. However this mechanism is somewhat fragile.
if (${WITH_MPI})
add_definitions(-DWITH_MPI)
endif (${WITH_MPI})

add_library(mf mf.c)
target_link_libraries(mf balancing utils)
add_executable(mf_scan ${EXCLUDE_FROM_RELEASE_BUILD} mf_scan.c)
target_link_libraries(mf_scan mf)

include_directories(${CADO_NFS_SOURCE_DIR}/linalg)
add_executable(mf_bal mf_bal.c ../rowset_heap.cpp)
target_link_libraries(mf_bal mf m)


configure_file(bwc.pl ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl ESCAPE_QUOTES @ONLY)

# test_hello does not compile. harmless
#add_executable(test_hello ${EXCLUDE_FROM_RELEASE_BUILD} test_hello.c parallelizing_info.c)
#target_link_libraries(test_hello utils)
#set_target_properties(test_hello PROPERTIES
#    COMPILE_FLAGS "--mpi"
#    LINK_FLAGS "--mpi")

add_executable(random ${EXCLUDE_FROM_RELEASE_BUILD} random.c)

# balancing.c itself is in fact standalone, and lightweight. We prefer to
# avoid having the whole thing entangled with too many dependencies all
# around, since all that split.c actually needs is very, very little.
add_executable(split split.c bw-common.c balancing.c)
#add_executable(shell
#    shell.c parallelizing_info.c bw-common.c bw-common-mpi.c)
#set_target_properties(shell
#    PROPERTIES COMPILE_FLAGS "--mpi" LINK_FLAGS "--mpi")

add_executable(acollect acollect.c bw-common.c)

install(TARGETS split acollect RUNTIME DESTINATION bin/linalg/bwc)

target_link_libraries(random utils m)
target_link_libraries(split utils)
# target_link_libraries(shell utils ${pthread_libs})
target_link_libraries(acollect utils)

set(balancing_sources
    balancing.c
    balancing_workhorse.c
    balancing_file_source.c
    balancing_mpi_source.c)

if (HAVE_CURL)
set(balancing_sources ${balancing_sources} balancing_rollbuf.c balancing_curl_source.c)
endif (HAVE_CURL)

add_library(balancing ${balancing_sources})
set_target_properties(balancing PROPERTIES
        COMPILE_FLAGS "--mpi"
        LINK_FLAGS "--mpi")

if (HAVE_CURL)
target_link_libraries(balancing curl)
endif (HAVE_CURL)

set(bwc_base_sources
    parallelizing_info.c
    intersections.cpp
    xvectors.c
    bw-common.c
    bw-common-mpi.c
    rusage.c
    async.c
    debug.c
    rolling.c
    worker-threads.c)

add_library(bwc_base ${bwc_base_sources})

target_link_libraries(bwc_base balancing utils ${pthread_libs})
set_target_properties(bwc_base PROPERTIES
        COMPILE_FLAGS "--mpi"
        LINK_FLAGS "--mpi")

add_executable(mf_twistvec ${EXCLUDE_FROM_RELEASE_BUILD} mf_twistvec.c)
target_link_libraries(mf_twistvec bwc_base)
# install(TARGETS mf_twistvec RUNTIME DESTINATION bin/linalg/bwc)

add_executable(mf_untwistvec mf_untwistvec.c)
target_link_libraries(mf_untwistvec bwc_base)
install(TARGETS mf_untwistvec RUNTIME DESTINATION bin/linalg/bwc)

set(all_variants)
macro(configure_bwc_variant v)
    set(${v}_matmul_sources 
        matmul-common.c matmul-basic.c matmul-sliced.cpp
        matmul-threaded.c matmul-bucket.cpp matmul.c)

    string(COMPARE EQUAL "${v}" "u64" doing_u64)
    if(doing_u64)
        if(HAVE_GCC_STYLE_AMD64_ASM)
            set(${v}_matmul_sources ${${v}_matmul_sources}
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
        endif(HAVE_GCC_STYLE_AMD64_ASM)
    endif(doing_u64)

    add_library(${v}_matmul ${${v}_matmul_sources})
    set_target_properties(${v}_matmul PROPERTIES
            COMPILE_FLAGS "--mpi -DSELECT_ABASE_${v}"
            LINK_FLAGS "--mpi")
    target_link_libraries(${v}_matmul bwc_base)
    target_link_libraries(${v}_matmul mat_io)


    add_library(${v}_bwc abase-${v}.c matmul_top.c)
    set_target_properties(${v}_bwc PROPERTIES
            COMPILE_FLAGS "--mpi -DSELECT_ABASE_${v}"
            LINK_FLAGS "--mpi")
    target_link_libraries(${v}_bwc bwc_base ${v}_matmul utils pthread)

    set(all_variants ${all_variants} ${v})
endmacro(configure_bwc_variant)

configure_bwc_variant(u64)
# configure_bwc_variant(u64k)
configure_bwc_variant(u64n)
# configure_bwc_variant(u128)

set(build_objs build.c)
set(bench_objs bench.c)
set(prep_objs prep.c xymats.c ../gauss.c)
set_source_files_properties (../gauss.c COMPILE_FLAGS "-DMULTI_ROW=3 -DNO_MAIN")
set(krylov_objs krylov.c xymats.c xdotprod.c)
set(dispatch_objs dispatch.c xymats.c xdotprod.c)
set(mksol_objs mksol.c xymats.c xdotprod.c)
set(gather_objs gather.c)
set(secure_objs secure.c)

macro(configure_program prg)
    foreach(v ${ARGN})
        add_executable(${v}_${prg} ${${prg}_objs})
        install(TARGETS ${v}_${prg} RUNTIME DESTINATION bin/linalg/bwc)
        target_link_libraries(${v}_${prg} ${v}_bwc)
        set_target_properties(${v}_${prg} PROPERTIES
                LINKER_LANGUAGE CXX
                COMPILE_FLAGS "-DSELECT_ABASE_${v} --mpi"
                LINK_FLAGS "--mpi")
    endforeach(v)
endmacro(configure_program)

set(used_variants       u64)
#set(used_variants ${all_variants})

configure_program(krylov ${used_variants})
configure_program(dispatch ${used_variants})
configure_program(mksol  ${used_variants})
configure_program(build  ${used_variants})
configure_program(bench  ${used_variants})
configure_program(prep   u64n)
configure_program(gather u64n)
configure_program(secure u64)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl DESTINATION bin/linalg/bwc)

# Also lingen.

link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)
link_directories (${CADO_NFS_BINARY_DIR}/gf2x/gf2x-fft/.libs)

include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/gf2x-fft)
add_executable(lingen lingen_binary.cpp)
install(TARGETS lingen RUNTIME DESTINATION bin/linalg/bwc)
target_link_libraries(lingen bwc_base utils gf2x-fft gf2x cantor128)
add_dependencies(lingen gf2x-build)

add_executable(short_matmul ${EXCLUDE_FROM_RELEASE_BUILD} short_matmul.c)

add_executable(cleanup cleanup.c ../blockmatrix.c ../bit_matrices.c
    ../gauss.c)
target_link_libraries(cleanup utils m)

include(CMakeLists-nodist.txt OPTIONAL)
