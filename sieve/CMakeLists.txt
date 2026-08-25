# Some compilation defaults specific to this directory

# The following is needed for some SSE stuff
# set(CMAKE_C_FLAGS   "-flax-vector-conversions ${CMAKE_C_FLAGS}")

include_directories(.)
#include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR}/ecm)

# ugly (copied from linalg/bwc/CMakeLists.txt...)
add_executable(makefb makefb.cpp)
cado_apply_pch(makefb)
target_link_libraries(makefb PRIVATE utils)

add_library(trialdiv OBJECT trialdiv.cpp)
cado_apply_pch(trialdiv)

add_library(las_core_0 OBJECT las-side-config.cpp)
target_link_libraries(las_core_0 PUBLIC utils)
cado_apply_pch(las_core_0)

add_executable(freerel freerel.cpp)
target_link_libraries(freerel PRIVATE las_core_0)
mark_targets_as_openmp(freerel)
cado_apply_pch(freerel)


# las_core_a has just very very basic types and code that can be used
# (and are actually used) by some unit tests.
add_library(las_core_a OBJECT las-qlattice.cpp)
target_link_libraries(las_core_a PUBLIC utils)
cado_apply_pch(las_core_a)

# las_core_b consists of all the object files that are necessary to build
# a las_info object. It is probably more than what would be strictly
# necessary, ideally. Everything extra is burden that we needlessly
# embark in dupsup.

# las_core_b0 is slighlty intermediary.
add_library(las_core_b0 OBJECT
        las-siever-config.cpp
        las-norms.cpp
        logapprox.cpp
        las-coordinates.cpp
        relation_cache.cpp
)
cado_apply_pch(las_core_b0)


# note that las_code_b does not always pull las_core_b0.
add_library(las_core_b OBJECT
        las-info.cpp
        las-special-q-task-collection.cpp
        las-special-q-task-tree.cpp
        las-memory.cpp
        las-parallel.cpp
        las-sieve-shared-data.cpp
        las-choose-sieve-area.cpp
        las-dlog-base.cpp
        las-cofactor.cpp
        las-threads.cpp
        las-auxiliary-data.cpp
        las-report-stats.cpp
        las-trialdiv.cpp
        las-bkmult.cpp
        las-todo-list.cpp
        las-galois.cpp
        las-globals.cpp
        las-threads-work-data.cpp
        las-dumpfile.cpp
        las-duplicate.cpp
        las-divide-primes.cpp
        las-sieve2357.cpp
)
cado_apply_pch(las_core_b)

add_library(las_core_c OBJECT las-arith.cpp)
target_link_libraries(las_core_c PUBLIC las_core_b)
cado_apply_pch(las_core_c)

set(las_multi_a_sources
        fb.cpp
        bucket.cpp
        las-config.cpp
        las-fill-in-buckets.cpp
        las-smallsieve.cpp
        siqs-smallsieve.cpp
        las-cofac-standalone.cpp
        las-detached-cofac.cpp
        las-process-bucket-region.cpp
        las-output.cpp
        las-multiobj-globals.cpp
        las-unsieve.cpp
)
if(HAVE_SSE2)
    list(APPEND las_multi_a_sources las-unsieve-sse2.cpp)
endif()
if(HAVE_AVX2)
    list(APPEND las_multi_a_sources las-unsieve-avx2.cpp)
endif()
add_library(las_multi_a_prod OBJECT ${las_multi_a_sources})
cado_apply_pch(las_multi_a_prod)
add_library(las_multi_a_debug OBJECT ${las_multi_a_sources})
cado_apply_pch(las_multi_a_debug)
set_property(TARGET las_multi_a_debug PROPERTY COMPILE_DEFINITIONS TRACE_K=1)
add_library(las_multi_a_descent OBJECT ${las_multi_a_sources})
cado_apply_pch(las_multi_a_descent)
set_property(TARGET las_multi_a_descent PROPERTY COMPILE_DEFINITIONS DLP_DESCENT SUPPORT_LARGE_Q)

add_library(siqs_multi_a_prod OBJECT ${las_multi_a_sources})
cado_apply_pch(siqs_multi_a_prod)
set_property(TARGET siqs_multi_a_prod PROPERTY COMPILE_DEFINITIONS SUPPORT_LARGE_Q=1 SIQS_SIEVE=1)
add_library(siqs_multi_a_debug OBJECT ${las_multi_a_sources})
cado_apply_pch(siqs_multi_a_debug)
set_property(TARGET siqs_multi_a_debug PROPERTY COMPILE_DEFINITIONS SUPPORT_LARGE_Q=1 SIQS_SIEVE=1 TRACE_K=1)

add_executable(las
    las-where-am-i-prod.cpp
    las.cpp
    )
target_link_libraries(las PRIVATE
    las_multi_a_prod las_core_c las_core_b las_core_b0 las_core_a las_core_0 trialdiv facul utils)
cado_apply_pch(las)

add_executable(las_tracek
    las-where-am-i-debug.cpp
    las.cpp
    )
target_link_libraries(las_tracek PRIVATE
    las_multi_a_debug las_core_c las_core_b las_core_b0 las_core_a las_core_0 trialdiv facul utils)
cado_apply_pch(las_tracek)
set_property(TARGET las_tracek PROPERTY COMPILE_DEFINITIONS TRACE_K=1)

add_executable(siqs las-where-am-i-prod.cpp las.cpp)
target_link_libraries(siqs PRIVATE
    siqs_multi_a_prod las_core_c las_core_b las_core_b0 las_core_a las_core_0 trialdiv facul utils)
cado_apply_pch(siqs)
set_property(TARGET siqs PROPERTY COMPILE_DEFINITIONS SUPPORT_LARGE_Q=1 SIQS_SIEVE=1)

add_executable(siqs_tracek las-where-am-i-debug.cpp las.cpp)
target_link_libraries(siqs_tracek PRIVATE
    siqs_multi_a_debug las_core_c las_core_b las_core_b0 las_core_a las_core_0 trialdiv facul utils)
cado_apply_pch(siqs_tracek)
set_property(TARGET siqs_tracek PROPERTY COMPILE_DEFINITIONS SUPPORT_LARGE_Q=1 SIQS_SIEVE=1 TRACE_K=1)

add_executable(las_descent las-where-am-i-prod.cpp las.cpp)
target_link_libraries(las_descent PRIVATE
    las_multi_a_descent las_core_c las_core_b las_core_b0 las_core_a las_core_0 trialdiv facul utils)
cado_apply_pch(las_descent)
set_property(TARGET las_descent PROPERTY COMPILE_DEFINITIONS DLP_DESCENT SUPPORT_LARGE_Q)

add_subdirectory(ecm)
add_subdirectory(strategies)

# pulling las_core_b is overkill, here
add_executable(dupsup dupsup.cpp las-where-am-i-prod.cpp)
cado_apply_pch(dupsup)
target_link_libraries(dupsup PRIVATE
    las_multi_a_descent las_core_b las_core_b0 las_core_a las_core_0 trialdiv facul utils)

add_executable(fake_rels EXCLUDE_FROM_ALL fake_rels.cpp)
cado_apply_pch(fake_rels)
target_link_libraries(fake_rels PRIVATE utils_with_io)

install(TARGETS makefb las las_descent freerel RUNTIME DESTINATION ${LIBSUFFIX}/sieve) 

include(CMakeLists-nodist.txt OPTIONAL)
