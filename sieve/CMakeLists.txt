# Some compilation defaults specific to this directory

# The following is needed for some SSE stuff
# set (CMAKE_C_FLAGS   "-flax-vector-conversions ${CMAKE_C_FLAGS}")

include_directories (.)
#include_directories (${CMAKE_CURRENT_BINARY_DIR})
link_directories (${CMAKE_CURRENT_BINARY_DIR}/ecm)

# ugly (copied from linalg/bwc/CMakeLists.txt...)
add_executable (makefb makefb.c)
target_link_libraries (makefb utils)

add_library (trialdiv OBJECT trialdiv.cpp)

add_library(las_core_0 OBJECT
        las-side-config.cpp)

add_executable(freerel freerel.cpp $<TARGET_OBJECTS:las_core_0>)
target_link_libraries(freerel utils)

# las_core_a has just very very basic types and code that can be used
# (and are actually used) by some unit tests.
add_library(las_core_a OBJECT
        las-qlattice.cpp
        )

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
)

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

add_library(las_core_c OBJECT
        las-arith.cpp
)

set(las_multi_a_sources
        fb.cpp
	bucket.cpp
        las-config.cpp
	las-fill-in-buckets.cpp
	las-smallsieve.cpp
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
add_library(las_multi_a_prod OBJECT ${las_multi_a_sources})
add_library(las_multi_a_debug OBJECT ${las_multi_a_sources})
set_property(TARGET las_multi_a_debug PROPERTY COMPILE_DEFINITIONS TRACE_K=1)
add_library(las_multi_a_descent OBJECT ${las_multi_a_sources})
set_property(TARGET las_multi_a_descent PROPERTY COMPILE_DEFINITIONS DLP_DESCENT SUPPORT_LARGE_Q)

add_executable(las
    las-where-am-i-prod.cpp
    las.cpp
    $<TARGET_OBJECTS:las_multi_a_prod>
    $<TARGET_OBJECTS:las_core_c>
    $<TARGET_OBJECTS:las_core_b>
    $<TARGET_OBJECTS:las_core_b0>
    $<TARGET_OBJECTS:las_core_a>
    $<TARGET_OBJECTS:las_core_0>
    $<TARGET_OBJECTS:trialdiv>
    )
target_link_libraries (las facul utils)

add_executable(las_tracek
    las-where-am-i-debug.cpp
    las.cpp
    $<TARGET_OBJECTS:las_multi_a_debug>
    $<TARGET_OBJECTS:las_core_c>
    $<TARGET_OBJECTS:las_core_b>
    $<TARGET_OBJECTS:las_core_b0>
    $<TARGET_OBJECTS:las_core_a>
    $<TARGET_OBJECTS:las_core_0>
    $<TARGET_OBJECTS:trialdiv>
    )
target_link_libraries (las_tracek facul utils)
set_property(TARGET las_tracek PROPERTY COMPILE_DEFINITIONS TRACE_K=1)

add_executable(las_descent
    las-where-am-i-prod.cpp
    las.cpp
    $<TARGET_OBJECTS:las_multi_a_descent>
    $<TARGET_OBJECTS:las_core_c>
    $<TARGET_OBJECTS:las_core_b>
    $<TARGET_OBJECTS:las_core_b0>
    $<TARGET_OBJECTS:las_core_a>
    $<TARGET_OBJECTS:las_core_0>
    $<TARGET_OBJECTS:trialdiv>
    )
target_link_libraries (las_descent facul utils)
set_property(TARGET las_descent PROPERTY COMPILE_DEFINITIONS DLP_DESCENT SUPPORT_LARGE_Q)

add_subdirectory(ecm)
add_subdirectory(strategies)

add_executable(dupsup dupsup.cpp
    las-where-am-i-prod.cpp
    $<TARGET_OBJECTS:las_multi_a_descent>
    $<TARGET_OBJECTS:las_core_b>        # overkill
    $<TARGET_OBJECTS:las_core_b0>
    $<TARGET_OBJECTS:las_core_a>
    $<TARGET_OBJECTS:las_core_0>
    $<TARGET_OBJECTS:trialdiv>
    )
target_link_libraries(dupsup facul utils)

add_executable(fake_rels EXCLUDE_FROM_ALL fake_rels.cpp)
target_link_libraries(fake_rels utils_with_io)

install(TARGETS makefb las las_descent freerel RUNTIME DESTINATION ${LIBSUFFIX}/sieve) 

include(CMakeLists-nodist.txt OPTIONAL)
