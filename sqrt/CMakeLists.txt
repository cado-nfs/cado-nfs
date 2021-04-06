
add_executable(sqrt sqrt.cpp)
target_link_libraries(sqrt utils_with_io ${pthread_libs})

CADO_MARK_OPENMP_CAPABLE_TARGET(sqrt)

add_executable(checkdep ${EXCLUDE_FROM_RELEASE_BUILD} checkdep.c)
target_link_libraries(checkdep utils)

install(TARGETS sqrt RUNTIME DESTINATION ${LIBSUFFIX}/sqrt)

# crtalgsqrt requires both threads and mpi.

if (HAVE_CABSL AND HAVE_CLOG)
add_executable(crtalgsqrt ${EXCLUDE_FROM_RELEASE_BUILD} crtalgsqrt.c powers_of_p.cpp knapsack.c)
target_link_libraries(crtalgsqrt utils ${pthread_libs})

CADO_MARK_MPI_CAPABLE_TARGET(crtalgsqrt)
endif()

CADO_MARK_TARGET_USING_OPENMP_CODE_INDIRECTLY(crtalgsqrt)

include(CMakeLists-nodist.txt OPTIONAL)
