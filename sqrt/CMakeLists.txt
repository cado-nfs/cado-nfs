
if (${WITH_MPI})
add_definitions(-DWITH_MPI)
endif()

add_executable(sqrt sqrt.cpp)
target_link_libraries(sqrt utils_with_io ${pthread_libs})
mark_targets_as_openmp(sqrt)

add_executable(checkdep ${EXCLUDE_FROM_RELEASE_BUILD} checkdep.c)
target_link_libraries(checkdep utils)

install(TARGETS sqrt RUNTIME DESTINATION ${LIBSUFFIX}/sqrt)

if (HAVE_CABSL AND HAVE_CLOG)
add_executable(crtalgsqrt ${EXCLUDE_FROM_RELEASE_BUILD} crtalgsqrt.c powers_of_p.cpp knapsack.c)
target_link_libraries(crtalgsqrt utils ${pthread_libs})

mark_targets_as_mpi(crtalgsqrt)
endif()

include(CMakeLists-nodist.txt OPTIONAL)
