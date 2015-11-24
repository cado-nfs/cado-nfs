
if (${WITH_MPI})
add_definitions(-DWITH_MPI)
endif()

add_executable(sqrt sqrt.c)
target_link_libraries(sqrt utils_with_io ${pthread_libs})

add_executable(checkdep ${EXCLUDE_FROM_RELEASE_BUILD} checkdep.c)
target_link_libraries(checkdep utils)

install(TARGETS sqrt RUNTIME DESTINATION ${LIBSUFFIX}/sqrt)

# crtalgsqrt requires both threads and mpi.
if(${WITH_PTHREADS})
add_definitions(-DWITH_PTHREADS)

if (HAVE_CABSL)
add_executable(crtalgsqrt ${EXCLUDE_FROM_RELEASE_BUILD} crtalgsqrt.c powers_of_p.cpp
    polyroots.c knapsack.c)
target_link_libraries(crtalgsqrt utils ${pthread_libs})

set_target_properties(crtalgsqrt PROPERTIES
    COMPILE_FLAGS "--mpi"
    LINK_FLAGS "--mpi")
endif()
endif()

include(CMakeLists-nodist.txt OPTIONAL)
