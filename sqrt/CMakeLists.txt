
if (${WITH_MPI})
add_definitions(-DWITH_MPI)
endif (${WITH_MPI})

add_executable(sqrt sqrt.c)
target_link_libraries(sqrt utils m)

install(TARGETS sqrt RUNTIME DESTINATION bin/sqrt)

add_executable(checkdep ${EXCLUDE_FROM_RELEASE_BUILD} checkdep.c)
target_link_libraries(checkdep utils)

install(TARGETS sqrt RUNTIME DESTINATION bin/sqrt)

# crtalgsqrt requires both threads and mpi.
if(${WITH_PTHREADS})
add_definitions(-DWITH_PTHREADS)

if (HAVE_CABSL)
add_executable(crtalgsqrt ${EXCLUDE_FROM_RELEASE_BUILD} crtalgsqrt.c powers_of_p.cpp
    polyroots.c knapsack.c)
target_link_libraries(crtalgsqrt utils m ${pthread_libs})

set_target_properties(crtalgsqrt PROPERTIES
    COMPILE_FLAGS "--mpi"
    LINK_FLAGS "--mpi")
endif(HAVE_CABSL)
endif(${WITH_PTHREADS})

# add_executable(dispatch dispatch.cpp)
# 
# set_target_properties(dispatch PROPERTIES
#     COMPILE_FLAGS "--mpi"
#     LINK_FLAGS "--mpi")
# 

include(CMakeLists-nodist.txt OPTIONAL)
