
if (${WITH_MPI})
add_definitions(-DWITH_MPI)
endif()

add_executable(sqrt sqrt.cpp)
target_link_libraries(sqrt utils_with_io)
mark_targets_as_openmp(sqrt)
install(TARGETS sqrt RUNTIME DESTINATION ${LIBSUFFIX}/sqrt)

add_executable(mgmsqrt mgmsqrt.cpp)
target_link_libraries(mgmsqrt utils_with_io)
mark_targets_as_openmp(mgmsqrt)
install(TARGETS mgmsqrt RUNTIME DESTINATION ${LIBSUFFIX}/sqrt)


add_executable(crtalgsqrt ${EXCLUDE_FROM_RELEASE_BUILD} crtalgsqrt.cpp
    powers_of_p.cpp abfiles.cpp knapsack.c sqrt_wq.cpp
    sqrt_cachefiles.cpp)
target_link_libraries(crtalgsqrt utils)

mark_targets_as_mpi(crtalgsqrt)

include(CMakeLists-nodist.txt OPTIONAL)
