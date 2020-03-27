
add_subdirectory(bblas)

include_directories(bwc) # for worker-threads.[ch]

include_directories(bblas)

# the bitlinalg library will almost surely go away at some point, or
# possibly be incorporated in bblas
add_library(bitlinalg blockmatrix.cpp)
target_link_libraries(bitlinalg bblas)
add_dependencies(bitlinalg gf2x-build)


add_executable(characters characters.cpp bwc/worker-threads.c)
target_link_libraries(characters bitlinalg utils_with_io pthread)
install(TARGETS characters RUNTIME DESTINATION ${LIBSUFFIX}/linalg)


include(CMakeLists-nodist.txt OPTIONAL)
