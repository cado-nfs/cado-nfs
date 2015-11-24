
include_directories(bwc) # for worker-threads.[ch] and matops

include_directories(${PROJECT_SOURCE_DIR}/gf2x/fft)
include_directories(${PROJECT_SOURCE_DIR}/gf2x/fft/mpfq)
add_library(bitlinalg blockmatrix.c matops.cpp gauss.c)
add_dependencies(bitlinalg gf2x-build)


add_executable(characters characters.c bwc/worker-threads.c)
target_link_libraries(characters bitlinalg utils_with_io pthread)
install(TARGETS characters RUNTIME DESTINATION ${LIBSUFFIX}/linalg)


include(CMakeLists-nodist.txt OPTIONAL)
