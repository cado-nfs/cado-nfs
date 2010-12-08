
# Some compilation defaults specific to this directory
# The following flags are specific to gcc, and fix a regression in gcc 4.3.0

set_property(DIRECTORY APPEND PROPERTY COMPILE_FLAGS fno-ivopts -fno-regmove -fno-tree-sra -fno-tree-ter)


# CFLAGS=-O0 -g
# CFLAGS+=-DPROFILE

# if you compile with GCC 4.2 or later, you can also try -march=native
# instead of the default -mtune=generic (cf man gcc and
# http://gcc.gnu.org/ml/gcc-help/2008-07/msg00060.html)

# using -DSSE_NORM_INIT seems to be slower (about 3% with gcc 4.3.0)
# CFLAGS+=-DSSE_NORM_INIT
# -fkeep-inline-functions does not yield any speedup, is it useful?
# CFLAGS+=-fkeep-inline-functions

include_directories (ecm)
include_directories (${CMAKE_CURRENT_BINARY_DIR})
link_directories (${CMAKE_CURRENT_BINARY_DIR}/ecm)

add_executable (makefb makefb.c)
target_link_libraries (makefb utils m)

add_executable (las las.c fb.c trialdiv.c)
target_link_libraries (las facul utils m ${pthread_libs})

#add_executable (las-copper las-copper.c fb.c trialdiv.c)
#target_link_libraries (las-copper facul utils m pthread)

add_executable(trialdiv EXCLUDE_FROM_ALL trialdiv.c)
target_link_libraries (trialdiv utils)
set_property(TARGET trialdiv APPEND PROPERTY COMPILE_FLAGS -DTESTDRIVE)

add_executable(freerel freerel.c fb.c)
target_link_libraries(freerel utils m)

add_subdirectory(ecm)

install(TARGETS makefb las freerel RUNTIME DESTINATION bin/sieve)

include(CMakeLists-nodist.txt OPTIONAL)

