
set(CMAKE_REQUIRED_LIBRARIES m)

include(${PROJECT_SOURCE_DIR}/config/search_for_function.cmake)
search_for_function(log2 HAVE_LOG2)
search_for_function(exp2 HAVE_EXP2)

search_for_function(cabsl HAVE_CABSL)

# clog is c99, but freebsd-11.1 (released in 2017) does not have it.
# https://www.freebsd.org/cgi/man.cgi?query=complex&manpath=FreeBSD+11.1-RELEASE
search_for_function(clog HAVE_CLOG)

set(math_libs ${CMAKE_REQUIRED_LIBRARIES_EXTRA})
