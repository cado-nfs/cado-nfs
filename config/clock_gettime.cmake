include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES)
set(OLD_CMAKE_REQUIRED_QUIET "${CMAKE_REQUIRED_QUIET}")
set(CMAKE_REQUIRED_QUIET 1)
foreach(clock CLOCK_MONOTONIC CLOCK_MONOTONIC_RAW
        CLOCK_THREAD_CPUTIME_ID)
    message(STATUS "Testing whether ${clock} can be used")
    CHECK_CXX_SOURCE_COMPILES(
"#define _POSIX_C_SOURCE 200112L
#include <time.h>
int main ()
{
  struct timespec ts[1];
  clock_gettime (${clock}, ts);
  return 0;
}" HAVE_${clock})
    if(HAVE_${clock})
    message(STATUS "Testing whether ${clock} can be used -- Success")
    else()
    message(STATUS "Testing whether ${clock} can be used -- Failed")
    endif()
endforeach()

set(CMAKE_REQUIRED_QUIET "${OLD_CMAKE_REQUIRED_QUIET}")
