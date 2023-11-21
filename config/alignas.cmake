
# alignas is c++11, but g++ 4.7 does not implement it.
# There are also "CMake compile features" which test for things like the
# alignas keyword, but they require CMake 3.1.

include(CheckCXXSourceCompiles)

message(STATUS "Testing whether C++11 alignas() can be used")

set(check_alignas_code "
#include <stdint.h>
int main()
{
    alignas(8) uint32_t x = 0;
    return x;
}
")

set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES)

set(OLD_CMAKE_REQUIRED_QUIET "${CMAKE_REQUIRED_QUIET}")
set(CMAKE_REQUIRED_QUIET 1)
CHECK_CXX_SOURCE_COMPILES("${check_alignas_code}" HAVE_ALIGNAS)
set(CMAKE_REQUIRED_QUIET "${OLD_CMAKE_REQUIRED_QUIET}")

if(HAVE_ALIGNAS)
    message(STATUS "Testing whether C++11 alignas() can be used -- yes")
else()
    message(STATUS "Testing whether C++11 alignas() can be used -- no")
endif()

