# We need to know what exact type uint64_t is typedef'd to. The reason is
# that some templates would otherwise resolve ambiguously if we define
# both.

include(CheckCXXSourceCompiles)

macro(testcode_uint64_is type)
    set(test_code "
#include <type_traits>
#include <stdint.h>
#include <stdio.h>
int main()
{
    static_assert(std::is_same<uint64_t, ${type}>::value, \"not this type\");
    return 0;
}
"
)
endmacro()

set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES)
testcode_uint64_is("unsigned long")
CHECK_CXX_SOURCE_COMPILES("${test_code}" UINT64_T_IS_UNSIGNED_LONG)
testcode_uint64_is("unsigned long long")
CHECK_CXX_SOURCE_COMPILES("${test_code}" UINT64_T_IS_UNSIGNED_LONG_LONG)

