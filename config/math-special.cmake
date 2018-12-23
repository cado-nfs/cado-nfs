include(CheckCXXSourceCompiles)

# I was a bit enthusiastic. These features are only in C++17 mainline,
# and optionally in boost.math.
#
# http://en.cppreference.com/w/cpp/experimental/special_math

set(check_math_special_functions "
#define __STDCPP_MATH_SPEC_FUNCS__ 201003L
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1       /* for expint() */
#include <cmath>

double foo(size_t p0, size_t p1) {
    return std::expint(log(p1)) - std::expint(log(p0));
}

int main() {
    return foo(17, 42) >= 3 ? 0 : 1;
}

")

set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES)
CHECK_CXX_SOURCE_COMPILES("${check_math_special_functions}" HAVE_STDCPP_MATH_SPEC_FUNCS)

