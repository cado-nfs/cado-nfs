#define __STDCPP_MATH_SPEC_FUNCS__ 201003L
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1       /* for expint() */
#include <cmath>

double foo(size_t p0, size_t p1) {
    return std::expint(log(p1)) - std::expint(log(p0));
}

int main() {
    return foo(17, 42) >= 3 ? 0 : 1;
}

