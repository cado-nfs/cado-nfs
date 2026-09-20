#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <exception>

#include "fmt/base.h"

#include "cado_main.hpp"

namespace cado {

int main_wrapper(int (*f)(int, char const **), int argc, char const * argv[])
{
    try {
        return f(argc, argv);
    } catch (std::exception const & e) {
        fmt::print(stderr, "{}\n", e.what());
        return EXIT_FAILURE;
    }
}

}
