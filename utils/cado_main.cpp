#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>

#include "fmt/base.h"

#include "cado_main.hpp"

namespace cado {

int main_wrapper(int (*f)(int, char const **), int argc, char const * argv[])
{
    try {
        return f(argc, argv);
    } catch (std::exception const & e) {
        /* So that the complaint comes after whatever the program printed
         * before it gave up, on both the C and the C++ stream. */
        fflush(stdout);
        std::cout.flush();
        /* Several of the messages we throw already end in a newline --
         * pl.fail() appends one, for instance -- and a good many begin
         * with "Error: " of their own. Print what we were given, and add
         * only what is missing. */
        std::string msg = e.what();
        if (!msg.ends_with('\n'))
            msg += '\n';
        fmt::print(stderr, "{}", msg);
        return EXIT_FAILURE;
    }
}

}
