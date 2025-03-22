#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include "decomp.hpp"
#include "utils_cxx.hpp"

std::ostream & operator<<(std::ostream & o, decomp const & D)
{
    return o << "[ " << join(D, " ") << " ]: " << D.nb_elem;
}
