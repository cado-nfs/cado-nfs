#include "cado.h"
#include <cstdio>
#include <array>
#include "facul_strategies.hpp"

int main(int argc, char * argv[])
{
    if (argc != 2)
        exit(EXIT_FAILURE);
    std::array<unsigned long, 2> B { 1ul<<20, 1ul<<20 };
    std::array<unsigned int, 2> lpb { 25, 25 };
    std::array<unsigned int, 2> mfb { 50, 50 };
    FILE * f = fopen(argv[1], "r");
    ASSERT_ALWAYS(f != NULL);
    auto F = facul_strategies(B, lpb, mfb, true, f);
    fclose(f);

}
