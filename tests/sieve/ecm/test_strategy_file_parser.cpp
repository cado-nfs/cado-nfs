#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>
#include <array>
#include "ecm/facul_strategies.hpp"

int main(int argc, char const * argv[])
{
    if (argc > 2)
        exit(EXIT_FAILURE);

    std::vector<unsigned long> const B { 1ul<<20, 1ul<<20 };
    std::vector<unsigned int> const lpb { 25, 25 };
    std::vector<unsigned int> const mfb { 50, 50 };

    FILE * f = stdin;
    if (argc == 2) {
        f = fopen(argv[1], "r");
        ASSERT_ALWAYS(f != NULL);
    }
    auto F = facul_strategies(B, lpb, mfb, true, f, 0);
    if (argc == 2) {
        fclose(f);
    }

    F.print(stdout);
}
