#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "ecm/facul_strategies.hpp"

int main(int argc, char const * argv[])
{
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " nsides [file]" << "\n";
        exit(EXIT_FAILURE);
    }

    int parsed = std::stoi(std::string(argv[1]));
    unsigned int nsides = parsed > 0 ? parsed : 0U;

    if (nsides == 0) {
        std::cerr << "Error for nsides\n";
        exit(EXIT_FAILURE);
    }

    std::vector<unsigned long> const B(nsides, 1ul<<20);
    std::vector<unsigned int> const lpb(nsides, 25);
    std::vector<unsigned int> const mfb(nsides, 50);

    FILE * f = stdin; /* read from stdin if no file is provided */
    if (argc == 3) {
        f = fopen(argv[2], "r");
        ASSERT_ALWAYS(f != NULL);
    }
    auto F = facul_strategies(B, lpb, mfb, true, f, 0);
    if (argc == 3) {
        fclose(f);
    }

    F.print(stdout);
}
