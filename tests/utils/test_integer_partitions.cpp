#include "cado.h" // IWYU pragma: keep

#include "fmt/base.h"

#include "integer_partitions.hpp"
#include "macros.h"

int main()
{
    const unsigned int a000041[] = {
        1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627, 792, 1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565, 5604, 6842, 8349, 10143, 12310, 14883, 17977, 21637, 26015, 31185, 37338, 44583, 53174, 63261, 75175, 89134, 105558, 124754, 147273, 173525
    };

    for(unsigned int i = 1 ; i < 20 ; i++) {
        unsigned int s = 0;
        for(auto const & p MAYBE_UNUSED : integer_partitions(i))
            s++;

        fmt::print("{} {} {}\n", i, a000041[i], s);
        ASSERT_ALWAYS(a000041[i] == s);
    }
    return 0;
}
