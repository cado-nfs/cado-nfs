#include "cado.h" // IWYU pragma: keep

#include "fmt/base.h"

#include "integer_partitions.hpp"
#include "macros.h"

int main()
{
    const unsigned int a000041[] = {
        1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627, 792, 1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565, 5604, 6842, 8349, 10143, 12310, 14883, 17977, 21637, 26015, 31185, 37338, 44583, 53174, 63261, 75175, 89134, 105558, 124754, 147273, 173525
    };

    for(unsigned int i = 0 ; i < 20 ; i++) {
        unsigned int s = 0;
        for(auto const & p MAYBE_UNUSED : integer_partitions(i))
            s++;

        fmt::print("{} {} {}\n", i, a000041[i], s);
        ASSERT_ALWAYS(a000041[i] == s);
    }

    const unsigned int a008483[] = {
        1, 0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 6, 9, 10, 13, 17, 21, 25, 33, 39, 49, 60, 73, 88, 110, 130, 158, 191, 230, 273, 331, 391, 468, 556, 660, 779, 927, 1087, 1284, 1510, 1775, 2075, 2438, 2842, 3323, 3872, 4510
    };
    for(unsigned int i = 0 ; i < 20 ; i++) {
        unsigned int s = 0;
        for(auto const & p MAYBE_UNUSED : integer_partitions(i, 3))
            s++;

        fmt::print("{} [3] {} {}\n", i, a008483[i], s);
        ASSERT_ALWAYS(a008483[i] == s);
    }

    const unsigned int a026812[] = {
        0, 0, 0, 0, 0, 0, 1, 1, 2, 3, 5, 7, 11, 14, 20, 26, 35, 44, 58, 71, 90, 110, 136, 163, 199, 235, 282, 331, 391, 454, 532, 612, 709, 811, 931, 1057, 1206, 1360, 1540, 1729,
    };
    for(unsigned int i = 0 ; i < 20 ; i++) {
        unsigned int s = 0;
        for(auto const & p MAYBE_UNUSED : integer_partitions_in_k_parts(i, 6))
            s++;

        fmt::print("{} [k=6] {} {}\n", i, a026812[i], s);
        ASSERT_ALWAYS(a026812[i] == s);
    }



    for(unsigned int min = 1 ; min < 8 ; min++) {
        for(unsigned int i = 0 ; i < 20 ; i++) {
            unsigned int s = 0;
            for(auto const & p MAYBE_UNUSED : integer_partitions(i, min))
                s++;

            fmt::print("{} [{}] {}\n", i, min, s);
        }
    }
    for(unsigned int k = 1 ; k < 8 ; k++) {
        for(unsigned int i = 0 ; i < 40 ; i++) {
            unsigned int s = 0;
            for(auto const & p MAYBE_UNUSED : integer_partitions_in_k_parts(i, k))
                s++;

            fmt::print("{} [k={}] {}\n", i, k, s);
        }
    }
    return 0;
}
