/* gen_decomp mfb fbb computes approximations of all decompositions of
 * mfb-bit integers into primes larger than fbb. The result is given for both methods in decomp.cpp, together with the ratio of what they predict.
 * Example:

$ test_gen_decomp 60 800000
20 40 3.113e+14 3.070e+14 0.986
21 39 4.581e+14 4.564e+14 0.996
22 38 4.485e+14 4.469e+14 0.996
23 37 4.405e+14 4.389e+14 0.996
24 36 4.337e+14 4.321e+14 0.996
25 35 4.282e+14 4.266e+14 0.996
26 34 4.237e+14 4.221e+14 0.996
27 33 4.203e+14 4.187e+14 0.996
28 32 4.179e+14 4.163e+14 0.996
29 31 4.165e+14 4.149e+14 0.996
30 30 2.080e+14 2.072e+14 0.996
20 41 1.007e+14 9.742e+13 0.967
21 40 5.827e+14 5.649e+14 0.970
22 39 5.699e+14 5.598e+14 0.982
23 38 5.589e+14 5.491e+14 0.982
24 37 5.496e+14 5.400e+14 0.983
25 36 5.420e+14 5.325e+14 0.983
26 35 5.358e+14 5.264e+14 0.982
27 34 5.309e+14 5.216e+14 0.982
28 33 5.273e+14 5.180e+14 0.982
29 32 5.249e+14 5.157e+14 0.982
30 31 5.237e+14 5.145e+14 0.982
20 20 20 9.775e+11 9.474e+11 0.969
20 20 21 3.682e+12 3.535e+12 0.960
20 21 21 6.840e+11 5.800e+11 0.848

*/

#include "cado.h" // IWYU pragma: keep

#include <climits>
#include <cstdlib>

#include "gen_decomp.hpp"
#include "macros.h"

int main(int argc, char const * argv[])
{
    ASSERT_ALWAYS(argc == 3);
    char * p;
    long mfb = strtol(argv[1], &p, 0);
    ASSERT_ALWAYS(*p == '\0');
    ASSERT_ALWAYS(mfb <= INT_MAX);
    unsigned long lim = strtoul(argv[2], &p, 0);
    ASSERT_ALWAYS(*p == '\0');

    generate_all_decomp_compare((int)mfb, lim);

    return EXIT_SUCCESS;
}
