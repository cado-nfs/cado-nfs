/* Speed and correctness test of invmod_redc_32() */

#include "cado.h" // IWYU pragma: keep
#include <cinttypes>               // for PRId64, PRIu64
#include <cstdint>                 // for uint32_t, uint64_t
#include <cstdio>                  // for fprintf
#include <cstdlib>                 // for exit
#include <cmath>                   // for sqrt
#include "las-arith.hpp"           // for invmod_redc_32
#include "gcd.h"
#include "tests_common.h"
#include "misc.h"

// coverity[root_function]
int main(int argc, const char **argv) {
    unsigned long N = 1000000;
    int test_correctness = 1, test_timing = 1;

    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER | PARSE_CHECK | PARSE_TIME);
    tests_common_get_iter(&N);
    tests_common_get_check_and_time(&test_correctness, &test_timing);

    if (test_timing) {
        uint32_t fake_sum = 0, to_invert = 0;
        volatile uint32_t m = 4294967291;
        uint32_t invm = -invmod_po2(m);

        for (unsigned int i = 0; i < N; i++) {
            to_invert += 1073741827;
            fake_sum += invmod_redc_32(to_invert, m, invm);
        }
        volatile uint32_t fake_sum_vol = fake_sum;
        if (fake_sum_vol) {}
    }

    if (test_correctness) {

        /* a few quick checks */
        ASSERT_ALWAYS(invmod_redc_32(1, 2315500393, 1575713575) != 0);

        uint32_t to_invert = 0;
        uint32_t m = (uint32_t) u64_random(state) | 1;
        uint32_t invm = -invmod_po2(m);
        uint32_t to_add = gmp_urandomm_ui(state, m);

        while (gcd_uint64(to_add, m) > 1) {
            to_add++;
            if (to_add == m)
                to_add = 0;
        }

        unsigned long sqrtN = sqrt(N) + 1;
        for (unsigned int i = 0; i < sqrtN; i++) {
            for (unsigned int j = 0; j < sqrtN; j++) {
                unsigned long t;
                ularith_addmod_ul_ul(&t, to_invert, to_add, m);
                to_invert = t;
                uint32_t inverse = invmod_redc_32(to_invert, m, invm);
                if (inverse == 0 || inverse == UINT32_MAX) {
                    /* Compute GCD. Must be > 1, otherwise error */
                    if (gcd_uint64(to_invert, m) == 1) {
                        fprintf (stderr, "Error: invmod_redc_32(%" PRIu32 ", %" PRIu32
                                 ") = %" PRIu32 " but GCD is 1\n", to_invert, m, inverse);
                        exit(EXIT_FAILURE);
                    }
                } else {
                    uint32_t one = mulmodredc_u32<true>(to_invert, inverse, m, invm);
                    if (one != 1) {
                        fprintf (stderr, "Error: invmod_redc_32(%" PRIu32 ", %" PRIu32
                                 ") = %" PRIu32 " is wrong.\n", to_invert, m, inverse);
                        fprintf (stderr, "mulmodredc_u32(%" PRIu32 ", %" PRIu32 ", %"
                                 PRIu32 ", %" PRIu32 ") = %" PRIu32 ".\n",
                                 to_invert, inverse, m, invm, one);
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }
    exit(EXIT_SUCCESS);
}
