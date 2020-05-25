#include "cado.h" // IWYU pragma: keep
#include "las-report-stats.hpp"
#include "verbose.h"    // verbose_output_print
#include "macros.h"


/* declared in las.cpp */
extern int trialdiv_first_side;

static double ratio(double a, unsigned long b)
{
    return b ? a/b : 0;
}

void las_report::display_survivor_counters() const
{
    auto const& S(survivors);
    verbose_output_print(0, 2, "# survivors before_sieve: %lu\n", S.before_sieve);
    verbose_output_print(0, 2, "# survivors after_sieve: %lu (ratio %.2e)\n", S.after_sieve, ratio(S.after_sieve, S.before_sieve));
    verbose_output_print(0, 2, "# survivors not_both_even: %lu\n", S.not_both_even);
    verbose_output_print(0, 2, "# survivors not_both_multiples_of_p: %lu\n", S.not_both_multiples_of_p);
    unsigned long s = S.not_both_multiples_of_p;
    for(int pside = 0 ; pside < 2 ; pside++) {
        int side = trialdiv_first_side ^ pside;
        unsigned long sx = S.trial_divided_on_side[side];
        ASSERT_ALWAYS(s == sx || sx == 0);
        if (s && !sx) {
            verbose_output_print(0, 2, "# trial_division skipped on side %d: %lu survivors kept\n", side, s);
        } else {
            ASSERT_ALWAYS(s == sx);
            sx = S.check_leftover_norm_on_side[side];
            verbose_output_print(0, 2, "# survivors trial_divided_on_side[%d]: %lu\n", side, sx);
            verbose_output_print(0, 2, "# survivors check_leftover_norm_on_side[%d]: %lu (%.1f%%)\n", side, sx, 100 * ratio(sx, s));
            s = sx;
        }
    }
    ASSERT_ALWAYS(S.enter_cofactoring == s);
    verbose_output_print(0, 2, "# survivors enter_cofactoring: %lu\n", S.enter_cofactoring);
    verbose_output_print(0, 2, "# survivors cofactored: %lu (%.1f%%)\n", S.cofactored, 100.0 * ratio(S.cofactored, S.enter_cofactoring));
    verbose_output_print(0, 2, "# survivors smooth: %lu\n", S.smooth);
}

