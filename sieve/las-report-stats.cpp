#include "cado.h" // IWYU pragma: keep
#include "las-report-stats.hpp"
#include "verbose.h"
#include "macros.h"
#include "utils_cxx.hpp"


/* declared in las.cpp */
extern int trialdiv_first_side;

void las_report::display_survivor_counters() const
{
    auto const& S(survivors);
    verbose_fmt_print(0, 2, "# survivors before_sieve: {}\n", S.before_sieve);
    verbose_fmt_print(0, 2, "# survivors after_sieve: {} (ratio {:.2e})\n", S.after_sieve, double_ratio(S.after_sieve, S.before_sieve));
    verbose_fmt_print(0, 2, "# survivors not_both_even: {}\n", S.not_both_even);
    verbose_fmt_print(0, 2, "# survivors not_both_multiples_of_p: {}\n", S.not_both_multiples_of_p);
    unsigned long s = S.not_both_multiples_of_p;
    for(int pside = 0 ; pside < 2 ; pside++) {
        int const side = trialdiv_first_side ^ pside;
        unsigned long sx = S.trial_divided_on_side[side];
        ASSERT_ALWAYS(s == sx || sx == 0);
        if (s && !sx) {
            verbose_fmt_print(0, 2, "# trial_division skipped on side {}: {} survivors kept\n", side, s);
        } else {
            ASSERT_ALWAYS(s == sx);
            sx = S.check_leftover_norm_on_side[side];
            verbose_fmt_print(0, 2, "# survivors trial_divided_on_side[{}]: {}\n", side, sx);
            verbose_fmt_print(0, 2, "# survivors check_leftover_norm_on_side[{}]: {} ({:.1f}%)\n", side, sx, 100 * double_ratio(sx, s));
            s = sx;
        }
    }
    ASSERT_ALWAYS(S.enter_cofactoring == s);

    verbose_fmt_print(0, 2, "# survivors enter_cofactoring: {}\n", S.enter_cofactoring);
    verbose_fmt_print(0, 2, "# survivors cofactored: {} ({:.1f}%)\n", S.cofactored, 100.0 * double_ratio(S.cofactored, S.enter_cofactoring));
    verbose_fmt_print(0, 2, "# survivors smooth: {}\n", S.smooth);
}

