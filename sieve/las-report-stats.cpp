#include "cado.h"
#include <stdlib.h>
#include "las-report-stats.hpp"
#include "verbose.h"

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
    verbose_output_print(0, 2, "# survivors trial_divided_on_side[0]: %lu\n", S.trial_divided_on_side[0]);
    verbose_output_print(0, 2, "# survivors check_leftover_norm_on_side[0]: %lu (%.1f%%)\n", S.check_leftover_norm_on_side[0], 100 * ratio(S.check_leftover_norm_on_side[0], S.trial_divided_on_side[0]));
    verbose_output_print(0, 2, "# survivors trial_divided_on_side[1]: %lu\n", S.trial_divided_on_side[1]);
    verbose_output_print(0, 2, "# survivors check_leftover_norm_on_side[1]: %lu (%.1f%%)\n", S.check_leftover_norm_on_side[1], 100 * ratio(S.check_leftover_norm_on_side[1], S.trial_divided_on_side[1]));
    verbose_output_print(0, 2, "# survivors enter_cofactoring: %lu\n", S.enter_cofactoring);
    verbose_output_print(0, 2, "# survivors cofactored: %lu (%.1f%%)\n", S.cofactored, 100.0 * ratio(S.cofactored, S.enter_cofactoring));
    verbose_output_print(0, 2, "# survivors smooth: %lu\n", S.smooth);
}

