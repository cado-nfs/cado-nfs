#include "cado.h" // IWYU pragma: keep

#include <stdio.h>

#include "las-config.h"
#include "macros.h"
#include "verbose.h" // verbose_output_print

int las_production_mode = 0;

int LOG_BUCKET_REGION = 16;
int LOG_BUCKET_REGIONS[FB_MAX_PARTS];
size_t BUCKET_REGION;
size_t BUCKET_REGIONS[FB_MAX_PARTS];

int NB_DEVIATIONS_BUCKET_REGIONS = 3;

void set_LOG_BUCKET_REGION()
{
    if (LOG_BUCKET_REGION > 16) {
        fprintf(stderr, "This binary only supports -B up to -B 16\n");
        ASSERT_ALWAYS(0);
        // we need to fix bucket_update_size_per_level<1>::type if we
        // want to explore larger B's.
    }

    BUCKET_REGION = ((size_t)1) << LOG_BUCKET_REGION;
    LOG_BUCKET_REGIONS[0] = -1;
    LOG_BUCKET_REGIONS[1] = LOG_BUCKET_REGION;
    LOG_BUCKET_REGIONS[2] = LOG_BUCKET_REGIONS[1] + LOG_NB_BUCKETS_2;
    LOG_BUCKET_REGIONS[3] = LOG_BUCKET_REGIONS[2] + LOG_NB_BUCKETS_3;
    BUCKET_REGIONS[0] = 0;
    BUCKET_REGIONS[1] = 1 << LOG_BUCKET_REGIONS[1];
    BUCKET_REGIONS[2] = 1 << LOG_BUCKET_REGIONS[2];
    BUCKET_REGIONS[3] = 1 << LOG_BUCKET_REGIONS[3];
}

void las_display_config_flags()
{
    verbose_output_print(0, 1, "# las flags:");
#ifdef SAFE_BUCKETS
    verbose_output_print(0, 1, " SAFE_BUCKETS");
#endif
#ifdef BUCKET_SIEVE_POWERS
    verbose_output_print(0, 1, " BUCKET_SIEVE_POWERS");
#endif
#ifdef PROFILE
    verbose_output_print(0, 1, " PROFILE");
#endif
#ifdef UNSIEVE_NOT_COPRIME
    verbose_output_print(0, 1, " UNSIEVE_NOT_COPRIME");
#endif
#ifdef WANT_ASSERT_EXPENSIVE
    verbose_output_print(0, 1, " WANT_ASSERT_EXPENSIVE");
#endif
#ifdef TRACE_K
    verbose_output_print(0, 1, " TRACE_K");
#endif
#ifdef TRACK_CODE_PATH
    verbose_output_print(0, 1, " TRACK_CODE_PATH");
#endif
#ifdef SUPPORT_LARGE_Q
    verbose_output_print(0, 1, " SUPPORT_LARGE_Q");
#endif
#ifdef SKIP_GCD3
    verbose_output_print(0, 1, " SKIP_GCD3");
#endif
    verbose_output_print(0, 1, " LOGNORM_GUARD_BITS=%1.2f",
                         (double)LOGNORM_GUARD_BITS);
    verbose_output_print(0, 1, "\n");
} /* }}} */
