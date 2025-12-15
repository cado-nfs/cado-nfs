#include "cado.h" // IWYU pragma: keep

#include <cstdio>

#include "fmt/format.h"

#include "las-config.hpp"
#include "macros.h"
#include "verbose.h"
#include "utils_cxx.hpp"

int las_production_mode = 0;

int LOG_BUCKET_REGION = 16;
int LOG_BUCKET_REGION_step = 8;

int LOG_BUCKET_REGIONS[MAX_TOPLEVEL + 1];
size_t BUCKET_REGION;
size_t BUCKET_REGIONS[MAX_TOPLEVEL + 1];

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
    BUCKET_REGIONS[0] = 0;

    LOG_BUCKET_REGIONS[1] = LOG_BUCKET_REGION;
    BUCKET_REGIONS[1] = 1 << LOG_BUCKET_REGIONS[1];

#if MAX_TOPLEVEL >= 2
    LOG_BUCKET_REGIONS[2] = LOG_BUCKET_REGIONS[1] + LOG_BUCKET_REGION_step;
    BUCKET_REGIONS[2] = BUCKET_REGIONS[1] << LOG_BUCKET_REGION_step;
#endif
#if MAX_TOPLEVEL >= 3
    LOG_BUCKET_REGIONS[3] = LOG_BUCKET_REGIONS[2] + LOG_BUCKET_REGION_step;
    BUCKET_REGIONS[3] = BUCKET_REGIONS[2] << LOG_BUCKET_REGION_step;
#endif
}

void las_display_config_flags()
{
    std::vector<std::string> flags {
#ifdef BUCKET_SIEVE_POWERS
        "BUCKET_SIEVE_POWERS",
#endif
#ifdef PROFILE
        "PROFILE",
#endif
#ifdef WANT_ASSERT_EXPENSIVE
        "WANT_ASSERT_EXPENSIVE",
#endif
#ifdef TRACE_K
        "TRACE_K",
#endif
#ifdef TRACK_CODE_PATH
        "TRACK_CODE_PATH",
#endif
#ifdef SUPPORT_LARGE_Q
        "SUPPORT_LARGE_Q",
#endif
        fmt::format("LOGNORM_GUARD_BITS={:1.2f}",
                (double) LOGNORM_GUARD_BITS)
    };


    verbose_fmt_print(0, 1, "# las flags: {}\n", join(flags, " "));
} /* }}} */
