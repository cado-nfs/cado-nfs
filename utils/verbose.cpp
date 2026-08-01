#include "cado.h" // IWYU pragma: keep

#include <cstdarg>
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <cstdlib>

#include <vector>
#include <utility>

#include <thread>
#include <mutex>
#include <condition_variable>

#include "verbose.hpp"
#include "portability.h" // strdup // IWYU pragma: keep
#include "macros.h"
#include "params.hpp"

#define G(X) CADO_VERBOSE_PRINT_ ## X
#define F(X) (UINT64_C(1) << G(X))

/* Mutex for verbose_output_*() functions */
// NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables)
static std::mutex io_mutex;
static std::condition_variable io_cond;
static bool batch_locked = false;
static std::thread::id batch_owner;

static uint64_t verbose_flag_word;
// NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables)

struct known_verbose_flag {
    unsigned int index;
    const char * name;
    int def;
};

static const std::vector<known_verbose_flag> verbose_flag_list = {
    {G(CMDLINE),                  "print-cmdline", 1 },
    {G(MODIFIED_FILES),           "print-modified-files", 1 },
    {G(COMPILATION_INFO),         "print-compilation-info", 1 },
    {G(BWC_DISPATCH_SLAVES),      "bwc-dispatch-slaves", 0 },
    {G(BWC_DISPATCH_MASTER),      "bwc-dispatch-master", 0 },
    {G(BWC_TIMING_GRIDS),         "bwc-timing-grids", 1 },
    {G(BWC_ITERATION_TIMINGS),    "bwc-iteration-timings", 1 },
    {G(BWC_CACHE_BUILD),          "bwc-cache-build", 0 },
    {G(BWC_DISPATCH_OUTER),       "bwc-dispatch-outer", 0 },
    {G(BWC_CPUBINDING),           "bwc-cpubinding", 1 },
    {G(BWC_CACHE_MAJOR_INFO),     "bwc-cache-major-info", 0 },
    {G(BWC_LOADING_MKSOL_FILES),  "bwc-loading-mksol-files", 1 },
};

const struct {
    const char * name;
    uint64_t mask;
} verbose_flag_groups[] = {
    { .name = "all-cmdline",
        .mask = F(CMDLINE) |
            F(MODIFIED_FILES) |
            F(COMPILATION_INFO) },
    { .name = "all-bwc-dispatch",
        .mask =  F(BWC_DISPATCH_SLAVES) |
            F(BWC_DISPATCH_MASTER) |
            F(BWC_DISPATCH_OUTER) |
            F(BWC_CACHE_BUILD) },
    { .name = "all-bwc-cache",
        .mask =       F(BWC_CACHE_MAJOR_INFO) |
            F(BWC_CACHE_BUILD) },
    { .name = "all-bwc-sub-timings", 
        .mask =       F(BWC_TIMING_GRIDS) |
            F(BWC_ITERATION_TIMINGS) },
};



/* This must be called in single-threaded context, preferably at program
 * start */
void verbose_interpret_parameters(cxx_param_list & pl)
{
    verbose_flag_word = ~0UL;

    /* mark these defaults. */
    for(auto const & [ index, name, def] : verbose_flag_list) {
        if (def == 0) {
            uint64_t mask = UINT64_C(1) << index;
            verbose_flag_word = verbose_flag_word & ~mask;
        }
    }

    const char * v = pl.lookup_old("verbose_flags");
    if (!v) return;

    char * w = strdup(v);
    char * p = w;
    char * q;
    for( ; *p != '\0' ; p = q) {
        q = strchr(p, ',');
        if (q) {
            *q++ = '\0';
        } else {
            q = p + strlen(p);
        }
        int enabled = 1;
        if (strncmp(p, "no-", 3) == 0) { enabled = 0; p += 3; }
        else if (strncmp(p, "no", 2) == 0) { enabled = 0; p += 2; }
        else if (*p == '^') { enabled = 0; p += 1; }
        else if (*p == '!') { enabled = 0; p += 1; }

        uint64_t mask = 0;
        for(auto const & [ index, name, def] : verbose_flag_list) {
            if (strcmp(p, name) == 0) {
                mask = UINT64_C(1) << index;
                break;
            }
        }
        for(size_t i = 0; i < sizeof(verbose_flag_groups) / sizeof(verbose_flag_groups[0]) ; i++) {
            if (strcmp(p, verbose_flag_groups[i].name) == 0) {
                mask = verbose_flag_groups[i].mask;
                break;
            }
        }
        if (!mask) {
            fprintf(stderr, "Verbose flag not recognized: %s\n", p);
            abort();
        }
        if (enabled) {
            verbose_flag_word |= mask;
        } else {
            verbose_flag_word &= ~mask;
        }
    }
    free(w);
}

void verbose_decl_usage(cxx_param_list & pl)
{
    pl.declare_usage("verbose_flags", "fine grained control on which messages get printed");
}

/* returns true if the following verbose flag is enabled */
int verbose_enabled(unsigned int flag) {
    return (verbose_flag_word & (UINT64_C(1) << flag)) != 0;
}

int verbose_vfprintf(FILE * f, int flag, const char * fmt, va_list ap)
{
    if (verbose_enabled(flag)) {
        return vfprintf(f, fmt, ap);
    }
    return 1;
}

int verbose_vprintf(int flag, const char * fmt, va_list ap)
{
    return verbose_vfprintf(stdout, flag, fmt, ap);
}
int verbose_fprintf(FILE * f, int flag, const char * fmt, ...)
{
    va_list ap;
    int rc;
    va_start(ap, fmt);
    rc = verbose_vfprintf(f, flag, fmt, ap);
    va_end(ap);
    return rc;
}
int verbose_printf(int flag, const char * fmt, ...)
{
    va_list ap;
    int rc;
    va_start(ap, fmt);
    rc = verbose_vprintf(flag, fmt, ap);
    va_end(ap);
    return rc;
}


/* Blocks until no other thread is in the monitor and
   no other thread holds a batch lock. A unique_lock object is returned,
   and its destruction in the parent scope will release the lock.
   */
static auto
monitor()
{
    std::unique_lock u(io_mutex);
    while (batch_locked && batch_owner != std::this_thread::get_id())
        io_cond.wait(u);
    return u;
}

int
verbose_output_start_batch()
{
    auto foo = monitor();
    batch_locked = true;
    batch_owner = std::this_thread::get_id();
    return 0;
}

int
verbose_output_end_batch()
{
    auto foo = monitor();
    ASSERT_ALWAYS(batch_locked);
    ASSERT_ALWAYS(batch_owner == std::this_thread::get_id());
    batch_locked = false;
    batch_owner = {};
    io_cond.notify_all();
    return 0;
}

/* Print formatted output to each output attached to this channel whose
   verbosity is at least the "verbosity" parameter.
   The "func" function is called with format string "fmt" and variable
   parameter list "va" for each output.
   If any output operation returns with an error, no further output is
   performed, and the error code of the failed operation is returned.
   Otherwise returns the return code of the last output operation.
   If no outputs are attached to this channel, returns 0. */
static int
vfprint_output(std::vector<std::pair<FILE *, int>> & output, const int verbosity,
               vfprintf_func_t func, const char * const fmt, va_list va)
{
    int rc = 0;
    /* For each output attached to this channel */
    for (auto [ F, v ] : output) {
        /* print string if output verbosity is at least "verbosity" */
        if (v >= verbosity) {
            va_list va_copied;
            va_copy(va_copied, va);
            rc = func(F, fmt, va_copied);
            va_end(va_copied);
            if (rc < 0)
                return rc;
        }
    }
    return rc;
}

/* we need this to avoid a static destruction order fiasco */
static std::vector<
        std::vector<
            std::pair<FILE *, int>
        >
        > & get_verbose_channel_outputs()
{
    static auto * res = new std::vector<
            std::vector<
                std::pair<FILE *, int>
            >
            >();
    return *res;
}

int
verbose_output_init(const size_t nr_channels)
{
    auto foo = monitor();
    get_verbose_channel_outputs().assign(nr_channels, {});
    return 1;
}

int
verbose_output_clear()
{
    auto foo = monitor();
    get_verbose_channel_outputs().clear();
    return 1;
}

int
verbose_output_add(const size_t channel, FILE * const out, const int verbose)
{
    auto foo = monitor();
    get_verbose_channel_outputs()[channel].emplace_back(out, verbose);
    return 0;
}

int
verbose_output_print(const size_t channel, const int verbose,
                     const char * const fmt, ...)
{
    va_list ap;
    int rc = 0;

    auto foo = monitor();
    va_start(ap, fmt);
    if (get_verbose_channel_outputs().empty()) {
        /* Default behaviour: print to stdout or stderr */
        ASSERT_ALWAYS(channel < 2);
        if (verbose <= 1) {
            FILE *out = (channel == 0) ? stdout : stderr;
            rc = vfprintf(out, fmt, ap);
        }
    } else {
        ASSERT_ALWAYS(channel < get_verbose_channel_outputs().size());
        rc = vfprint_output(get_verbose_channel_outputs()[channel], verbose, &vfprintf,
                            fmt, ap);
    }
    va_end(ap);
    return rc;
}

bool verbose_would_print(const size_t channel, const int verbosity)
{
    for (auto [ F, v ] : get_verbose_channel_outputs()[channel]) {
        /* print string if output verbosity is at least "verbosity" */
        if (v >= verbosity)
            return true;
    }
    return false;
}

void verbose_output_flush(const size_t channel, const int verbose)
{
    if (get_verbose_channel_outputs().empty()) {
        if (verbose > 1)
            return;
        FILE *out = (channel == 0) ? stdout : stderr;
        fflush(out);
    } else {
        for (auto [ F, v ] : get_verbose_channel_outputs()[channel]) {
            if (v >= verbose)
                fflush(F);
        }
    }
}

FILE *
verbose_output_get(const size_t channel, const int verbose, size_t index)
{
    auto foo = monitor();

    FILE *output = nullptr;
    if (get_verbose_channel_outputs().empty()) {
        /* Default behaviour: channel 0 has stdout, channel 1 has stderr,
           each with verbosity 1. */
        ASSERT_ALWAYS(channel < 2);
        if (verbose <= 1) {
            output = (channel == 0) ? stdout : stderr;
        }
    } else {
        /* Iterate through all the outputs for this channel */
        for (auto [ F, v ] : get_verbose_channel_outputs()[channel]) {
            /* Count those outputs that have verbosity at least "verbose" */
            if (v >= verbose && index-- == 0) {
                /* If that's the index-th output with enough verbosity,
                   return it. */
                return F;
            }
        }
    }

    return output;
}

int
verbose_output_vfprint(const size_t channel, const int verbose,
                       vfprintf_func_t func, const char * const fmt, ...)
{
    va_list ap;
    int rc = 0;

    auto foo = monitor();
    va_start(ap, fmt);
    if (get_verbose_channel_outputs().empty()) {
        /* Default behaviour: print to stdout or stderr */
        ASSERT_ALWAYS(channel < 2);
        if (verbose <= 1) {
            FILE *out = (channel == 0) ? stdout : stderr;
            rc = func(out, fmt, ap);
        }
    } else {
        rc = vfprint_output(get_verbose_channel_outputs()[channel], verbose, func, fmt,
                            ap);
    }
    va_end(ap);
    return rc;
}
