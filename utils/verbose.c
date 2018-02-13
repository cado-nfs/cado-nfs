#include "cado.h"
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <pthread.h>
#include <inttypes.h>
#include "portability.h" /* For vasprintf */
#include "verbose.h"

#define G(X) CADO_VERBOSE_PRINT_ ## X
#define F(X) (UINT64_C(1) << G(X))

/* Mutex for verbose_output_*() functions */
static pthread_mutex_t io_mutex[1] = {PTHREAD_MUTEX_INITIALIZER};
static pthread_cond_t io_cond[1] = {PTHREAD_COND_INITIALIZER};
static int batch_locked = 0;
static pthread_t batch_owner;

const char * verbose_flag_list[] = 
{
    [G(CMDLINE)]                = "print-cmdline",
    [G(MODIFIED_FILES)]         = "print-modified-files",
    [G(COMPILATION_INFO)]       = "print-compilation-info",
    [G(BWC_DISPATCH_SLAVES)]    = "bwc-dispatch-slaves",
    [G(BWC_DISPATCH_MASTER)]    = "bwc-dispatch-master",
    [G(BWC_TIMING_GRIDS)]       = "bwc-timing-grids",
    [G(BWC_ITERATION_TIMINGS)]  = "bwc-iteration-timings",
    [G(BWC_CACHE_BUILD)]        = "bwc-cache-build",
    [G(BWC_DISPATCH_OUTER)]     = "bwc-dispatch-outer",
    [G(BWC_CPUBINDING)]         = "bwc-cpubinding",
    [G(BWC_CACHE_MAJOR_INFO)]   = "bwc-cache-major-info",
    [G(BWC_LOADING_MKSOL_FILES)]= "bwc-loading-mksol-files",
};

struct {
    const char * name;
    uint64_t mask;
} verbose_flag_groups[] = {
    { "all-cmdline",
            F(CMDLINE) |
            F(MODIFIED_FILES) |
            F(COMPILATION_INFO) },
    { "all-bwc-dispatch",
            F(BWC_DISPATCH_SLAVES) |
            F(BWC_DISPATCH_MASTER) |
            F(BWC_DISPATCH_OUTER) |
            F(BWC_CACHE_BUILD) },
    { "all-bwc-sub-timings", 
            F(BWC_TIMING_GRIDS) |
            F(BWC_ITERATION_TIMINGS) },
};


uint64_t verbose_flag_word;


/* This must be called in single-threaded context, preferably at program
 * start */
void verbose_interpret_parameters(param_list pl)
{
    verbose_flag_word = ~0UL;

    const char * v = param_list_lookup_string(pl, "verbose_flags");
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
        for(size_t i = 0 ; i < sizeof(verbose_flag_list) / sizeof(verbose_flag_list[0]) ; i++) {
            if (strcmp(p, verbose_flag_list[i]) == 0) {
                mask = UINT64_C(1) << (int) i;
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

void verbose_decl_usage(param_list pl)
{
    param_list_decl_usage(pl, "verbose_flags", "fine grained control on which messages get printed");
}

/* returns true if the following verbose flag is enabled */
int verbose_enabled(int flag) {
    return verbose_flag_word & (UINT64_C(1) << flag);
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
   no other thread holds a batch lock */
static int
monitor_enter()
{
    if (pthread_mutex_lock(io_mutex) != 0)
        return 1;
    while (batch_locked && !pthread_equal(batch_owner, pthread_self())) {
        /* Queue this thread as waiting for the condition variable, release
           the mutex and put thread to sleep */
        if (pthread_cond_wait(io_cond, io_mutex) != 0)
          return 1; /* Should we unlock first? */
    }
    /* Now we own the mutex and no other thread holds the batch lock */
    return 0;
}

static int
monitor_leave()
{
    return pthread_mutex_unlock(io_mutex);
}

int
verbose_output_start_batch()
{
    if (monitor_enter() != 0)
        return 1;
    batch_locked = 1;
    batch_owner = pthread_self();
    if (monitor_leave() != 0)
        return 1;
    return 0;
}

int
verbose_output_end_batch()
{
    if (monitor_enter() != 0)
        return 1;
    ASSERT_ALWAYS(batch_locked);
    ASSERT_ALWAYS(pthread_equal(batch_owner, pthread_self()));
    batch_locked = 0;
    batch_owner = 0;
    if (pthread_cond_broadcast(io_cond) != 0)
        return 1;
    if (monitor_leave() != 0)
        return 1;
    return 0;
}

struct outputs_s {
    size_t nr_outputs;
    int *verbosity;
    FILE **outputs;
};

static void
init_output(struct outputs_s * const output)
{
    output->nr_outputs = 0;
    output->verbosity = NULL;
    output->outputs = NULL;
}

static void
clear_output(struct outputs_s * const output)
{
    free (output->outputs);
    free (output->verbosity);
    output->nr_outputs = 0;
    output->outputs = NULL;
    output->verbosity = NULL;
}

static int
add_output(struct outputs_s *output, FILE * const out, const int verbosity)
{
    const size_t new_nr = output->nr_outputs + 1;

    FILE ** const new_output = (FILE **) realloc(output->outputs, new_nr * sizeof(FILE *));
    if (new_output == NULL)
        return 0;
    output->outputs = new_output;

    int * const new_verbosity = (int *) realloc(output->verbosity, new_nr * sizeof(int));
    if (new_verbosity == NULL)
        return 0;
    output->verbosity = new_verbosity;

    output->nr_outputs = new_nr;
    output->outputs[new_nr - 1] = out;
    output->verbosity[new_nr - 1] = verbosity;
    return 1;
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
vfprint_output(const struct outputs_s * const output, const int verbosity,
               vfprintf_func_t func, const char * const fmt, va_list va)
{
    int rc = 0;
    /* For each output attached to this channel */
    for (size_t i = 0; i < output->nr_outputs; i++) {
        /* print string if output verbosity is at least "verbosity" */
        if (output->verbosity[i] >= verbosity) {
            va_list va_copied;
            va_copy(va_copied, va);
            rc = func(output->outputs[i], fmt, va_copied);
            va_end(va_copied);
            if (rc < 0)
                return rc;
        }
    }
    return rc;
}

/* Static variables, the poor man's Singleton. */
static size_t _nr_channels = 0;
static struct outputs_s *_channel_outputs = NULL;

int
verbose_output_init(const size_t nr_channels)
{
    if (monitor_enter() != 0)
        return 0;
    _channel_outputs = (struct outputs_s *) malloc(nr_channels * sizeof(struct outputs_s));
    if (_channel_outputs == NULL) {
        pthread_mutex_unlock(io_mutex);
        return 0;
    }
    _nr_channels = nr_channels;
    for (size_t i = 0; i < nr_channels; i++)
        init_output(&_channel_outputs[i]);
    if (monitor_leave() != 0)
        return 0;
    return 1;
}

int
verbose_output_clear()
{
    if (monitor_enter() != 0)
        return 0;
    for (size_t i = 0; i < _nr_channels; i++)
        clear_output(&_channel_outputs[i]);
    free(_channel_outputs);
    _nr_channels = 0;
    _channel_outputs = NULL;
    if (monitor_leave() != 0)
        return 0;
    return 1;
}

int
verbose_output_add(const size_t channel, FILE * const out, const int verbose)
{
    if (monitor_enter() != 0)
        return 0;
    ASSERT_ALWAYS(channel < _nr_channels);
    int rc = add_output(&_channel_outputs[channel], out, verbose);
    if (monitor_leave() != 0)
        return 0;
    return rc;
}

int
verbose_output_print(const size_t channel, const int verbose,
                     const char * const fmt, ...)
{
    va_list ap;
    int rc = 0;

    if (monitor_enter() != 0)
        return -1;
    va_start(ap, fmt);
    if (_channel_outputs == NULL) {
        /* Default behaviour: print to stdout or stderr */
        ASSERT_ALWAYS(channel < 2);
        if (verbose <= 1) {
            FILE *out = (channel == 0) ? stdout : stderr;
            rc = vfprintf(out, fmt, ap);
        }
    } else {
        ASSERT_ALWAYS(channel < _nr_channels);
        rc = vfprint_output(&_channel_outputs[channel], verbose, &vfprintf,
                            fmt, ap);
    }
    va_end(ap);
    if (monitor_leave() != 0)
        return -1;
    return rc;
}

FILE *
verbose_output_get(const size_t channel, const int verbose, const size_t index)
{
    if (monitor_enter() != 0)
        return NULL;

    FILE *output = NULL;
    if (_channel_outputs == NULL) {
        /* Default behaviour: channel 0 has stdout, channel 1 has stderr,
           each with verbosity 1. */
        ASSERT_ALWAYS(channel < 2);
        if (verbose <= 1) {
            output = (channel == 0) ? stdout : stderr;
        }
    } else {
        ASSERT_ALWAYS(channel < _nr_channels);
        struct outputs_s * const chan = &_channel_outputs[channel];
        size_t j = 0;
        /* Iterate through all the outputs for this channel */
        for (size_t i = 0; i < chan->nr_outputs; i++) {
            /* Count those outputs that have verbosity at least "verbose" */
            if (chan->verbosity[i] >= verbose) {
                /* If that's the index-th output with enough verbosity,
                   return it. */
                if (index == j) {
                    output = chan->outputs[j];
                    break;
                }
                j++;
            }
        }
    }

    if (monitor_leave() != 0)
        return NULL;
    return output;
}

/* On Windows, the format strings for PRId64 and PRIu64 are "I64d" and "I64u",
   resp., but GMP's printf does not understand those. This is the source of
   hard-to-debug crashes, so we try to protect against this error here.
   We don't include gmp.h here and we don't require that this module be linked
   against GMP, so we try to guess whether the format string is intended for
   GMP's printf by looking for the %Zd conversion, rather than, say, comparing
   the function pointer to GMP's gmp_vfprintf() function. */
static void
die_on_MINGW_PRI64(const char * const fmt MAYBE_UNUSED)
{
#ifdef HAVE_MINGW
  if (strstr(fmt, "%Zd") == NULL)
    return;
  ASSERT_ALWAYS(strstr(fmt, "%" PRId64) == NULL);
  ASSERT_ALWAYS(strstr(fmt, "%" PRIu64) == NULL);
#endif
}


int
verbose_output_vfprint(const size_t channel, const int verbose,
                       vfprintf_func_t func, const char * const fmt, ...)
{
    va_list ap;
    int rc = 0;

    die_on_MINGW_PRI64(fmt);

    if (monitor_enter() != 0)
        return -1;
    va_start(ap, fmt);
    if (_channel_outputs == NULL) {
        /* Default behaviour: print to stdout or stderr */
        ASSERT_ALWAYS(channel < 2);
        if (verbose <= 1) {
            FILE *out = (channel == 0) ? stdout : stderr;
            rc = func(out, fmt, ap);
        }
    } else {
        ASSERT_ALWAYS(channel < _nr_channels);
        rc = vfprint_output(&_channel_outputs[channel], verbose, func, fmt,
                            ap);
    }
    va_end(ap);
    if (monitor_leave() != 0)
        return -1;
    return rc;
}
