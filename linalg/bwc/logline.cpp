#include "cado.h" // IWYU pragma: keep
#include <cstdint>      // for SIZE_MAX
#include <cstring>      // for memset, strlen
#include <cstdio>       // fprintf // IWYU pragma: keep
#include <cstdarg>      // va_list // IWYU pragma: keep
#include <cstdlib>
#include <string>        // for string, basic_string
#include <vector>
#include "logline.h"
#include "memusage.h"   // Memusage2
#include "params.h"     // param_list_parse_*
#include "select_mpi.h"
#include "timing.h"     // seconds
#include "portability.h" // asprintf // IWYU pragma: keep
#include "misc.h"       // size_disp_fine
#include "macros.h"

/* This is intended to provide progress info, and optionally also more
 * detailed progress info when needed. Detailed progress info is assumed
 * to possibly take several lines, while terse progress info is expected
 * to take only one line.
 */

/* Log messages are provided in the form [[x, "text"]], with x an integer
 * indicating a level of detail. 0 is the lowest level of detail, and
 * higher means more detailed (and printed less often).
 *
 * Whenever we output an info which is more detailed than what we had
 * before, and with no finishing newline, then the unfinished previous
 * line is terminated with ..., and the terminating messaged will be
 * issued alone later on.
 *
 * An exception is for level 0 (introduced by logline_begin()), for which
 * the "title" is then echoed when finishing the printing.
 */
struct logline {
    /* Only log messages prefixed with <x> with x <= max will be printed
     */
    int max = 0;
    /* Title string introducing this log set */
    std::string header;
    /* Did last print contain a trailing newline ? */
    bool eol = false;
    /* Last printed detail level */
    int lastlevel = 0;

    FILE * f = NULL;

    /* Number of newlines seen (used to tell whether the title has to be
     * echoed again) */
    int nnl = 0;

    double start = 0;

    std::vector<std::string> prefixes;
};

double logline_serialize();
void logline_unserialize(double);

static struct logline * current;

size_t logline_thresholds[10] = {
    /* This list should be increasing */
    SIZE_MAX,
    SIZE_MAX,
    SIZE_MAX,
    SIZE_MAX,
    SIZE_MAX,
    SIZE_MAX,
    SIZE_MAX,
    SIZE_MAX,
    SIZE_MAX,
    SIZE_MAX,
};

int logline_timings = 1;

int logline_print_all_mpi_nodes = 0;

static double start_time = -1;

void logline_init_timer()
{
    start_time = wct_seconds();
}

int logline_report_wct = 0;

static double logline_timer()
{
    return logline_report_wct ? wct_seconds() : seconds();
}

double logline_serialize()
{
    double tt = wct_seconds() - start_time;
    return tt;
}

void logline_unserialize(double tt)
{
    start_time = wct_seconds() - tt;
}


void logline_decl_usage(param_list_ptr pl)
{
    param_list_decl_usage(pl, "logline_threshold",
            "print log lines of verbosity level i only for sizes greater than i-th item in this comma-separated list");
    param_list_decl_usage(pl, "logline_timings",
            "print timings with log lines");
    param_list_decl_usage(pl, "logline_report_wct",
            "print wct taken by each step marked by log lines");
    param_list_decl_usage(pl, "logline_print_all_mpi_nodes",
            "enable logline printing on all MPI nodes");

}

int logline_interpret_parameters(param_list_ptr pl)
{
    int thr[10];
    int n = param_list_parse_int_list(pl, "logline_threshold", thr, 10, ",");
    for(int i = 0 ; i < n ; i++) {
        logline_thresholds[i] = thr[i];
    }
    param_list_parse_int(pl, "logline_timings", &logline_timings);
    param_list_parse_int(pl, "logline_report_wct", &logline_report_wct);
    param_list_parse_int(pl, "logline_print_all_mpi_nodes", &logline_print_all_mpi_nodes);
    return 0;
}

static void logline_puts_raw(int level, const char * s)
{
    if (!current) return;
    if (level > current->max) return;
    if (level != current->lastlevel && !current->eol) {

        fputs(" ...\n", current->f);
        current->nnl += (current->eol = 1);
    }
    for( ; current->prefixes.size() > (unsigned int) level ; current->prefixes.pop_back()) ;
    if (current->eol) {
        if (logline_timings) {
            char buf1[16];
            char buf2[16];
            size_disp_fine(1024UL * Memusage2(), buf1, 10000.0);
            size_disp_fine(1024UL * PeakMemusage(), buf2, 10000.0);
            fprintf(current->f, "[%.3f %.3f %s %s] ", wct_seconds() - start_time, seconds(), buf1, buf2);
        }
        for(auto const & c : current->prefixes) {
            fputs(c.c_str(), current->f);
            fputc(' ', current->f);
        }
    } else if (level == current->lastlevel) {
        fputc(' ', current->f);
    }
    fputs(s, current->f);
    current->prefixes.push_back(s);
    size_t n = strlen(s);
    current->nnl += (current->eol = s[n-1] == '\n');
    current->lastlevel = level;
}


int logline_begin(FILE * f, size_t size, const char * fmt, ...)
{
    va_list ap;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank && !logline_print_all_mpi_nodes) return 0;
    ASSERT_ALWAYS(!current);
    if (size < logline_thresholds[0]) return 0;
    int level;
    for(level = 0 ; level < 10 && size >= logline_thresholds[level] ; level++);
    if (level == 0) return 0;
    level--;
    current = new logline;
    current->max = level;
    current->lastlevel = 0;
    current->f = f;
    current->start = logline_timer();
    char * tmp;
    va_start(ap, fmt);
    int rc = vasprintf(&(tmp), fmt, ap);
    ASSERT_ALWAYS(rc >= 0);
    current->header = tmp;
    free(tmp);
    current->eol = true;
    logline_puts_raw(0, current->header.c_str());
    va_end(ap);
    return 1;
}

int logline_end(double * rr, const char * fmt, ...)
{
    if (!current) return 0;
    va_list ap;
    va_start(ap, fmt);
    char * text2;
    double tt = logline_timer() - current->start;
    if (fmt) {
        char * text;
        int rc;
        rc = vasprintf(&text, fmt, ap);
        ASSERT_ALWAYS(rc >= 0);
        rc = asprintf(&text2, "%s [%.2f]\n", text, tt);
        ASSERT_ALWAYS(rc >= 0);
        free(text);
    } else {
        int rc;
        rc = asprintf(&text2, "[%.2f]\n", tt);
        ASSERT_ALWAYS(rc >= 0);
    }
    if (current->nnl)
        logline_puts_raw(0, current->header.c_str());
    logline_puts_raw(0, text2);
    free(text2);
    delete current;
    current = NULL;
    va_end(ap);
    if (rr) *rr += tt;
    return 1;
}


int logline_vprintf(int level, const char * fmt, va_list ap)
{
    if (!current) return 0;
    char * text;
    int rc = vasprintf(&text, fmt, ap);
    ASSERT_ALWAYS(rc >= 0);
    logline_puts_raw(level, text);
    free(text);
    return 1;
}

int logline_printf(int level, const char * fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    logline_vprintf(level, fmt, ap);
    va_end(ap);
    return 1;
}


