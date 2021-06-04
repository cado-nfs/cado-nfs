#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

/* specifically for las-output.cpp ; I'd say it's a bug, we should move
 * it to the independent tier */

#include <cstring>       // strerror (in DIE_ERRNO_DIAG) // IWYU pragma: keep
#include <cstdio>          // for NULL, stderr, fflush, fopen, fprintf, setvbuf
#include <cstdlib>         // for exit, EXIT_FAILURE
#include "gzip.h"       // fopen_maybe_compressed
#include "las-config.h"    // for las_display_config_flags
#include "las-output.hpp"  // for las_output, NR_CHANNELS, TRACE_CHANNEL
#include "macros.h"        // for ASSERT_ALWAYS, DIE_ERRNO_DIAG
#include "verbose.h"    // verbose_output_print
#include "params.h"


/* {{{ las_verbose things */
static void las_verbose_enter(cxx_param_list & pl, FILE * output, int verbose)
{
    verbose_interpret_parameters(pl);
    verbose_output_init(NR_CHANNELS);
    verbose_output_add(0, output, verbose + 1);
    verbose_output_add(1, stderr, 1);
    /* Channel 2 is for statistics. We always print them to las' normal output */
    verbose_output_add(2, output, 1);
    if (param_list_parse_switch(pl, "-stats-stderr")) {
        /* If we should also print stats to stderr, add stderr to channel 2 */
        verbose_output_add(2, stderr, 1);
    }
#ifdef TRACE_K
    const char *trace_file_name = param_list_lookup_string(pl, "traceout");
    FILE *trace_file = stderr;
    if (trace_file_name != NULL) {
        trace_file = fopen(trace_file_name, "w");
        DIE_ERRNO_DIAG(trace_file == NULL, "fopen(%s)", trace_file_name);
    }
    verbose_output_add(TRACE_CHANNEL, trace_file, 1);
#endif
}

static void las_verbose_leave()
{
    verbose_output_clear();
}
/* }}} */

/*{{{ stuff related to las output: -out, -stats-stderr, and so on. */

void las_output::fflush()
{
    verbose_output_start_batch();
    ::fflush(output);
    verbose_output_end_batch();
}

void las_output::set(cxx_param_list & pl)
{
    ASSERT_ALWAYS(output == NULL);
    output = stdout;
    outputname = param_list_lookup_string(pl, "out");
    if (outputname) {
	if (!(output = fopen_maybe_compressed(outputname, "w"))) {
	    fprintf(stderr, "Could not open %s for writing\n", outputname);
	    exit(EXIT_FAILURE);
	}
    }
    verbose = param_list_parse_switch(pl, "-v");
    setvbuf(output, NULL, _IOLBF, 0);      /* mingw has no setlinebuf */
    las_verbose_enter(pl, output, verbose);

    param_list_print_command_line(output, pl);
    las_display_config_flags();
}

/* This cannot be a dtor because las_output is a global, and it holds
 * references to pl which is local -- Thus we want to control the exact
 * time where fclose is called.
 */
void las_output::release()
{
    if (outputname)
        fclose_maybe_compressed(output, outputname);
    las_verbose_leave();
    outputname = NULL;
}


void las_output::configure_switches(cxx_param_list & pl)
{
    param_list_configure_switch(pl, "-stats-stderr", NULL);
    param_list_configure_switch(pl, "-v", NULL);
}

void las_output::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "v",    "verbose mode, also prints sieve-area checksums");
    param_list_decl_usage(pl, "out",  "filename where relations are written, instead of stdout");
#ifdef TRACE_K
    param_list_decl_usage(pl, "traceout", "Output file for trace output, default: stderr");
#endif
    param_list_decl_usage(pl, "stats-stderr", "print stats to stderr in addition to stdout/out file");
}
/*}}}*/

