#ifndef LAS_OUTPUT_HPP_
#define LAS_OUTPUT_HPP_

#include <cstdio>
struct cxx_param_list;

struct las_output {
    int verbose = 1;
    /* outputname is owned by pl, and output by the libc. We want to
     * control when these get released, with the release() function */
    FILE *output = NULL;
    const char * outputname = NULL; /* keep track of whether it's gzipped or not */
    void set(cxx_param_list & pl);
    void release();
    void fflush();
    static void declare_usage(cxx_param_list & pl);
    static void configure_aliases(cxx_param_list &) { }
    static void configure_switches(cxx_param_list & pl);
};

/* used in verbose_output_print.
 *
 * This is a potential bug, as the channels are tied globally to FILE*
 * pointers. This should be fixed. (_channel_outputs in verbose.c)
 */
enum {
  output,
  ERROR_CHANNEL,
  STATS_CHANNEL,
  TRACE_CHANNEL,
  NR_CHANNELS /* This must be the last element of the enum */
};

#endif	/* LAS_OUTPUT_HPP_ */
