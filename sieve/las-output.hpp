#ifndef CADO_LAS_OUTPUT_HPP
#define CADO_LAS_OUTPUT_HPP

#include <cstdio>
#include <string>

struct cxx_param_list;

struct las_output {
    int verbose = 1;
    FILE *output = nullptr;
    std::string outputname;
    explicit las_output(cxx_param_list & pl);
    void fflush() const;
    static void declare_usage(cxx_param_list & pl);
    static void configure_aliases(cxx_param_list &) { }
    static void configure_switches(cxx_param_list & pl);
    ~las_output();
    las_output(las_output const &) = delete;
    las_output(las_output &&) = delete;
    las_output& operator=(las_output const &) = delete;
    las_output& operator=(las_output &&) = delete;
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

#endif	/* CADO_LAS_OUTPUT_HPP */
