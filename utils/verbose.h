#ifndef VERBOSE_H_
#define VERBOSE_H_

#include <stdarg.h>
#include <stdio.h>
#include "params.h"

#ifdef __cplusplus
extern "C" {
#endif

/* This provides fine-grain control over which messages we want. In order
 * to add a new flag, you must define its name as a preprocessor integer
 * macro here, and define a corresponding text string in verbose.c
 */

#define CADO_VERBOSE_PRINT_CMDLINE                      0
#define CADO_VERBOSE_PRINT_MODIFIED_FILES               1
#define CADO_VERBOSE_PRINT_COMPILATION_INFO             2
#define CADO_VERBOSE_PRINT_BWC_DISPATCH_SLAVES          3
#define CADO_VERBOSE_PRINT_BWC_DISPATCH_MASTER          4
#define CADO_VERBOSE_PRINT_BWC_TIMING_GRIDS             5
#define CADO_VERBOSE_PRINT_BWC_ITERATION_TIMINGS        6
#define CADO_VERBOSE_PRINT_BWC_CACHE_BUILD              7
#define CADO_VERBOSE_PRINT_BWC_DISPATCH_OUTER           8
#define CADO_VERBOSE_PRINT_BWC_CPUBINDING               9
#define CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO         10
#define CADO_VERBOSE_PRINT_BWC_LOADING_MKSOL_FILES      11

typedef int (*vfprintf_func_t)(FILE *, const char *, va_list);

/* This must be called in single-threaded context, preferably at program
 * start */
extern void verbose_interpret_parameters(param_list pl);

extern int verbose_enabled(int flag);
extern int verbose_vfprintf(FILE * f, int flag, const char * fmt, va_list ap);
extern int verbose_vprintf(int flag, const char * fmt, va_list ap);
extern int verbose_fprintf(FILE * f, int flag, const char * fmt, ...);
extern int verbose_printf(int flag, const char * fmt, ...);
extern void verbose_decl_usage(param_list pl);


/*
  The program can initialise zero or more "channels".
  To each "channel", zero or more "outputs" (FILE handles) can be attached,
  each with an output verbosity value.
  Text can be printed to a channel, with a verbosity v, then the text is
  sent to each output attached to this channel for which the output verbosity
  is at least v.
  The default behaviour, before calling the verbose_output_init() function or
  after calling verbose_output_clear(), has 2 channels:
  channel 0 with 1 output, going to stdout with verbosity 1, and
  channel 1 with 1 output, going to stderr with verbosity 1.
  All functions are mutex-protected, i.e., the output system forms a
  "monitor."
*/

/* Lock the output system for exclusive use for a while.
   Useful for, e.g., printing lines with multiple calls, or multiple
   consecutive lines, through multiple verbose_output_*() calls. */
int verbose_output_start_batch();

/* Unlock the output system again. Only the locking thread may call it. */
int verbose_output_end_batch();

/* Init nr_channels channels, with no outputs attached
   Returns 1 on success, and 0 on error. */
int verbose_output_init(size_t nr_channels);

/* Reset channels back to the 2 default channels.
   Returns 1 on success, and 0 on error. */
int verbose_output_clear();

/* Add an output FILE handle to a channel with the given verbosity.
   Returns 1 on success, and 0 on error. */
int verbose_output_add(size_t channel, FILE * out, int verbose);

/* Print formatted output (using libc's printf() formatting) to all outputs
   attached to the given channel where the output's verbosity is at least
   the given verbosity parameter.
   Returns the return value of the last printf() call on success,
   or 0 if nothing gets printed. Returns a negative value on error. */
int verbose_output_print(size_t channel, int verbose, const char *fmt, ...) ATTR_PRINTF (3, 4);

/* Get the index-th FILE handle attached to a channel, counting only those
   outputs whose output verbosity is at least "verbose."
   If there are fewer than index such handles, returns NULL.

   This is rather hackish, but we need it to be able to pass a FILE handle
   to functions that, e.g., print complex data to a stream.

   The idiom to use to print to every attached output would be something like

   FILE *out;
   for (size_t i=0; (out=verbose_output_get(c, v, i)) != NULL; i++)
       print_something(out, ...);
*/
FILE *verbose_output_get(size_t channel, int verbose, size_t index);

/* Print to every attached output, using a print function with a vfprintf()-
   like interface. E.g.,

   verbose_output_vfprint(0, 1, gmp_vfprintf, "%Zd\n", bigint);

   Note that GMP requires stdarg.h to be included BEFORE gmp.h to be able
   to declare prototypes for functions that take va_list.
   
   WARNING: do not print PRI?64 format strings with gmp_vfprintf(), this
   crashes on MinGW! */
int verbose_output_vfprint(size_t channel, int verbose, vfprintf_func_t func,
                           const char * fmt, ...);

#ifdef __cplusplus
}
#endif

#endif	/* VERBOSE_H_ */
