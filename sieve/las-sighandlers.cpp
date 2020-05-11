#include "cado.h" // IWYU pragma: keep

#include <cstdio>     // for fprintf, stderr
#include <cstring>    // for strsignal

#ifdef HAVE_CXXABI_H
/* We use that to demangle C++ names */
#include <cxxabi.h> // IWYU pragma: keep
#endif
#ifdef HAVE_GLIBC
#include <execinfo.h>                    // for backtrace, backtrace_symbols
#include <csignal>                      // for signal, raise, SIGABRT, SIGSEGV
#endif
#include "utils.h"      // IWYU pragma: keep

#ifdef HAVE_GLIBC
static void signal_handling (int signum)/*{{{*/
{
   fprintf (stderr, "*** Error: caught signal \"%s\"\n", strsignal (signum));

   int sz = 100, i;
   void *buffer [sz];
   char** text;

   sz = backtrace (buffer, sz);
   text = backtrace_symbols (buffer, sz);

   fprintf(stderr, "======= Backtrace: =========\n");
   for (i = 0; i < sz; i++)
       fprintf (stderr, "%s\n", text [i]);

   signal (signum, SIG_DFL);
   raise (signum);
}/*}}}*/
#endif

void las_sighandlers_install()
{
#ifdef HAVE_GLIBC
    signal (SIGABRT, signal_handling);
    signal (SIGSEGV, signal_handling);
#else
    verbose_output_print(0, 0, "# Cannot catch signals, lack glibc support\n");
#endif
}

