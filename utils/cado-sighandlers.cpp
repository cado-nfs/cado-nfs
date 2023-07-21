#include "cado.h" // IWYU pragma: keep

#include <cstdio>     // for fprintf, stderr
#include <cstring>    // for strsignal

#ifdef HAVE_EXECINFO
#include <execinfo.h>                    // for backtrace, backtrace_symbols
#include <csignal>                      // for signal, raise, SIGABRT, SIGSEGV
#else
#include "verbose.h"    // verbose_output_print
#endif
#include "cado-sighandlers.h"

#ifdef HAVE_EXECINFO
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

void cado_sighandlers_install()
{
#ifdef HAVE_EXECINFO
    signal (SIGABRT, signal_handling);
    signal (SIGSEGV, signal_handling);
#else
    verbose_output_print(0, 0, "# Cannot catch signals in an interesting way, lack library support\n");
#endif
}

