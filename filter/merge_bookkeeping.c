#include "cado.h"
#include <stdio.h>
#include "timing.h"  // seconds
#include "merge_bookkeeping.h"

int merge_pass = 0;

#ifdef BIG_BROTHER
unsigned char *touched_columns = NULL;
#endif

/* define CANCEL to count column cancellations */
// #define CANCEL
#ifdef CANCEL
unsigned long cancel_rows = 0;
unsigned long cancel_cols[CANCEL_MAX] = {0,};
#endif

double cpu_t[8] = {0};
double wct_t[8] = {0};

int merge_verbose = 0; /* verbosity level */


void
print_timings (char *s, double cpu, double wct)
{
  printf ("%s %.1fs (cpu), %.1fs (wct) [cpu/wct=%.1f]\n",
	  s, cpu, wct, cpu / wct);
  fflush (stdout);
}

