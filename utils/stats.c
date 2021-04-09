#include "cado.h" // IWYU pragma: keep
#include <math.h>       // INFINITY
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <pthread.h>
#include "macros.h"
#include "stats.h"
#include "timing.h"     // wct_seconds

void
stats_init (stats_data_t r, FILE *f, uint64_t *followed_var,
            uint8_t max_log_report, const char *verb, const char *name,
            const char *outofname, const char *abbrv)
{
  r->followed_var = followed_var;
  r->last_report = 0;
  r->t0 = wct_seconds();
  r->out = f;
  r->log_report = MIN_LOG_REPORT;
  if (max_log_report >= MIN_LOG_REPORT)
    r->max_log_report = max_log_report;
  else
    r->max_log_report = MIN_LOG_REPORT;
  r->name = name;
  r->abbrv = abbrv;
  r->verb = verb;
  r->outofname = outofname;
  r->end_done = 0;
}

void stats_clear(stats_data_t S MAYBE_UNUSED)
{
    /* no-op, in fact */
}

/* Return 1 if more than 2^r->log_report relations were read since last
   progress report.
   Otherwise return 0 */
int
stats_test_progress (stats_data_t r)
{
  uint64_t v;
  uint64_t w;
  int l;
#pragma omp atomic read
  v = *(r->followed_var);
#pragma omp atomic read
  l = r->log_report;
#pragma omp atomic read
  w = r->last_report;
  if ((v >> l) != w)
    return 1;
  else
    return 0;
}

/* Print a line of the form: 
 *    Read 42 relations in 1.4s -- 30.0 rels/s
 * Prepend a "Done: " if end is non-zero.
 * Add "-- xy.z MB/s " in the middle, if nByte > 0, with xy.z being
 * nByte/time spent.
 * Add "(out of <outof> ssss) " in the middle, if outof > 0, with ssss being
 * r->outofname.
 */
void
stats_print_progress (stats_data_t r, uint64_t i, uint64_t outof, size_t nByte,
                      int end)
{
  char MBpart[32] = "";
  char outofpart[64] = "";
  double t, dt, speed;

  if (!end && !stats_test_progress(r)) {
      return;
  }

  static pthread_mutex_t mm = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_lock(&mm);
  if (!(end && !r->end_done) && !stats_test_progress(r)) {
      pthread_mutex_unlock(&mm);
      return;
  }
  if (end)
      r->end_done = 1;

  t = wct_seconds();
  dt = t - r->t0;
  speed = dt > 0.01 ? i/dt : INFINITY;
  if (nByte > 0)
  {
    double mb_s = dt > 0.01 ? (nByte/dt * 1.0e-6) : INFINITY;
    snprintf (MBpart, 32, "-- %.1f MB/s ", mb_s);
  }
  if (outof > 0)
    snprintf (outofpart, 64, "(out of %" PRIu64 " %s) ", outof, r->outofname);
  const char * prefix = (end) ? "# Done: " : "# ";
  fprintf(r->out, "%s%s %" PRIu64 " %s %sin %.1fs %s-- %.1f %s/s\n",
          prefix, r->verb, i, r->name, outofpart, dt, MBpart, speed, r->abbrv);
  fflush(r->out);

  /* We only want to ensure that the reads are always consistent. We're
   * the only place doing writes anyway.
   */
  int l;
  int ml = r->max_log_report;
  uint64_t v;
#pragma omp atomic read
  l = r->log_report;
#pragma omp atomic update
  r->log_report += l < ml;
  l += l < ml;
#pragma omp atomic read
  v = *(r->followed_var);
#pragma omp atomic write
  r->last_report = v >> l;

  pthread_mutex_unlock(&mm);
}
