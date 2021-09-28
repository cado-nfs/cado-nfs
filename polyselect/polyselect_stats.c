#include "cado.h"
#include <math.h>
#include "polyselect_stats.h"
#include "timing.h"
#include "misc.h"
#include "memusage.h"

void polyselect_stats_init(polyselect_stats_ptr stats)
{
  memset(stats, 0, sizeof(polyselect_stats));
  polyselect_data_series_init(stats->raw_lognorm);
  polyselect_data_series_init(stats->opt_lognorm);
  polyselect_data_series_init(stats->exp_E);
  polyselect_data_series_init(stats->beta);
  polyselect_data_series_init(stats->eta);
  polyselect_data_series_init(stats->best_exp_E_Weibull);
  polyselect_data_series_init(stats->raw_proj_alpha);
  polyselect_data_series_init(stats->opt_proj_alpha);
  gmp_randinit_default(stats->rstate);
  stats->best_exp_E_Weibull->rstate = stats->rstate;
  stats->st0 = seconds();
  stats->wct0 = wct_seconds();
}

void polyselect_stats_clear(polyselect_stats_ptr stats)
{
  polyselect_data_series_clear(stats->raw_lognorm);
  polyselect_data_series_clear(stats->opt_lognorm);
  polyselect_data_series_clear(stats->exp_E);
  polyselect_data_series_clear(stats->beta);
  polyselect_data_series_clear(stats->eta);
  polyselect_data_series_clear(stats->best_exp_E_Weibull);
  polyselect_data_series_clear(stats->raw_proj_alpha);
  polyselect_data_series_clear(stats->opt_proj_alpha);
  gmp_randclear(stats->rstate);
  memset(stats, 0, sizeof(polyselect_stats));
}

void polyselect_stats_setup_keep_best(polyselect_stats_ptr stats, int keep)
{
  stats->keep = keep;
  stats->best_opt_logmu = (double *) malloc(keep * sizeof(double));
  ASSERT_ALWAYS(stats->best_opt_logmu != NULL);
  stats->best_exp_E = (double *) malloc(keep * sizeof(double));
  ASSERT_ALWAYS(stats->best_exp_E != NULL);
  for (int i = 0; i < keep; i++)
    {
      stats->best_opt_logmu[i] = NAN;	/* best logmu after size optimization */
      stats->best_exp_E[i] = NAN;	/* best logmu after rootsieve */
    }
}

void polyselect_stats_display_final(polyselect_stats_ptr stats, int verbose)
{
  /* finishing up statistics */
  if (verbose >= 0)
    {
      printf("# Stat: potential collisions=%1.2f (%1.2e/s)\n",
	     stats->potential_collisions, 1000.0 * stats->potential_collisions
	     / (double) milliseconds());
      if (stats->collisions > 0)
	{
          char tmp[128] = { '\0' };

          ASSERT_ALWAYS(stats->collisions == stats->raw_lognorm->size);
          ASSERT_ALWAYS(stats->collisions == stats->raw_proj_alpha->size);

          polyselect_data_series_snprintf_summary(tmp, sizeof(tmp), stats->raw_lognorm);
	  printf ("# Stat: raw lognorm%s\n", tmp);

          polyselect_data_series_snprintf_summary(tmp, sizeof(tmp), stats->raw_proj_alpha);
	  printf ("# Stat: raw proj. alpha%s\n", tmp);
	  printf
	      ("# Stat: discarded %lu polynomials because f[d]*f[d-2] > 0\n",
	       stats->discarded1);
	  if (stats->collisions_good > 0)
	    {
              ASSERT_ALWAYS(stats->collisions_good == stats->opt_lognorm->size);
              ASSERT_ALWAYS(stats->collisions_good == stats->opt_proj_alpha->size);
              ASSERT_ALWAYS(stats->collisions_good == stats->exp_E->size);
              polyselect_data_series_snprintf_summary(tmp, sizeof(tmp), stats->opt_lognorm);
	      printf ("# Stat: optimized lognorm%s\n", tmp);
              polyselect_data_series_snprintf_summary(tmp, sizeof(tmp), stats->opt_proj_alpha);
	      printf ("# Stat: opt proj. alpha%s\n", tmp);
	      /* the exp_E statistics can be used as follows: if the mean is
	         m and the standard deviation s, then assuming a normal
	         distribution, the minimum order statistic for K polynomials
	         is given by:
	         m - s*(sqrt(2log(K))-(log(log(K))+1.377)/(2*sqrt(2log(K))))
	       */
	      printf ("# Stat: discarded %lu polynomials because f[d-1]*f[d-3] > 0\n",
		   stats->discarded2);
              polyselect_data_series_snprintf_summary(tmp, sizeof(tmp), stats->exp_E);
	      printf ("# Stat: exp_E%s\n", tmp);
	    }
	}
    }

  printf
      ("# Stat: tried %lu ad-value(s), %lu size-optimized polynomials, kept %lu\n",
       stats->number_of_ad_values, stats->opt_found, stats->collisions_good);

  /* print best keep values of logmu */
  if (stats->collisions_good > 0)
    {
      printf("# Stat: best exp_E after size optimization:");
      for (int i = 0; i < stats->keep; i++)
	if (!isnan(stats->best_exp_E[i]))
	  printf(" %1.2f", stats->best_exp_E[i]);
      printf("\n");
    }

  /* print total time (this gets parsed by the scripts) */
  printf("# Stat: total phase took %.2fs\n", seconds() - stats->st0);

  /* If we can't time this because of multithreading and missing
   * RUSAGE_THREAD, it gets set to -1, so don't print it in that case.
   */
  if (stats->optimize_time >= 0) {
    printf("# Stat: size-optimization took %.2fs\n", stats->optimize_time);
  }

  {
    char buf[16];
    printf("# Stat: peak mem usage %s\n",
	   size_disp(PeakMemusage() << 10, buf));
  }
}

/* typically used to accumulate thread-level stats into global stats.
 * Don't use without proper locking!
 *
 * Regarding polyselect: prefer polyselect_main_data_commit_stats
 */
void polyselect_stats_accumulate(polyselect_stats_ptr to, polyselect_stats_srcptr from)
{
    to->number_of_ad_values += from->number_of_ad_values;
    to->tot_found += from->tot_found;
    to->potential_collisions += from->potential_collisions;
    to->discarded1 += from->discarded1;
    to->collisions += from->collisions;
    to->discarded2 += from->discarded2;
    to->collisions_good += from->collisions_good;
    to->opt_found += from->opt_found;
    to->optimize_time += from->optimize_time;

    polyselect_data_series_combine(to->raw_lognorm, from->raw_lognorm);
    polyselect_data_series_combine(to->opt_lognorm, from->opt_lognorm);
    polyselect_data_series_combine(to->exp_E, from->exp_E);
    polyselect_data_series_combine(to->beta, from->beta);
    polyselect_data_series_combine(to->eta, from->eta);
    polyselect_data_series_combine(to->best_exp_E_Weibull, from->best_exp_E_Weibull);
    polyselect_data_series_combine(to->raw_proj_alpha, from->raw_proj_alpha);
    polyselect_data_series_combine(to->opt_proj_alpha, from->opt_proj_alpha);
}

