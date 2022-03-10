#ifndef POLYSELECT_STATS_H_
#define POLYSELECT_STATS_H_

#include "polyselect_data_series.h"
#include "polyselect_priority_queue.h"

#ifdef __cplusplus
extern "C" {
#endif

struct polyselect_stats_s {
    unsigned long number_of_ad_values;    /* number of ad values that are
                                             represented by this stats
                                             object */

    int tot_found;		/* total number of polynomials */

    double potential_collisions;
    unsigned long discarded1;	/* f[d] * f[d-2] > 0 */
    unsigned long collisions;
    unsigned long discarded2;	/* f[d-1] * f[d-3] > 0 */
    unsigned long collisions_good;
    unsigned long opt_found;	/* number of size-optimized polynomials */
    double optimize_time;

    /* we've still got to see which are really updated when */
    polyselect_data_series_t raw_lognorm;
    polyselect_data_series_t opt_lognorm;
    polyselect_data_series_t exp_E;
    polyselect_data_series_t beta;
    polyselect_data_series_t eta;
    /* 8184e1a188 introduced a seconday "best_exp_E" series, which is a
     * name clash with the one above. We don't really know if these
     * should be the same thing or not.
     */
    polyselect_data_series_t best_exp_E_Weibull;
    polyselect_data_series_t raw_proj_alpha;
    polyselect_data_series_t opt_proj_alpha;

    /* currently used only for best_exp_E with weibull_moments2
     */
    gmp_randstate_t rstate;

    /* record the start times */
    double st0, wct0;

    polyselect_priority_queue_t best_opt_logmu;
    polyselect_priority_queue_t best_exp_E;
};

typedef struct polyselect_stats_s polyselect_stats[1];
typedef struct polyselect_stats_s * polyselect_stats_ptr;
typedef const struct polyselect_stats_s * polyselect_stats_srcptr;

extern void polyselect_stats_init(polyselect_stats_ptr stats, size_t keep);
extern void polyselect_stats_clear(polyselect_stats_ptr stats);
extern void polyselect_stats_reset(polyselect_stats_ptr stats);

extern void polyselect_stats_update_keep(polyselect_stats_ptr stats, size_t keep);

extern void polyselect_stats_accumulate(polyselect_stats_ptr to, polyselect_stats_srcptr from);

extern void polyselect_stats_display_final(polyselect_stats_ptr stats, int verbose);

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_STATS_H_ */
