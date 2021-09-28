#ifndef POLYSELECT_MAIN_DATA_H_
#define POLYSELECT_MAIN_DATA_H_

#include <gmp.h>
#include <stdint.h>
#include <pthread.h>

#include "cado_poly.h"
#include "polyselect_data_series.h"
#include "polyselect_stats.h"
#include "polyselect_poly_header.h"
#include "polyselect_qroots.h"
#include "params.h"

#ifdef __cplusplus
extern "C" {
#endif


struct polyselect_main_data_s {
    uint32_t *Primes;
    unsigned long lenPrimes;
    unsigned long incr;

    mpz_t N;
    unsigned int d;

    unsigned long P;

    mpz_t admin, admax;
    cado_poly best_poly, curr_poly;

    double maxtime;
    double target_E;		/* target E-value, 0.0 if not given */

    unsigned int sopt_effort;

    unsigned long nq;

    /* These stats are the same as the local stats, and they are merged
     * by the polyselect_main_data_commit_stats call */
    polyselect_stats stats;

    pthread_mutex_t stats_lock;

    int verbose;
};

typedef struct polyselect_main_data_s polyselect_main_data[1];
typedef struct polyselect_main_data_s * polyselect_main_data_ptr;
typedef const struct polyselect_main_data_s * polyselect_main_data_srcptr;

extern void polyselect_main_data_init_defaults(polyselect_main_data_ptr main);
extern void polyselect_main_data_clear(polyselect_main_data_ptr main);

extern int64_t polyselect_main_data_get_M(polyselect_main_data_srcptr main);

extern size_t polyselect_main_data_expected_number_of_pairs(polyselect_main_data_srcptr main);

/* the number of expected collisions is 8*lenPrimes^2/2/(2P)^2 */
extern double polyselect_main_data_expected_collisions(polyselect_main_data_srcptr main);

extern int polyselect_main_data_check_parameters(polyselect_main_data_srcptr main, mpz_srcptr m0, double q);

extern unsigned long
find_suitable_lq(polyselect_poly_header_srcptr header,
		 polyselect_qroots_srcptr SQ_R,
                 unsigned long *k,
                 polyselect_main_data_ptr main);

extern void polyselect_main_data_commit_stats(polyselect_main_data_ptr main, polyselect_stats_ptr stats);


extern void polyselect_main_data_prepare_primes(polyselect_main_data_ptr main);

extern unsigned long polyselect_main_data_number_of_ad_tasks(polyselect_main_data_ptr main);

extern void polyselect_main_data_parse_Nd(polyselect_main_data_ptr main, param_list_ptr pl);

extern void polyselect_main_data_parse_ad_range(polyselect_main_data_ptr main, param_list_ptr pl);

extern void polyselect_main_data_parse_maxtime_or_target(polyselect_main_data_ptr main, param_list_ptr pl);

extern void polyselect_main_data_parse_P(polyselect_main_data_ptr main, param_list_ptr pl);

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_MAIN_DATA_H_ */
