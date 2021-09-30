#ifndef POLYSELECT_MAIN_QUEUE_H_
#define POLYSELECT_MAIN_QUEUE_H_

#include <gmp.h>
#include <stdint.h>

#include "cado_poly.h"
#include "polyselect_main_data.h"
#include "polyselect_locals.h"

/* A few configuration flags first, which affect some specific points.
 */
/* This is used in the collisions calls */
#define INIT_FACTOR 8UL

/* number of special (q, r) per batch */
#define BATCH_SIZE 20

#define DEFAULT_NQ 1000		/* default max num of nq considered for each ad */

#define DEFAULT_INCR 60 /* we want a positive integer with many divisors,
                           other values are 210, 2310, 30030, 510510, 9699690,
                           223092870 */

#define DEFAULT_POLYSELECT_KEEP 10			/* number of best raw polynomials kept */

#ifdef __cplusplus
extern "C" {
#endif


/* This is one the functions through which a raw polynomial pair,
 * found from a collision in the tables, undergoes further processing.
 *
 * Note that the stats structure to be used is _not_ protected by locks,
 * so that it's much better if it's a thread-local stats object. It will
 * be pushed to the global stats with polyselect_main_data_commit_stats
 */
extern int optimize_raw_poly(mpz_poly_ptr f, mpz_poly_ptr g,
        polyselect_main_data_srcptr,
        polyselect_stats_ptr stats);

extern void polyselect_fprintf_poly_pair(FILE * fp, mpz_srcptr N,                    mpz_poly_srcptr f, mpz_poly_srcptr g, int raw);

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_MAIN_QUEUE_H_ */
