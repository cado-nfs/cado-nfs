#ifndef POLYSELECT_MAIN_QUEUE_H_
#define POLYSELECT_MAIN_QUEUE_H_

#include <gmp.h>
#include <stdint.h>

#include "cado_poly.h"
#include "polyselect_main_data.h"

/* XXX This catch-all header is almost empty and obsolete to a large
 * extent now. Most of its content has gone to other places of the code.
 */

/* A few configuration flags first, which affect some specific points.
 */

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

extern void polyselect_fprintf_poly_pair(FILE * fp, mpz_srcptr N,                    mpz_poly_srcptr f, mpz_poly_srcptr g, int raw);

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_MAIN_QUEUE_H_ */
