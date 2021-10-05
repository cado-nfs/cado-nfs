#ifndef POLYSELECT_LOCALS_H_
#define POLYSELECT_LOCALS_H_

#include <gmp.h>

#include "polyselect_main_data.h"
#include "polyselect_poly_header.h"
#include "polyselect_proots.h"
#include "polyselect_stats.h"

#ifdef __cplusplus
extern "C" {
#endif



struct polyselect_thread_locals_s {
    int idx;
    mpz_t ad;
    polyselect_stats stats;
    polyselect_main_data_ptr main;
    
    /* This is just a copy (const pointer) of stats->rstate. Since we
     * have a random state there, let's just use it...
     */
    gmp_randstate_ptr rstate;

    polyselect_poly_header_t header;
    polyselect_proots_t R;
};

typedef struct polyselect_thread_locals_s polyselect_thread_locals[1];
typedef struct polyselect_thread_locals_s * polyselect_thread_locals_ptr;
typedef const struct polyselect_thread_locals_s * polyselect_thread_locals_srcptr;


extern void polyselect_thread_locals_init(polyselect_thread_locals_ptr loc, polyselect_main_data_ptr main, int idx);
extern void polyselect_thread_locals_clear(polyselect_thread_locals_ptr loc);


#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_LOCALS_H_ */
