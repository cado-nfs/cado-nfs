#ifndef POLYSELECT_LOCALS_H_
#define POLYSELECT_LOCALS_H_

#include <gmp.h>

#include "polyselect_main_data.h"
#include "polyselect_poly_header.h"
#include "polyselect_proots.h"
#include "polyselect_stats.h"
#include "dllist.h"

#ifdef __cplusplus
extern "C" {
#endif



struct polyselect_thread_locals_s {
    int idx;
    mpz_t ad;
    polyselect_stats stats;

    /* This pointer is const on purpose. We do _not_ want concurrent
     * write access to this. In particular, access to the global stats is
     * meant to be policed, and we mostly rely on
     * polyselect_main_data_commit_stats to do the write access (which
     * are lock-protected).
     */
    polyselect_main_data_srcptr main;
    
    /* This is just a copy (const pointer) of stats->rstate. Since we
     * have a random state there, let's just use it...
     */
    gmp_randstate_ptr rstate;

    /* used by several routines. Actually a replica of data which exists
     * elsewhere, mostly */
    polyselect_poly_header_t header;

    /* not entirely clear that this deserves here. I think that it could
     * probably live on the stack somehere in the collisions_ functions
     * now */
    polyselect_proots_t R;

    /* Are there any asynchronous tasks that this local thread in
     * particular has set aside, and that should be picked up for further
     * processing ?
     */
    struct dllist_head async_jobs;
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
