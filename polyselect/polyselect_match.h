#ifndef CADO_POLYSELECT_MATCH_H
#define CADO_POLYSELECT_MATCH_H

#include <stdint.h>
#include "polyselect_poly_header.h"
#include "dllist.h"
#include "polyselect_thread.h"

#ifdef __cplusplus
extern "C" {
#endif

struct polyselect_match_info_s {
    polyselect_poly_header_t header;
    unsigned long p1;
    unsigned long p2;
    int64_t i;
    uint64_t q;
    mpz_t rq;
    struct dllist_head queue;
};
typedef struct polyselect_match_info_s polyselect_match_info_t[1];
typedef struct polyselect_match_info_s * polyselect_match_info_ptr;
typedef const struct polyselect_match_info_s * polyselect_match_info_srcptr;

extern void polyselect_match_info_clear(polyselect_match_info_ptr task);
extern void polyselect_match_info_init_trivial(polyselect_match_info_ptr task);
extern void polyselect_match_info_init(polyselect_match_info_ptr job, unsigned long p1, unsigned long p2, const int64_t i, uint64_t q, mpz_srcptr rq, polyselect_thread_srcptr loc);

extern void polyselect_match_info_set(polyselect_match_info_ptr job, unsigned long p1, unsigned long p2, const int64_t i, uint64_t q, mpz_srcptr rq, polyselect_thread_srcptr loc);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_POLYSELECT_MATCH_H */
