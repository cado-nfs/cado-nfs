#ifndef CADO_ROPT_SINGLE_SUBLATTICE_PRIORITY_QUEUE_H
#define CADO_ROPT_SINGLE_SUBLATTICE_PRIORITY_QUEUE_H

#include <stddef.h>
#include "gmp.h"

/* This represents a single sublattice modulo one p^e. Oddly enough, p is
 * not stored in the structure.
 *
 * The single_sublattice_info type is exposed here because we find it
 * convenient to have it available in the caller code as well. It's not
 * just an implementation detail.
 */
struct single_sublattice_info {
    /* val is the average p-valuation for u,v congruent to these values
     * mod p^e */
    float val;
    unsigned int u, v;
    char e;
#ifdef __cplusplus
    /* Now this is actually implementation details... */
    inline bool operator<(single_sublattice_info const & o) const {
        if (val < o.val) return true;
        if (val > o.val) return false;
        if (u < o.u) return true;
        if (u > o.u) return false;
        if (v < o.v) return true;
        if (v > o.v) return false;
        if (e < o.e) return true;
        if (e > o.e) return false;
        return false;
    }
#endif
};

/**
 * Priority queue for sublattices over a single p^e.
 * 
 * used in ropt_stage1.c
 *
 */
struct single_sublattice_priority_queue_s {
    void * impl;
};
typedef struct single_sublattice_priority_queue_s single_sublattice_priority_queue[1];
typedef struct single_sublattice_priority_queue_s * single_sublattice_priority_queue_ptr;
typedef const struct single_sublattice_priority_queue_s * single_sublattice_priority_queue_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

void single_sublattice_priority_queue_init(single_sublattice_priority_queue_ptr q,
                         size_t max_len);


void single_sublattice_priority_queue_push(single_sublattice_priority_queue_ptr q, const struct single_sublattice_info * z);

void single_sublattice_priority_queue_clear (single_sublattice_priority_queue_ptr q);
size_t single_sublattice_priority_queue_size (single_sublattice_priority_queue_srcptr q);
int single_sublattice_priority_queue_empty (single_sublattice_priority_queue_srcptr q);

void single_sublattice_priority_queue_pop_all(single_sublattice_priority_queue_srcptr q, struct single_sublattice_info * out);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_ROPT_SINGLE_SUBLATTICE_PRIORITY_QUEUE_H */
