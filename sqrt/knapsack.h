#ifndef CADO_KNAPSACK_H
#define CADO_KNAPSACK_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*knapsack_object_callback_t)(void *, unsigned long v, int64_t x);

struct knapsack_object_s {
    const int64_t * tab;
    unsigned int stride;
    unsigned int offset;
    unsigned int nelems;
    int64_t bound;
    knapsack_object_callback_t cb;
    void * cb_arg;
};

typedef struct knapsack_object_s knapsack_object[1];
typedef struct knapsack_object_s * knapsack_object_ptr;
typedef const struct knapsack_object_s * knapsack_object_srcptr;

/* This interface solves a modular knapsack with the standard O(2^(n/2))
 * time, O(2^(n/2)) memory algorithm (one may do better both in time and
 * memory).
 *
 * Solutions must have all coefficients in [-1,+1], and the combination
 * must fit within the interval [-bound,+bound]
 *
 * The stride value exists because we also allow the possibility for the
 * knapsack to be multi-dimensional.
 *
 * there's a trivial factor of two to be gained by constraining one
 * coordinate -- however I keep the stuff this way for a while, because
 * wraparounds are sometimes tricky, and may lead to missing solutions.
 */
extern int knapsack_solve(knapsack_object_ptr ks);

/* these are almost empty placeholders. the caller must proved the tab
 * pointer himself */
void knapsack_object_init(knapsack_object_ptr);
void knapsack_object_clear(knapsack_object_ptr);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_KNAPSACK_H */
