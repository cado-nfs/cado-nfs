#ifndef CADO_ROPT_SUBLATTICE_PRIORITY_QUEUE_HPP
#define CADO_ROPT_SUBLATTICE_PRIORITY_QUEUE_HPP

#include <stddef.h>
#include "gmp.h"

/**
 * Priority queue for sublattices over a product of p^e.
 * 
 * used in ropt_stage1.c
 *
 */
struct sublattice_priority_queue_s {
    void * impl;
};
typedef struct sublattice_priority_queue_s sublattice_priority_queue[1];
typedef struct sublattice_priority_queue_s * sublattice_priority_queue_ptr;
typedef const struct sublattice_priority_queue_s * sublattice_priority_queue_srcptr;

void sublattice_priority_queue_init(sublattice_priority_queue_ptr q,
                         size_t max_len);

void sublattice_priority_queue_push(sublattice_priority_queue_ptr q,
                            mpz_srcptr u,
                            mpz_srcptr v,
                            mpz_srcptr mod, 
                            float val );

void sublattice_priority_queue_clear (sublattice_priority_queue_ptr q);
size_t sublattice_priority_queue_size (sublattice_priority_queue_srcptr q);
int sublattice_priority_queue_empty (sublattice_priority_queue_srcptr q);

void sublattice_priority_queue_do(sublattice_priority_queue_srcptr q, void (*f)(mpz_srcptr u, mpz_srcptr v, mpz_srcptr modulus, void * arg), void *arg);

#endif	/* CADO_ROPT_SUBLATTICE_PRIORITY_QUEUE_HPP */
