#ifndef MPFQ_FAKE_HPP_
#define MPFQ_FAKE_HPP_

#include <gmp.h>
#include "gmp_aux.h"

#if defined(MPFQ_FAKE_HPP_) && defined(MPFQ_LAYER_H_)
#error "mpfq_layer.h and mpfq_fake.hpp are incomaptible"
#endif

typedef void * abdst_field;
typedef const void * absrc_field;
typedef void * abfield;
typedef unsigned long abelt[1];
typedef unsigned long * abdst_elt;
typedef const unsigned long * absrc_elt;
typedef unsigned long * abvec;
typedef unsigned long * abdst_vec;
typedef const unsigned long * absrc_vec;
#define abinit(ab, x)
#define abclear(ab, x)
#define abcxx_out(ab, os, x) ((os) << (*(x)))
#define abfscan(ab, f, x) (1)
static inline void abfield_init(abdst_field) {}
static inline void abfield_clear(abdst_field) {}
#define MPFQ_PRIME_MPZ 0
static inline void abfield_specify(abdst_field, ...) {}
static inline mpz_srcptr abfield_characteristic_srcptr(absrc_field) { return NULL; }
static inline size_t abvec_elt_stride(abdst_field, size_t k) { return (k / 64) * sizeof(mp_limb_t); }
static inline bool abis_zero(abdst_field, absrc_elt x) { return *x==0; }
static inline void abrandom(abdst_field, abdst_elt x, gmp_randstate_ptr rstate) { *x = gmp_urandomb_ui(rstate, 1); }

#endif	/* MPFQ_FAKE_HPP_ */
