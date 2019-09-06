#ifndef MPFQ_FAKE_HPP_
#define MPFQ_FAKE_HPP_

#if defined(MPFQ_FAKE_HPP_) && defined(MPFQ_LAYER_H_)
#error "mpfq_layer.h and mpfq_fake.hpp are incomaptible"
#endif

typedef void * abdst_field;
typedef const void * absrc_field;
typedef void * abfield;
typedef unsigned long * abdst_elt;
typedef const unsigned long * absrc_elt;
typedef unsigned long * abvec;
typedef unsigned long * abdst_vec;
typedef const unsigned long * absrc_vec;
static inline void abfield_init(abdst_field) {}
static inline void abfield_clear(abdst_field) {}
#define MPFQ_PRIME_MPZ 0
static inline void abfield_specify(abdst_field, ...) {}
static inline mpz_srcptr abfield_characteristic_srcptr(absrc_field) { return NULL; }

#endif	/* MPFQ_FAKE_HPP_ */
