#ifndef ARITH_GENERIC_HPP_
#define ARITH_GENERIC_HPP_

#include <stddef.h>
#include <string>
#include <ostream>
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "arith-concrete-base.hpp"

/* The goal of this file is to provide a level of abstraction in order to
 * call the field-specific routines via virtual bases.
 *
 * We need to decide whether we want a virtual hierarchy with only
 * (simili-) static functions, operating on objects with
 * instance-specific meaning. The other option would be to have the
 * function layer on the objects themselves. In both cases, we would have
 * several fragile points.
 */

/* matmul_init is one of the main code multiplexing points.  The
 * instantiation of the class below is a derived class that also inherits
 * from a concrete class. Internally, the low-level implementation will
 * only see the concrete class (with some wild type-punning).
 */

struct arith_generic {
    // virtual arith_vbase& operator=(arith_vbase const &) = 0;
    // add
    // sub
    // addmul_ui
    // submul_ui
    // reduce
    
    struct elt : public arith_concrete_base::elt {};

    virtual void vec_add_and_reduce(elt * dst, elt const * b, size_t n) const = 0;
    virtual void vec_set(elt *, elt const *, size_t) const = 0;
    virtual void vec_neg(elt *, elt const *, size_t) const = 0;
    virtual void vec_set_zero(elt *, size_t) const = 0;
    virtual bool vec_is_zero(elt const *, size_t) const = 0;
    virtual void vec_set_random(elt *, size_t, gmp_randstate_ptr) const = 0;
    virtual void vec_add_dotprod(elt &, elt const *, elt const *, size_t) const = 0;
    virtual void vec_addmul_and_reduce(elt *, elt const *, elt const &, size_t) const = 0;
    virtual int vec_cmp(elt const * a, elt const * b, size_t k) const = 0;
    virtual int cmp(elt const & a, elt const & b) const = 0;
    virtual int cmp(elt const & a, unsigned long b) const = 0;

    /* We use elt * as a vector type, but it should be understood as
     * something opaque, really.
     */
    virtual elt * vec_subvec(elt *, size_t) const = 0;
    virtual elt const * vec_subvec(elt const *, size_t) const = 0;

    inline elt & vec_item(elt * p, size_t k) { return *vec_subvec(p, k); }
    inline elt const & vec_item(elt const * p, size_t k) { return *vec_subvec(p, k); }
    virtual size_t vec_elt_stride(size_t) const = 0;
    inline size_t elt_stride() { return vec_elt_stride(1); }
    virtual bool is_zero(elt const &) const = 0;
    virtual void simd_set_ui_at(elt &, size_t, int) const = 0;
    virtual void simd_add_ui_at(elt &, size_t, int) const = 0;
    virtual int simd_hamming_weight(elt const &) const = 0;
    virtual void vec_simd_set_ui_at(elt * p, size_t k, int v) const = 0;
    virtual void vec_simd_add_ui_at(elt * p, size_t k, int v) const = 0;
    virtual int vec_simd_hamming_weight(elt const * p, size_t n) const = 0;
    virtual int vec_simd_find_first_set(elt &, elt const * p, size_t n) const = 0;
    virtual std::ostream& cxx_out(std::ostream&, elt const &) const = 0;
    virtual elt * alloc(size_t = 1, size_t = 64) const = 0;
    virtual elt * realloc(elt *, size_t, size_t, size_t = 64) const = 0;
    virtual void free(elt *) const = 0;
    virtual void set(elt &, elt const &) const = 0;
    virtual void neg(elt &, elt const &) const = 0;
    virtual void set(elt &, cxx_mpz const &) const = 0;
    virtual void add_and_reduce(elt &, elt const &) const = 0;
    virtual void sub_and_reduce(elt &, elt const &) const = 0;
    virtual void reduce(elt &) const = 0;

    virtual std::string impl_name() const = 0;
    virtual size_t simd_groupsize() const = 0;
    virtual bool is_characteristic_two() const = 0;
    virtual mpz_srcptr characteristic() const = 0;
    virtual arith_concrete_base const * concrete() const = 0;
    virtual arith_concrete_base * concrete() = 0;

    virtual ~arith_generic() = default;

    static arith_generic * instance(mpz_srcptr p, int simd_groupsize);
};

/* the counterpart would be arith-hard, probably (and arith-modp is
 * "hard").
 *
 * we need, for the binary case: add set set_zero vec_set_zero vec_init,
 * and the actual type for contiguous storage.
 */

#endif	/* ARITH_GENERIC_HPP_ */
