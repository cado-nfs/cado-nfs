#ifndef ARITH_CROSS_HPP_
#define ARITH_CROSS_HPP_

#include "arith-generic.hpp"

/* This structure provides abstract access to operations that deal with
 * two distinct SIMD contexts (same _field_, but different SIMD group
 * sizes). This is relevant mostly for GF(2), and then the code
 * eventually falls back to routines that are defined in bblas
 */
struct arith_cross_generic {
    static arith_cross_generic * instance(arith_generic * A0, arith_generic * A1);

    /*
     * This computes w = transpose(u)*v, following the conventions
     * detailed here.
     *
     * Let g0 = A0->simd_groupsize() and g1 = A1->simd_groupsize()
     *
     * w is a g0*g1 matrix (a sequence of g0 elements wrt A1)
     * u is a n*g0 matrix (a sequence of n elements wrt A0)
     * v is a n*g1 matrix (a sequence of n elements wrt A1)
     */
    virtual void add_dotprod(arith_generic::elt * w, arith_generic::elt const * u, arith_generic::elt const * v, unsigned int n) const = 0;

    /*
     * This computes w += u*v, following the conventions detailed here.
     *
     * Let g0 = A0->simd_groupsize() and g1 = A1->simd_groupsize()
     *
     * w is an n*g1 matrix (a sequence of n elements wrt A1)
     * u is an n*g0 matrix (a sequence of n elements wrt A0)
     * v is a g0*g1 matrix (a sequence of g0 elements wrt A1).
     */
    virtual void addmul_tiny(arith_generic::elt * w, arith_generic::elt const * u, arith_generic::elt const * v, unsigned int n) const = 0;

    /*
     * This computes w = transpose(u), following the conventions
     * detailed here.
     *
     * Let g0 = A0->simd_groupsize() and g1 = A1->simd_groupsize()
     *
     * u is a g0*g1 matrix (a sequence of g0 elements wrt A1)
     * w is a g1*g0 matrix (a sequence of g1 elements wrt A0)
     */
    virtual void transpose(arith_generic::elt * w, arith_generic::elt const * u) const = 0;

    virtual ~arith_cross_generic() = default;
};


#endif	/* ARITH_CROSS_HPP_ */
