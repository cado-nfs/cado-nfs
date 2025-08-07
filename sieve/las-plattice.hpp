#ifndef CADO_LAS_PLATTICE_HPP
#define CADO_LAS_PLATTICE_HPP

#include <cstdint>

#include <limits>
#include <vector>

#include "macros.h"

#include "fb-types.hpp"
#include "gcd.h"        // for gcd_ul
#include "las-config.h" // for NOPROFILE_INLINE
#include "misc.h"       // for UMAX

// scan-headers: stop here

/*******************************************************************/
/********        Walking in p-lattices                    **********/

/*  p-lattice stuff */

#if 0 // this comment belongs somewhere else
// Compute the root r describing the lattice inside the q-lattice
// corresponding to the factor base prime (p,R).
// Formula: r = - (a1-R*b1)/(a0-R*b0) mod p
// Assumes p < 2^32

/* General version of the lattice transform function. Allows projective
   roots in input and output, and handles prime powers.
   In input, if the root is projective, say s/t (mod p) and t is
   non-invertible (mod p), then we expect R = p + (t/s mod p).
   On output, if the root is projective, say u/v (mod p) and v is
   non-invertible (mod p), then return value r = p + (v/u mod p).
   So projective roots are stored as their reciprocal, and have p added
   to signal the fact that it's a reciprocal value.
*/
#endif

/*
 * Algorithm by Franke and Kleinjung for lattice sieving of largish
 * primes.
 *
 * see also plattice.sage
 */

// Proposition 1 of [FrKl05]:
// Compute a basis <(i0, j0), (i1, j1)> of the p-lattice
// inside the q-lattice, such that
//    j0, j1 > 0
//    -I < i0 <= 0 <= i1 < I
//    i1-i0 >= I
//
// Sizes:
//    p is less than 32 bits and I fits easily in 32 bits.
//    So, i0 and j0 fit easily in 32 bits, since they are less than I
//    Now, i1 and j1 are also bounded by p, so 32 bits is enough
//    However: a and c can be as large as p*I (not both ?).
//    We still store them in 32 bits, since if they are larger, it means
//    that as soon as they are added to the offset for S, the index will
//    be out of range for S and the loop stops. Hence, this is safe to
//    replace a and c by a large value within 32 bits, when they are
//    larger than 32 bits.
//    Except that the choice of this ``large value'' requires some
//    caution. We need a value which can be used either for a or c, or
//    both, so that adding a, or c, or both, to a value within [0,IJ[ is
//    guaranteed to exceed IJ, but can't wrap around. Up to I=15, it's
//    rather easy. With the rescaling of J, at worst we may have IJ
//    within the range [2^29,2^30[. Thus if a and c are set to 2^30-1,
//    a.k.a INT32_MAX/2, then adding either, or the sum of both, to a
//    valid value x is guaranteed to be at most 3*2^30, which fits within
//    32 bits.
//    For I=16, it's much, much harder. Unless somebody comes up with a
//    nice idea, I see no way to avoid 64-bit arithmetic (which has some
//    cost, sure, but not terribly expensive). For consistency, we make
//    all data types for x, a, and c 64-bit in this case, and choose the
//    overflow constants as UINT32_MAX.

typedef uint64_t plattice_x_t;

class plattice_info; // IWYU pragma: keep

struct plattice_info_dense_t;

class plattice_enumerator;

/* Information on a p-lattice: lattice basis coordinates, and auxiliary
   values that can be derived directly from them, such as the constants
   used by the Franke-Kleinjung algorithm */
class plattice_info
{
    // vectors are (-mi0,j0) (the "warp" vector)
    //         and (i1,j1) (the "step" vector)
    friend struct plattice_info_dense_t;
    friend class plattice_enumerator;

  protected:
    uint32_t mi0; /* This encodes i0 = -mi0 */
    uint32_t j0;
    uint32_t i1;
    uint32_t j1;

  public:
    bool operator==(plattice_info const & o)
    {
        return mi0 == o.mi0 && i1 == o.i1 && j0 == o.j0 && j1 == o.j1;
    }
    bool operator!=(plattice_info const & o) { return !operator==(o); }

  public:
    int32_t get_i0() const { return -(int32_t)mi0; }
    int32_t get_i1() const { return i1; }
    uint32_t get_j0() const { return j0; }
    uint32_t get_j1() const { return j1; }

    /* the rule for advancing from (i,j) to the next is
     * (here imin = -I/2 )
     *
     * add the warp vector if (i - imin) >= mi0
     * add the step vector if (i - imin) < I - i1
     * add both in between.
     *
     * or, equivalenty:
     * add the warp vector if (i - imin) >= I - i1
     * add the step vector if (i - imin) < mi0
     *
     * Now in terms of x enconding:
     *
     * (i,j) is encoded as x = (i + I/2) + (j << logI)
     *
      if (i >= bound_warp)
        x += inc_warp;
      if (i < bound_step)
        x += inc_step;
     */
  protected:
    uint32_t get_bound_step(int const logI MAYBE_UNUSED) const
    {
        /* Note that if we're a vertical line, we have mi0 = 0, and
         * therefore we'll never have (i - imin) < mi0 and we'll
         * never add the step vector */
        return mi0;
    }

    uint32_t get_bound_warp(int const logI) const
    {
        if (i1 >> logI)
            return 0;
        return (1U << logI) - get_i1();
    }

    plattice_x_t get_inc_warp(int const logI) const
    {
        /* the warp step will add a negative integer to the i
         * coordinate, but that will often compensate with what is
         * present in x. At least, after adding warp or warp + step,
         * we won't get a carry.
         *
         */
        return ((uint64_t)get_j0() << logI) + (int64_t)get_i0();
    }

    plattice_x_t get_inc_step(int const logI) const
    {
        /* The case when i1 >= I can happen for vertical lines. This
         * implies mi0=0, and therefore we'll never ever trigger the
         * addition of the step vector. Bottom line: it's fine if we
         * put junk in here. We may even use that as a marker.
         */
        if (i1 >> logI)
            return UINT64_MAX;
        return ((uint64_t)get_j1() << logI) + (int64_t)get_i1();
    }

  public:
    uint32_t determinant() const { return mi0 * j1 + j0 * i1; };

    /* it seems that we no longer discard any vector, in fact */
    bool is_discarded() const { return mi0 == 0 && j0 == 0; }

    bool is_vertical_line(int const logI) const
    {
        // there are several possible shapes of vertical lines. the
        // most classical is
        // [ 0  1 ]
        // [ p  0 ] with p >= I ; IOW, root=0
        //
        // but more generally, any lattice that comes out of
        // reduce_plattice as:
        //
        // [ 0  p/g ]
        // [ g  *   ] with g >= I
        //
        // in this case, the "step" vector (g, *) will never be used,
        // and we'll only use the "warp" vector (0, p/g)
        return i1 >> logI; // the post-conditions imply mi0=0
        // another way to write it:
        // return get_inc_step() == UINT64_MAX;
    }

    bool is_projective_like(int const logI) const
    {
        // a projective root has a step vector of the form (*, 0).
        // here, we slightly extend this definition to also recognize
        // as "projective-like" a lattice whose step vector is (<I, 0)
        return j1 == 0 && !(i1 >> logI);
        // another way to write it:
        // return !(get_inc_step() >> logI);
    }

  protected:
    plattice_info()
        : mi0(0)
        , j0(0)
        , i1(0)
        , j1(0)
    {
    }

    void reduce_with_vertical_vector(uint32_t I)
    {
        /* At this point, (mi0,j0) represents itself, i.e. a vector with
         * two positive coordinates.
         */
        const uint32_t a = (I - 1) % mi0;
        uint32_t mi2 = (I - 1) - a - i1;
        const uint32_t j2 = j1;
        if (mi0 + mi2 < I)
            mi2 += mi0;
        i1 = mi0;
        j1 = j0;
        mi0 = mi2;
        j0 = j2;
    }

    void initial_basis(unsigned long const q, unsigned long const r, bool proj)
    {
        if (!proj) {
            mi0 = q;
            j0 = 0;
            i1 = r;
            j1 = 1;
        } else {
            const unsigned long gs = r;
            unsigned long t; // , h;
            const unsigned long g = xgcd_ul(&t, gs, q);
            mi0 = q / g;
            j0 = 0;
            i1 = t;
            j1 = g;
        }
    }

    /* If a lattice satisfies this condition, then its reduced form is
     * given by calling reduce_with_vertical_vector(I). Note that
     * calling reduce() works in all cases, though.
     *
     * Note also that this check is only valid for the initial basis
     * lattice.
     */
    bool needs_special_treatment(uint32_t I) const
    {
        return (i1 == 0 || (j1 > 1 && mi0 < I));
    }

#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
#include "las-reduce-plattice-production-asm.hpp"
#else
#include "las-reduce-plattice-simplistic.hpp"
#endif /* HAVE_GCC_STYLE_AMD64_INLINE_ASM */

    void reduce(uint32_t I)
    {
        /* We tried many different versions. See
         * https://gitlab.inria.fr/cado-nfs/cado-nfs/-/merge_requests/43
         * as well as tests/sieve/test-reduce-plattice.cpp and
         * the various contender implementations in
         * tests/sieve/reduce-plattice/
         */
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
        reduce_plattice_asm(I);
#else
        reduce_plattice_simplistic(I);
#endif
    }

  public:
    bool check_pre_conditions(uint32_t I) const
    {
        if (!(j0 == 0))
            return false;
        if (!(mi0 > 0))
            return false;
        if (!(j1 > 0))
            return false;
        // if (!(i1 >= 0)) return false;
        if (!(mi0 * j1 >= I))
            return false;
        return true;
    }

    bool check_post_conditions(uint32_t I) const
    {
        if (!(mi0 < I))
            return false;
        if (!(j0 > 0))
            return false;
        // since we want to handle the projective case, j1 may be zero
        // if (!(j1 >= 0)) return false;
        // this is an empty condition because of the type
        // if (!(0 <= i1)) return false;
        // This assertion is possibly violated for vertical lattices
        // i1 < I
        if (!(i1 < I || mi0 == 0))
            return false;
        if (!(i1 + mi0 >= I))
            return false;
        return true;
    }

    plattice_info(unsigned long const q, unsigned long const r, bool proj,
                  int const logI)
    {
        initial_basis(q, r, proj);
        /* At this point, (mi0,j0) represents itself, i.e. a vector with
         * two positive coordinates.
         * Note that j0==0
         */
        uint32_t I = UINT32_C(1) << logI;
        reduce(I);
    }
};

/* Like plattice_info, but remembers the offset of the factor base entry
   that generated each lattice basis. This offset becomes the "prime hint"
   in the bucket updates. */
struct plattice_sieve_entry : public plattice_info {
    slice_offset_t hint;
    plattice_sieve_entry(fbprime_t const p, fbroot_t const r, bool proj,
                         int const logI, slice_offset_t const hint)
        : plattice_info(p, r, proj, logI)
        , hint(hint) {};
};

// Compute 1/x mod m
// For m=2,3,6, this is trivial. Otherwise, has to be implemented
static inline uint32_t invmod(uint32_t x, uint32_t m)
{
    if (m == 2 || m == 3) {
        return x;
    } else if (m == 6) {
        ASSERT((x == 1) || (x == 5));
        return x;
    } else {
        ASSERT_ALWAYS(0);
    }
}

/* This is (currently) common across all instantiations of
 * plattice_enumerator. However if someday plattice_x_t changes type to
 * vary depending on the level (which, in my reckoning, is now doable),
 * then we should probably just move this function to the
 * plattice_enumerator class anyway (or, equivalently, make it a template
 * based on plattice_x_t).
 */
struct plattice_enumerator_base {
    static plattice_x_t starting_point(plattice_info const & pli,
                                       int const logI, sublat_t const & sublat)
    {
        int64_t I = int64_t(1) << logI;
        uint32_t m = sublat.m;

        // first FK vector a
        int64_t i0 = pli.get_i0();
        int64_t j0 = pli.get_j0();
        // second FK vector b
        int64_t i1 = pli.get_i1();
        int64_t j1 = pli.get_j1();

        // FIXME: should have a better understanding of those cases
        // (and not only in the sublat case, to be honest)
        // Right now, we desactivate them, by putting a starting point
        // above the limit.
        if ((i1 == 0) || (j1 == 0)) {
            return plattice_x_t(UMAX(plattice_x_t));
        }

        // If j0 or j1 reaches 20 bits, then it was saturated during the
        // dense storage. Let's skip this prime for which the FK basis is
        // very skewed. (only a few hits are missed)
        if ((j0 == ((1 << 20) - 1)) || (j1 == ((1 << 20) - 1))) {
            return plattice_x_t(UMAX(plattice_x_t));
        }

        // Look for alpha and beta such that
        //   alpha*a + beta*b == (sublat.i0,sublat.j0) mod m
        // This is a 2x2 system of determinant p, coprime to m.
        int64_t det =
            pli.determinant(); // det() returns the opposite of what we want
        det = (-det) % m;
        if (det < 0)
            det += m;
        det = invmod(det, m);
        int64_t al = (j1 * sublat.i0 - i1 * sublat.j0) % m;
        int64_t be = (-j0 * sublat.i0 + i0 * sublat.j0) % m;
        al = (al * det) % m;
        be = (be * det) % m;
        if (al < 0)
            al += m;
        if (be < 0)
            be += m;

        // Now, compute this potential starting point:
        int64_t ii = (al * i0 + be * i1 - sublat.i0) / m; // exact divisions
        int64_t jj = (al * j0 + be * j1 - sublat.j0) / m;

        // But here, ii might be beyond the bounds. So, let's fix.
        // It should be enough to subtract one of the FK vectors.
        // Note that a is the vector with negative abscissa.
        if (ii < -I / 2) {
            ASSERT(ii - i0 >= 0);
            ii -= i0;
            jj -= j0;
        } else if (ii > I / 2 - 1) {
            ASSERT(ii - i1 <= 0);
            ii -= i1;
            jj -= j1;
        }
        ASSERT((ii >= -I / 2) && (ii < I / 2));

        // But now, jj might be negative! So let's start the FK walk until we
        // go positive.
        while (jj < 0) {
            int64_t aux = ii;
            if (aux >= I / 2 - i1) {
                ii += i0;
                jj += j0;
            }
            if (aux < -I / 2 - i0) {
                ii += i1;
                jj += j1;
            }
        }

        // Now, (ii,jj) is the starting point we are looking for. Let's
        // convert it to the plattice_x_t type.
        plattice_x_t res = (ii + I / 2) + (jj << logI);
        return res;
    }
};

/* Class for enumerating lattice points with the Franke-Kleinjung algorithm */
class plattice_enumerator : public plattice_enumerator_base
{
  protected:
    // Maybe at some point, the plattice_x_t type could be templated in
    // order to have it 32 bits for non-top levels.
    plattice_x_t inc_warp, inc_step;
    uint32_t bound_step, bound_warp;
    slice_offset_t hint;
    plattice_x_t x;

  public:
    plattice_x_t get_inc_step() const { return inc_step; }

    struct fence {
        uint32_t maskI;
        plattice_x_t even_mask;
        plattice_x_t end;
        fence(int const logI, int J)
            : maskI((1U << logI) - 1U)
            , even_mask((plattice_x_t(1) << logI) | plattice_x_t(1))
            , end(plattice_x_t(J) << logI)
        {
        }
        fence(int const logI, int J, plattice_x_t cap)
            : fence(logI, J)
        {
            end = std::min(end, cap);
        }
    };

    plattice_enumerator(plattice_info const & basis, slice_offset_t const hint,
                        int const logI, sublat_t const & sublat)
        : inc_warp(basis.get_inc_warp(logI))
        , inc_step(basis.get_inc_step(logI))
        , bound_step(basis.get_bound_step(logI))
        , bound_warp(basis.get_bound_warp(logI))
        , hint(hint)
    {
        if (!sublat.m)
            x = plattice_x_t(1) << (logI - 1);
        else {
            x = plattice_enumerator_base::starting_point(basis, logI, sublat);
        }
    }

    plattice_enumerator(plattice_info const & basis, slice_offset_t const hint,
                        int const logI)
        : inc_warp(basis.get_inc_warp(logI))
        , inc_step(basis.get_inc_step(logI))
        , bound_step(basis.get_bound_step(logI))
        , bound_warp(basis.get_bound_warp(logI))
        , hint(hint)
        , x(plattice_x_t(1) << (logI - 1))
    {
    }

    // plattice_enumerator(const plattice_enumerator&) = default;

    /* This function is quite critical */
    void next(fence const & F)
    {
        const uint32_t i = x & F.maskI;
        if (i >= bound_warp)
            x += inc_warp;
        if (i < bound_step)
            x += inc_step;
    }

    /* Currently merely checks that not both are even */
    bool probably_coprime(fence const & F) const
    {
        return (x & F.even_mask) != 0;
    }

    bool done(fence const & F) const { return x >= F.end; }

    void advance_to_next_area(fence const & F)
    {
        if (x != std::numeric_limits<plattice_x_t>::max())
            x -= F.end;
    }

    int advance_to_end_of_row_or_smallest_region(uint32_t mask)
    {
        /* This function shouldn't be critical. We want the last matching
         * position on the row where x stands, or in the lowest-level
         * region if that happens to span less than a row.  This depends
         * on the "step" vector (i1,j1).  Note that we need this only for
         * projective primes (and powers), so that the corresponding
         * plattice has i1>0 and j1==0
         *
         * (for projective non-powers, we should have (i1,j1)=(1,0),
         * inc_step=1, and bound_warp=I-1. However we might have
         * something different for projective powers.)
         *
         * This returns the number of updates (not counting the current
         * x) that are skipped over by this update.
         */
#ifdef BUCKET_SIEVE_POWERS
        if (inc_step > 1) {
            int n = (mask - (x & mask)) / inc_step;
            x += n * inc_step;
            return n;
        }
#endif
        /* This should be the dominant case, even when bucket-sieving
         * powers. Hence the special-case above, so that we don't do a
         * division by 1 for all projective primes.
         */
        int n = (mask - (x & mask));
        x |= mask;
        return n;
    }
    plattice_x_t get_x() const { return x; }
    void set_x(plattice_x_t xx) { x = xx; }
    slice_offset_t get_hint() const { return hint; }

    void finish() { x = std::numeric_limits<plattice_x_t>::max(); }
    bool is_projective_like(int const logI) const
    {
        // see plattice_info::is_projective_like : step vector (<I, 0)
        return !(inc_step >> logI);
    }

    bool is_vertical_line(int const logI MAYBE_UNUSED) const
    {
        return inc_step == UINT64_MAX;
    }
};

/* Also enumerates lattice points, but probably_coprime() does a full gcd()
   to ensure that points are really coprime. Very slow. Not used, just there
   for experimenting.*/
class plattice_enumerator_coprime : public plattice_enumerator
{
    unsigned long u, v;
    typedef plattice_enumerator super;
    typedef typename super::fence fence;

  public:
    plattice_enumerator_coprime(plattice_info const & basis,
                                slice_offset_t const hint, int const logI,
                                sublat_t const & sublat)
        : plattice_enumerator(basis, hint, logI, sublat)
        , u(0)
        , v(0)
    {
    }
    void next(fence const & F)
    {
        uint32_t i = super::x & F.maskI;
        if (i >= super::bound_warp) {
            super::x += super::inc_warp;
            u++;
        }
        if (i < super::bound_step) {
            super::x += super::inc_step;
            v++;
        }
    }
    bool probably_coprime(plattice_enumerator::fence const & F) const
    {
        return (super::x & F.even_mask) != 0 && gcd_ul(u, v) == 1;
    }
};

class plattices_vector_t : public std::vector<plattice_enumerator>
{
    /* The index here is the global index, across all fb parts */
    slice_index_t index = -1;
    double weight = 0;

  public:
    plattices_vector_t() = default;
    plattices_vector_t(slice_index_t index, double weight)
        : index(index)
        , weight(weight)
    {
    }
    /* returns a global index */
    slice_index_t get_index() const { return index; };
    slice_index_t get_weight() const { return weight; };
};

/* Dense version of plattice_info and friends for long-term storage in
 * sublat mode. */

struct plattice_info_dense_t {
    uint32_t pack[3];
    // This pack of 96 bits is enough to contain
    //   mi0, i1, j0, j1
    // as 20-bit unsigned integers and
    //   hint
    // as a 16-bit integer.
    //
    // Note that mi0 and i1 are less than I, so this is ok, but
    // j0 and j1 could be larger. However, this is for very skewed
    // plattices, and we lose only a few hits by skipping those primes.
    // So we saturate them at 2^20-1 for later detection.
    //
    // uint16_t hint; // FIXME: this could be recovered for free...

    plattice_info_dense_t(plattice_info const & pli, uint16_t _hint)
    {
        uint32_t mi0;
        uint32_t i1;
        uint32_t j0;
        uint32_t j1;
        uint16_t hint;
        hint = _hint;
        // Handle orthogonal lattices (proj and r=0 cases)
        if (pli.i1 == 1 && pli.j1 == 0) {
            i1 = 1;
            j1 = 0;
            mi0 = UMAX(uint32_t);
            j0 = pli.j0;
        } else if (pli.i1 == 0 && pli.j1 == 1) {
            i1 = 0;
            j1 = 1;
            mi0 = UMAX(uint32_t);
            j0 = pli.j0;
        } else {
            // generic case: true FK-basis
            mi0 = pli.mi0;
            j0 = pli.j0;
            i1 = pli.i1;
            j1 = pli.j1;
        }
        constexpr uint32_t mask8 = (1 << 8) - 1;
        constexpr uint32_t mask16 = (1 << 16) - 1;
        constexpr uint32_t mask20 = (1 << 20) - 1;

        // Saturate skewed lattices, for later detection and skipping.
        j0 = std::min(j0, mask20);
        j1 = std::min(j1, mask20);

        pack[0] = (mi0 & mask20) | (i1 << 20);
        pack[1] = ((i1 >> 12) & mask8) | ((j0 & mask20) << 8) | (j1 << 28);
        pack[2] = ((j1 >> 4) & mask16) | (hint << 16);
    }

    plattice_info unpack(int const logI) const
    {
        plattice_info pli;
        constexpr uint32_t mask8 = (1 << 8) - 1;
        constexpr uint32_t mask16 = (1 << 16) - 1;
        constexpr uint32_t mask20 = (1 << 20) - 1;
        pli.mi0 = pack[0] & mask20;
        pli.i1 = (pack[0] >> 20) | ((pack[1] & mask8) << 12);
        pli.j0 = (pack[1] >> 8) & mask20;
        pli.j1 = (pack[1] >> 28) | ((pack[2] & mask16) << 4);

        // Orthogonal bases
        if (pli.i1 == 1 && pli.j1 == 0) {
            pli.mi0 = ((int32_t)1 << logI) - 1;
        } else if (pli.i1 == 0 && pli.j1 == 1) {
            pli.j0 = ((int32_t)1 << logI) + 1;
        }
        return pli;
    }

    uint16_t get_hint() const { return pack[2] >> 16; }
};

class plattices_dense_vector_t : public std::vector<plattice_info_dense_t>
{
  public:
#if !GNUC_VERSION(4, 7, 2)
    /* Apparently there's a bug in gcc 4.7.2, which does not recognize
     * that a vector of plattices_dense_vector_t is ok with a move ctor
     * and no copy ctor. I'm not entirely sure why, and I don't really
     * want to investigate. Theoretically, defining ctors as below should
     * be fine.
     */
    plattices_dense_vector_t(plattices_dense_vector_t const &) = delete;
    plattices_dense_vector_t(plattices_dense_vector_t &&) = default;
    plattices_dense_vector_t & operator=(plattices_dense_vector_t &&) = default;
    plattices_dense_vector_t() = default;
#endif
};

// std::vector<plattices_dense_vector_t> is for remembering the FK
// basis in sublat mode, between two different congruences of (i,j) mod
// m. For simplicity, we remember them only for the toplevel.
template <int> struct precomp_plattice_dense_t {
    typedef std::vector<plattices_dense_vector_t> type;
};

#endif /* CADO_LAS_PLATTICE_HPP */
