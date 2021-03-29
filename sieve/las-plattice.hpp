#ifndef LAS_PLATTICE_H
#define LAS_PLATTICE_H

#include <cstdint>
#include <vector>

#include "macros.h"

#include "fb-types.h"
#include "las-config.h"  // for NOPROFILE_INLINE
#include "misc.h"       // for UMAX
#include "gcd.h"       // for gcd_ul

/*******************************************************************/
/********        Walking in p-lattices                    **********/

/*  p-lattice stuff */


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


/*
 * Algorithm by Franke and Kleinjung for lattice sieving of largish
 * primes.
 *
 * see also plattice.sage
 */

// Proposition 1 of [FrKl05]:
// Compute a basis <(a0, a1), (b0, b1)> of the p-lattice
// inside the q-lattice, such that
//    a1, b1 > 0
//    -I < a0 <= 0 <= b0 < I
//    b0-a0 >= I
//
// Sizes:
//    p is less than 32 bits and I fits easily in 32 bits.
//    So, a0 and a1 fit easily in 32 bits, since they are less than I
//    Now, b0 and b1 are also bounded by p, so 32 bits is enough
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
//

// OBSOLETE: the large value is now UINT64MAX>>2. So, plattice_x_t must
// be 64 bits. We could deal until logI <= 31 when logJ = logI -1.


// Return value:
// * non-zero if everything worked ok
// * zero when the algorithm failed. This can happen when p is a prime power,
//   and g, gcd(p,r) >= I, since then the subtractive Euclidean algorithm will
//   yield (a0=g, b0=0) at some point --- or the converse --- and the loop
//   while (|a0| >= I) a0 += b0 will loop forever.
//
// Note that on a c166 example, this code alone accounts for almost 20%
// of the computation time.


typedef uint64_t plattice_x_t;

struct plattice_info_t;  // IWYU pragma: keep
static int reduce_plattice (plattice_info_t *, fbprime_t, fbroot_t, uint32_t);


/* Information on a p-lattice: lattice basis coordinates, and auxiliary
   values that can be derived directly from them, such as the constants
   used by the Franke-Kleinjung algorithm */
struct plattice_info_t {
  static const int shift = 0;

  // vectors are a=(a0,a1) and b=(b0,b1)
  int32_t a0;
  uint32_t a1;
  int32_t b0;
  uint32_t b1;

  plattice_info_t(const fbprime_t p, const fbroot_t r, const bool proj, const int logI) {
    /* placate compiler. Sure, we shouldn't need it really, but
     * admittedly the control flow through reduce_plattice is contrived.
     */
    a0 = a1 = b0 = b1 = 0;
    if (UNLIKELY(proj && r == 0)) {
      /* This lattice basis might work in principle, but it generates hits in
         all locations i=1, ..., I/2-1, j = 0, all of which are useless except
         i=1, j=0.  */
      // Note: if J>p, we are missing some hits, here.
      a1 = p;
      a0 = -((int32_t)1 << logI) + 1;
      b1 = 0;
      b0 = 1;
    } else if (UNLIKELY(!proj && r == 0)) {
      a1 = ((int32_t)1 << logI) + 1;
      a0 = -p; /* Thus bound0 = p > I, inc_c is always added */
      b1 = 1;
      b0 = 0; /* Thus bound1 = I - b0 = I, inc_a is never added */
    } else {
      ASSERT_ALWAYS(!proj);
      /* One essential condition of proposition 1 in FrKl05 is that p >= I */
      ASSERT(p >> logI);
      int rc = reduce_plattice (this, p, r, 1U << logI);
      /* These are the constraints that reduce_plattice should meet */
      ASSERT(-(INT32_C(1) << logI) < a0);
      ASSERT(a0 <= 0);
      ASSERT(b0 >= 0);
      ASSERT(b0 < (INT32_C(1) << logI));
      ASSERT((b0-a0) >= (INT32_C(1) << logI));
      ASSERT(a1 > 0);
      ASSERT(b1 > 0);
      if (UNLIKELY(rc == 0)) {
        /* gcd(r, p) > I. We currently can't handle this case. Set everything
           to zero to signal to calling code that this is not a valid basis.
           Is there a better way? Exception, maybe? */
        a0 = a1 = b0 = b1 = 0;
      }
    }
  }

  plattice_info_t(int32_t aa0, uint32_t aa1, int32_t bb0, uint32_t bb1) {
      a0 = aa0;
      a1 = aa1;
      b0 = bb0;
      b1 = bb1;
  }

  /* Return the four coordinates, but multiplied by 2 if we use mod 2 classes */
  int32_t get_a0() const {return a0 << shift;}
  uint32_t get_a1() const {return a1 << shift;}
  int32_t get_b0() const {return b0 << shift;}
  uint32_t get_b1() const {return b1 << shift;}

  uint32_t get_bound0(const int logI MAYBE_UNUSED) const {
      return -get_a0();
  }

  uint32_t get_bound1(const int logI) const {
      return (1U << logI) - get_b0();
  }

  plattice_x_t get_inc_a(const int logI) const {
      return ((uint64_t)get_a1() << logI) + (int64_t)get_a0();
  }

  plattice_x_t get_inc_c(const int logI) const {
    return ((uint64_t)get_b1() << logI) + (int64_t)get_b0();
  }

  uint32_t det() const {return (-a0)*b1 + a1*b0;};
};

// Original version of reduce_plattice, always with a division for each step.
// This version is pretty fast on new powerful processors, but slow on others.
// I keep it because in fews years, if the int32_t division is faster, it's
// possible this code might be the fastest.
NOPROFILE_INLINE int
reduce_plattice (plattice_info_t *pli, const fbprime_t p, const fbroot_t r, uint32_t I)
{
  int32_t a0 = -((int32_t) p), b0 = (int32_t) r, a1 = 0, b1 = 1, k;
  const int32_t hI = (int32_t) I;
  const int32_t mhI = -hI;
  while (LIKELY(b0 >= hI)) {
    k = a0 / b0; a0 %= b0; a1 -= k * b1;
    if (UNLIKELY(a0 > mhI)) break;
    k = b0 / a0; b0 %= a0; b1 -= k * a1;
    /* We may conceivably unroll a bit more, or a bit less, here. Just
     * tuck in as many copies of the following block as you wish. */
#if 1
    if (UNLIKELY(b0 < hI )) break;
    k = a0 / b0; a0 %= b0; a1 -= k * b1;
    if (UNLIKELY(a0 > mhI)) break;
    k = b0 / a0; b0 %= a0; b1 -= k * a1;
#endif
  }
  k = b0 - hI - a0;
  if (b0 > -a0) {
    if (UNLIKELY(!a0)) return 0;
    k /= a0; b0 -= k * a0; b1 -= k * a1;
  } else {
    if (UNLIKELY(!b0)) return 0;
    k /= b0; a0 += k * b0; a1 += k * b1;
  }
  pli->a0 = (int32_t) a0; pli->a1 = (uint32_t) a1; pli->b0 = (int32_t) b0; pli->b1 = (uint32_t) b1;
  return 1;
}


/* Like plattice_info_t, but remembers the offset of the factor base entry
   that generated each lattice basis. This offset becomes the "prime hint"
   in the bucket updates. */
struct plattice_sieve_entry : public plattice_info_t {
  slice_offset_t hint;
  plattice_sieve_entry(const fbprime_t p, const fbroot_t r, const bool proj, const int logI, const slice_offset_t hint)
     : plattice_info_t(p, r, proj, logI), hint(hint) {};
};

// Compute 1/x mod m
// For m=2,3,6, this is trivial. Otherwise, has to be implemented
static inline uint32_t invmod(uint32_t x, uint32_t m) {
    if (m==2 || m==3) {
        return x;
    } else if (m==6) {
        ASSERT((x==1) || (x==5));
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
static plattice_x_t starting_point(const plattice_info_t &pli,
        const int logI, const sublat_t &sublat) {
    int64_t I = int64_t(1)<<logI;
    uint32_t m  = sublat.m;
    uint32_t i0 = sublat.i0;
    uint32_t j0 = sublat.j0;

    // first FK vector a
    int64_t a0 = pli.get_a0();
    int64_t a1 = pli.get_a1();
    // second FK vector b
    int64_t b0 = pli.get_b0();
    int64_t b1 = pli.get_b1();

    // FIXME: should have a better understanding of those cases
    // (and not only in the sublat case, to be honest)
    // Right now, we desactivate them, but putting a starting point
    // higher than the limit.
    if ((b0 == 0) || (b1 == 0)) {
        return plattice_x_t(UMAX(plattice_x_t));
    }

    // If a1 or b1 reaches 20 bits, then it was saturated during the
    // dense storage. Let's skip this prime for which the FK basis is
    // very skewed. (only a few hits are missed)
    if ((a1 == ((1<<20)-1)) || (b1 == ((1<<20)-1))) {
        return plattice_x_t(UMAX(plattice_x_t));
    }

    // Look for alpha and beta such that
    //   alpha*a + beta*b == (i0,j0) mod m
    // This is a 2x2 system of determinant p, coprime to m.
    int64_t det = pli.det(); // det() returns the opposite of what we want
    det = (-det) % m;
    if (det < 0)
        det += m;
    det = invmod(det, m);
    int64_t al = ( b1*i0 - b0*j0) % m;
    int64_t be = (-a1*i0 + a0*j0) % m;
    al = (al*det) % m;
    be = (be*det) % m;
    if (al < 0)
        al += m;
    if (be < 0)
        be += m;

    // Now, compute this potential starting point:
    int64_t ii = (al*a0 + be*b0 - i0) / m; // exact divisions
    int64_t jj = (al*a1 + be*b1 - j0) / m;

    // But here, ii might be beyond the bounds. So, let's fix.
    // It should be enough to subtract one of the FK vectors.
    // Note that a is the vector with negative abscissa.
    if (ii < -I/2) {
        ASSERT(ii - a0 >= 0);
        ii -= a0;
        jj -= a1;
    } else if (ii > I/2-1) {
        ASSERT(ii - b0 <= 0);
        ii -= b0;
        jj -= b1;
    }
    ASSERT ((ii >= -I/2) && (ii < I/2));

    // But now, jj might be negative! So let's start the FK walk until we
    // go positive.
    while (jj < 0) {
        int64_t aux = ii;
        if (aux >= I/2 - b0) {
            ii += a0;
            jj += a1;
        }
        if (aux < -I/2 - a0) {
            ii += b0;
            jj += b1;
        }
    }

    // Now, (ii,jj) is the starting point we are looking for. Let's
    // convert it to the plattice_x_t type.
    plattice_x_t res = (ii+I/2) + (jj<<logI);
    return res;
}
};

/* Class for enumerating lattice points with the Franke-Kleinjung algorithm */
template<int LEVEL>
class plattice_enumerator : public plattice_enumerator_base {
protected:
    // Maybe at some point, the plattice_x_t type could be templated in
    // order to have it 32 bits for non-top levels.
    plattice_x_t inc_a, inc_c;
    uint32_t bound0, bound1;
    slice_offset_t hint;
    plattice_x_t x;

public:

    struct fence {
        uint32_t maskI;
        plattice_x_t even_mask;
        plattice_x_t end;
        fence(const int logI, int J) {
            maskI = (1U << logI) - 1U;
            even_mask = (plattice_x_t(1) << logI) | plattice_x_t(1);
            end = plattice_x_t(J) << logI;
        }
        fence(const int logI, int J, plattice_x_t cap) : fence(logI, J) {
            if (end >= cap) end = cap;
        }
    };

    plattice_enumerator(const plattice_info_t &basis,
            const slice_offset_t hint, const int logI, const sublat_t &sublat)
        : hint(hint)
    {
        inc_a = basis.get_inc_a(logI);
        inc_c = basis.get_inc_c(logI);
        bound0 = basis.get_bound0(logI);
        bound1 = basis.get_bound1(logI);
        if (!sublat.m) 
            x = plattice_x_t(1) << (logI-1);
        else {
            x = plattice_enumerator_base::starting_point(basis, logI, sublat);
        }
    }

    plattice_enumerator(const plattice_info_t &basis,
            const slice_offset_t hint, const int logI)
        : hint(hint)
    {
        inc_a = basis.get_inc_a(logI);
        inc_c = basis.get_inc_c(logI);
        bound0 = basis.get_bound0(logI);
        bound1 = basis.get_bound1(logI);
        x = plattice_x_t(1) << (logI-1);
    }


    plattice_enumerator(const plattice_enumerator&) = default;

    /* This function is quite critical */
    void next(fence const & F) {
      uint32_t i = x & F.maskI;
      if (i >= bound1)
        x += inc_a;
      if (i < bound0)
        x += inc_c;
    }

    /* Currently merely checks that not both are even */
    bool probably_coprime(fence const & F) const {return (x & F.even_mask) != 0;}

    inline bool done(fence const & F) { return x >= F.end; }
    void advance_to_next_area(fence const & F) { x -= F.end; }

    inline void advance_to_end_of_projective_first_line(fence const & F)
    {
        /* This function is not critical at all. We want the last
         * matching position on the line (i,0). This depends on the
         * (b0,b1) vector, and works **ONLY** in the projective case, and
         * **ONLY** while we're on the first line !
         *
         * for projective non-powers, we should have (b0,b1)=(1,0),
         * inc_c=1, and bound1=I-1. However we might have something
         * different for projective powers. Presently, powers are not
         * bucket-sieved anyway, so there's little point in bothering.
         * (see "Shall we enable bucket-sieving for powers" in
         * las-fill-in-buckets.cpp)
         */
        x = F.maskI;
        ASSERT(inc_c == 1);
        ASSERT(bound1 == F.maskI);
    }
    plattice_x_t get_x() const {return x;}
    void set_x(plattice_x_t xx) {x = xx;}
    plattice_x_t get_bound1() const {return bound1;}
    plattice_x_t get_inc_c() const {return inc_c;}
    slice_offset_t get_hint() const {return hint;}
};

/* Also enumerates lattice points, but probably_coprime() does a full gcd()
   to ensure that points are really coprime. Very slow. Not used, just there
   for experimenting.*/
template<int LEVEL>
class plattice_enumerator_coprime : public plattice_enumerator<LEVEL> {
  unsigned long u, v;
  typedef plattice_enumerator<LEVEL> super;
  typedef typename super::fence fence;
public:
  plattice_enumerator_coprime(const plattice_info_t &basis,
          const slice_offset_t hint, const int logI, const sublat_t &sublat)
    : plattice_enumerator<LEVEL>(basis, hint, logI, sublat), u(0), v(0) {}
  void next(fence const & F) {
    uint32_t i = super::x & F.maskI;
    if (i >= super::bound1) { super::x += super::inc_a; u++; }
    if (i < super::bound0) { super::x += super::inc_c; v++; }
  }
  bool probably_coprime(typename plattice_enumerator<LEVEL>::fence const & F) const {
    return (super::x & F.even_mask) != 0 && gcd_ul(u, v) == 1;
  }
};

template<int LEVEL>
class plattices_vector_t:
        public std::vector<plattice_enumerator<LEVEL>> {
    /* The index here is the global index, across all fb parts */
    slice_index_t index;
    double weight;
public:
    plattices_vector_t() = default;
    plattices_vector_t(slice_index_t index, double weight) : index(index), weight(weight) {}
    /* returns a global index */
    slice_index_t get_index() const {return index;};
    slice_index_t get_weight() const {return weight;};
};

/* Dense version of plattice_info_t and friends for long-term storage in
 * sublat mode. */

/* This can now be a template. We don't use it, so far.
 */
template<int /* LEVEL */>
struct plattice_info_dense_t {
    uint32_t pack[3];
    // This pack of 96 bits is enough to contain
    //   minus_a0, b0, a1, b1
    // as 20-bit unsigned integers and
    //   hint
    // as a 16-bit integer.
    //
    // Note that minus_a0 and b0 are less than I, so this is ok, but
    // a1 and b1 could be larger. However, this is for very skewed
    // plattices, and we lose only a few hits by skipping those primes.
    // So we saturate them at 2^20-1 for later detection.
    //
    // uint16_t hint; // FIXME: this could be recovered for free...

    plattice_info_dense_t(const plattice_info_t & pli, uint16_t _hint) {
        uint32_t minus_a0;
        uint32_t b0; 
        uint32_t a1;
        uint32_t b1;
        uint16_t hint;
        hint = _hint;
        // Handle orthogonal lattices (proj and r=0 cases)
        if (pli.b0 == 1 && pli.b1 == 0) {
            b0 = 1;
            b1 = 0;
            minus_a0 = UMAX(uint32_t);
            a1 = pli.a1;
        } else if (pli.b0 == 0 && pli.b1 == 1) {
            b0 = 0;
            b1 = 1;
            minus_a0 = UMAX(uint32_t);
            a1 = pli.a1;
        } else {
            // generic case: true FK-basis
            ASSERT(pli.b0 >= 0);
            ASSERT(pli.a0 <= 0);
            minus_a0 = -pli.a0;
            a1 = pli.a1;
            b0 = pli.b0;
            b1 = pli.b1;
        }
        uint32_t mask8  = (1<<8)-1;
        uint32_t mask16 = (1<<16)-1;
        uint32_t mask20 = (1<<20)-1;

        // Saturate skewed lattices, for later detection and skipping.
        if (a1 > mask20)
            a1 = mask20;
        if (b1 > mask20)
            b1 = mask20;

        pack[0] = (minus_a0 & mask20) | (b0 << 20);
        pack[1] = ((b0 >> 12) & mask8 ) | ((a1 & mask20) << 8) | (b1 << 28);
        pack[2] = ((b1 >> 4) & mask16) | (hint << 16);
    }

    plattice_info_t unpack(const int logI) const {
        uint32_t minus_a0;
        uint32_t b0; 
        uint32_t a1;
        uint32_t b1;
        
        uint32_t mask8 = (1<<8)-1;
        uint32_t mask16 = (1<<16)-1;
        uint32_t mask20 = (1<<20)-1;
        minus_a0 = pack[0] & mask20;
        b0 = (pack[0] >> 20) | ((pack[1] & mask8) << 12);
        a1 = (pack[1] >> 8) & mask20;
        b1 = (pack[1] >> 28) | ((pack[2] & mask16) << 4);

        plattice_info_t pli(-int32_t(minus_a0), uint32_t(a1),
        int32_t(b0), uint32_t(b1));
        // Orthogonal bases
        if (pli.b0 == 1 && pli.b1 == 0) {
            pli.a0 = -((int32_t)1 << logI) + 1;
        } else if (pli.b0 == 0 && pli.b1 == 1) {
            pli.a1 = ((int32_t)1 << logI) + 1;
        }
        return pli;
    }

    uint16_t get_hint() const {
        return pack[2] >> 16;
    }

};

template<int LEVEL>
class plattices_dense_vector_t:
        public std::vector<plattice_info_dense_t<LEVEL>> {
public:
#if !GNUC_VERSION(4,7,2)
    /* Apparently there's a bug in gcc 4.7.2, which does not recognize
     * that a vector of plattices_dense_vector_t is ok with a move ctor
     * and no copy ctor. I'm not entirely sure why, and I don't really
     * want to investigate. Theoretically, defining ctors as below should
     * be fine.
     */
    plattices_dense_vector_t(plattices_dense_vector_t const&) = delete;
    plattices_dense_vector_t(plattices_dense_vector_t&&) = default;
    plattices_dense_vector_t& operator=(plattices_dense_vector_t&&) = default;
    plattices_dense_vector_t() = default;
#endif
};

// std::vector<plattices_dense_vector_t<LEVEL>> is for remembering the FK
// basis in sublat mode, between two different congruences of (i,j) mod
// m. For simplicity, we remember them only for the toplevel.
template<int LEVEL>
struct precomp_plattice_dense_t {
    typedef std::vector<plattices_dense_vector_t<LEVEL>> type;
};


/* All of these are defined in las-plattice.cpp */
extern template class plattice_enumerator<1>;
extern template class plattice_enumerator_coprime<1>;
extern template class plattices_vector_t<1>;
extern template struct plattice_info_dense_t<1>;
extern template class plattices_dense_vector_t<1>;
extern template struct precomp_plattice_dense_t<1>;
extern template class plattice_enumerator<2>;
extern template class plattice_enumerator_coprime<2>;
extern template class plattices_vector_t<2>;
extern template struct plattice_info_dense_t<2>;
extern template class plattices_dense_vector_t<2>;
extern template struct precomp_plattice_dense_t<2>;
extern template class plattice_enumerator<3>;
extern template class plattice_enumerator_coprime<3>;
extern template class plattices_vector_t<3>;
extern template struct plattice_info_dense_t<3>;
extern template class plattices_dense_vector_t<3>;
extern template struct precomp_plattice_dense_t<3>;



#endif
