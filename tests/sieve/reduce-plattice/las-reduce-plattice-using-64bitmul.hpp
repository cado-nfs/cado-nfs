#ifndef LAS_REDUCE_PLATTICE_USING_64BITMUL_HPP_
#define LAS_REDUCE_PLATTICE_USING_64BITMUL_HPP_

/* This is example code that is meant to go __into__ a class which has
 * several data members and function members defined. For maintenance,
 * all these pieces of code are put in separate files, but these must be
 * #included from within the class itself.
 */
#ifndef LAS_PLATTICE_H
#error "This prologue is only here to please clangd !!!"
#include <cstdint>
#include "macros.h"
struct mock_plattice_info {
        uint32_t mi0;  /* This encodes i0 = -mi0 */
        uint32_t j0;
        uint32_t i1;
        uint32_t j1;
        void reduce_with_vertical_vector(uint32_t I);
        bool needs_special_treatment(uint32_t I) const;
#endif

/* The idea is tempting, but the cost of the full-width 64-bit
 * multiplication is really killing us.
 *
 * Also note that there are several nasty corner cases here and there
 * that force us to keep track of mi0 and i1 almost all the time
 * anyway.
 *
 * This code has many asserts that check that the 64-bit
 * representatives are kept in sync with the (-mi0,j0) and (i1,j1)
 * vectors. Only the debug code checks that. The non-debug code skips
 * this bookkeeping.
 */
void using_64bit_mul(uint32_t I) {
    if (needs_special_treatment(I)) {
        reduce_with_vertical_vector(I);
        return;
    }

    uint64_t Ix = uint64_t(I) << 32;
    constexpr uint64_t M = uint64_t(-1) << 32;
    constexpr uint64_t L MAYBE_UNUSED = uint32_t(-1);
    uint64_t ij0 = (uint64_t(-mi0) << 32) | j0;
    uint64_t ij1 = (uint64_t(i1) << 32) | j1;
    uint64_t t;
#define ASSERT_CONSISTENCY() do {					\
        ASSERT_PLATTICE(mi0 == -((uint32_t)(ij0>>32)));		\
        ASSERT_PLATTICE(i1 == ij1 >> 32);				\
        ASSERT_PLATTICE(j0 == (ij0 & L));				\
        ASSERT_PLATTICE(j1 == (ij1 & L));				\
} while (0)

    for( ;; ) {
        ASSERT_PLATTICE(j0 <= j1);
        ASSERT_CONSISTENCY();
        t = ij1 - Ix;
        i1 = ij1 >> 32;
        ASSERT_PLATTICE((i1 < I) == (t >= ij1));
        mi0 = -((uint32_t)(ij0>>32));
        if (t >= ij1) { // i1 < I) {
            j1 = ij1;
            if (!i1) {
                j0 = ij1 - ij0;
                reduce_with_vertical_vector(I);
                return;
            }
            ASSERT_PLATTICE(mi0 + i1 >= I);

            int a = (mi0 + i1 - I) / i1;
            ij0 += a * ij1;
#ifndef NDEBUG
            mi0 -= a * i1;
            j0  += a * j1;
#endif
            ASSERT_CONSISTENCY();
            mi0 = -((uint32_t)(ij0>>32));
            j0 = ij0;
            return;
        }
        if (mi0 < i1 * 3) {
            {
                ij0 += ij1;
                mi0 -= i1; // see below
#ifndef NDEBUG
                j0 += j1;
#endif
            }
            /* what we really want to test is mi0 >= i1. However the
             * j coordinates come into play when comparing -ij0 with
             * ij1. 
             */
            if (mi0 >= i1) { // if (-ij0 >= ij1) {
                ij0 += ij1;
#ifndef NDEBUG
                mi0 -= i1; j0 += j1;
#endif
            }
        } else
        {
            int k = mi0 / i1;
            // int k = (-(ij0 & M)) / (ij1 & M);
            ASSERT_PLATTICE(k);
            ij0 += k * ij1;
#ifndef NDEBUG
            mi0 -= k * i1;
            j0 += k * j1;
#endif
            ASSERT_CONSISTENCY();
        }

        ASSERT_PLATTICE(j1 <= j0);

        ASSERT_CONSISTENCY();
        t = (ij0&M) + Ix;
        ASSERT_PLATTICE((mi0 < I) == ((-(ij0&M)) < Ix));
        ASSERT_PLATTICE((mi0 < I) == (t != 0 && t <= Ix));
        mi0 = -((uint32_t)(ij0>>32));
        i1 = ij1 >> 32;
        if (t != 0 && t <= Ix) { // mi0 < I) {
            if (!mi0) {
                mi0 = i1;
                j0 = ij1;
                j1 = ij0;
                i1 = 0;
                reduce_with_vertical_vector(I);
                return;
            }
            int a = (mi0 + i1 - I) / mi0;
            ij1 += a * ij0;
#ifndef NDEBUG
            i1 -= a * mi0;
            j1 += a * j0;
#endif
            ASSERT_CONSISTENCY();
            i1 = ij1 >> 32;
            j1 = ij1;
            j0 = ij0;
            return;
        }
        if (i1 < 3 * mi0) {
            {
                ij1 += ij0;
                i1 -= mi0; // see below
#ifndef NDEBUG
                j1 += j0;
#endif
            }
            // same remark as above
            if (i1 >= mi0) { // if (ij1 >= -ij0) {
                ij1 += ij0;
#ifndef NDEBUG
                i1 -= mi0; j1 += j0;
#endif
            }
        } else
        {
            int k = i1 / mi0;
            ASSERT(k);
            ij1 += k * ij0;
#ifndef NDEBUG
            i1 -= k * mi0; j1 += k * j0;
#endif
            ASSERT_CONSISTENCY();
        }
    }
#undef ASSERT_CONSISTENCY
}

#ifndef LAS_PLATTICE_H
#error "This epilogue is only here to please clangd !!!"
};
#endif

#endif	/* LAS_REDUCE_PLATTICE_USING_64BITMUL_HPP_ */
