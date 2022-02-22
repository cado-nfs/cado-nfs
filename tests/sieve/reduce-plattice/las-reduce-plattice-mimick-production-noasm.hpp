#ifndef LAS_REDUCE_PLATTICE_MIMICK_PRODUCTION_NOASM_HPP_
#define LAS_REDUCE_PLATTICE_MIMICK_PRODUCTION_NOASM_HPP_

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

/* This is a C version of the production asm code. In fact, we don't
 * use this code when asm fails us, because it's very badly dealt
 * with by the compiler.
 */
void mimick_production_noasm(uint32_t I) {
    if (needs_special_treatment(I)) {
        reduce_with_vertical_vector(I);
        return;
    }

    uint64_t u0 = (((uint64_t)-j0)<<32) | (uint64_t) mi0;
    uint64_t u1 = (((uint64_t) j1)<<32) | (uint64_t) i1;
    int branch;
    for( ;; ) {
        branch = 0;
        if (((uint32_t) u1) < I) break;
        {
            uint64_t saved = u0;
            u0 -= u1;
            if (((uint32_t) saved) < ((uint32_t) u1)) { u0 = saved; }
        }
        {
            uint64_t saved = u0;
            u0 -= u1;
            if (((uint32_t) saved) < ((uint32_t) u1)) { u0 = saved; }
            else {
                /* do the division */
                int k = ((uint32_t) u0) / ((uint32_t) i1);
                /* the new mi0 is accessible as the remainder of the
                 * division, so we can assume that the compiler will
                 * do the right thing.
                 */
                uint32_t new_mi0 = (((uint32_t) u0) % ((uint32_t) i1));
                uint32_t old_mj0 = u0 >> 32;
                uint32_t old_j1 = u1 >> 32;
                uint32_t new_mj0 = old_mj0 - k * old_j1;
                u0 = (((uint64_t) new_mj0) << 32) + new_mi0;
            }
        }
        branch = I;
        if (((uint32_t) u0) < I) break;
        {
            uint64_t saved = u1;
            u1 -= u0;
            if (((uint32_t) saved) < ((uint32_t) u0)) { u1 = saved; }
        }
        {
            uint64_t saved = u1;
            u1 -= u0;
            if (((uint32_t) saved) < ((uint32_t) u0)) { u1 = saved; }
            else {
                /* do the division */
                int k = ((uint32_t) u1) / ((uint32_t) u0);
                /* the new i1 is accessible as the remainder of the
                 * division, so we can assume that the compiler will
                 * do the right thing.
                 */
                uint32_t new_i1 = (((uint32_t) u1) % ((uint32_t) u0));
                uint32_t old_j1 = u1 >> 32;
                uint32_t old_mj0 = u0 >> 32;
                uint32_t new_j1 = old_j1 - k * old_mj0;
                u1 = (((uint64_t) new_j1) << 32) + new_i1;
            }
        }
    }
    mi0 = u0;
    i1 = u1;
    j0 = -(int32_t)(u0>>32);
    j1 = (int32_t)(u1>>32);
    if (branch) {
        if (mi0 == 0) {
            mi0 = i1;
            i1 = j0 ; j0 = j1 ; j1 = i1;
            i1 = 0;
            reduce_with_vertical_vector(I);
            return;
        }
        ASSERT(mi0 + i1 >= I);
        int a = (mi0 + i1 - I) / mi0;
        i1 -= a * mi0;
        j1 += a * j0;
    } else {
        if (i1 == 0) {
            // Lo=matrix([ (mi0, j1-j0), (i1, j1)])
            j0 = j1 - j0;
            reduce_with_vertical_vector(I);
            return;
        }
        ASSERT(mi0 + i1 >= I);
        int a = (mi0 + i1 - I) / i1;
        mi0 -= a * i1;
        j0  += a * j1;
    }
}

#ifndef LAS_PLATTICE_H
#error "This epilogue is only here to please clangd !!!"
};
#endif

#endif	/* LAS_REDUCE_PLATTICE_MIMICK_PRODUCTION_NOASM_HPP_ */
