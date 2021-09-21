#ifndef LAS_REDUCE_PLATTICE_SWAPPING_LOOP_HPP_
#define LAS_REDUCE_PLATTICE_SWAPPING_LOOP_HPP_

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

void swapping_loop(uint32_t I)
{
    if (needs_special_treatment(I)) {
        reduce_with_vertical_vector(I);
        return;
    }

    /* This is the main reduce_plattice loop */
    int flip;
    for(flip = 0 ; i1 >= I; flip ^= 1 ) {
        /* do partial unrolling for the frequent case where the
         * quotient is either 1 or 2.
         * this has a significant overall impact
         */
#if 0
        if (mi0 < i1 * 3) {
            { mi0 -= i1; j0 += j1; }
            if (mi0 >= i1) { mi0 -= i1; j0 += j1; }
        } else
#endif
        {
            int k = mi0 / i1; mi0 -= k * i1; j0 += k * j1;
        }
        std::swap(mi0, i1);
        std::swap(j0, j1);
    }
    /* an "UNLIKELY" macro here actually has an adverse
     * effect...  */
    if (i1 == 0) {
        if (!flip)
            // Lo=matrix([ (mi0, j1-j0), (i1, j1)])
            j0 = j1 - j0;
        reduce_with_vertical_vector(I);
        return;
    } else {
        int a = (mi0 + i1 - I) / i1;
        mi0 -= a * i1;
        j0  += a * j1;
    }
    if (flip) {
        std::swap(mi0, i1);
        std::swap(j0, j1);
    }
}

#ifndef LAS_PLATTICE_H
#error "This epilogue is only here to please clangd !!!"
};
#endif

#endif	/* LAS_REDUCE_PLATTICE_SWAPPING_LOOP_HPP_ */
