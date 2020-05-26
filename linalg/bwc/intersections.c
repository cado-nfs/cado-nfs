#include "cado.h" // IWYU pragma: keep
#include "macros.h"
#include "intersections.h"

unsigned int intersect_two_intervals(unsigned int * offset_a, unsigned int * offset_b, unsigned int a0, unsigned int a1, unsigned int b0, unsigned int b1)
{
    if (a0 >= b1 || b0 >= a1) {
        *offset_a = 0;
        *offset_b = 0;
        return 0;
    }

    int d0 = (int) (b0 - a0);
    int d1 = (int) (b1 - a1);

    unsigned int u0 = 0;
    unsigned int u1 = a1 - a0;
    unsigned int v0 = 0;
    unsigned int v1 = b1 - b0;

    if (d0 > 0) {
        u0 += d0;
        u1 -= d0;
    } else {
        v0 += -d0;
        v1 -= -d0;
    }

    if (d1 < 0) {
        u1 -= -d1;
    } else {
        v1 -= d1;
    }

    /* Our data will be of interest for indices in the
     * range [u0..u0+u1[. It goes to range [v0..v0+v1[
     */

    ASSERT_ALWAYS(u1 == v1);

    *offset_a = u0;
    *offset_b = v0;
    return u1;
}
