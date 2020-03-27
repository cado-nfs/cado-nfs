#include "cado.h"
#include "bblas.hpp"

/* at least theoretically, we should be able to use the following unicode
 * characters to obtain more condensed printing.
 *
    U+2588	█	Full block
    U+2599	▙	Quadrant upper left and lower left and lower right
    U+259B	▛	Quadrant upper left and upper right and lower left
    U+259C	▜	Quadrant upper left and upper right and lower right
    U+259F	▟	Quadrant upper right and lower left and lower right
    U+2580	▀	Upper half block
    U+2584	▄	Lower half block
    U+258C	▌	Left half block
    U+2590	▐	Right half block
    U+259A	▚	Quadrant upper left and lower right
    U+259E	▞	Quadrant upper right and lower left
    U+2596	▖	Quadrant lower left
    U+2597	▗	Quadrant lower right
    U+2598	▘	Quadrant upper left
    U+259D	▝	Quadrant upper right
 *
 * But this should go in a separate python script, naturally.
 *
 */
void pmat_6464(mat64 m)
{
    for(int i = 0; i < 64 ; i++) {
        uint64_t mask=1;
        for(int j = 0 ; j < 64 ; j++, mask<<=1) {
            putchar((m[i]&mask) ?'1':'0');
        }
        putchar('\n');
    }
}

void pmat_mn(mat64 * m, int rb, int cb)
{
    ASSERT_ALWAYS(rb%64 == 0); rb /= 64;
    ASSERT_ALWAYS(cb%64 == 0); cb /= 64;
    for(int ib = 0 ; ib < rb ; ib++) {
        mat64 * mr = m + ib * cb;
        for(int i = 0; i < 64 ; i++) {
            for(int jb = 0 ; jb < cb ; jb++) {
                uint64_t mask=1;
                for(int j = 0 ; j < 64 ; j++, mask<<=1) {
                    putchar((mr[jb][i]&mask) ?'1':'0');
                }
            }
            putchar('\n');
        }
    }
}
