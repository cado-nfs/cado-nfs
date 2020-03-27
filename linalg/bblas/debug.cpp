#include "cado.h"
#include "bblas.hpp"

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
