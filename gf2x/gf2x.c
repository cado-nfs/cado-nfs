/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2014, 2015
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of either:
    - If the archive contains a file named toom-gpl.c (not a trivial
    placeholder), the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.
    - If the archive contains a file named toom-gpl.c which is a trivial
    placeholder, the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.
   
   You should have received a copy of the GNU General Public License as
   well as the GNU Lesser General Public License along with this program;
   see the files COPYING and COPYING.LIB.  If not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301, USA.
*/

#include "gf2x/gf2x-config.h"

#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "gf2x.h"
#include "gf2x/gf2x-impl.h"

/* We _do_ need the inlines here. */
#include "gf2x/gf2x-small.h"

#include "gf2x/gf2x-thresholds.h"

int gf2x_lib_version_code = GF2X_VERSION_CODE;

void gf2x_mul_basecase_inner(unsigned long * c, const unsigned long * a,
			 long na, const unsigned long * b, long nb)
{
    assert(c != a);
    assert(c != b);

    if (na == nb) {
	switch (na) {
	case 0:		// This can occur in call from KarMul return;
	case 1: gf2x_mul1(c, a[0], b[0]); return;
	case 2: gf2x_mul2(c, a, b); return;
	case 3: gf2x_mul3(c, a, b); return;
	case 4: gf2x_mul4(c, a, b); return;
	case 5: gf2x_mul5(c, a, b); return;
	case 6: gf2x_mul6(c, a, b); return;
	case 7: gf2x_mul7(c, a, b); return;
	case 8: gf2x_mul8(c, a, b); return;
	case 9: gf2x_mul9(c, a, b); return;
	default:
#if GF2X_GNUC_VERSION_ATLEAST(4,5,0)
            __builtin_unreachable();
#endif
            return;
	}
    } else if (na < nb) {
        /* FIXME -- this does not seem efficient */
        long i;
        /* No need to clear c[], it's filled up progressively */
        c[nb] = gf2x_mul_1_n(c, b, nb, a[0]);
        for (i = 1; i < na; i++) {
            c[nb + i] = gf2x_addmul_1_n(c + i, c + i, b, nb, a[i]);
        }
    } else {
        gf2x_mul_basecase_inner(c, b, nb, a, na);
    }
}

#ifdef GF2X_MUL_FFT_TABLE
int64_t T_FFT_TAB[][2] = GF2X_MUL_FFT_TABLE;
#endif

/* This is the toplevel multiplication routine. It handles the temporary
 * storage if necessary.
 */

void gf2x_mul_pool_init(gf2x_mul_pool_t p)
{
    memset(p, 0, sizeof(gf2x_mul_pool_t));
}

void gf2x_mul_pool_clear(gf2x_mul_pool_t p)
{
    free(p->stk);
    p->stk_size = 0;
}

int gf2x_mul(unsigned long * c,
        const unsigned long * a, unsigned long sa,
        const unsigned long * b, unsigned long sb)
{
    return gf2x_mul_r(c, a, sa, b, sb, NULL);
}

int gf2x_mul_r(unsigned long * c,
        const unsigned long * a, unsigned long sa,
        const unsigned long * b, unsigned long sb, gf2x_mul_pool_t pool)
{
    int rc = 0;
    unsigned long sc = sa + sb;
    /* As a starting guess, assume that c may alias a or b */
    unsigned long * dst = c;
    gf2x_mul_pool_t xxpool;
    struct gf2x_mul_pool_s * xpool = NULL;

    if (sa > sb) {
        return gf2x_mul_r(c, b, sb, a, sa, pool);
    }
    // now sa <= sb (note: sa and sb are interchanged in Toom3uMul etc)

    if (sa < GF2X_MUL_KARA_THRESHOLD) {
        /* This calls the hand-crafted code if sa == sb */
        if (dst == a || dst == b) {
            dst = malloc(sc * sizeof(unsigned long));
            if (dst == NULL) {
                rc = GF2X_ERROR_OUT_OF_MEMORY;
                goto end_of_gf2x_mul_r;
            }
        }
        gf2x_mul_basecase(dst, a, sa, b, sb);
        goto end_of_gf2x_mul_r;
    }

    /* This ugly cpp block entirely disables the FFT if it has
     * not yet been tuned */
#ifdef GF2X_MUL_FFT_TABLE
    long ix, K, sab = sc / 2;
    long max_ix = sizeof(T_FFT_TAB)/sizeof(T_FFT_TAB[0]);
    for (ix = 0; ix + 1 < max_ix && T_FFT_TAB[ix + 1][0] <= sab; ix++);
    /* now T_FFT_TAB[ix][0] <= sab < T_FFT_TAB[ix+1][0] */
    K = T_FFT_TAB[ix][1];

    /* gf2x_mul_fft can handle unbalanced operands if not too
     * small: return the result in {dst, sa+sb} */

    /* Note that dst aliasing a or b works with fft. */

    if (sc >= GF2X_TERNARY_FFT_MINIMUM_SIZE && K && K != 1) {
        rc = gf2x_mul_fft(dst, a, sa, b, sb, K);
        goto end_of_gf2x_mul_r;
    }
#endif

    unsigned long sp, sp2;
    sp = gf2x_toomspace(sa); // Space for balanced TC routines

    if (sa != sb)
      {
        sp2 = gf2x_toomuspace(2 * sa); // Space for unbalanced TC routines
        if (sp < sp2)
          sp = sp2; /* worst-case required */
        sp2 = 2 * sa + gf2x_toomspace(sa); // Space for unbalanced TC routines w/ lazy cut
        if (sp < sp2)
          sp = sp2; /* worst-case required */
      }

    if (pool == NULL) {
        gf2x_mul_pool_init(xxpool);
        xpool = xxpool;
    } else {
        xpool = pool;
    }

    if (xpool->stk_size < sp) {
        void * p = realloc(xpool->stk, sp * sizeof(unsigned long));
        if (p == NULL) {
            rc = GF2X_ERROR_OUT_OF_MEMORY;
            goto end_of_gf2x_mul_r;
        }
        xpool->stk = p;
        xpool->stk_size = sp;
    }

    /* None of the alternatives below supports c aliasing a or b */
    if (dst == a || dst == b) {
        dst = malloc(sc * sizeof(unsigned long));
        if (dst == NULL) {
            rc = GF2X_ERROR_OUT_OF_MEMORY;
            goto end_of_gf2x_mul_r;
        }
    }

    if (sa == sb) {
        // Avoid copy in common case
        gf2x_mul_toom(dst, a, b, sa, xpool->stk);
#if GPL_CODE_PRESENT
    } else if ((sa == (sb + 1) / 2) && gf2x_best_utoom(sb)) {
        // Another common case
        // due to GCD algorithm
        gf2x_mul_tc3u(dst, b, sb, a, xpool->stk);
#endif  /* GPL_CODE_PRESENT */
    } else {
        unsigned long *v = xpool->stk + gf2x_toomspace(sa);

        unsigned long i, j;

        memset(dst, 0, sc * sizeof(unsigned long));

        unsigned long * ptr = dst;

        for (;;) {
            if (sa == 0)
                break;

            if (sa == 1) {
                ptr[sb] ^= gf2x_addmul_1_n(ptr, ptr, b, sb, a[0]);
                break;
            }

            // finally: the general case
            for (i = 0; i + sa <= sb; i += sa) {
                // Generic (balanced) Toom-Cook mult.
                gf2x_mul_toom(v, a, b + i, sa, xpool->stk);
                for (j = 0; j < 2 * sa; j++)
                    ptr[i + j] ^= v[j];
            }

            {
                const unsigned long *t;
                unsigned long st;

                /* Swap a and b, and go for the next spin */
                t = a;
                st = sa;
                a = b + i;
                sa = sb - i;
                b = t;
                sb = st;
            }
            ptr = ptr + i;
        }
    }
end_of_gf2x_mul_r:
    if (pool == NULL && xpool != NULL) {
        gf2x_mul_pool_clear(xxpool);
    }
    if (dst && dst != c) {
        /* Then we have allocated a temp buffer */
        memcpy(c, dst, sc * sizeof(unsigned long));
        free(dst);
    }
    return rc;
}

/* vim: set sw=4 sta et: */
