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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <time.h>

#include "gf2x/gf2x-impl.h"
#include "gf2x/gf2x-thresholds.h"

#include "test-tools.h"
#include "check-and-bench-common.h"

static long init_extra_arg = ENGINE_EXTRA_ARG_DEFAULT;

// cputime in millisec.
static int cputime()
{
  return (int) (clock () / ((double) CLOCKS_PER_SEC / 1000.0));
}

int time_total;
int time_dft;
int time_conv;
int time_ift;

/* using the selected FFT engine, multiply two n by n matrices of
 * polynomials with Fl and Gl words in each entry, respectively */
long ENGINE_mul(unsigned long ** H, unsigned long ** F, size_t Fl, unsigned long ** G, size_t Gl, int n)
{
    int t;
    ENGINE_info_t order;

    t=cputime(); time_total -= t;
    ENGINE_init(order, Fl * GF2X_WORDSIZE, Gl * GF2X_WORDSIZE, init_extra_arg);

    ENGINE_ptr f = ENGINE_alloc(order, n*n);

    t=cputime(); time_dft -= t;
    for(int i = 0 ; i < n ; i++)
        for(int j = 0 ; j < n ; j++)
            ENGINE_dft(order, ENGINE_get(order, f, i*n+j), F[i*n+j], Fl * GF2X_WORDSIZE);
    t=cputime(); time_dft += t;

    ENGINE_ptr g = ENGINE_alloc(order, n*n);

    t=cputime(); time_dft -= t;
    for(int i = 0 ; i < n ; i++)
        for(int j = 0 ; j < n ; j++)
            ENGINE_dft(order, ENGINE_get(order, g, i*n+j), G[i*n+j], Gl * GF2X_WORDSIZE);
    t=cputime(); time_dft += t;

    ENGINE_ptr h = ENGINE_alloc(order, n*n);

    t=cputime(); time_conv -= t;
    for(int i = 0 ; i < n ; i++)
        for(int j = 0 ; j < n ; j++) {
            ENGINE_srcptr * ff = malloc(n * sizeof(ENGINE_srcptr));
            ENGINE_srcptr * gg = malloc(n * sizeof(ENGINE_srcptr));
            for(int k = 0 ; k < n ; k++) {
                ff[k] = ENGINE_get_const(order, (ENGINE_srcptr) f, i*n+k);
                gg[k] = ENGINE_get_const(order, (ENGINE_srcptr) g, k*n+j);
            }
            ENGINE_zero(order, ENGINE_get(order, h, i*n+j), 1);
            ENGINE_addcompose_n(order, ENGINE_get(order, h, i*n+j), ff, gg, n);
            free(ff);
            free(gg);
        }
    t=cputime(); time_conv += t;

    ENGINE_free(order, g, n*n);
    ENGINE_free(order, f, n*n);

    t=cputime(); time_ift -= t;
    for(int i = 0 ; i < n ; i++)
        for(int j = 0 ; j < n ; j++)
            ENGINE_ift(order, H[i*n+j], (Fl+Gl) * GF2X_WORDSIZE - 1,  ENGINE_get(order, h, i*n+j));
    t=cputime(); time_ift += t;

    long res = ENGINE_recoverorder(order);

    ENGINE_free(order, h, n*n);
    ENGINE_clear(order);

    t=cputime(); time_total += t;

    return res;
}

void set_extra_arg_from_n32bitwords(size_t n32 GF2X_MAYBE_UNUSED, long supplied)
{
    if (supplied) {
        init_extra_arg = supplied;
        return;
    }
#ifdef GF2X_MUL_FFT_TABLE
#ifdef ENGINE_TERNARY
    int nwords = n32;
#if GF2X_WORDSIZE == 64
    nwords = nwords / 2 + (nwords & 1);
#endif
    int64_t T_FFT_TAB[][2] = GF2X_MUL_FFT_TABLE;
    long max_ix = sizeof(T_FFT_TAB)/sizeof(T_FFT_TAB[0]);
    int j;
#ifdef ARTIFICIAL_NON_SPLIT_VERSION
    /* remove completely the FFT2 setups */
    j = 0;
    for(int i = 0 ; i < max_ix ; i++) {
        if (T_FFT_TAB[i][1] > 0) {
            T_FFT_TAB[j][0] = T_FFT_TAB[i][0];
            T_FFT_TAB[j][1] = T_FFT_TAB[i][1];
            j++;
        }
    }
    max_ix = j;
#endif
    /* artificial: use FFT always */
    assert(T_FFT_TAB[0][0] == 1);
    if (T_FFT_TAB[0][1] == 1) T_FFT_TAB[0][1] = 3;
    /* also for intermediary steps */
    j = 0;
    for(int i = 0 ; i < max_ix ; i++) {
        if (T_FFT_TAB[i][1] != 1) {
            T_FFT_TAB[j][0] = T_FFT_TAB[i][0];
            T_FFT_TAB[j][1] = T_FFT_TAB[i][1];
            j++;
        }
    }
    max_ix = j;
    int i;
    for (i = 0; i + 1 < max_ix && T_FFT_TAB[i + 1][0] <= nwords; i++);
    /* now T_FFT_TAB[i][0] <= nwords < T_FFT_TAB[i+1][0] */
    init_extra_arg = T_FFT_TAB[i][1];
#endif  /* ENGINE_TERNARY */
#endif  /* GF2X_MUL_FFT_TABLE */
}
