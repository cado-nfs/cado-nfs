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

#define _GNU_SOURCE
#include "gf2x.h"
#include "gf2x/gf2x-config.h"
#include "gf2x/gf2x-impl.h"

#include "gf2x-fake-fft.h"
#include "macros.h"


#define ULONG_BITS      GF2X_WORDSIZE

/* Number of words holding B bits ; better naming sought. */
#ifndef iceildiv
/* unfortunately this fails miserably if x+y-1 overflows */
#define iceildiv(x,y)	(((x)+(y)-1)/(y))
#endif
#define BITS_TO_WORDS(B,W)      iceildiv((B),(W))

/* nF, nG : number of coefficients */
void gf2x_fake_fft_init(gf2x_fake_fft_info_ptr p, size_t nF, size_t nG, ...)
{
    p->n1 = nF;
    p->n2 = nG;
    p->n3 = nF + nG - 1;
    size_t nc = 0;
    if (p->n1 > nc) nc = p->n1;
    if (p->n2 > nc) nc = p->n2;
    p->size = 2 * BITS_TO_WORDS(nc, ULONG_BITS);
}

/* n is a number of coefficients ! */
void gf2x_fake_fft_dft(gf2x_fake_fft_info_srcptr p GF2X_MAYBE_UNUSED, gf2x_fake_fft_ptr dst, const unsigned long * src, size_t n) {
    ASSERT(n <= p->n1 || n <= p->n2);
    size_t s = BITS_TO_WORDS(n, ULONG_BITS);
    memcpy(dst, src, s * sizeof(unsigned long));
    if (n % ULONG_BITS) {
        /* Just as we are computing this assertion, we could easily mask out
         * the bits ourselves. However, our interface mandates that the high
         * bits be cleared in any case. So make sure we properly enforce this
         * constraint.
         */
        ASSERT((src[s-1] & ~((1UL << (n % ULONG_BITS)) - 1)) == 0);
    }
    memset(dst + s, 0, (p->size - s) * sizeof(unsigned long));
}

/* n is a number of coefficients ! -- requiring less than the compute size is
 * okay, as long as it's understood as a means for truncating the data. So we
 * don't do checking for zero high bits.
 */
void gf2x_fake_fft_ift(gf2x_fake_fft_info_srcptr p GF2X_MAYBE_UNUSED, unsigned long * dst, size_t n, gf2x_fake_fft_srcptr src) {
    ASSERT(n <= p->n3);
    size_t t = BITS_TO_WORDS(n, ULONG_BITS);
    memcpy(dst, src, t * sizeof(unsigned long));
}
void gf2x_fake_fft_compose(gf2x_fake_fft_info_srcptr p GF2X_MAYBE_UNUSED, gf2x_fake_fft_ptr dst, gf2x_fake_fft_srcptr s1, gf2x_fake_fft_srcptr s2) {
    size_t n1 = BITS_TO_WORDS(p->n1, ULONG_BITS);
    size_t n2 = BITS_TO_WORDS(p->n2, ULONG_BITS);
    gf2x_mul(dst, s1, n1, s2, n2);
}
void gf2x_fake_fft_addcompose_n(gf2x_fake_fft_info_srcptr p GF2X_MAYBE_UNUSED, gf2x_fake_fft_ptr dst, gf2x_fake_fft_srcptr * s1, gf2x_fake_fft_srcptr * s2, size_t n) {
    size_t n1 = BITS_TO_WORDS(p->n1, ULONG_BITS);
    size_t n2 = BITS_TO_WORDS(p->n2, ULONG_BITS);
    unsigned long * h = malloc(p->size * sizeof(unsigned long));
    /* lacking addmul in gf2x, we do some extra allocation */
    memset(h, 0, p->size * sizeof(unsigned long));
    for(size_t k = 0 ; k < n ; k++) {
        gf2x_mul(h, s1[k], n1, s2[k], n2);
        for(unsigned int i = 0 ; i < p->size ; i++) {
            dst[i] ^= h[i];
        }
    }
    free(h);
}
void gf2x_fake_fft_addcompose(gf2x_fake_fft_info_srcptr p GF2X_MAYBE_UNUSED, gf2x_fake_fft_ptr dst, gf2x_fake_fft_srcptr s1, gf2x_fake_fft_srcptr s2) {
    gf2x_fake_fft_addcompose_n(p, dst, &s1, &s2, 1);
}
void gf2x_fake_fft_add(gf2x_fake_fft_info_srcptr p GF2X_MAYBE_UNUSED, gf2x_fake_fft_ptr dst, gf2x_fake_fft_srcptr s1, gf2x_fake_fft_srcptr s2) {
    size_t i;
    for(i = 0 ; i < p->size ; i++) {
        dst[i] = s1[i] ^ s2[i];
    }
}

void gf2x_fake_fft_cpy(gf2x_fake_fft_info_srcptr p GF2X_MAYBE_UNUSED, gf2x_fake_fft_ptr dst, gf2x_fake_fft_srcptr s) {
    memcpy(dst, s, (p->size)*sizeof(unsigned long));
}

size_t gf2x_fake_fft_size(gf2x_fake_fft_info_srcptr p) {
    return p->size;
}

void gf2x_fake_fft_init_similar(gf2x_fake_fft_info_ptr o, size_t bits_a, size_t bits_b, gf2x_fake_fft_info_srcptr other GF2X_MAYBE_UNUSED)
{
    gf2x_fake_fft_init(o, bits_a, bits_b);
}

int gf2x_fake_fft_compatible(gf2x_fake_fft_info_srcptr o1 GF2X_MAYBE_UNUSED, gf2x_fake_fft_info_srcptr o2 GF2X_MAYBE_UNUSED)
{
    return 1;
}

/* vim: set sw=4 sta et: */
