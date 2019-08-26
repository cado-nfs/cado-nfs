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

#define _XOPEN_SOURCE   600
#define _POSIX_C_SOURCE   600
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <assert.h>
#ifdef  _POSIX_C_SOURCE
#include <sys/types.h>
#include <unistd.h>
#endif

#include "gf2x.h"

#ifdef GF2X_SOURCE_TREE
#include "fft/gf2x-fft.h"
#else
#include "gf2x/gf2x-fft.h"
#endif

/* Handy, and does not require libm */
#ifndef iceildiv
/* C99 defines division as truncating, i.e. rounding to 0. To get ceil in all
   cases, we have to distinguish by sign of arguments. Note that (afaik) the
   sign and truncation direction of the quotient with any negative operads was
   undefined in C89. */
/* iceildiv requires *unsigned* operands in spite of the name suggesting
   (signed) integer type. For negative operands, the result is wrong. */
#define iceildiv(x,y) ((x) == 0 ? 0 : ((x)-1)/(y)+1)
#endif

#ifndef ULONG_BITS
#define	ULONG_BITS	((int) (sizeof(unsigned long) * CHAR_BIT))
#endif

/* Checks the different FFT codes, and see whether they give any
 * consistent results.
 *
 * We compute f1g1+f2g2 for several sizes of f1, f2, using all the
 * available algorithms.
 *
 *
 * This tests _only_ the functionality in the libgf2x-fft that is
 * eventually shipped.
 *
 * Note that this part of the library is still clearly development work.
 *
 * When this fail, we print the polynomials computed as lists of 64-bit
 * ints, to be parsed with magma as follows:
 *
 * f2:=Polynomial(GF(2),&cat[Intseq(x,2,64):x in f2]);
 *
 * F:=func<f2|Polynomial(GF(2),&cat[Intseq(x,2,64):x in f2])>;
 * f2:=F(f2);
 *
 */

void usage()
{
    fprintf(stderr, "Usage: ./gf2x_fft_check <max size>\n");
}

void display(const char * text, unsigned long * f, size_t nf)
{
    printf("%s:=[", text);
    for(size_t i = 0 ; i < nf ; i++) {
        if (i) printf(",");
        printf("%lu",f[i]);
    }
    printf("];\n");
}

/* returns an unsigned long with long strings of zeroes and ones */
unsigned long rand2_ulong()
{
    /* for 64-bit ulongs, we need 1 + nrounds * 6 bits as given by
     * random(), which is ok until nchunks == 5.
     *
     * Note that we don't care much about randomness in the strict sense
     * (in case anybody still wonders).
     */
    unsigned int r = random();
    unsigned long bit = r & 1UL;
    unsigned long res = -r;
    unsigned int k = 0;
    for( ; (1 << k) < ULONG_BITS ; k++);
    r >>= 1;
    for(unsigned int i = 0 ; i < 4 ; i++) {
        unsigned int chunksize = 1 + (r & (ULONG_BITS - 1));
        bit = !bit;
        res += bit;
        r >>= k;
        res <<= chunksize;
        res -= bit;
    }
    return res;
}

int docheck(size_t nf, size_t ng, int nrep)
{
    // doing f1g1+f2g2
    unsigned long * f1, * f2, * g1, * g2;
    unsigned long * gf2x_fake_fft_h;
    unsigned long * gf2x_cantor_fft_h;
    unsigned long * gf2x_ternary_fft_h;
    unsigned long * gf2x_fake_fft_hx;
    unsigned long * gf2x_cantor_fft_hx;
    unsigned long * gf2x_ternary_fft_hx;

    unsigned int nf1 = nf;
    unsigned int ng1 = ng;
    unsigned int nf2 = nf;
    unsigned int ng2 = ng;

    unsigned int nh = nf + ng - 1;

    size_t nwf1 = iceildiv(nf1, ULONG_BITS);
    size_t nwg1 = iceildiv(ng1, ULONG_BITS);
    size_t nwf2 = iceildiv(nf2, ULONG_BITS);
    size_t nwg2 = iceildiv(ng2, ULONG_BITS);

    size_t nwh = iceildiv(nh, ULONG_BITS);

    f1 = malloc(nwf1 * sizeof(unsigned long));
    g1 = malloc(nwg1 * sizeof(unsigned long));
    f2 = malloc(nwf2 * sizeof(unsigned long));
    g2 = malloc(nwg2 * sizeof(unsigned long));

    gf2x_fake_fft_h = malloc(nwh * sizeof(unsigned long));
    gf2x_cantor_fft_h = malloc(nwh * sizeof(unsigned long));
    gf2x_ternary_fft_h = malloc(nwh * sizeof(unsigned long));

    gf2x_fake_fft_hx = malloc(nwh * sizeof(unsigned long));
    gf2x_cantor_fft_hx = malloc(nwh * sizeof(unsigned long));
    gf2x_ternary_fft_hx = malloc(nwh * sizeof(unsigned long));

    /*
    for(size_t i = 0 ; i < nwf1 ; i++) printf("%016lx", f1[i]); printf("\n");
    for(size_t i = 0 ; i < nwg1 ; i++) printf("%016lx", g1[i]); printf("\n");
    for(size_t i = 0 ; i < nwf2 ; i++) printf("%016lx", f2[i]); printf("\n");
    for(size_t i = 0 ; i < nwg2 ; i++) printf("%016lx", g2[i]); printf("\n");
    */

#define ALLOC1(E,x) E ## _t * E ## _ ## x = E ## _alloc(E, 1)

#define SETUP(E, nf, ng)               				\
    E ## _info_t E;						\
    E ## _init(E, nf, ng);					\
    ALLOC1(E, tf1); ALLOC1(E, tg1);				\
    ALLOC1(E, tf2); ALLOC1(E, tg2);				\
    E ## _ptr E ## _tfs[2] = { E ## _tf1, E ## _tf2, };           \
    E ## _ptr E ## _tgs[2] = { E ## _tg1, E ## _tg2, };           \
    ALLOC1(E, th);                                              \
    ALLOC1(E, thx);

#define SETUP3(E, nf, ng, K)           				\
    E ## _info_t E;						\
    E ## _init(E, nf, ng, K);					\
    ALLOC1(E, tf1); ALLOC1(E, tg1);				\
    ALLOC1(E, tf2); ALLOC1(E, tg2);				\
    E ## _ptr E ## _tfs[2] = { E ## _tf1, E ## _tf2, };           \
    E ## _ptr E ## _tgs[2] = { E ## _tg1, E ## _tg2, };           \
    ALLOC1(E, th);                                              \
    ALLOC1(E, thx);

#define DO_zero_th(E) E ## _zero(E, E ## _th, 1)
#define DO_zero_thx(E) E ## _zero(E, E ## _thx, 1)
#define DO_dft_f1(E) E ## _dft(E, E ## _tf1, f1, nf1)
#define DO_dft_g1(E) E ## _dft(E, E ## _tg1, g1, ng1)
#define DO_dft_f2(E) E ## _dft(E, E ## _tf2, f2, nf2)
#define DO_dft_g2(E) E ## _dft(E, E ## _tg2, g2, ng2)
#define DO_comp1(E)  E ## _addcompose(E, E ## _th, E ## _tf1, E ## _tg1)
#define DO_comp2(E)  E ## _addcompose(E, E ## _th, E ## _tf2, E ## _tg2)
#define DO_compx(E)  E ## _addcompose_n(E, E ## _thx, (E ## _srcptr *) E ## _tfs, (E ## _srcptr *) E ## _tgs, 2)
#define DO_ift_h(E)  E ## _ift(E, E ## _h, nh, E ## _th)
#define DO_ift_hx(E) E ## _ift(E, E ## _hx, nh, E ## _thx)

#define EXTRA_DISPLAY(E)  do {       					\
    display(#E "_tf1", (unsigned long *) E ## _tf1,                     \
            E ## _size(E) * sizeof(E ## _t) / sizeof(unsigned long));   \
    display(#E "_tg1", (unsigned long *) E ## _tg1,                     \
            E ## _size(E) * sizeof(E ## _t) / sizeof(unsigned long));	\
    display(#E "_tf2", (unsigned long *) E ## _tf2,                     \
            E ## _size(E) * sizeof(E ## _t) / sizeof(unsigned long));	\
    display(#E "_tg2", (unsigned long *) E ## _tg2,                     \
            E ## _size(E) * sizeof(E ## _t) / sizeof(unsigned long));	\
    display(#E "_th",  (unsigned long *) E ## _th,                      \
            E ## _size(E) * sizeof(E ## _t) / sizeof(unsigned long));	\
    display(#E "_thx",  (unsigned long *) E ## _thx,                    \
            E ## _size(E) * sizeof(E ## _t) / sizeof(unsigned long));   \
} while (0)

#define FREE1(E, x) E ## _free(E, E ## _ ## x, 1)

#define LEAVE(E)        						\
    FREE1(E, tf1); FREE1(E, tg1);			        	\
    FREE1(E, tf2); FREE1(E, tg2);		        		\
    FREE1(E, th);							\
    FREE1(E, thx);							\
    E ## _clear(E);

    unsigned int seed = random();
    srandom(seed);
    SETUP(gf2x_fake_fft, nf, ng)
    SETUP(gf2x_cantor_fft, nf, ng)
    SETUP3(gf2x_ternary_fft, nf, ng, 81)
    for(int i = 0 ; i < nrep ; i++) {
        for(size_t i = 0 ; i < nwf1 ; i++) f1[i] = rand2_ulong();
        for(size_t i = 0 ; i < nwg1 ; i++) g1[i] = rand2_ulong();
        for(size_t i = 0 ; i < nwf2 ; i++) f2[i] = rand2_ulong();
        for(size_t i = 0 ; i < nwg2 ; i++) g2[i] = rand2_ulong();
        if (nf1 % ULONG_BITS) f1[nwf1-1] &= (1UL << (nf1 % ULONG_BITS)) - 1;
        if (ng1 % ULONG_BITS) g1[nwg1-1] &= (1UL << (ng1 % ULONG_BITS)) - 1;
        if (nf2 % ULONG_BITS) f2[nwf2-1] &= (1UL << (nf2 % ULONG_BITS)) - 1;
        if (ng2 % ULONG_BITS) g2[nwg2-1] &= (1UL << (ng2 % ULONG_BITS)) - 1;

        DO_dft_f1(gf2x_fake_fft);
        DO_dft_g1(gf2x_fake_fft);
        DO_dft_f2(gf2x_fake_fft);
        DO_dft_g2(gf2x_fake_fft);
        DO_zero_th(gf2x_fake_fft);
        DO_comp1(gf2x_fake_fft);
        DO_comp2(gf2x_fake_fft);
        DO_ift_h(gf2x_fake_fft);

        DO_zero_thx(gf2x_fake_fft);
        DO_compx(gf2x_fake_fft);
        DO_ift_hx(gf2x_fake_fft);

        if (memcmp(gf2x_fake_fft_h, gf2x_fake_fft_hx, nwh * sizeof(unsigned long)) != 0) {
            fprintf(stderr, "fake != fakex for %zu*%zu ; seed=%u\n",
                    nf,ng,seed);
            printf("w:=%d;\n", (int) ULONG_BITS);
            display("f1",f1,nwf1);
            display("g1",g1,nwg1);
            display("f2",f2,nwf2);
            display("g2",g2,nwg2);
            display("gf2x_fake_fft_h",gf2x_fake_fft_h,nwh);
            display("gf2x_fake_fft_hx",gf2x_fake_fft_hx,nwh);
            EXTRA_DISPLAY(gf2x_fake_fft);
            fflush(stdout);
            fflush(stderr);
            abort();
        }


        DO_dft_f1(gf2x_cantor_fft);
        DO_dft_g1(gf2x_cantor_fft);
        DO_dft_f2(gf2x_cantor_fft);
        DO_dft_g2(gf2x_cantor_fft);
        DO_zero_th(gf2x_cantor_fft);
        DO_comp1(gf2x_cantor_fft);
        DO_comp2(gf2x_cantor_fft);
        DO_ift_h(gf2x_cantor_fft);

        DO_zero_thx(gf2x_cantor_fft);
        DO_compx(gf2x_cantor_fft);
        DO_ift_hx(gf2x_cantor_fft);

        if (memcmp(gf2x_cantor_fft_h, gf2x_cantor_fft_hx, nwh * sizeof(unsigned long)) != 0) {
            fprintf(stderr, "cantor != cantorx for %zu*%zu ; seed=%u\n",
                    nf,ng,seed);
            printf("w:=%d;\n", (int) ULONG_BITS);
            display("f1",f1,nwf1);
            display("g1",g1,nwg1);
            display("f2",f2,nwf2);
            display("g2",g2,nwg2);
            display("gf2x_cantor_fft_h",gf2x_cantor_fft_h,nwh);
            display("gf2x_cantor_fft_hx",gf2x_cantor_fft_hx,nwh);
            EXTRA_DISPLAY(gf2x_cantor_fft);
            fflush(stdout);
            fflush(stderr);
            abort();
        }

        if (memcmp(gf2x_fake_fft_h, gf2x_cantor_fft_h, nwh * sizeof(unsigned long)) != 0) {
            fprintf(stderr, "fake != cantor for %zu*%zu ; seed=%u\n",
                    nf,ng,seed);
            printf("w:=%d;\n", (int) ULONG_BITS);
            display("f1",f1,nwf1);
            display("g1",g1,nwg1);
            display("f2",f2,nwf2);
            display("g2",g2,nwg2);
            display("gf2x_fake_fft_h",gf2x_fake_fft_h,nwh);
            display("gf2x_cantor_fft_h",gf2x_cantor_fft_h,nwh);
            EXTRA_DISPLAY(gf2x_fake_fft);
            EXTRA_DISPLAY(gf2x_cantor_fft);
            fflush(stdout);
            fflush(stderr);
            abort();
        }

        DO_dft_f1(gf2x_ternary_fft);
        DO_dft_g1(gf2x_ternary_fft);
        DO_dft_f2(gf2x_ternary_fft);
        DO_dft_g2(gf2x_ternary_fft);
        DO_zero_th(gf2x_ternary_fft);
        DO_comp1(gf2x_ternary_fft);
        DO_comp2(gf2x_ternary_fft);
        DO_ift_h(gf2x_ternary_fft);

        DO_zero_thx(gf2x_ternary_fft);
        DO_compx(gf2x_ternary_fft);
        DO_ift_hx(gf2x_ternary_fft);

        if (memcmp(gf2x_ternary_fft_h, gf2x_ternary_fft_hx, nwh * sizeof(unsigned long)) != 0) {
            fprintf(stderr, "gf2x_ternary_fft != gf2x_ternary_fftx for %zu*%zu ; seed=%u\n",
                    nf,ng,seed);
            printf("w:=%d;\n", (int) ULONG_BITS);
            display("f1",f1,nwf1);
            display("g1",g1,nwg1);
            display("f2",f2,nwf2);
            display("g2",g2,nwg2);
            display("gf2x_ternary_fft_h",gf2x_ternary_fft_h,nwh);
            display("gf2x_ternary_fft_hx",gf2x_ternary_fft_hx,nwh);
            EXTRA_DISPLAY(gf2x_ternary_fft);
            fflush(stdout);
            fflush(stderr);
            abort();
        }

        if (memcmp(gf2x_fake_fft_h, gf2x_ternary_fft_h, nwh * sizeof(unsigned long)) != 0) {
            fprintf(stderr, "fake != gf2x_ternary_fft for %zu*%zu ; seed=%u\n",
                    nf,ng,seed);
            printf("w:=%d;\n", (int) ULONG_BITS);
            display("f1",f1,nwf1);
            display("g1",g1,nwg1);
            display("f2",f2,nwf2);
            display("g2",g2,nwg2);
            display("gf2x_fake_fft_h",gf2x_fake_fft_h,nwh);
            display("gf2x_ternary_fft_h",gf2x_ternary_fft_h,nwh);
            EXTRA_DISPLAY(gf2x_fake_fft);
            EXTRA_DISPLAY(gf2x_ternary_fft);
            fflush(stdout);
            fflush(stderr);
            abort();
        }
    }
    LEAVE(gf2x_fake_fft)
    LEAVE(gf2x_cantor_fft)
    LEAVE(gf2x_ternary_fft)

    free(gf2x_fake_fft_h);
    free(gf2x_cantor_fft_h);
    free(gf2x_ternary_fft_h);
    free(g2);
    free(f2);
    free(g1);
    free(f1);

    return 0;
}


int main(int argc, char * argv[])
{
    if (argc < 1 || argc > 4)
        usage();
    unsigned int max_nbits_a = 1000;
    if (argc >= 2) max_nbits_a = atoi(argv[1]);
    unsigned int max_nbits_b = max_nbits_a;
    unsigned int seed;
    if (argc >= 3) 
        max_nbits_b = atoi(argv[2]);
    if (argc >= 4) {
        seed = atoi(argv[3]);
    } else {
        /* make this completely deterministic in debug mode */
#ifdef NDEBUG
        srandom(time(NULL));
#ifdef  _POSIX_C_SOURCE
        for(int i = 0 ; i < getpid() % 1009 ; i++) random();
#endif
#endif
        seed = random();
    }
    printf("// Seed is %u\n", seed);
    srandom(seed);

    docheck(max_nbits_a,max_nbits_b,10);

    for(unsigned int i = 0 ; i < 200 ; i++) {
        unsigned int na = 100 + (random() % max_nbits_a);
        unsigned int nb = 100 + (random() % max_nbits_b);
        docheck(na, nb,10);
        fputc('.',stderr); fflush(stderr);
    }
    printf("\n");
    for(unsigned int i = 0 ; i < max_nbits_a ; ) {
        // test some fancy sizes as well.
        i += 2 * ULONG_BITS / 3 + (random() % 3);
        docheck(i,i * max_nbits_b / max_nbits_a,10);
        fputc('.',stderr); fflush(stderr);
    }
    printf("\n");
    return 0;
}

