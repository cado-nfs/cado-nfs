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

/* 20080122 -- shortened the input file somewhat. Old disabled code has
 * been deleted in this operation.
 *
 *
 * XXX The inverse transform operation should be unified between the
 * truncated and non-truncated case.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

#include "gf2x/gf2x-config.h"
#include "gf2x/gf2x-impl.h"

#include "gf2x-cantor-fft.h"

/* The following flags affect the behaviour of the program */

#define xxxCANTOR_GM            /* Use Gao-Mateer recursion */
#define xxxCANTOR_GM_TRUNCATE   /* Use truncated variant */

#include "mpfq/mpfq_name_K.h"

/* It's a placeholder, really. After macro expansion, nobody really uses
 * this */
Kfield K;


/*
#define Kdst_elt        unsigned long *
#define Ksrc_elt        const unsigned long *
*/

// Kelt * fbase;
// size_t * findex;

// Some constants related to Cantor's algorithm
// Beta_i are such that Beta_{i-1} = Beta_i^2 + Beta_i
// and { Beta_i : i in [1..128] } is a GF(2)-basis of GF(2^128).
// The following is (the begining of) a solution.
// NB: the defining polynomial for GF(2^128) is x^128 + x^7 + x^2 + x + 1
//
//
/* Magma code:
 _<x>:=PolynomialRing(GF(2));
 K<z>:=ext<GF(2)|x^128 + x^7 + x^2 + x + 1>;
 beta:=K!1;
 for i in [1..31] do
      ti:=Polynomial([beta,1,1]);
      beta:=[r[1]:r in Roots(ti) | Eltseq(r[1])[1] eq 0][1];
      print Intseq(Seqint(ChangeUniverse(Eltseq(beta),Integers()),2),2^32);
 end for;

 // same deal for base field 2^64
 _<x>:=PolynomialRing(GF(2));
 K<z>:=ext<GF(2)|x^64 + x^4 + x^3 + x + 1>;
 beta:=K!1;
 for i in [1..31] do
      ti:=Polynomial([beta,1,1]);
      beta:=[r[1]:r in Roots(ti) | Eltseq(r[1])[1] eq 0][1];
      print Intseq(Seqint(ChangeUniverse(Eltseq(beta),Integers()),2),2^32);
 end for;

*/

#if CANTOR_BASE_FIELD_SIZE == 128
#if GF2X_WORDSIZE == 32
#define BETA(x,y,z,t) { x, y, z, t, }
#else
#define BETA(x,y,z,t) { y << 32 | x, t << 32 | z }
#endif
static const Kelt Betai[32] = {
    BETA(1UL,0UL,0UL,0UL),
    BETA(2763131656UL, 1735044255UL, 4101184249UL, 693813425UL),
    BETA(465473314UL, 4279278019UL, 3443545289UL, 2719171618UL),
    BETA(45213768UL, 1404590948UL, 3790498114UL, 2035930892UL),
    BETA(1766741572UL, 3401373924UL, 3912612066UL, 150385427UL),
    BETA(1642619116UL, 1296609638UL, 891822488UL, 1237922326UL),
    BETA(1432407564UL, 2705320108UL, 2845284892UL, 87917909UL),
    BETA(3169845528UL, 1555107770UL, 977112198UL, 1297233226UL),
    BETA(3766302226UL, 3784070584UL, 2372619472UL, 751494893UL),
    BETA(3319219312UL, 1485529668UL, 1919364198UL, 554536894UL),
    BETA(2221231524UL, 3494370010UL, 2246133939UL, 3255900204UL),
    BETA(4004748792UL, 4275287155UL, 1207038226UL, 273485854UL),
    BETA(4134924972UL, 3206127284UL, 3660212267UL, 3465895584UL),
    BETA(1810844824UL, 432364429UL, 3275102095UL, 3685546794UL),
    BETA(2239846096UL, 4115864360UL, 3774983060UL, 1484196804UL),
    BETA(572310756UL, 2341975643UL, 1205907074UL, 1233095387UL),
    BETA(374654328UL, 2701982826UL, 4048512218UL, 95057245UL),
    BETA(1069419604UL, 3135871566UL, 603196539UL, 2195404101UL),
    BETA(2819376738UL, 2390789826UL, 3206934123UL, 3126567229UL),
    BETA(2667030844UL, 1546355395UL, 894286960UL, 806061570UL),
    BETA(3589162292UL, 2284726595UL, 3987361705UL, 3261030005UL),
    BETA(488495394UL, 2513641547UL, 3724856370UL, 956902965UL),
    BETA(894705622UL, 3412990576UL, 2985357841UL, 2327012336UL),
    BETA(691430644UL, 2704377464UL, 1647445035UL, 3189688948UL),
    BETA(3013060094UL, 3867471049UL, 66167789UL, 3014841169UL),
    BETA(1101590000UL, 1839657742UL, 3478054432UL, 884001900UL),
    BETA(494046998UL, 2899707219UL, 2156601667UL, 2796513881UL),
    BETA(3398486548UL, 57111011UL, 1428951442UL, 894796186UL),
    BETA(3950556970UL, 3334574894UL, 3366886047UL, 3938361620UL),
    BETA(999564750UL, 4233655517UL, 706258300UL, 428364188UL),
    BETA(1913007172UL, 1894536580UL, 871798750UL, 1621883580UL),
    BETA(2252056148UL, 2850320680UL, 3455532848UL, 234262497UL),
};
#elif CANTOR_BASE_FIELD_SIZE == 64
#if GF2X_WORDSIZE == 32
#define BETA(x,y) { x, y, }
#else
#define BETA(x,y) { y << 32 | x }
#endif
static const Kelt Betai[32] = {
    BETA(1UL,0UL),
    BETA(663411714UL, 432617119UL),
    BETA(1868560276UL, 2709645270UL),
    BETA(3463994632UL, 1572356951UL),
    BETA(4123644106UL, 3111375974UL),
    BETA(2199521118UL, 1377486985UL),
    BETA(3720914518UL, 54323390UL),
    BETA(129568784UL, 3045354574UL),
    BETA(3011130092UL, 1082636491UL),
    BETA(2349523324UL, 13685640UL),
    BETA(595534370UL, 2950016112UL),
    BETA(3426396920UL, 3824529305UL),
    BETA(1779128504UL, 1511895985UL),
    BETA(4269728328UL, 160726100UL),
    BETA(6878532UL, 244231950UL),
    BETA(189484890UL, 2728465232UL),
    BETA(3603827632UL, 3916404741UL),
    BETA(1553915212UL, 1397321611UL),
    BETA(4104923806UL, 2905758902UL),
    BETA(3622692902UL, 1250807691UL),
    BETA(3057861136UL, 2750006610UL),
    BETA(3384573648UL, 1604922073UL),
    BETA(3447141598UL, 2743728400UL),
    BETA(3393587672UL, 3910595551UL),
    BETA(1272090518UL, 4170545265UL),
    BETA(595211090UL, 1447942745UL),
    BETA(3344624952UL, 2905848940UL),
    BETA(3440271156UL, 1470304335UL),
    BETA(3445621658UL, 2968847076UL),
    BETA(2448727948UL, 3947030431UL),
    BETA(1957951568UL, 1424342839UL),
    BETA(3960779986UL, 3060294886UL),
};
#else
#error  "Define CANTOR_BASE_FIELD_SIZE to 64 or 128"
#endif

// The S_i polynomials of Cantor's algorithm: linearized polynomials such
// that S_i | S_{i+1}.
static const size_t ind_S[31][15] = {
    {},                         // S0
    {1},                        // S1
    {1},                        // S2
    {4, 2, 1},                  // S3
    {1},                        // S4
    {16, 2, 1},                 // S5
    {16, 4, 1},                 // S6
    {64, 32, 16, 8, 4, 2, 1},   // S7
    {1},                        // S8
    {256, 2, 1},                // S9
    {256, 4, 1},                // S10
    {1024, 512, 256, 8, 4, 2, 1},       // S11
    {256, 16, 1},               // S12
    {4096, 512, 256, 32, 16, 2, 1},     // S13
    {4096, 1024, 256, 64, 16, 4, 1},    // S14
    {16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1},     // S15
    {1},                        // S16
    {65536, 2, 1},              // S17
    {65536, 4, 1},              // S18
    {262144, 131072, 65536, 8, 4, 2, 1},        // S19
    {65536, 16, 1},             // S20
    {1048576, 131072, 65536, 32, 16, 2, 1},     // S21
    {1048576, 262144, 65536, 64, 16, 4, 1},     // S22
    {4194304, 2097152, 1048576, 524288, 262144, 131072, 65536, 128, 64, 32, 16, 8, 4, 2, 1},    // S23
    {65536, 256, 1},            // S24
    {16777216, 131072, 65536, 512, 256, 2, 1},  // S25
    {16777216, 262144, 65536, 1024, 256, 4, 1}, // S26
    {67108864, 33554432, 16777216, 524288, 262144, 131072, 65536, 2048, 1024, 512, 256, 8, 4, 2, 1},    // S27
    {16777216, 1048576, 65536, 4096, 256, 16, 1},       // S28
    {268435456, 33554432, 16777216, 2097152, 1048576, 131072, 65536, 8192, 4096, 512, 256, 32, 16, 2, 1},       // S29
    {268435456, 67108864, 16777216, 4194304, 1048576, 262144, 65536, 16384, 4096, 1024, 256, 64, 16, 4, 1},     // S30
};

// The number of terms in S_i.
static const unsigned int ind_number[31] =
    { 0, 1, 1, 3, 1, 3, 3, 7, 1, 3, 3, 7, 3, 7, 7, 15,
    1, 3, 3, 7, 3, 7, 7, 15, 3, 7, 7, 15, 7, 15, 15
};

// At index i=i0+2*i1+4*i2+8*i3..., one computes 
//     allBetai[i] = i0*Betai[1] + i1*Betai[2] + i2*Betai[3] + ...
static inline void allBetai(Kdst_elt x, size_t i)
{
    size_t j;
    Kset_zero(x);
    j = 0;
    for( ; i ; i>>=1, j++) {
        if (!(i & 1)) continue;
        Kadd(x, x, Betai[j]);
    }
}

static inline size_t transform_datasize(gf2x_cantor_fft_info_srcptr p)
{
    size_t n = p->n;
#ifdef WITHOUT_CANTOR_TRUNCATION
    n = 1 << p->k;
#endif
    return n;
}

#ifdef  CANTOR_GM
#if 0
/* These versions would make sense if Kelt were defined to be _v2di *,
 * which it isn't */
/* TODO: anyway we would need to redefine all this to use __m128i's */
static inline void allBetai1(Kdst_elt x, size_t i)
{
    size_t j;
    Kset_zero(x);
    j = 0;
    __v2di m[2] = { (__v2di) {0,0}, (__v2di) { -1L, -1L } };
    for( ; i ; i>>=1, j++) {
            * (__v2di*)  x ^= m[i&1] & * (__v2di*) Betai[j];
    }
}

// compute omega_i and omega_{2i}
static inline void allBetai2(Kdst_elt x, Kdst_elt y, size_t i)
{
    size_t j;
    Kset_zero(x);
    Kset_zero(y);
    j = 0;
    __v2di m[2] = { (__v2di) {0,0}, (__v2di) { -1L, -1L } };
    for( ; i ; i>>=1, j++) {
            * (__v2di*)  x ^= m[i&1] & * (__v2di*) Betai[j];
            * (__v2di*)  y ^= m[i&1] & * (__v2di*) Betai[j+1];
    }
}
static inline void allBetai4(Kdst_elt x, Kdst_elt y, Kdst_elt z, Kdst_elt t, size_t i)
{
    size_t j;
    Kset_zero(x);
    Kset_zero(y);
    Kset_zero(z);
    Kset_zero(t);
    j = 0;
    __v2di m[2] = { (__v2di) {0,0}, (__v2di) { -1L, -1L } };
    for( ; i ; i>>=1, j++) {
            * (__v2di*) x ^= m[i&1] & * (__v2di*) Betai[j];
            * (__v2di*) y ^= m[i&1] & * (__v2di*) Betai[j+1];
            * (__v2di*) z ^= m[i&1] & * (__v2di*) Betai[j+2];
            * (__v2di*) t ^= m[i&1] & * (__v2di*) Betai[j+3];
    }
}
#else
static inline void allBetai1(Kdst_elt x, size_t i)
{
    size_t j;
    Kset_zero(x);
    j = 0;
    for( ; i ; i>>=1, j++) {
        if (!(i & 1)) continue;
        Kadd(x, x, Betai[j]);
    }
}
static inline void allBetai2(Kdst_elt x, Kdst_elt y, size_t i)
{
    size_t j;
    Kset_zero(x);
    Kset_zero(y);
    j = 0;
    for( ; i ; i>>=1, j++) {
        if (!(i & 1)) continue;
        Kadd(x, x, Betai[j]);
        Kadd(y, y, Betai[j+1]);
    }
}
static inline void allBetai4(Kdst_elt x, Kdst_elt y, Kdst_elt z, Kdst_elt t, size_t i)
{
    size_t j;
    Kset_zero(x);
    Kset_zero(y);
    Kset_zero(z);
    Kset_zero(t);
    j = 0;
    for( ; i ; i>>=1, j++) {
        if (!(i & 1)) continue;
        Kadd(x, x, Betai[j]);
        Kadd(y, y, Betai[j+1]);
        Kadd(z, z, Betai[j+2]);
        Kadd(t, t, Betai[j+3]);
    }
}
#endif  /* v2di or not v2di */
#endif  /* CANTOR_GM */

#ifdef  CANTOR_GM
static inline void expand(Kelt * f, unsigned int u, unsigned int t)
{
    // f has at most 2^u coefficients. Split in pieces of 2^t
    // coefficients.
    if (u == t) return;
    size_t K1 = 1UL << (u-1);
    size_t K0 = 1UL << (u-1-t);
    for (size_t i = 2 * K1 - 1; i >= K1; --i) {
        Kadd(f[i+K0-K1], f[i+K0-K1], f[i]); 
    }
    expand(f, u-1, t);
    expand(f+K1, u-1, t);
}
static inline void expand_trunc(Kelt * f, unsigned int u, unsigned int t, size_t n)
{
    // f has at most n coefficients. Split in pieces of 2^t
    // coefficients.
    size_t K1;
    for( ; n <= (K1 = (1UL << (u-1))) ; u--);
    if (u <= t) return;
    size_t K0 = 1UL << (u-1-t);
    for (size_t i = n - 1; i >= K1; --i) {
        Kadd(f[i+K0-K1], f[i+K0-K1], f[i]); 
    }
    expand(f, u-1, t);
    expand_trunc(f+K1, u-1, t, n-K1);
}

#if 0
static inline void expand_unroll(Kelt * f, unsigned int u, unsigned int t)
{
    // f has at most 2^u coefficients. Split in pieces of 2^t
    // coefficients.
    if (u == t) return;
    size_t K1 = 1UL << (u-1);
    size_t K0 = 1UL << (u-1-t);
    for (size_t i = 2 * K1 - 1; i >= K1; --i) {
        Kadd(f[i+K0-K1], f[i+K0-K1], f[i]); 
        fprintf(stderr, "Kadd(f[%ld],f[%ld],f[%ld]);\n",
                findex[i+K0-K1+f-fbase],
                findex[i+K0-K1+f-fbase],
                findex[i+f-fbase]);
    }
    expand_unroll(f, u-1, t);
    expand_unroll(f+K1, u-1, t);
}
#endif

void transpose_inplace(Kelt * f, unsigned int t)
{
#define Kswap(a,b) do { Kelt x; Kset(x,a); Kset(a,b); Kset(b,x); } while (0)
    // transpose
    for(size_t i = 0 ; i < (1UL<<t) ; i++) {
        for(size_t j = i+1 ; j < (1UL<<t) ; j++) {
            Kswap(f[(i<<t)+j],f[(j<<t)+i]);
        }
    }
}

// assuming f has 1UL<<t1 rows of 1UL<<t2 values, transpose into g. Set striding
// in g to be 1UL<<t1 values, for 1UL<<t2 rows.
//
// doing this in place would be somewhat tricky, but maybe reachable
// assuming that we stick to the case of powers of 2 for the dimensions.
void transpose_outofplace(Kelt * g, Kelt * f, unsigned int t1, unsigned int t2)
{
    for(size_t j = 0 ; j < (1UL<<t2) ; j++) {
        for(size_t i = 0 ; i < (1UL<<t1) ; i++) {
            Kset(g[i],f[i<<t2]);
        }
        g += 1UL << t1;
        f++;
    }
}

#if 0
void gm_trick_unroll(unsigned int two_t, Kelt * f, unsigned int j, size_t m, size_t d)
{
    assert((two_t & (two_t-1)) == 0);
    if (two_t == 1) {
        if (j) {
            Kelt z;
            allBetai(z,2*j);
            Kmul(z,z,f[1]);
            Kadd(f[0],z,f[0]);
            fprintf(stderr, "allBetai(z,%d*j+%d);\n",2*m,2*d);
            fprintf(stderr, "Kmul(z,z,f[%ld]);\n",findex[1 + f - fbase]);
            fprintf(stderr, "Kadd(f[%ld],z,f[%ld]);\n",findex[f - fbase], findex[f-fbase]);
        }
        Kadd(f[1],f[1],f[0]);
        fprintf(stderr, "Kadd(f[%ld],f[%ld],f[%ld]);\n",findex[f-fbase+1],findex[f-fbase+1], findex[f-fbase]);
        return;
    }

    unsigned int t = two_t >> 1;
    size_t tau = 1UL << t;
    expand_unroll(f,2*t,t);

    for(size_t i = 0 ; i < tau ; i++) {
        for(size_t j = i+1 ; j < tau ; j++) {
            Kswap(f[(i<<t)+j],f[(j<<t)+i]);
#if 1
            size_t idx = findex[(i<<t)+j+f-fbase];
            findex[(i<<t)+j+f-fbase] = findex[(j<<t)+i+f-fbase];
            findex[(j<<t)+i+f-fbase] = idx;
#else
            fprintf(stderr, "Kswap(f[%ld],f[%ld]);\n",
                    (i<<t)+j+f-fbase,
                    (j<<t)+i+f-fbase);
#endif
        }
    }
    
    // evaluate ; we can use f as a buffer, now.
    for(size_t l = 0 ; l < tau ; l++) {
        gm_trick_unroll(t, f + (l << t), j, m, d);
    }

    for(size_t i = 0 ; i < tau ; i++) {
        for(size_t j = i+1 ; j < tau ; j++) {
            Kswap(f[(i<<t)+j],f[(j<<t)+i]);
#if 1
            size_t idx = findex[(i<<t)+j+f-fbase];
            findex[(i<<t)+j+f-fbase] = findex[(j<<t)+i+f-fbase];
            findex[(j<<t)+i+f-fbase] = idx;
#else
            fprintf(stderr, "Kswap(f[%ld],f[%ld]);\n",
                    (i<<t)+j+f-fbase,
                    (j<<t)+i+f-fbase);
#endif
        }
    }
    
    // evaluate again
    for(size_t i = 0 ; i < tau ; i++) {
        gm_trick_unroll(t, f + (i << t), (j << t) + i, m << t, (d << t) + i);
    }
}
#endif

static inline void gm_trick2(Kdst_elt f0, Kdst_elt f1, Kdst_elt f2, Kdst_elt f3, Ksrc_elt w2, Ksrc_elt w4)
{
    Kelt z;
    // Kelt w2,w4;
    Kadd(f2,f2,f3);
    Kadd(f1,f1,f2);
    // allBetai2(w2,w4,2*j);
    Kmul(z,w2,f2); Kadd(f0,z,f0); Kadd(f2,f2,f0);
    Kmul(z,w2,f3); Kadd(f1,z,f1); Kadd(f3,f3,f1);
    Kmul(z,w4,f1); Kadd(f0,z,f0); Kadd(f1,f1,f0);
    Kadd(z,w4,Betai[1]);Kmul(z,z,f3); Kadd(f2,z,f2); Kadd(f3,f3,f2);
}
static inline void gm_trick2z(Kdst_elt f0, Kdst_elt f1, Kdst_elt f2, Kdst_elt f3)
{
    Kelt z;
    Kadd(z,f2,f3);                      
    Kadd(f3,f1,f2);                     
    Kadd(f2,f0,z);                      
    Kadd(f1,f1,f2);                     
    Kmul(z,Betai[1],f3);
    Kadd(f2,z,f2);
    Kadd(f3,f3,f2);
}

void gm_trick4(Kelt * f, size_t j)
{
    Kelt f0; Kelt f1; Kelt f2; Kelt f3; 
    Kelt f4; Kelt f5; Kelt f6; Kelt f7; 
    Kelt f8; Kelt f9; Kelt fa; Kelt fb;
    Kelt fc; Kelt fd; Kelt fe; Kelt ff;

    Kset(f0,f[0]);	Kset(f1,f[1]);	Kset(f2,f[2]);	Kset(f3,f[3]);
    Kset(f4,f[4]);	Kset(f5,f[5]);	Kset(f6,f[6]);	Kset(f7,f[7]);
    Kset(f8,f[8]);	Kset(f9,f[9]);	Kset(fa,f[10]); Kset(fb,f[11]);
    Kset(fc,f[12]);	Kset(fd,f[13]); Kset(fe,f[14]); Kset(ff,f[15]);

    Kadd(f9,f9,ff);	Kadd(f8,f8,fe); Kadd(f7,f7,fd);
    Kadd(f6,f6,fc);	Kadd(f5,f5,fb); Kadd(f4,f4,fa);
    Kadd(f3,f3,f9);	Kadd(f2,f2,f8);	Kadd(f4,f4,f7);
    Kadd(f3,f3,f6);	Kadd(f2,f2,f5);	Kadd(f1,f1,f4);
    Kadd(fc,fc,ff);	Kadd(fb,fb,fe);	Kadd(fa,fa,fd);
    Kadd(f9,f9,fc);

    Kelt w2;
    Kelt w4;
    Kelt w8;
    Kelt w16;
    allBetai4(w2,w4,w8,w16,j<<1);

    gm_trick2(f0,f4,f8,fc,w2,w4);
    gm_trick2(f1,f5,f9,fd,w2,w4);
    gm_trick2(f2,f6,fa,fe,w2,w4);
    gm_trick2(f3,f7,fb,ff,w2,w4);
    gm_trick2(f0,f1,f2,f3,w8,w16);
    Kadd(w2,w8,Betai[1]);
    Kadd(w4,w16,Betai[2]);
    gm_trick2(f4,f5,f6,f7,w2,w4);
    Kadd(w8,w8,Betai[2]);
    Kadd(w16,w16,Betai[3]);
    gm_trick2(f8,f9,fa,fb,w8,w16);
    Kadd(w2,w2,Betai[2]);
    Kadd(w4,w4,Betai[3]);
    gm_trick2(fc,fd,fe,ff,w2,w4);
    Kset(f[0],f0);	Kset(f[1],f1);	Kset(f[2],f2);	Kset(f[3],f3);
    Kset(f[4],f4);	Kset(f[5],f5);	Kset(f[6],f6);	Kset(f[7],f7);
    Kset(f[8],f8);	Kset(f[9],f9);	Kset(f[10],fa); Kset(f[11],fb);
    Kset(f[12],fc);	Kset(f[13],fd); Kset(f[14],fe); Kset(f[15],ff);
}
void gm_trick4z(Kelt * f)
{
    Kelt f0; Kelt f1; Kelt f2; Kelt f3; 
    Kelt f4; Kelt f5; Kelt f6; Kelt f7; 
    Kelt f8; Kelt f9; Kelt fa; Kelt fb;
    Kelt fc; Kelt fd; Kelt fe; Kelt ff;

    Kset(f0,f[0]);	Kset(f1,f[1]);	Kset(f2,f[2]);	Kset(f3,f[3]);
    Kset(f4,f[4]);	Kset(f5,f[5]);	Kset(f6,f[6]);	Kset(f7,f[7]);
    Kset(f8,f[8]);	Kset(f9,f[9]);	Kset(fa,f[10]); Kset(fb,f[11]);
    Kset(fc,f[12]);	Kset(fd,f[13]); Kset(fe,f[14]); Kset(ff,f[15]);

    Kadd(f9,f9,ff);	Kadd(f8,f8,fe); Kadd(f7,f7,fd);
    Kadd(f6,f6,fc);	Kadd(f5,f5,fb); Kadd(f4,f4,fa);
    Kadd(f3,f3,f9);	Kadd(f2,f2,f8);	Kadd(f4,f4,f7);
    Kadd(f3,f3,f6);	Kadd(f2,f2,f5);	Kadd(f1,f1,f4);
    Kadd(fc,fc,ff);	Kadd(fb,fb,fe);	Kadd(fa,fa,fd);
    Kadd(f9,f9,fc);

    Kelt w2;
    Kelt w4;

    gm_trick2z(f0,f4,f8,fc);
    gm_trick2z(f1,f5,f9,fd);
    gm_trick2z(f2,f6,fa,fe);
    gm_trick2z(f3,f7,fb,ff);
    gm_trick2z(f0,f1,f2,f3);
    gm_trick2(f4,f5,f6,f7,Betai[1],Betai[2]);
    gm_trick2(f8,f9,fa,fb,Betai[2],Betai[3]);
    Kadd(w2,Betai[1],Betai[2]);
    Kadd(w4,Betai[2],Betai[3]);
    gm_trick2(fc,fd,fe,ff,w2,w4);

    Kset(f[0],f0);	Kset(f[1],f1);	Kset(f[2],f2);	Kset(f[3],f3);
    Kset(f[4],f4);	Kset(f[5],f5);	Kset(f[6],f6);	Kset(f[7],f7);
    Kset(f[8],f8);	Kset(f[9],f9);	Kset(f[10],fa); Kset(f[11],fb);
    Kset(f[12],fc);	Kset(f[13],fd); Kset(f[14],fe); Kset(f[15],ff);
}


void gm_trick(unsigned int two_t, Kelt * f, size_t j)
{
    assert((two_t & (two_t-1)) == 0);
    if (two_t == 4) {
        if (j) gm_trick4(f,j); else gm_trick4z(f);
        return;
    } else if (two_t == 2) {
        if (j) {
            Kelt w2,w4;
            allBetai2(w2,w4,2*j);
            gm_trick2(f[0],f[1],f[2],f[3],w2,w4);
        } else {
            gm_trick2z(f[0],f[1],f[2],f[3]);
        }
        return;
    } else if (two_t == 1) {
        if (j) {
            Kelt z;
            allBetai(z,2*j);
            Kmul(z,z,f[1]);
            Kadd(f[0],z,f[0]);
        }
        Kadd(f[1],f[1],f[0]);
        return;
    }
#if 0
    if (two_t == 4) {
        fprintf(stderr, "// Unrolling gm_trick(%d,f,%d)\n",two_t, j);
        fbase = f;
        /* beware ; on 32-bit machines, two_t might be equal to 32... */
        findex = malloc(sizeof(long int) << two_t);
        for(int i = 0 ; i < (1UL << two_t) ; i++) {
            findex[i] = i;
        }
        gm_trick_unroll(two_t, f, j, 1, 0);
        free(findex);
        findex = NULL;
        fbase = NULL;
        fprintf(stderr, "// Done unrolling\n\n");
        return;
    }
#endif

    unsigned int t = two_t >> 1;
    size_t tau = 1UL << t;

    expand(f,2*t,t);

    transpose_inplace(f,t);
    
    // evaluate ; we can use f as a buffer, now.
    for(size_t l = 0 ; l < tau ; l++) {
        gm_trick(t, f + (l << t), j);
    }
    transpose_inplace(f,t);

    // evaluate again
    for(size_t i = 0 ; i < tau ; i++) {
        gm_trick(t, f + (i << t), (j << t) + i);
    }
}

void multievaluateKrec(Kelt * f, unsigned int i, size_t rep_beta);

/* This is the same gao-mateer trick. It reduces a length 2^two_t
 * polynomial to 2^t polynomials of length 2^t. Recursive work for the
 * inner evaluation is handled
 * by multievaluateKrec (which in turn may call us). The outer
 * evaluations are also processed the same way, and also handed over to
 * multievaluateKrec
 */
void multievaluateKrecGM(Kelt * f, unsigned int two_t, size_t j)
{
    assert((two_t & (two_t-1)) == 0);
    unsigned int t = two_t >> 1;
    size_t tau = 1UL << t;
    expand(f,2*t,t);
    transpose_inplace(f,t);
    // evaluate ; we can use f as a buffer, now.
    for(size_t l = 0 ; l < tau ; l++) {
        multievaluateKrec(f + (l << t), t, j);
    }
    transpose_inplace(f,t);
    // evaluate again
    for(size_t i = 0 ; i < tau ; i++) {
        multievaluateKrec(f + (i << t), t, (j << t) + i);
    }
}
#endif

#ifdef  CANTOR_GM_TRUNCATE
// compute only n' = 1UL<<k values -- f must not have more than n coefficients !
//
// n' is larger than or equal to n.
//
// f and buf must hold space for 1UL << k coefficients.
//
// The trailing values in f must be equal to zero.
void gm_trick_trunc(unsigned int two_t, Kelt * f, Kelt * buf, size_t j, unsigned int k, size_t n)
{
    assert((two_t & (two_t-1)) == 0);

    if (n <= 1)
        return;
    if (two_t == 1) {
        if (j) {
            Kelt z;
            allBetai(z,2*j);
            Kmul(z,z,f[1]);
            Kadd(f[0],z,f[0]);
        }
        Kadd(f[1],f[1],f[0]);
        return;
    }
    assert(n <= (1UL << k));
    if (k <= (two_t>>1)) {
        gm_trick_trunc(two_t>>1, f, buf, j, k, n);
        return;
    }

    unsigned int t = two_t >> 1;

    assert(k > t);
    assert(k <= two_t);
    assert(n > (1UL << t));
    /* beware ; on 32-bit machines, two_t might be equal to 32... */
    assert(two_t == GF2X_WORDSIZE || n <= (1UL << two_t));

    size_t tau = 1UL << t;
    // long cn = (n+tau-1)>>t;      // Ceiling(n/tau);
    size_t cn = 1UL << (k-t);

    expand_trunc(f,2*t,t,n);
    memset(f + n, 0, ((1UL<<k)-n) * sizeof(Kelt));

    // need 1UL<<k in buf
    transpose_outofplace(buf,f,k-t,t);
    
    // evaluate ; we can use f as a buffer, now.
    Kelt * h = buf;
    for(size_t l = 0 ; l < tau ; l++) {
        // we're going to need 1UL << (k-t) in f, no big deal.
        gm_trick_trunc(t, h, f, j, k-t, cn);
        h += cn;
    }
    transpose_outofplace(f,buf,t,k-t);

    // evaluate again. Here we can8t truncate.
    for(size_t i = 0 ; i < cn ; i++) {
        gm_trick(t, f + (i << t), (j << t) + i);
    }
}
#endif


/* 20080122 -- hard-coded versions of reduceSi for small sizes have been
 * tried, but do not seem to pay off. Look for the svn data to
 * investigate again the potential benefit.
 */

// (sk+beta)*(sk+beta+1) s(sk+beta)=s(sk)+s(beta)=s_(k+1) + s(beta)
// The generic reduction function.
static inline void reduceSi(unsigned int k, Kelt * f, Kelt beta)
{
    size_t i, K1;
    unsigned int j;
    Kelt coeff;
    Ksrc_elt fi;

    K1 = 1UL << k;
    // put rest mod Sk+beta in low part of f and quo in high part.
    for (i = 2 * K1 - 1; i >= K1; --i) {
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            size_t index = ind_S[k][j] + i - K1;
            Kadd(f[index], f[index], fi);
        }
        Kmul(coeff, fi, beta);
        Kadd(f[i - K1], f[i - K1], coeff);
    }

    // The other rest is just the sum of low and hi parts
    for (i = 0; i < K1; ++i)
        Kadd(f[i + K1], f[i + K1], f[i]);
}

#ifdef  CANTOR_GM
// same as multieval, but only to chunks of size kappa
void reduce_top_cantor(unsigned int k, unsigned int kappa, Kelt * f, size_t j)
{
    if (k == kappa)
        return;
    // fprintf(stderr, "top-reduce %d %d\n",k,kappa);
    j *= 2;
    Kelt x;
    allBetai(x, j);
    reduceSi(k-1, f, x);
    reduce_top_cantor(k-1,kappa,f,j);
    reduce_top_cantor(k-1,kappa,f+(1L<<(k-1)),j+1);
}
#endif


// The generic reduction function.
// truncated version
static inline void reduceSiNobeta_trunc(unsigned int k, Kelt * f, size_t L)
{
    size_t j;
    size_t i, K1;
    Ksrc_elt fi;

    K1 = 1UL << k;
    // put rest mod Sk+beta in low part of f and quo in high part.
    for (i = L - 1; i >= K1; --i) {
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            long index = ind_S[k][j] + i - K1;
            Kadd(f[index], f[index], fi);
        }
    }

    // The other rest is just the sum of low and hi parts
    for (i = 0; i < K1; ++i)
        Kadd(f[i + K1], f[i + K1], f[i]);
}

// A generic interpolation function.
// TODO: unroll small sizes, like for reduction.
static inline void interpolateSi(unsigned int k, Kelt * f, Kelt beta)
{
    unsigned int j;
    size_t i, K1 = 1UL << k;

    for (i = 0; i < K1; ++i)
        Kadd(f[i + K1], f[i + K1], f[i]);

    for (i = K1; i < K1 << 1; ++i) {
        Kelt coeff;
        Ksrc_elt fi;
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            size_t index = ind_S[k][j] + i - K1;
            Kadd(f[index], f[index], fi);
        }
        Kmul(coeff, fi, beta);
        Kadd(f[i - K1], f[i - K1], coeff);
    }
}

static inline void interpolateSiNobeta(unsigned int k, Kelt * f)
{
    unsigned int j;
    size_t i, K1 = 1UL << k;

    for (i = 0; i < K1; ++i)
        Kadd(f[i + K1], f[i + K1], f[i]);

    for (i = K1; i < K1 << 1; ++i) {
        Ksrc_elt fi;
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            size_t index = ind_S[k][j] + i - K1;
            Kadd(f[index], f[index], fi);
        }
    }
}

// f has length 2^i
void multievaluateKrec(Kelt * f, unsigned int i, size_t rep_beta)
{
    Kelt beta;

#if 0
#ifdef  CANTOR_GM
    assert(i < 32);
    if ((i&-i)==i) {
        static int n[17]={0,};
        n[i]++;
        if ((n[i]&-n[i])==n[i]) {
            fprintf(stderr, "call for i==%d (%d-th)\n",i,n[i]);
        }
        if (i==4) {
            gm_trick4(f, rep_beta);
        } else {
            multievaluateKrecGM(f, i, rep_beta);
        }
        return;
    }
#endif
#endif

    rep_beta <<= 1;
    allBetai(beta, rep_beta);
    reduceSi(i - 1, f, beta);

    if (i == 2) {
        rep_beta <<= 1;
        allBetai(beta, rep_beta);
        reduceSi(0, f, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(0, f + 2, beta);
        return;
    } else if (i == 3) {
        rep_beta <<= 1;
        allBetai(beta, rep_beta);
        reduceSi(1, f, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(1, f + 4, beta);

        rep_beta <<= 1;
        allBetai(beta, rep_beta);
        reduceSi(0, f, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(0, f + 2, beta);
        Kadd(beta, beta, Betai[2]);
        reduceSi(0, f + 6, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(0, f + 4, beta);
        return;
    }
    if (i > 1) {
        {
            multievaluateKrec(f, i - 1, rep_beta);
            multievaluateKrec(f + (1UL << (i - 1)), i - 1, rep_beta + 1);
        }
    }
}

#ifdef  CANTOR_GM
// f must have 2^k coeffs exactly
void multievaluateGM(Kelt * f, unsigned int k, size_t length GF2X_MAYBE_UNUSED)
{
    unsigned int t = 1;
    unsigned int two_t;
    for( ; (two_t = t << 1) < k ; t = two_t) ;
    assert(t < k && k <= two_t);
#ifdef CANTOR_GM_TRUNCATE
    Kelt * buf = malloc(sizeof(Kelt) << k);
    gm_trick_trunc(two_t,f,buf,0,k,length);
    free(buf);
#else
    reduce_top_cantor(k,t,f,0);
    size_t K1 = 1UL << k;
    size_t tau = 1UL << t;
    size_t j = 0;
    for(size_t u = 0 ; u < K1 ; u += tau) {
        // treat coeffs [u..u+tau[ = [u..u+2^t[
        gm_trick(t, f + u, j);
        j++;
    }
#endif
}
#endif

// evaluate f, having _length_ coefficients, at the 2^i roots following
// rep_beta
void multievaluateKrec_trunc(Kelt * f, unsigned int i, size_t rep_beta, size_t length)
{
    Kelt beta;

    assert(i >= 1);

    rep_beta <<= 1;
    allBetai(beta, rep_beta);
    reduceSi(i - 1, f, beta);

    if (i == 2) {
        rep_beta <<= 1;
        allBetai(beta, rep_beta);
        reduceSi(0, f, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(0, f + 2, beta);
        return;
    }
    if (i == 3) {
        rep_beta <<= 1;
        allBetai(beta, rep_beta);
        reduceSi(1, f, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(1, f + 4, beta);

        rep_beta <<= 1;
        allBetai(beta, rep_beta);
        reduceSi(0, f, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(0, f + 2, beta);
        Kadd(beta, beta, Betai[2]);
        reduceSi(0, f + 6, beta);
        Kadd(beta, beta, Betai[1]);
        reduceSi(0, f + 4, beta);
        return;
    }
    if (i > 1) {
        if (length < (1UL << (i - 1))) {
            multievaluateKrec_trunc(f, i - 1, rep_beta, length);
        } else if (length > (1UL << (i - 1))) {
            multievaluateKrec(f, i - 1, rep_beta);
            multievaluateKrec_trunc(f + (1UL << (i - 1)), i - 1, rep_beta + 1,
                                    length - (1UL << (i - 1)));
        } else { // length == 1UL<<(i-1)
            multievaluateKrec(f, i - 1, rep_beta);
        }
    }
}

// evaluate f, having length coefficients, at all the roots of the polynomial
// s_k
void multievaluateKnew_trunc(Kelt * f, unsigned int k, size_t length)
{
    assert(length <= (1UL << k));
    if (k == 0) {
        // s_0 is x. We're evaluating a constant polynomial at a point, which
        // turns out to be 0, but that hardly matters. So there's actually
        // _nothing_ at all to do.
        return;
    } else if (k == 1) {
        // we need to evaluate at 0 and 1. So we must basically add f[0] to
        // f[1], unless there is no f1.
        if (length <= 1)
            return;
        else
            Kadd(f[1], f[0], f[1]);
        return;
    }

    reduceSiNobeta_trunc(k - 1, f, length);

    multievaluateKrec(f, k - 1, 0);
    multievaluateKrec_trunc(f + (1UL << (k - 1)), k - 1, 1,
                            length - (1UL << (k - 1)));
}

// Interpolation with subproduct tree.
void interpolateK(Kelt * f, unsigned int k)
{
    size_t i, j;
    Kelt beta;

    for (i = 0; i < k - 1; ++i) {
        for (j = 0; j < (1UL << (k - 1 - i)); ++j) {
            size_t index = j * (1UL << (i + 1));
            allBetai(beta, 2 * j);
            interpolateSi(i, f + index, beta);
        }
    }
    interpolateSiNobeta(k - 1, f);
}


// Quotrem Si. f is possibly larger than 2 times 1UL<<k.
// Degree of f is given by L (with L > 1UL<<k)
// Inplace computation: low part of f receive rem, high part gets quo.
static inline void quotremSi(unsigned int k, Kelt * f, size_t L, Kelt beta)
{
    size_t j;
    size_t i, K1;
    Kelt coeff;
    Ksrc_elt fi;

    K1 = 1UL << k;
    // put rest mod Sk+beta in low part of f and quo in high part.
    for (i = L - 1; i >= K1; --i) {
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            size_t index = ind_S[k][j] + i - K1;
            Kadd(f[index], f[index], fi);
        }
        Kmul(coeff, fi, beta);
        Kadd(f[i - K1], f[i - K1], coeff);
    }
}

// Mul Si. f is possibly smaller than half of 1UL<<k.
// Degree of f is given by L (with L<= (1UL<<k))
// Inplace computation: Top part of the result is stored in f, and low
// part of the result gets added to what is at the right of f.
static inline void mulSi(unsigned int k, Kelt * f, size_t L, Kelt beta)
{
    unsigned int j;
    size_t i, K1 = 1UL << k;

    for (i = 0; i < L; ++i) {
        Kelt coeff;
        Ksrc_elt fi;
        fi = f[i];
        for (j = 0; j < ind_number[k]; ++j) {
            size_t index = ind_S[k][j] + i - K1;
            Kadd(f[index], f[index], fi);
        }
        Kmul(coeff, fi, beta);
        Kadd(f[i - K1], f[i - K1], coeff);
    }
}

#if GNUC_VERSION_ATLEAST(3,4,0)
#define gf2x_clzl(x)         __builtin_clzl(x)
#define gf2x_ctzl(x)         __builtin_ctzl(x)
#else
/* provide slow fallbacks */
static inline int gf2x_clzl(unsigned long x)
{
        static const int t[4] = { 2, 1, 0, 0 };
        int a = 0;
        int res;
#if (GF2X_WORDSIZE == 64)
        if (x >> 32) { a += 32; x >>= 32; }
#endif  
        if (x >> 16) { a += 16; x >>= 16; }
        if (x >>  8) { a +=  8; x >>=  8; }
        if (x >>  4) { a +=  4; x >>=  4; }
        if (x >>  2) { a +=  2; x >>=  2; }
        res = GF2X_WORDSIZE - 2 - a + t[x];
        return res;
}

/* the following code is correct because if x = 0...0abc10...0, then
   -x = ~x + 1, where ~x = 1...1(1-a)(1-b)(1-c)01...1, thus
   -x = 1...1(1-a)(1-b)(1-c)10...0, and x & (-x) = 0...000010...0 */
static inline int gf2x_ctzl(unsigned long x)
{
  ASSERT(GF2X_WORDSIZE == sizeof(unsigned long) * CHAR_BIT);
  return (GF2X_WORDSIZE - 1) - gf2x_clzl(x & - x);
}
#endif

// Reduce f (with 1UL<<k coeffs) modulo the product of s_i corresponding
// to length.
static inline void reduceModTrunc(Kelt * f, unsigned int k, size_t length)
{
    unsigned int i;
    size_t len, ii;
    len = length;
    Kelt *pf = f;
    Kelt beta;
    size_t offset_f;

    // go down, computing quotrems
    offset_f = 0;
    assert(k);
    for (i = k ; i-- > 0 ; ) {
        if (len >> i) {
            allBetai(beta, offset_f >> i);
            quotremSi(i, f + offset_f, (1UL << k) - offset_f, beta);
            offset_f += 1UL << i;
            len -= 1UL << i;
        }
    }

    // go up, reconstructing general rem
    len = length;
    i = gf2x_ctzl(len);
    len >>= i;
    pf = f + length;
    ii = 1UL << i;
    i++;
    len >>= 1;
    for (; i <= k - 1; ++i) {
        if (len & 1) {
            allBetai(beta, (length - (1UL << i) - ii) >> i);
            mulSi(i, pf - ii, ii, beta);
            ii += 1UL << i;
        }
        len >>= 1;
    }
}


// Interpolation with subproduct tree.
// This is a truncated version, where we reconstruct length coefficients,
// from length values. (length <= 1UL<<k)
// Assume that non-computed values are set to zero.
void interpolateK_trunc(Kelt * f, unsigned int k, size_t length)
{
    unsigned int i;
    size_t j;
    Kelt beta;

    if (k == 0) {
        return;
    } else if (k == 1) {
        // we know f(0) and f(1). They are related to the calues we want by:
        // f(0) = f[0]
        // f(1) = f[0] ^ f[1]
        // so f[0] stays unchanged, and f[1] gets xored w/ f[0].
        if (length <= 1)
            return;
        Kadd(f[1],f[1],f[0]);
        return;
    }
    for (i = 0; i < k - 1; ++i) {
        size_t maxj = 1UL << (k - 1 - i);
        for (j = 0; j < maxj; ++j) {
            size_t index = j * (1UL << (i + 1));
            if (length <= (index + (1UL << i))) {
                break;
            }
            allBetai(beta, 2 * j);
            interpolateSi(i, f + index, beta);
        }
    }
    if (length > (1UL << (k - 1)))
        interpolateSiNobeta(k - 1, f);

    reduceModTrunc(f, k, length);
}

#if (GF2X_WORDSIZE == 64)
#if CANTOR_BASE_FIELD_SIZE == 128
void decomposeK(Kelt * f, const unsigned long * F, size_t Fl, int k)
{
    assert(Fl <= (1UL << k));
    for (size_t i = 0; i < Fl; ++i) {
        /* Set low part of 128 bits from 64-bits in F[i] */
        f[i][0] = F[i];
        f[i][1] = 0;
    }
    memset(f + Fl, 0, ((1UL << k) - Fl) * sizeof(Kelt));
}

void recomposeK(unsigned long * F, Kelt * f, size_t Fl, int k GF2X_MAYBE_UNUSED)
{
    assert(Fl <= (1UL << k));
    F[0] = f[0][0];
    for (size_t i = 1; i < Fl ; ++i)
        F[i] = f[i][0] ^ f[i - 1][1];
}
#elif CANTOR_BASE_FIELD_SIZE == 64
void decomposeK(Kelt * f, const unsigned long * F, size_t Fl, int k)
{
    size_t i;
    /* We're computing a DFT of length 2^k, so we can accomodate 2^k*32
     * bits in the coefficients. We need to make Fl*64 bits fit. Hence
     * the following assertion.
     */
    assert(Fl <= (1UL << (k-1)));
    for (i = 0; i < 2*Fl ; i+=2) {
        /* Set low part of 64 bits from low 32 bits in F[i/2] */
        f[i][0] = F[i/2] & ((1UL<<32)-1);
        /* Set low part of 64 bits from low 32 bits in F[i/2] */
        f[i+1][0] = F[i/2] >> 32;
    }
    memset(f + i, 0, ((1UL << k) - 2*Fl) * sizeof(Kelt));
}

void recomposeK(unsigned long * F, Kelt * f, size_t Fl, int k GF2X_MAYBE_UNUSED)
{
    size_t i;
    assert(Fl <= (1UL << (k-1)));
    i=0;
        F[i] = f[2*i][0] ^ (f[2*i+1][0] << 32) ;
    for (i = 1; i < Fl ;i++) {
        F[i] = f[2*i][0] ^ (f[2*i+1][0] << 32) ^ (f[2*i-1][0] >> 32);
    }
}
#else
#error  "Define CANTOR_BASE_FIELD_SIZE to 64 or 128"
#endif

#elif (GF2X_WORDSIZE == 32)
#if CANTOR_BASE_FIELD_SIZE == 128
void decomposeK(Kelt * f, const unsigned long * F, size_t Fl, int k)
{
    size_t i;
    assert(Fl <= (1UL << (k+1)));
    for (i = 0; i < Fl / 2 ; ++i) {
        /* Set low part of 128 bits from 64-bits in F[i] */
        f[i][0] = F[2*i];
        f[i][1] = F[2*i + 1];
        f[i][2] = 0;
        f[i][3] = 0;
    }
    if (Fl & 1) {
        /* Fix the tail, too */
        f[i][0] = F[2*i];
        f[i][1] = 0;
        f[i][2] = 0;
        f[i][3] = 0;
        i++;
    }
    memset(f + i, 0, ((1UL << k) - i) * sizeof(Kelt));
}

void recomposeK(unsigned long * F, Kelt * f, size_t Fl, unsigned int k GF2X_MAYBE_UNUSED)
{
    size_t i;

    assert(Fl <= (1UL << (k+1)));
    F[0] = f[0][0];
    F[1] = f[0][1];
    for (i = 2; i + 1 < Fl; i += 2) {
        F[i] = f[i/2][0] ^ f[i/2 - 1][2];
        F[i+1] = f[i/2][1] ^ f[i/2 - 1][3];
    }
    if (Fl & 1) {
            F[i] = f[i/2][0] ^ f[i/2 - 1][2];
    }
}
#elif CANTOR_BASE_FIELD_SIZE == 64
void decomposeK(Kelt * f, const unsigned long * F, size_t Fl, int k)
{
    /* We're computing a DFT of length 2^k, so we can accomodate 2^k*32
     * bits in the coefficients. We need to make Fl*32 bits fit. Hence
     * the following assertion.
     */
    assert(Fl <= (1UL << k));
    for (size_t i = 0; i < Fl; ++i) {
        /* Set low part of 64 bits from 32-bits in F[i] */
        f[i][0] = F[i];
        f[i][1] = 0;
    }
    memset(f + Fl, 0, ((1UL << k) - Fl) * sizeof(Kelt));
}

void recomposeK(unsigned long * F, Kelt * f, size_t Fl, int k GF2X_MAYBE_UNUSED)
{
    assert(Fl <= (1UL << k));
    F[0] = f[0][0];
    for (size_t i = 1; i < Fl ; ++i)
        F[i] = f[i][0] ^ f[i - 1][1];
}
#else
#error  "Define CANTOR_BASE_FIELD_SIZE to 64 or 128"
#endif
#else
#error "define GF2X_WORDSIZE"
#endif

/* nF is a number of coefficients == number of bits ; a.k.a. degree + 1 */
void gf2x_cantor_fft_init(gf2x_cantor_fft_info_t p, size_t nF, size_t nG, ...)
{
    unsigned int k;
    size_t Hl;
    size_t n;

    /* Since internally we're working with 64-bit data, then it's really
     * a hard 64 here, not GF2X_WORDSIZE : We're just deciding on the order
     * of things.
     */
    const int keep_bits_per_coeff = CANTOR_BASE_FIELD_SIZE / 2;
    size_t Fl = (nF + keep_bits_per_coeff - 1) / keep_bits_per_coeff;
    size_t Gl = (nG + keep_bits_per_coeff - 1) / keep_bits_per_coeff;

    n = Hl = Fl + Gl;               // nb of Kelt of the result, with padding.
    for(k = 1; (1UL << k) < n ; k++) ;
    /* We used to refuse k < 2 here. Now we've got sufficient provision
     * in here to accomodate for the case where k==1, so a safe fallback
     * works. */
    p->k = k;
    p->n = n;

}

void gf2x_cantor_fft_clear(gf2x_cantor_fft_info_t p)
{
    memset(p, 0, sizeof(gf2x_cantor_fft_info_t));
}

/* nF is a number of coefficients */
void gf2x_cantor_fft_dft(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_ptr x, const unsigned long * F, size_t nF)
{
    size_t Fl = (nF + GF2X_WORDSIZE - 1) / GF2X_WORDSIZE;
    if (nF % GF2X_WORDSIZE) {
        /* Just as we are computing this assertion, we could easily mask out
         * the bits ourselves. However, our interface mandates that the high
         * bits be cleared in any case. So make sure we properly enforce this
         * constraint.
         */
        assert((F[Fl-1] & ~((1UL << (nF % GF2X_WORDSIZE)) - 1)) == 0);
    }
#ifdef WITHOUT_CANTOR_TRUNCATION
    memset(x, 0, sizeof(Kelt) << p->k);
#endif
    decomposeK(x,F,Fl,p->k);
#ifdef  CANTOR_GM
    multievaluateGM(x, p->k, p->n);
#else
#ifdef WITHOUT_CANTOR_TRUNCATION
    multievaluateKrec(x, p->k, 0);
#else
    multievaluateKnew_trunc(x, p->k, p->n);
#endif
#endif
}


void gf2x_cantor_fft_compose(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_ptr y, gf2x_cantor_fft_srcptr x1, gf2x_cantor_fft_srcptr x2)
{
    for (size_t j = 0; j < transform_datasize(p) ; j++) {
        Kmul(y[j], x1[j], x2[j]);
    }
}

void gf2x_cantor_fft_addcompose_n(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_ptr y, gf2x_cantor_fft_srcptr * x1, gf2x_cantor_fft_srcptr * x2, size_t n)
{
    Kelt er;
#if 0
    /* strategy 1: run all pointers simultaneously, and run many
     * simultaneous accumulators */
#define ACCUMULATE      256
    size_t j;
    for (j = 0; j + ACCUMULATE - 1 < transform_datasize(p) ; j+= ACCUMULATE) {
        Kelt_ur e[ACCUMULATE], s[ACCUMULATE];
        for(int t = 0 ; t < ACCUMULATE ; t++)
            Kelt_ur_set_zero(s[t]);
        for(size_t k = 0 ; k < n ; k++) {
            for(int t = 0 ; t < ACCUMULATE ; t++) {
                Kmul_ur(e[t], x1[k][j+t], x2[k][j+t]);
                Kelt_ur_add(s[t], s[t], e[t]);
            }
        }
        for(int t = 0 ; t < ACCUMULATE ; t++) {
            Kreduce(er, s[t]);
            Kadd(y[j+t], y[j+t], er);
        }
    }
    for (; j < transform_datasize(p) ; j++) {
        Kelt_ur e, s;
        Kelt_ur_set_zero(s);
        for(size_t k = 0 ; k < n ; k++) {
            Kmul_ur(e, x1[k][j], x2[k][j]);
            Kelt_ur_add(s, s, e);
        }
        Kreduce(er, s);
        Kadd(y[j], y[j], er);
    }
#else
    /* strategy 2: run only a fixed number of pointers at the same time.
     * We do more reductions, but control how fewer we do. 
     */
#define SIMULTANEOUS    8
    size_t k;
    for(k = 0 ; k + SIMULTANEOUS - 1 < n ; k += SIMULTANEOUS) {
        for(size_t j = 0 ; j < transform_datasize(p) ; j++) {
            Kelt_ur s;
            Kelt_ur_set_zero(s);
            for(size_t r = 0 ; r < SIMULTANEOUS ; r++) {
                Kelt_ur e;
                Kmul_ur(e, x1[r][j], x2[r][j]);
                Kelt_ur_add(s, s, e);
            }
            Kreduce(er, s);
            Kadd(y[j], y[j], er);
        }
        x1 += SIMULTANEOUS;
        x2 += SIMULTANEOUS;
    }
    for( ; k < n ; k ++) {
        for(size_t j = 0 ; j < transform_datasize(p) ; j++) {
            Kelt_ur s;
            Kelt_ur_set_zero(s);
            {
                size_t r = 0;
                Kelt_ur e;
                Kmul_ur(e, x1[r][j], x2[r][j]);
                Kelt_ur_add(s, s, e);
            }
            Kreduce(er, s);
            Kadd(y[j], y[j], er);
        }
        x1 ++;
        x2 ++;
    }
#endif
}

void gf2x_cantor_fft_addcompose(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_ptr y, gf2x_cantor_fft_srcptr x1, gf2x_cantor_fft_srcptr x2)
{
    gf2x_cantor_fft_addcompose_n(p, y, &x1, &x2, 1);
}

void gf2x_cantor_fft_add(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_ptr y, gf2x_cantor_fft_srcptr x1, gf2x_cantor_fft_srcptr x2)
{
    for (size_t j = 0; j < transform_datasize(p) ; j++) {
        Kadd(y[j], x1[j], x2[j]);
    }
}

void gf2x_cantor_fft_cpy(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_ptr y, gf2x_cantor_fft_srcptr x)
{
    memcpy(y, x, transform_datasize(p)*sizeof(Kelt));
}

size_t gf2x_cantor_fft_size(const gf2x_cantor_fft_info_t p)
{
    // the _size() is a number of Kelt of the result.
    // n normally, except when not truncating...
    return transform_datasize(p);
}

/* nH is a number of coefficients */
/* this destroys the input */
void gf2x_cantor_fft_ift(
        const gf2x_cantor_fft_info_t p,
        unsigned long * H,
        size_t nH,
        gf2x_cantor_fft_ptr h)
{
    size_t Hl = (nH + GF2X_WORDSIZE - 1) / GF2X_WORDSIZE;

#ifdef WITHOUT_CANTOR_TRUNCATION
    interpolateK(h, p->k);
#else
    // fill in with zeros to facilitate interpolation
    memset(h + p->n, 0, ((1UL << p->k) - p->n) * sizeof(Kelt));
    if (p->n & (p->n - 1)) {
        /* n is not a power of 2 */
        interpolateK_trunc(h, p->k, p->n);
    } else {
        interpolateK(h, p->k);
    }
#endif
    recomposeK(H, h, Hl, p->k);
}


gf2x_cantor_fft_ptr gf2x_cantor_fft_alloc(const gf2x_cantor_fft_info_t p, size_t n)
{
    return (Kelt *) malloc((n+(n << p->k)) * sizeof(Kelt));
}
void gf2x_cantor_fft_free(
        const gf2x_cantor_fft_info_t p GF2X_MAYBE_UNUSED,
        gf2x_cantor_fft_ptr x,
        size_t n GF2X_MAYBE_UNUSED)
{
    free(x);
}
gf2x_cantor_fft_srcptr gf2x_cantor_fft_get_const(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_srcptr x, size_t k)
{
    return x + (k << p->k) + k;
}
gf2x_cantor_fft_ptr gf2x_cantor_fft_get(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_ptr x, size_t k)
{
    return x + (k << p->k) + k;
}
void gf2x_cantor_fft_zero(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_ptr x, size_t n)
{
    memset(x, 0, ((n << p->k) + n) * sizeof(Kelt));
}
void gf2x_cantor_fft_init_similar(gf2x_cantor_fft_info_ptr o, size_t bits_a, size_t bits_b, gf2x_cantor_fft_info_srcptr other GF2X_MAYBE_UNUSED)
{
    gf2x_cantor_fft_init(o, bits_a, bits_b);
}
int gf2x_cantor_fft_compatible(gf2x_cantor_fft_info_srcptr o1, gf2x_cantor_fft_info_srcptr o2)
{
    return o1->k == o2->k;
}

/* vim: set sw=4 sta et: */
