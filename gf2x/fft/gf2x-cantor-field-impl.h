/* This file is part of the gf2x library.

   Copyright 2022, 2023
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
#ifndef GF2X_CANTOR_FIELD_IMPL_H_
#define GF2X_CANTOR_FIELD_IMPL_H_

#include <string.h>

#include "fft/gf2x-cantor-fft.h"
#include "fft/gf2x-fft-impl-utils.h"
#include "gf2x.h"
#include "gf2x/gf2x-impl.h"
#include "gf2x/gf2x-small.h"

#ifndef GF2X_CANTOR_BASE_FIELD_SIZE
#error "GF2X_CANTOR_BASE_FIELD_SIZE must be defined"
#endif

#ifndef GF2X_WORDSIZE
#error "GF2X_WORDSIZE must be defined"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Our goal is to define a few functions, inherited from mpfq in
 * their interface, to accomodate our finite field implementation
 * needs. These are pretty simple.
 *
 * We have two influential macros
 *   GF2X_CANTOR_BASE_FIELD_SIZE  (128 or 64)
 *   GF2X_WORDSIZE  (64 or 32)
 *
 * The interfaces that we need are
 *  elt
 *  dst_elt
 *  src_elt

 *  set
 *  swap

 *  set_zero
 *  is_zero
 *  add
 *  mul
 *  elt_ur
 *  mul_ur
 *  elt_ur_set_zero
 *  reduce
 *  elt_ur_add
 *
 */

typedef gf2x_cantor_fft_elt Kelt;
typedef gf2x_cantor_fft_dst_elt Kdst_elt;
typedef gf2x_cantor_fft_src_elt Ksrc_elt;

typedef unsigned long Kelt_ur[2 * GF2X_CANTOR_BASE_FIELD_SIZE / GF2X_WORDSIZE];
typedef unsigned long* Kdst_elt_ur;
typedef const unsigned long* Ksrc_elt_ur;

static inline void Kset(Kdst_elt r, Ksrc_elt s)
{
    if (r != s)
        memcpy(r, s, sizeof(Kelt));
}

static inline void Kset_zero(Kdst_elt r) { memset(r, 0, sizeof(Kelt)); }

static inline void Kelt_ur_set_zero(Kdst_elt_ur r)
{
    memset(r, 0, sizeof(Kelt_ur));
}

static inline int Kis_zero(Ksrc_elt r)
{
    for (unsigned int i = 0; i < sizeof(Kelt) / sizeof(r[0]); i++)
        if (r[i])
            return 0;
    return 1;
}

static inline void Kswap(Kdst_elt a, Kdst_elt b)
{
    Kelt x;
    Kset(x, a);
    Kset(a, b);
    Kset(b, x);
}

static inline void Kadd(Kdst_elt r, Ksrc_elt s1, Ksrc_elt s2)
{
    for (unsigned int i = 0; i < sizeof(Kelt) / sizeof(r[0]); i++)
        r[i] = s1[i] ^ s2[i];
}

static inline void Kelt_ur_add(Kdst_elt_ur r,
                               Ksrc_elt_ur s1,
                               Ksrc_elt_ur s2)
{
    for (unsigned int i = 0; i < sizeof(Kelt_ur) / sizeof(r[0]); i++)
        r[i] = s1[i] ^ s2[i];
}

static inline void Kmul_ur(Kdst_elt_ur r, Ksrc_elt s1, Ksrc_elt s2);
static inline void Kreduce(Kdst_elt, Kdst_elt_ur r);

static inline void Kmul(Kdst_elt r, Ksrc_elt s1, Ksrc_elt s2)
{
    Kelt_ur t;
    Kmul_ur(t, s1, s2);
    Kreduce(r, t);
}

static inline void Kmul_ur(Kdst_elt_ur t, Ksrc_elt s1, Ksrc_elt s2)
{
#if GF2X_CANTOR_BASE_FIELD_SIZE == 128 && GF2X_WORDSIZE == 64
    gf2x_mul2(t, s1, s2);
#elif GF2X_CANTOR_BASE_FIELD_SIZE == 128 && GF2X_WORDSIZE == 32
    gf2x_mul4(t, s1, s2);
#elif GF2X_CANTOR_BASE_FIELD_SIZE == 64 && GF2X_WORDSIZE == 64
    gf2x_mul1(t, *s1, *s2);
#elif GF2X_CANTOR_BASE_FIELD_SIZE == 64 && GF2X_WORDSIZE == 32
    gf2x_mul2(t, s1, s2);
#else
#error "Unsupported combination"
#endif
}

/* *Mpfq::gf2n::squaring::code_for_sqr_ur */

static inline void Kreduce(Kdst_elt r, Kdst_elt_ur t)
{
#if GF2X_CANTOR_BASE_FIELD_SIZE == 128 && GF2X_WORDSIZE == 64
    unsigned long z;
    z = t[3];
    /* t[2],t[1] += z * X^128 = z * (X^7 + X^2 + X + 1) */
    t[1] ^= (z << 7) ^ (z << 2) ^ (z << 1) ^ z;
    t[2] ^= (z >> 57) ^ (z >> 62) ^ (z >> 63);
    /* t[1],t[0] += z * X^128 = z * (X^7 + X^2 + X + 1) */
    z = t[2];
    r[0] = t[0] ^ (z << 7) ^ (z << 2) ^ (z << 1) ^ z;
    r[1] = t[1] ^ (z >> 57) ^ (z >> 62) ^ (z >> 63);
#elif GF2X_CANTOR_BASE_FIELD_SIZE == 64 && GF2X_WORDSIZE == 64
    unsigned long s[2];
    /* 63 excess bits */
    {
        unsigned long z;
        z = t[0];
        s[0] = z;
    }
    memset(s + 1, 0, 1 * sizeof(unsigned long));
    {
        unsigned long z;
        z = t[1];
        s[0] ^= z << 4;
        s[0] ^= z << 3;
        s[0] ^= z << 1;
        s[0] ^= z;
        z >>= 60;
        s[1] ^= z;
        z >>= 1;
        s[1] ^= z;
    }
    /* 3 excess bits */
    {
        unsigned long z;
        z = s[0];
        r[0] = z;
    }
    {
        unsigned long z;
        z = s[1];
        r[0] ^= z << 4;
        r[0] ^= z << 3;
        r[0] ^= z << 1;
        r[0] ^= z;
    }
#elif GF2X_CANTOR_BASE_FIELD_SIZE == 128 && GF2X_WORDSIZE == 32
    unsigned long s[5];
    /* 127 excess bits */
    {
        unsigned long z;
        z = t[0];
        s[0] = z;
        z = t[1];
        s[1] = z;
        z = t[2];
        s[2] = z;
        z = t[3];
        s[3] = z;
    }
    memset(s + 4, 0, 1 * sizeof(unsigned long));
    {
        unsigned long z;
        z = t[4];
        s[0] ^= z << 7;
        s[0] ^= z << 2;
        s[0] ^= z << 1;
        s[0] ^= z;
        z >>= 25;
        z ^= t[5] << 7;
        s[1] ^= z;
        z >>= 5;
        z ^= t[5] >> 25 << 27;
        s[1] ^= z;
        z >>= 1;
        z ^= t[5] >> 30 << 31;
        s[1] ^= z;
        z >>= 1;
        z ^= (t[5] & ~0x7fffffffUL);
        s[1] ^= z;
        z >>= 25;
        z ^= t[6] << 7;
        s[2] ^= z;
        z >>= 5;
        z ^= t[6] >> 25 << 27;
        s[2] ^= z;
        z >>= 1;
        z ^= t[6] >> 30 << 31;
        s[2] ^= z;
        z >>= 1;
        z ^= (t[6] & ~0x7fffffffUL);
        s[2] ^= z;
        z >>= 25;
        z ^= t[7] << 7;
        s[3] ^= z;
        z >>= 5;
        z ^= t[7] >> 25 << 27;
        s[3] ^= z;
        z >>= 1;
        z ^= t[7] >> 30 << 31;
        s[3] ^= z;
        z >>= 1;
        s[3] ^= z;
        z >>= 25;
        s[4] ^= z;
        z >>= 5;
        s[4] ^= z;
    }
    /* 6 excess bits */
    {
        unsigned long z;
        z = s[0];
        r[0] = z;
        z = s[1];
        r[1] = z;
        z = s[2];
        r[2] = z;
        z = s[3];
        r[3] = z;
    }
    {
        unsigned long z;
        z = s[4];
        r[0] ^= z << 7;
        r[0] ^= z << 2;
        r[0] ^= z << 1;
        r[0] ^= z;
    }
#elif GF2X_CANTOR_BASE_FIELD_SIZE == 64 && GF2X_WORDSIZE == 32
    unsigned long s[3];
    /* 63 excess bits */
    {
        unsigned long z;
        z = t[0];
        s[0] = z;
        z = t[1];
        s[1] = z;
    }
    memset(s + 2, 0, 1 * sizeof(unsigned long));
    {
        unsigned long z;
        z = t[2];
        s[0] ^= z << 4;
        s[0] ^= z << 3;
        s[0] ^= z << 1;
        s[0] ^= z;
        z >>= 28;
        z ^= t[3] << 4;
        s[1] ^= z;
        z >>= 1;
        z ^= t[3] >> 28 << 31;
        s[1] ^= z;
        z >>= 2;
        z ^= t[3] >> 29 << 30;
        s[1] ^= z;
        z >>= 1;
        s[1] ^= z;
        z >>= 28;
        s[2] ^= z;
        z >>= 1;
        s[2] ^= z;
    }
    /* 3 excess bits */
    {
        unsigned long z;
        z = s[0];
        r[0] = z;
        z = s[1];
        r[1] = z;
    }
    {
        unsigned long z;
        z = s[2];
        r[0] ^= z << 4;
        r[0] ^= z << 3;
        r[0] ^= z << 1;
        r[0] ^= z;
    }
#endif
}

#ifdef __cplusplus
}
#endif

#endif /* GF2X_CANTOR_FIELD_IMPL_H_ */
