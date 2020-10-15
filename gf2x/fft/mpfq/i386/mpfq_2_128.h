#ifndef MPFQ_2_128_H_
#define MPFQ_2_128_H_

/* MPFQ generated file -- do not edit */

#include "gf2x.h"
#include "gf2x/gf2x-small.h"

#include "mpfq.h"
#include "mpfq_gf2n_common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>

#include "assert.h"
#ifdef	MPFQ_LAST_GENERATED_TAG
#undef	MPFQ_LAST_GENERATED_TAG
#endif
#define MPFQ_LAST_GENERATED_TAG      2_128

/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::gf2n::field */
/* Automatically generated code for GF(2^128) */
/* Definition polynomial P = X^128 + X^7 + X^2 + X + 1 */
/* Active handler: Mpfq::gf2n::trivialities */
/* Active handler: Mpfq::gf2n::io */
/* Active handler: Mpfq::gf2n::linearops */
/* Active handler: Mpfq::gf2n::inversion */
/* Active handler: Mpfq::gf2n::reduction */
/* Active handler: Mpfq::gf2n::mul */
/* Active handler: Mpfq::defaults::poly */
/* Options used:{
   coeffs=[ 128, 7, 2, 1, 0, ],
   helper=/tmp/mpfq-cado/gf2n/helper/helper,
   n=128,
   no_gmp=1,
   output_path=i386,
   slice=4,
   table=/tmp/mpfq-cado/gf2x/wizard.table,
   tag=2_128,
   w=32,
   } */

typedef mpfq_2_field mpfq_2_128_field;
typedef mpfq_2_dst_field mpfq_2_128_dst_field;
typedef mpfq_2_src_field mpfq_2_128_src_field;

typedef unsigned long mpfq_2_128_elt[4];
typedef unsigned long * mpfq_2_128_dst_elt;
typedef const unsigned long * mpfq_2_128_src_elt;

typedef unsigned long mpfq_2_128_elt_ur[8];
typedef unsigned long * mpfq_2_128_dst_elt_ur;
typedef const unsigned long * mpfq_2_128_src_elt_ur;

typedef mpfq_2_128_elt * mpfq_2_128_vec;
typedef mpfq_2_128_elt * mpfq_2_128_dst_vec;
typedef mpfq_2_128_elt * mpfq_2_128_src_vec;

typedef struct {
  mpfq_2_128_vec c;
  unsigned int alloc;
  unsigned int size;
} mpfq_2_128_poly_struct;
typedef mpfq_2_128_poly_struct mpfq_2_128_poly [1];
typedef mpfq_2_128_poly_struct * mpfq_2_128_dst_poly;
typedef mpfq_2_128_poly_struct * mpfq_2_128_src_poly;
/* Extra types defined by implementation: */
typedef mpfq_2_128_elt_ur * mpfq_2_128_dst_vec_ur;
typedef mpfq_2_128_elt_ur * mpfq_2_128_vec_ur;
typedef mpfq_2_128_elt_ur * mpfq_2_128_src_vec_ur;

#ifdef  __cplusplus
extern "C" {
#endif

/* Elementary assignment functions */
static inline
void mpfq_2_128_set(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_set_zero(mpfq_2_128_dst_field, mpfq_2_128_dst_elt);

/* Comparison functions */
static inline
int mpfq_2_128_is_zero(mpfq_2_128_dst_field, mpfq_2_128_src_elt);

/* Arithmetic operations on elements */
static inline
void mpfq_2_128_add(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_mul(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_sqr(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt);

/* Operations involving unreduced elements */
static inline
void mpfq_2_128_elt_ur_set(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur, mpfq_2_128_src_elt_ur);
static inline
void mpfq_2_128_elt_ur_set_elt(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_elt_ur_set_zero(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur);
static inline
void mpfq_2_128_elt_ur_add(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur, mpfq_2_128_src_elt_ur, mpfq_2_128_src_elt_ur);
static inline
void mpfq_2_128_mul_ur(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur, mpfq_2_128_src_elt, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_sqr_ur(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_reduce(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_dst_elt_ur);
static inline
void mpfq_2_128_vec_set_zero(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, unsigned int);
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
/* *Mpfq::defaults::flatdata::code_for_set, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_set(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s)
{
    if (r != s) memcpy(r,s,sizeof(mpfq_2_128_elt));
}

/* *Mpfq::defaults::flatdata::code_for_set_zero, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_set_zero(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r)
{
    mpfq_2_128_vec_set_zero(K,(mpfq_2_128_dst_vec)r,1);
}

/* *Mpfq::defaults::flatdata::code_for_is_zero, Mpfq::gf2n::trivialities */
static inline
int mpfq_2_128_is_zero(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_elt r)
{
        unsigned int i;
        for(i = 0 ; i < sizeof(mpfq_2_128_elt)/sizeof(r[0]) ; i++) {
            if (r[i]) return 0;
        }
        return 1;
}

/* *Mpfq::gf2n::trivialities::code_for_add */
static inline
void mpfq_2_128_add(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s1, mpfq_2_128_src_elt s2)
{
    int i;
    for(i = 0 ; i < 4 ; i++)
        r[i] = s1[i] ^ s2[i];
}

/* *Mpfq::gf2n::trivialities::code_for_mul */
static inline
void mpfq_2_128_mul(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s1, mpfq_2_128_src_elt s2)
{
    mpfq_2_128_elt_ur t;
    mpfq_2_128_mul_ur(K, t, s1, s2);
    mpfq_2_128_reduce(K, r, t);
}

/* *Mpfq::gf2n::trivialities::code_for_sqr */
static inline
void mpfq_2_128_sqr(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s)
{
    mpfq_2_128_elt_ur t;
    mpfq_2_128_sqr_ur(K, t, s);
    mpfq_2_128_reduce(K, r, t);
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_elt_ur_set(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur r, mpfq_2_128_src_elt_ur s)
{
    if (r != s) memcpy(r,s,sizeof(mpfq_2_128_elt_ur));
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_elt, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_elt_ur_set_elt(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur r, mpfq_2_128_src_elt s)
{
    memset(r, 0, sizeof(mpfq_2_128_elt_ur)); memcpy(r,s,sizeof(mpfq_2_128_elt));
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_zero, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_elt_ur_set_zero(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur r)
{
    memset(r, 0, sizeof(mpfq_2_128_elt_ur));
}

/* *Mpfq::gf2n::trivialities::code_for_elt_ur_add */
static inline
void mpfq_2_128_elt_ur_add(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur r, mpfq_2_128_src_elt_ur s1, mpfq_2_128_src_elt_ur s2)
{
    int i;
    for(i = 0 ; i < 8 ; i++)
        r[i] = s1[i] ^ s2[i];
}

/* *Mpfq::gf2n::mul::code_for_mul_ur */
static inline
void mpfq_2_128_mul_ur(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur t, mpfq_2_128_src_elt s1, mpfq_2_128_src_elt s2)
{
    gf2x_mul4(t, s1, s2);
}

/* *Mpfq::gf2n::squaring::code_for_sqr_ur */
static inline
void mpfq_2_128_sqr_ur(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur t, mpfq_2_128_src_elt s)
{
    static const unsigned long g[16] = {
        0, 1, 4, 5, 16, 17, 20, 21,
        64, 65, 68, 69, 80, 81, 84, 85,
    };
    {
        unsigned long u;
        u = g[s[0]       & 15];
    t[0]  = u;
        u = g[s[0] >>  4 & 15];
    t[0] ^= u <<  8;
        u = g[s[0] >>  8 & 15];
    t[0] ^= u << 16;
        u = g[s[0] >> 12 & 15];
    t[0] ^= u << 24;
        u = g[s[0] >> 16 & 15];
    t[1]  = u;
        u = g[s[0] >> 20 & 15];
    t[1] ^= u <<  8;
        u = g[s[0] >> 24 & 15];
    t[1] ^= u << 16;
        u = g[s[0] >> 28 & 15];
    t[1] ^= u << 24;
        u = g[s[1]       & 15];
    t[2]  = u;
        u = g[s[1] >>  4 & 15];
    t[2] ^= u <<  8;
        u = g[s[1] >>  8 & 15];
    t[2] ^= u << 16;
        u = g[s[1] >> 12 & 15];
    t[2] ^= u << 24;
        u = g[s[1] >> 16 & 15];
    t[3]  = u;
        u = g[s[1] >> 20 & 15];
    t[3] ^= u <<  8;
        u = g[s[1] >> 24 & 15];
    t[3] ^= u << 16;
        u = g[s[1] >> 28 & 15];
    t[3] ^= u << 24;
        u = g[s[2]       & 15];
    t[4]  = u;
        u = g[s[2] >>  4 & 15];
    t[4] ^= u <<  8;
        u = g[s[2] >>  8 & 15];
    t[4] ^= u << 16;
        u = g[s[2] >> 12 & 15];
    t[4] ^= u << 24;
        u = g[s[2] >> 16 & 15];
    t[5]  = u;
        u = g[s[2] >> 20 & 15];
    t[5] ^= u <<  8;
        u = g[s[2] >> 24 & 15];
    t[5] ^= u << 16;
        u = g[s[2] >> 28 & 15];
    t[5] ^= u << 24;
        u = g[s[3]       & 15];
    t[6]  = u;
        u = g[s[3] >>  4 & 15];
    t[6] ^= u <<  8;
        u = g[s[3] >>  8 & 15];
    t[6] ^= u << 16;
        u = g[s[3] >> 12 & 15];
    t[6] ^= u << 24;
        u = g[s[3] >> 16 & 15];
    t[7]  = u;
        u = g[s[3] >> 20 & 15];
    t[7] ^= u <<  8;
        u = g[s[3] >> 24 & 15];
    t[7] ^= u << 16;
        u = g[s[3] >> 28 & 15];
    t[7] ^= u << 24;
    }
}

/* *Mpfq::gf2n::reduction::code_for_reduce */
static inline
void mpfq_2_128_reduce(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_dst_elt_ur t)
{
    {
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
            s[0]^= z <<  7;
            s[0]^= z <<  2;
            s[0]^= z <<  1;
            s[0]^= z;
            z >>= 25;
            z^= t[5] <<  7;
            s[1]^= z;
            z >>= 5;
            z^= t[5] >> 25 << 27;
            s[1]^= z;
            z >>= 1;
            z^= t[5] >> 30 << 31;
            s[1]^= z;
            z >>= 1;
            z^= (t[5] & ~0x7fffffffUL);
            s[1]^= z;
            z >>= 25;
            z^= t[6] <<  7;
            s[2]^= z;
            z >>= 5;
            z^= t[6] >> 25 << 27;
            s[2]^= z;
            z >>= 1;
            z^= t[6] >> 30 << 31;
            s[2]^= z;
            z >>= 1;
            z^= (t[6] & ~0x7fffffffUL);
            s[2]^= z;
            z >>= 25;
            z^= t[7] <<  7;
            s[3]^= z;
            z >>= 5;
            z^= t[7] >> 25 << 27;
            s[3]^= z;
            z >>= 1;
            z^= t[7] >> 30 << 31;
            s[3]^= z;
            z >>= 1;
            s[3]^= z;
            z >>= 25;
            s[4]^= z;
            z >>= 5;
            s[4]^= z;
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
            r[0]^= z <<  7;
            r[0]^= z <<  2;
            r[0]^= z <<  1;
            r[0]^= z;
        }
    }
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set_zero, Mpfq::defaults::flatdata, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_vec_set_zero(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec r, unsigned int n)
{
    memset(r, 0, n*sizeof(mpfq_2_128_elt));
}


#endif  /* MPFQ_2_128_H_ */

/* vim:set ft=cpp: */
