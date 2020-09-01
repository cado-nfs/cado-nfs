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
   output_path=x86_64,
   slice=4,
   table=/tmp/mpfq-cado/gf2x/wizard.table,
   tag=2_128,
   w=64,
   } */

typedef mpfq_2_field mpfq_2_128_field;
typedef mpfq_2_dst_field mpfq_2_128_dst_field;
typedef mpfq_2_src_field mpfq_2_128_src_field;

typedef unsigned long mpfq_2_128_elt[2];
typedef unsigned long * mpfq_2_128_dst_elt;
typedef const unsigned long * mpfq_2_128_src_elt;

typedef unsigned long mpfq_2_128_elt_ur[4];
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
typedef mpfq_2_128_elt_ur * mpfq_2_128_vec_ur;
typedef mpfq_2_128_elt_ur * mpfq_2_128_dst_vec_ur;
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
    for(i = 0 ; i < 2 ; i++)
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
    for(i = 0 ; i < 4 ; i++)
        r[i] = s1[i] ^ s2[i];
}

/* *Mpfq::gf2n::mul::code_for_mul_ur */
static inline
void mpfq_2_128_mul_ur(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur t, mpfq_2_128_src_elt s1, mpfq_2_128_src_elt s2)
{
    gf2x_mul2(t, s1, s2);
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
    t[0] ^= u << 32;
        u = g[s[0] >> 20 & 15];
    t[0] ^= u << 40;
        u = g[s[0] >> 24 & 15];
    t[0] ^= u << 48;
        u = g[s[0] >> 28 & 15];
    t[0] ^= u << 56;
        u = g[s[0] >> 32 & 15];
    t[1]  = u;
        u = g[s[0] >> 36 & 15];
    t[1] ^= u <<  8;
        u = g[s[0] >> 40 & 15];
    t[1] ^= u << 16;
        u = g[s[0] >> 44 & 15];
    t[1] ^= u << 24;
        u = g[s[0] >> 48 & 15];
    t[1] ^= u << 32;
        u = g[s[0] >> 52 & 15];
    t[1] ^= u << 40;
        u = g[s[0] >> 56 & 15];
    t[1] ^= u << 48;
        u = g[s[0] >> 60 & 15];
    t[1] ^= u << 56;
        u = g[s[1]       & 15];
    t[2]  = u;
        u = g[s[1] >>  4 & 15];
    t[2] ^= u <<  8;
        u = g[s[1] >>  8 & 15];
    t[2] ^= u << 16;
        u = g[s[1] >> 12 & 15];
    t[2] ^= u << 24;
        u = g[s[1] >> 16 & 15];
    t[2] ^= u << 32;
        u = g[s[1] >> 20 & 15];
    t[2] ^= u << 40;
        u = g[s[1] >> 24 & 15];
    t[2] ^= u << 48;
        u = g[s[1] >> 28 & 15];
    t[2] ^= u << 56;
        u = g[s[1] >> 32 & 15];
    t[3]  = u;
        u = g[s[1] >> 36 & 15];
    t[3] ^= u <<  8;
        u = g[s[1] >> 40 & 15];
    t[3] ^= u << 16;
        u = g[s[1] >> 44 & 15];
    t[3] ^= u << 24;
        u = g[s[1] >> 48 & 15];
    t[3] ^= u << 32;
        u = g[s[1] >> 52 & 15];
    t[3] ^= u << 40;
        u = g[s[1] >> 56 & 15];
    t[3] ^= u << 48;
        u = g[s[1] >> 60 & 15];
    t[3] ^= u << 56;
    }
}

/* *Mpfq::gf2n::reduction::code_for_reduce */
static inline
void mpfq_2_128_reduce(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_dst_elt_ur t)
{
    {
                unsigned long z;
                z = t[3];
                /* t[2],t[1] += z * X^128 = z * (X^7 + X^2 + X + 1) */
                t[1] ^= (z << 7) ^ (z << 2) ^ (z << 1) ^ z;
                t[2] ^= (z >> 57) ^ (z >> 62) ^ (z >> 63);
                /* t[1],t[0] += z * X^128 = z * (X^7 + X^2 + X + 1) */
                z = t[2];
                r[0] = t[0] ^ (z << 7) ^ (z << 2) ^ (z << 1) ^ z;
                r[1] = t[1] ^ (z >> 57) ^ (z >> 62) ^ (z >> 63);
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
