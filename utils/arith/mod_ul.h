/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. Due to inlining, this file must be included
   in the caller's source code with #include */

/* Naming convention: all function start with modul, for
   MODulus Unsigned Long, followed by underscore, functionality of function
  (add, mul, etc), and possibly underscore and specification of what argument
  types the function takes (_ul, etc). */

#ifndef CADO_MOD_UL_H
#define CADO_MOD_UL_H

#include "cado_config.h" // just because we're a header.
/**********************************************************************/
#include "macros.h"  // for MAYBE_UNUSED, ASSERT_EXPENSIVE, ASSERT, ASS...
#include "ularith.h" // for ularith_mul_ul_ul_2ul, ularith_div_2ul_ul_ul_r
#include <limits.h>  // for ULONG_BITS, ULONG_MAX
#include <stdint.h>  // for uint64_t, int64_t
#include <stdio.h>
#include <stdlib.h> // for size_t, llabs

/*********************************************************************/
/* Helper macros, see also ularith.h */

/* A macro for function renaming. All functions here start with modul_ */
#define MODUL_RENAME(x) modul_##x

#define MODUL_SIZE 1
#define MODUL_MAXBITS ULONG_BITS

typedef unsigned long residueul_t[MODUL_SIZE];
typedef unsigned long modintul_t[MODUL_SIZE];
typedef unsigned long modulusul_t[MODUL_SIZE];

/* ================= Functions that are part of the API ================= */

/* Some functions for integers of the same width as the modulus */

MAYBE_UNUSED static inline void modul_intinit(modintul_t r)
{
    r[0] = 0;
}

MAYBE_UNUSED static inline void modul_intclear(modintul_t r MAYBE_UNUSED)
{
}

MAYBE_UNUSED static inline void modul_intset(modintul_t r, modintul_t const s)
{
    r[0] = s[0];
}

MAYBE_UNUSED static inline void modul_intset_ul(modintul_t r,
                                                unsigned long const s)
{
    r[0] = s;
}

/* The two mod*_uls() functions import/export modint_t from/to an array of
   unsigned longs. For modul_intset_ul, the size of the array is passed as
   a parameter n. For mod_intget_uls(), the required array size can be
   determined via mod_intbits(); if the modint_t is zero, mod_intget_uls()
   writes 0 to the first output unsigned long. It returns the number of
   unsigned longs written. */
MAYBE_UNUSED static inline void
modul_intset_uls(modintul_t r, unsigned long const * s, size_t const n)
{
    ASSERT_ALWAYS(n <= MODUL_SIZE);
    if (n == 0)
        r[0] = 0;
    else
        r[0] = s[0];
}

MAYBE_UNUSED static inline unsigned long modul_intget_ul(modintul_t const s)
{
    return s[0];
}

MAYBE_UNUSED static inline size_t modul_intget_uls(unsigned long * r,
                                                   modintul_t const s)
{
    r[0] = s[0];
    return 1;
}

MAYBE_UNUSED static inline int modul_intequal(modintul_t const a,
                                              modintul_t const b)
{
    return (a[0] == b[0]);
}

MAYBE_UNUSED static inline int modul_intequal_ul(modintul_t const a,
                                                 unsigned long const b)
{
    return (a[0] == b);
}

MAYBE_UNUSED static inline int modul_intcmp(modintul_t const a,
                                            modintul_t const b)
{
    return (a[0] < b[0]) ? -1 : (a[0] == b[0]) ? 0 : 1;
}

MAYBE_UNUSED static inline int modul_intcmp_ul(modintul_t const a,
                                               unsigned long const b)
{
    return (a[0] < b) ? -1 : (a[0] == b) ? 0 : 1;
}

MAYBE_UNUSED static inline int modul_intcmp_uint64(modintul_t const a,
                                                   uint64_t const b)
{
    if (b > ULONG_MAX)
        return -1;
    return (a[0] < b) ? -1 : (a[0] == b) ? 0 : 1;
}

MAYBE_UNUSED static inline int modul_intfits_ul(modintul_t const a MAYBE_UNUSED)
{
    return 1;
}

MAYBE_UNUSED static inline void modul_intadd(modintul_t r, modintul_t const a,
                                             modintul_t const b)
{
    r[0] = a[0] + b[0];
}

MAYBE_UNUSED static inline void modul_intsub(modintul_t r, modintul_t const a,
                                             modintul_t const b)
{
    r[0] = a[0] - b[0];
}

MAYBE_UNUSED static inline void modul_intshr(modintul_t r, modintul_t const s,
                                             int const i)
{
    r[0] = s[0] >> i;
}

MAYBE_UNUSED static inline void modul_intshl(modintul_t r, modintul_t const s,
                                             int const i)
{
    r[0] = s[0] << i;
}

/* Returns the number of bits in a, that is, floor(log_2(n))+1.
   For n==0 returns 0. */
MAYBE_UNUSED static inline size_t modul_intbits(modintul_t const a)
{
    if (a[0] == 0)
        return 0;
    return ULONG_BITS - ularith_clz(a[0]);
}

/* r = n/d. We require d|n */
MAYBE_UNUSED static inline void
modul_intdivexact(modintul_t r, modintul_t const n, modintul_t const d)
{
    r[0] = n[0] / d[0];
}

/* r = n%d */
MAYBE_UNUSED static inline void modul_intmod(modintul_t r, modintul_t const n,
                                             modintul_t const d)
{
    r[0] = n[0] % d[0];
}

/* Functions for the modulus */

MAYBE_UNUSED static inline void modul_initmod_ul(modulusul_t m,
                                                 unsigned long const s)
{
    m[0] = s;
}

MAYBE_UNUSED static inline void modul_initmod_int(modulusul_t m,
                                                  modintul_t const s)
{
    m[0] = s[0];
}

MAYBE_UNUSED static inline unsigned long modul_getmod_ul(modulusul_t const m)
{
    return m[0];
}

MAYBE_UNUSED static inline void modul_getmod_int(modintul_t r,
                                                 modulusul_t const m)
{
    r[0] = m[0];
}

MAYBE_UNUSED static inline void modul_clearmod(modulusul_t m MAYBE_UNUSED)
{
}

/* Functions for residues */
static inline void modul_neg(residueul_t, residueul_t const, modulusul_t const);

/* Initialises a residue_t type and sets it to zero */
MAYBE_UNUSED static inline void modul_init(residueul_t r,
                                           modulusul_t const m MAYBE_UNUSED)
{
    r[0] = 0UL;
}

/* Initialises a residue_t type, but does not set it to zero. For fixed length
   residue_t types, that leaves nothing to do at all. */
MAYBE_UNUSED static inline void
modul_init_noset0(residueul_t r MAYBE_UNUSED, modulusul_t const m MAYBE_UNUSED)
{
}

MAYBE_UNUSED static inline void modul_clear(residueul_t r MAYBE_UNUSED,
                                            modulusul_t const m MAYBE_UNUSED)
{
}

MAYBE_UNUSED static inline void modul_set(residueul_t r, residueul_t const s,
                                          modulusul_t const m MAYBE_UNUSED)
{
    ASSERT_EXPENSIVE(s[0] < m[0]);
    r[0] = s[0];
}

MAYBE_UNUSED static inline void
modul_set_ul(residueul_t r, unsigned long const s, modulusul_t const m)
{
    r[0] = s % m[0];
}

/* Sets the residue_t to the class represented by the integer s. Assumes that
   s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED static inline void modul_set_ul_reduced(residueul_t r,
                                                     unsigned long const s,
                                                     modulusul_t const m
                                                     MAYBE_UNUSED)
{
    ASSERT(s < m[0]);
    r[0] = s;
}

MAYBE_UNUSED static inline void modul_set_int(residueul_t r, modintul_t const s,
                                              modulusul_t const m)
{
    r[0] = s[0] % m[0];
}

MAYBE_UNUSED static inline void modul_set_int_reduced(residueul_t r,
                                                      modintul_t const s,
                                                      modulusul_t const m
                                                      MAYBE_UNUSED)
{
    ASSERT(s[0] < m[0]);
    r[0] = s[0];
}

MAYBE_UNUSED static inline void
modul_set_uint64(residueul_t r, uint64_t const s, modulusul_t const m)
{
    r[0] = s % m[0];
}

MAYBE_UNUSED static inline void modul_set_int64(residueul_t r, int64_t const s,
                                                modulusul_t const m)
{
    r[0] = llabs(s) % m[0];
    if (s < 0)
        modul_neg(r, r, m);
}

/* These two are so trivial that we don't really require m in the
 * interface. For 1 we might, as the internal representation might
 * not use "1" for 1 (e.g. when using Montgomery's REDC.)
 * For interface homogeneity we make even modul_set0 take the m parameter.
 */
MAYBE_UNUSED static inline void modul_set0(residueul_t r,
                                           modulusul_t const m MAYBE_UNUSED)
{
    r[0] = 0UL;
}

MAYBE_UNUSED static inline void modul_set1(residueul_t r,
                                           modulusul_t const m MAYBE_UNUSED)
{
    r[0] = m[0] != 1UL;
}

/* Exchanges the values of the two arguments */

MAYBE_UNUSED static inline void modul_swap(residueul_t a, residueul_t b,
                                           modulusul_t const m MAYBE_UNUSED)
{
    unsigned long t;
    ASSERT_EXPENSIVE(a[0] < m[0] && b[0] < m[0]);
    t = a[0];
    a[0] = b[0];
    b[0] = t;
}

MAYBE_UNUSED static inline unsigned long
modul_get_ul(residueul_t const s, modulusul_t const m MAYBE_UNUSED)
{
    ASSERT_EXPENSIVE(s[0] < m[0]);
    return s[0];
}

MAYBE_UNUSED static inline void modul_get_int(modintul_t r, residueul_t const s,
                                              modulusul_t const m MAYBE_UNUSED)
{
    ASSERT_EXPENSIVE(s[0] < m[0]);
    r[0] = s[0];
}

MAYBE_UNUSED static inline int modul_equal(residueul_t const a,
                                           residueul_t const b,
                                           modulusul_t const m MAYBE_UNUSED)
{
    ASSERT_EXPENSIVE(a[0] < m[0] && b[0] < m[0]);
    return (a[0] == b[0]);
}

MAYBE_UNUSED static inline int modul_is0(residueul_t const a,
                                         modulusul_t const m MAYBE_UNUSED)
{
    ASSERT_EXPENSIVE(a[0] < m[0]);
    return (a[0] == 0UL);
}

MAYBE_UNUSED static inline int modul_is1(residueul_t const a,
                                         modulusul_t const m MAYBE_UNUSED)
{
    ASSERT_EXPENSIVE(a[0] < m[0]);
    return (a[0] == 1UL);
}

MAYBE_UNUSED static inline void modul_add(residueul_t r,
        residueul_t const a,
        residueul_t const b,
        modulusul_t const m)
{
#ifdef MODTRACE
    printf("modul_add: a = %lu, b = %lu", a[0], b[0]);
#endif

    ularith_addmod_ul_ul(r, a[0], b[0], m[0]);

#ifdef MODTRACE
    printf(", r = %lu\n", r[0]);
#endif
}

MAYBE_UNUSED static inline void modul_add1(residueul_t r, residueul_t const a,
                                           modulusul_t const m)
{
    ASSERT_EXPENSIVE(a[0] < m[0]);
    r[0] = a[0] + 1;
    if (r[0] == m[0])
        r[0] = 0;
}

MAYBE_UNUSED static inline void modul_add_ul(residueul_t r, residueul_t const a,
                                             unsigned long const b,
                                             modulusul_t const m)
{
    ularith_addmod_ul_ul(r, a[0], b % m[0], m[0]);
}

MAYBE_UNUSED static inline void modul_sub(residueul_t r, residueul_t const a,
                                          residueul_t const b,
                                          modulusul_t const m)
{
#ifdef MODTRACE
    printf("submod_ul: a = %lu, b = %lu", a[0], b[0]);
#endif

    ularith_submod_ul_ul(r, a[0], b[0], m[0]);

#ifdef MODTRACE
    printf(", r = %lu\n", r[0]);
#endif
}

MAYBE_UNUSED static inline void modul_sub_ul(residueul_t r, residueul_t const a,
                                             unsigned long const b,
                                             modulusul_t const m)
{
    ularith_submod_ul_ul(r, a[0], b % m[0], m[0]);
}

MAYBE_UNUSED static inline void modul_neg(residueul_t r, residueul_t const a,
                                          modulusul_t const m)
{
    ASSERT_EXPENSIVE(a[0] < m[0]);
    if (a[0] == 0UL)
        r[0] = a[0];
    else
        r[0] = m[0] - a[0];
}

MAYBE_UNUSED static inline void modul_mul(residueul_t r, residueul_t const a,
                                          residueul_t const b,
                                          modulusul_t const m)
{
    unsigned long t1, t2;
#if defined(MODTRACE)
    printf("mulmod_ul: a = %lu, b = %lu", a[0], b[0]);
#endif

    ASSERT_EXPENSIVE(a[0] < m[0] && b[0] < m[0]);
    ularith_mul_ul_ul_2ul(&t1, &t2, a[0], b[0]);
    ularith_div_2ul_ul_ul_r(r, t1, t2, m[0]);

#if defined(MODTRACE)
    printf(", r = %lu\n", r);
#endif
}

MAYBE_UNUSED static inline void modul_sqr(residueul_t r, residueul_t const a,
                                          modulusul_t const m)
{
    unsigned long t1, t2;

    ASSERT_EXPENSIVE(a[0] < m[0]);
    ularith_mul_ul_ul_2ul(&t1, &t2, a[0], a[0]);
    ularith_div_2ul_ul_ul_r(r, t1, t2, m[0]);
}

/* Computes (a * 2^wordsize) % m */
MAYBE_UNUSED static inline void
modul_tomontgomery(residueul_t r, residueul_t const a, modulusul_t const m)
{
    ASSERT_EXPENSIVE(a[0] < m[0]);
    ularith_div_2ul_ul_ul_r(r, 0UL, a[0], m[0]);
}

/* Computes (a / 2^wordsize) % m */
MAYBE_UNUSED static inline void modul_frommontgomery(residueul_t r,
                                                     residueul_t const a,
                                                     unsigned long const invm,
                                                     modulusul_t const m)
{
    unsigned long t_low, t_high;
    t_low = a[0] * invm;
    ularith_mul_ul_ul_2ul(&t_low, &t_high, t_low, m[0]);
    r[0] = t_high + (a[0] != 0UL ? 1UL : 0UL);
}

/* Computes (a / 2^wordsize) % m, but result can be r = m.
   Input a must not be equal 0 */
MAYBE_UNUSED static inline void modul_redcsemi_ul_not0(residueul_t r,
                                                       unsigned long const a,
                                                       unsigned long const invm,
                                                       modulusul_t const m)
{
    unsigned long t_low, t_high;

    ASSERT(a != 0);

    t_low = a * invm; /* t_low <= 2^w-1 */
    ularith_mul_ul_ul_2ul(&t_low, &t_high, t_low, m[0]);
    /* t_high:t_low <= (2^w-1) * m */
    r[0] = t_high + 1UL;
    /* (t_high+1):t_low <= 2^w + (2^w-1) * m  <= 2^w + 2^w*m - m
                      <= 2^w * (m + 1) - m */
    /* r <= floor ((2^w * (m + 1) - m) / 2^w) <= floor((m + 1) - m/2^w)
         <= m */
}

/* Computes ((a + b) / 2^wordsize) % m. a <= m is permissible */
MAYBE_UNUSED static inline void
modul_addredc_ul(residueul_t r, residueul_t const a, unsigned long const b,
                 unsigned long const invm, modulusul_t const m)
{
    unsigned long s_low, s_high, t_low, t_high;

    ASSERT_EXPENSIVE(a[0] <= m[0]);
    s_low = b;
    s_high = 0UL;
    ularith_add_ul_2ul(&s_low, &s_high, a[0]);

    t_low = s_low * invm;
    ularith_mul_ul_ul_2ul(&t_low, &t_high, t_low, m[0]);
    ASSERT_EXPENSIVE(s_low + t_low == 0UL);
    r[0] = t_high + s_high + (s_low != 0UL ? 1UL : 0UL);

    /* r = ((a+b) + (((a+b)%2^w * invm) % 2^w) * m) / 2^w  Use a<=m-1, b<=2^w-1
       r <= (m + 2^w - 1 + (2^w - 1) * m) / 2^w
          = (m - 1 + 2^w + m*2^w - m) / 2^w
          = (- 1 + 2^w + m2^w) / 2^w
          = m + 1 - 1/2^w
       r <= m, since r is an integer
    */
    if (r[0] == m[0])
        r[0] = 0UL;
}

/* Computes ((a + b) / 2^wordsize) % m, but result can be == m.
   a <= m is permissible */
MAYBE_UNUSED static inline void
modul_addredcsemi_ul(residueul_t r, residueul_t const a, unsigned long const b,
                     unsigned long const invm, modulusul_t const m)
{
    unsigned long s_low, s_high, t_low;

    ASSERT_EXPENSIVE(a[0] <= m[0]);
    s_low = b;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    {
        unsigned char sb;
        __asm__ __VOLATILE("addq %2, %0\n\t" /* cy * 2^w + s_low = a + b */
                           "setne %1\n\t"    /* if (s_low != 0) sb = 1 */
                           "adcb $0, %1\n"   /* sb += cy */
                           : "+&r"(s_low), "=qm"(sb)
                           : "rm"(a[0])
                           : "cc");
        s_high = sb;
    }
#elif defined(__i386__) && defined(__GNUC__)
    {
        unsigned char sb;
        __asm__ __VOLATILE("addl %2, %0\n\t"
                           "setne %1\n\t"
                           "adcb $0, %1\n"
                           : "+&r"(s_low), "=qm"(sb)
                           : "rm"(a[0])
                           : "cc");
        s_high = sb;
    }
#else
    s_high = 0UL;
    ularith_add_ul_2ul(&s_low, &s_high, a[0]);
    s_high += (s_low != 0UL ? 1UL : 0UL);
#endif

    t_low = s_low * invm;
    ularith_mul_ul_ul_2ul(&t_low, r, t_low, m[0]);
    ASSERT_EXPENSIVE(s_low + t_low == 0UL);
    r[0] += s_high;

    /* r = ((a+b) + (((a+b)%2^w * invm) % 2^w) * m) / 2^w
       r <= ((a+b) + (2^w - 1) * m) / 2^w
       r <= (m + 2^w-1 + m*2^w - m) / 2^w
       r <= (2^w -1 + p2^w) / 2^w
       r <= p + 1 - 1/2^w
       r <= p
    */
}

MAYBE_UNUSED static inline void
modul_mulredc(residueul_t r, residueul_t const a, residueul_t const b,
              unsigned long const invm, modulusul_t const m)
{
    unsigned long p_low, p_high, t_low, t_high;

    ASSERT_EXPENSIVE(m[0] % 2 != 0);
    ASSERT_EXPENSIVE(a[0] < m[0] && b[0] < m[0]);
#if defined(MODTRACE)
    printf("(%lu * %lu / 2^%ld) %% %lu", a[0], b[0], 8 * sizeof(unsigned long),
           m[0]);
#endif

    ularith_mul_ul_ul_2ul(&p_low, &p_high, a[0], b[0]);
    t_low = p_low * invm;
    ularith_mul_ul_ul_2ul(&t_low, &t_high, t_low, m[0]);
    /* Let w = 2^wordsize. We know (p_high * w + p_low) + (t_high * w + t_low)
       == 0 (mod w) so either p_low == t_low == 0, or p_low !=0 and t_low != 0.
       In the former case we want p_high + t_high + 1, in the latter
       p_high + t_high */
    /* Since a <= p-1 and b <= p-1, and p <= w-1, a*b <= w^2 - 4*w + 4, so
       adding 1 to p_high is safe */
#if 0
  /* S_lower? */
  ularith_add_ul_2ul (&p_low, &p_high, t_low);
#else
    p_high += (p_low != 0UL ? 1UL : 0UL);
#endif

    modul_add(r, &p_high, &t_high, m);

#if defined(MODTRACE)
    printf(" == %lu /* PARI */ \n", r[0]);
#endif
}

/* FIXME: check for overflow if b > m */
MAYBE_UNUSED static inline void
modul_mulredc_ul(residueul_t r, residueul_t const a, unsigned long const b,
                 unsigned long const invm, modulusul_t const m)
{
    unsigned long p_low, p_high, t_low, t_high;
    ASSERT_EXPENSIVE(m[0] % 2 != 0);
#if defined(MODTRACE)
    printf("(%lu * %lu / 2^%ld) %% %lu", a[0], b, 8 * sizeof(unsigned long),
           m[0]);
#endif

    ularith_mul_ul_ul_2ul(&p_low, &p_high, a[0], b);
    t_low = p_low * invm;
    ularith_mul_ul_ul_2ul(&t_low, &t_high, t_low, m[0]);
    p_high += (p_low != 0UL ? 1UL : 0UL);
    r[0] = (p_high >= m[0] - t_high) ? (p_high - (m[0] - t_high))
                                     : (p_high + t_high);

#if defined(MODTRACE)
    printf(" == %lu /* PARI */ \n", r[0]);
#endif
}

/* Computes (a * b + c)/ 2^wordsize % m. Requires that
   a * b + c < 2^wordsize * m */

MAYBE_UNUSED static inline void
modul_muladdredc_ul(residueul_t r, residueul_t const a, unsigned long const b,
                    unsigned long const c, unsigned long const invm,
                    modulusul_t const m)
{
    unsigned long p_low, p_high, t_low, t_high;
    ASSERT_EXPENSIVE(m[0] % 2 != 0);
#if defined(MODTRACE)
    printf("(%lu * %lu / 2^%ld) %% %lu", a[0], b, 8 * sizeof(unsigned long),
           m[0]);
#endif

    ularith_mul_ul_ul_2ul(&p_low, &p_high, a[0], b);
    ularith_add_ul_2ul(&p_low, &p_high, c);
    t_low = p_low * invm;
    ularith_mul_ul_ul_2ul(&t_low, &t_high, t_low, m[0]);
    p_high += (p_low != 0UL ? 1UL : 0UL);
    r[0] = (p_high >= m[0] - t_high) ? (p_high - (m[0] - t_high))
                                     : (p_high + t_high);

#if defined(MODTRACE)
    printf(" == %lu /* PARI */ \n", r[0]);
#endif
}

MAYBE_UNUSED static inline void modul_div2(residueul_t r, residueul_t const a,
                                           modulusul_t const m)
{
    r[0] = ularith_div2mod(a[0], m[0]);
}

MAYBE_UNUSED static inline int modul_next(residueul_t r, modulusul_t const m)
{
    return (++r[0] == m[0]);
}

MAYBE_UNUSED static inline int modul_finished(residueul_t const r,
                                              modulusul_t const m)
{
    return (r[0] == m[0]);
}

#ifdef __cplusplus
extern "C" {
#endif

/* prototypes of non-inline functions */
int modul_div3(residueul_t, residueul_t const, modulusul_t const);
int modul_div5(residueul_t, residueul_t const, modulusul_t const);
int modul_div7(residueul_t, residueul_t const, modulusul_t const);
int modul_div11(residueul_t, residueul_t const, modulusul_t const);
int modul_div13(residueul_t, residueul_t const, modulusul_t const);
void modul_gcd(modintul_t, residueul_t const, modulusul_t const);
void modul_pow_ul(residueul_t, residueul_t const, unsigned long,
                  modulusul_t const);
void modul_2pow_ul(residueul_t, unsigned long, modulusul_t const);
void modul_pow_mp(residueul_t, residueul_t const, unsigned long const *,
                  int, modulusul_t const);
void modul_2pow_mp(residueul_t, unsigned long const *, int,
                   modulusul_t const);
void modul_V_ul(residueul_t, residueul_t const, unsigned long,
                modulusul_t const);
void modul_V_mp(residueul_t, residueul_t const, unsigned long const *,
                int, modulusul_t);
int modul_sprp(residueul_t const, modulusul_t const);
int modul_sprp2(modulusul_t const);
int modul_isprime(modulusul_t const);
int modul_inv(residueul_t, residueul_t const, modulusul_t const);
int modul_inv_odd(residueul_t, residueul_t const, modulusul_t const);
int modul_inv_powerof2(residueul_t, residueul_t const, modulusul_t const);
int modul_batchinv(residueul_t *, residueul_t const *, size_t,
                   residueul_t const, modulusul_t const);
int modul_jacobi(residueul_t const, modulusul_t const);

#ifdef __cplusplus
}
#endif

#endif /* CADO_MOD_UL_H */
