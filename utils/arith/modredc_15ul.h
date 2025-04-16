/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. The modulus can be up to 1.5 unsigned longs 
   in size (meaning 48 bits if unsigned long has 32 bits, or 96 bits if 
   unsigned long has 64 bits). Residues are stored in Montgomery form,
   reduction after multiplication is done with REDC. Due to inlining, 
   this file must be included in the caller's source code with #include */

/* Naming convention: all function start with modredc15ul, for 
   MODulus REDC 1.5 Unsigned Longs, followed by underscore, functionality of 
   function (add, mul, etc), and possibly underscore and specification of 
   what argument types the function takes (_ul, etc). */

#ifndef CADO_MODREDC_15UL_H
#define CADO_MODREDC_15UL_H

/**********************************************************************/
#include "cado_config.h"  // just because we're a header.
#include <stdlib.h>       // for size_t, abort
#if defined(MODTRACE)
#include <stdio.h>
#endif
#include <limits.h>
#include <stdint.h>
#include "macros.h"
#include "ularith.h"

// as it turns out, modintredc15ul_t and modintredc2ul2_t are the same
// type, and there's no point in having diverging implementations.
#include "modredc_2ul2.h"

/*********************************************************************/
/* Helper macros */

/* A macro for function renaming. All functions here start with 
   modredc15ul_ */
#define MODREDC15UL_RENAME(x) modredc15ul_##x

#define MODREDC15UL_SIZE 2
#define MODREDC15UL_MINBITS ULONG_BITS
#define MODREDC15UL_MAXBITS (ULONG_BITS + ULONG_BITS/2)

typedef residueredc2ul2_t residueredc15ul_t;
typedef modintredc2ul2_t modintredc15ul_t;

/* We can't simply typedef, because of the nasty tricks played by
 * facul_doit_onefm, which is only resolved based on argument overloads.
 *
 * The stopgap solution is to redeclare a struct, which must be exactly
 * the same as modulusredc2ul2_t, of course. And then do wild casts.
 * Eeek...
 *
 * Pretty much everything here piggies-back onto modredc2ul2.
 */

// typedef modulusredc2ul2_t modulusredc15ul_t;
struct __modulusredc15ul_t { 
  modintredc15ul_t m;
  residueredc15ul_t one;
  unsigned long invm;
};
typedef struct __modulusredc15ul_t modulusredc15ul_t[1];

/* ==================== Functions used internally ==================== */

static inline void
modredc15ul_add (residueredc15ul_t r, const residueredc15ul_t a, 
		 const residueredc15ul_t b, const modulusredc15ul_t m);
static inline void
modredc15ul_get_int (modintredc15ul_t r, const residueredc15ul_t s, 
		     const modulusredc15ul_t m MAYBE_UNUSED);

MAYBE_UNUSED
static inline void
modredc15ul_tomontgomery (residueredc15ul_t r, const residueredc15ul_t s,
			  const modulusredc15ul_t m)
{
    modredc2ul2_tomontgomery(r, s, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_frommontgomery (residueredc15ul_t r, const residueredc15ul_t s,
			    const modulusredc15ul_t m)
{
    modredc2ul2_frommontgomery(r, s, (const struct __modulusredc2ul2_t *) m);
}

/* Do a one-word REDC, i.e., divide by 2^ULONG_BITS */
MAYBE_UNUSED
static inline void
modredc15ul_redc1 (residueredc15ul_t r, const residueredc15ul_t s,
		   const modulusredc15ul_t m)
{
    modredc2ul2_redc1(r, s, (const struct __modulusredc2ul2_t *) m);
}

/* ================= Functions that are part of the API ================= */

/* Some functions for integers of the same width as the modulus */
/* All functions here are just aliases for the 2ul2 functions.
 */

MAYBE_UNUSED
static inline void
modredc15ul_intinit (modintredc15ul_t r)
{
    modredc2ul2_intinit(r);
}


MAYBE_UNUSED
static inline void
modredc15ul_intclear (modintredc15ul_t r MAYBE_UNUSED)
{
    modredc2ul2_intclear(r);
}


MAYBE_UNUSED
static inline void
modredc15ul_intset (modintredc15ul_t r, const modintredc15ul_t s)
{
    modredc2ul2_intset(r, s);
}

MAYBE_UNUSED
static inline void
modredc15ul_intset_ul (modintredc15ul_t r, const unsigned long s)
{
    modredc2ul2_intset_ul(r, s);
}

MAYBE_UNUSED
static inline void
modredc15ul_intset_uls (modintredc15ul_t r, const unsigned long *s, 
                        const size_t n)
{
    modredc2ul2_intset_uls(r, s, n);
}

/* Get the least significant unsigned long of r */
MAYBE_UNUSED
static inline unsigned long 
modredc15ul_intget_ul (const modintredc15ul_t r)
{
    return modredc2ul2_intget_ul(r);
}

MAYBE_UNUSED
static inline size_t  
modredc15ul_intget_uls (unsigned long *r, const modintredc15ul_t s)
{
    return modredc2ul2_intget_uls(r, s);
}

MAYBE_UNUSED
static inline double
modredc15ul_intget_double (const modintredc15ul_t s)
{
    return modredc2ul2_intget_double(s);
}

MAYBE_UNUSED
static inline int
modredc15ul_intequal (const modintredc15ul_t a, const modintredc15ul_t b)
{
    return modredc2ul2_intequal(a, b);
}

MAYBE_UNUSED
static inline int
modredc15ul_intequal_ul (const modintredc15ul_t a, const unsigned long b)
{
    return modredc2ul2_intequal_ul(a, b);
}

/* Returns 1 if a < b, 0 otherwise */
MAYBE_UNUSED
static inline int
modredc15ul_intlt (const modintredc15ul_t a, const modintredc15ul_t b)
{
    return modredc2ul2_intlt(a, b);
}

MAYBE_UNUSED
static inline int
modredc15ul_intcmp (const modintredc15ul_t a, const modintredc15ul_t b)
{
    return modredc2ul2_intcmp(a, b);
}

MAYBE_UNUSED
static inline int
modredc15ul_intcmp_ul (const modintredc15ul_t a, const unsigned long b)
{
    return modredc2ul2_intcmp_ul(a, b);
}

MAYBE_UNUSED
static inline int
modredc15ul_intcmp_uint64 (const modintredc15ul_t a, const uint64_t b)
{
    return modredc2ul2_intcmp_uint64(a, b);
}

MAYBE_UNUSED
static inline int
modredc15ul_intfits_ul (const modintredc15ul_t a)
{
    return modredc2ul2_intfits_ul(a);
}

MAYBE_UNUSED
static inline void
modredc15ul_intadd (modintredc15ul_t r, const modintredc15ul_t a,
		    const modintredc15ul_t b)
{
    modredc2ul2_intadd(r, a, b);
}

MAYBE_UNUSED
static inline void
modredc15ul_intsub (modintredc15ul_t r, const modintredc15ul_t a,
		    const modintredc15ul_t b)
{
    modredc2ul2_intsub(r, a, b);
}

/* Returns the number of bits in a, that is, floor(log_2(a))+1. 
   For a == 0 returns 0. */
MAYBE_UNUSED
static inline size_t 
modredc15ul_intbits (const modintredc15ul_t a)
{
    return modredc2ul2_intbits(a);
}


/* r = trunc(s / 2^i) */
MAYBE_UNUSED
static inline void
modredc15ul_intshr (modintredc15ul_t r, const modintredc15ul_t s, const unsigned int i)
{
    modredc2ul2_intshr(r, s, i);
}


/* r = (s * 2^i) % (2^(2 * ULONG_BITS)) */
MAYBE_UNUSED
static inline void
modredc15ul_intshl (modintredc15ul_t r, const modintredc15ul_t s, const unsigned int i)
{
    modredc2ul2_intshl(r, s, i);
}


/* r = n/d. We require d|n */
MAYBE_UNUSED
static inline void
modredc15ul_intdivexact (modintredc15ul_t r, const modintredc15ul_t n,
                         const modintredc15ul_t d)
{
    return modredc2ul2_intdivexact(r, n, d);
}

MAYBE_UNUSED
static inline unsigned long
modredc15ul_intmod_ul (const modintredc15ul_t n, const unsigned long d)
{
    return modredc2ul2_intmod_ul(n, d);
}

/* r = n%d */
MAYBE_UNUSED
static inline void
modredc15ul_intmod (modintredc15ul_t r, const modintredc15ul_t n,
                    const modintredc15ul_t d)
{
    modredc2ul2_intmod(r, n, d);
}


/* Functions for the modulus */

/* Init the modulus from a modintredc15ul_t. */
MAYBE_UNUSED
static inline void
modredc15ul_initmod_int (modulusredc15ul_t m, const modintredc15ul_t s)
{
    modredc2ul2_initmod_int((struct __modulusredc2ul2_t *) m, s);
}


/* Returns the modulus as an modintredc15ul_t. */
MAYBE_UNUSED
static inline void
modredc15ul_getmod_int (modintredc15ul_t r, const modulusredc15ul_t m)
{
    modredc2ul2_getmod_int(r, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_clearmod (modulusredc15ul_t m)
{
    modredc2ul2_clearmod((struct __modulusredc2ul2_t *) m);
}


/* Functions for residues */

/* Initialises a residueredc15ul_t and sets it to zero */
MAYBE_UNUSED
static inline void
modredc15ul_init (residueredc15ul_t r, const modulusredc15ul_t m MAYBE_UNUSED)
{
    modredc2ul2_init(r, (const struct __modulusredc2ul2_t *) m);
}


/* Initialises a residueredc15ul_t, but does not set it to zero. For fixed 
   length residueredc15ul_t, that leaves nothing to do at all. */
MAYBE_UNUSED
static inline void
modredc15ul_init_noset0 (residueredc15ul_t r MAYBE_UNUSED, 
			 const modulusredc15ul_t m MAYBE_UNUSED)
{
    modredc2ul2_init_noset0(r, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_clear (residueredc15ul_t r MAYBE_UNUSED, 
		   const modulusredc15ul_t m MAYBE_UNUSED)
{
    modredc2ul2_clear(r, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_set (residueredc15ul_t r, const residueredc15ul_t s, 
		 const modulusredc15ul_t m MAYBE_UNUSED)
{
    modredc2ul2_set(r, s, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_set_ul (residueredc15ul_t r, const unsigned long s, 
		    const modulusredc15ul_t m)
{
    modredc2ul2_set_ul(r, s, (const struct __modulusredc2ul2_t *) m);
}


/* Sets the residueredc15ul_t to the class represented by the integer s. 
   Assumes that s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED
static inline void
modredc15ul_set_ul_reduced (residueredc15ul_t r, const unsigned long s, 
			    const modulusredc15ul_t m MAYBE_UNUSED)
{
    modredc2ul2_set_ul_reduced(r, s, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_set_int (residueredc15ul_t r, const modintredc15ul_t s, 
		     const modulusredc15ul_t m)
{
    modredc2ul2_set_int(r, s, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_set_int_reduced (residueredc15ul_t r, const modintredc15ul_t s, 
			     const modulusredc15ul_t m)
{
    modredc2ul2_set_int_reduced(r, s, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED 
static inline void 
modredc15ul_set0 (residueredc15ul_t r, const modulusredc15ul_t m MAYBE_UNUSED) 
{ 
    modredc2ul2_set0(r, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED 
static inline void 
modredc15ul_set1 (residueredc15ul_t r, const modulusredc15ul_t m) 
{ 
    modredc2ul2_set1(r, (const struct __modulusredc2ul2_t *) m);
}


/* Exchanges the values of the two arguments */

MAYBE_UNUSED
static inline void
modredc15ul_swap (residueredc15ul_t a, residueredc15ul_t b, 
		  const modulusredc15ul_t m MAYBE_UNUSED)
{
    modredc2ul2_swap(a, b, (const struct __modulusredc2ul2_t *) m);
}


/* Returns the least significant unsigned long of the residue. How to signal
   if the residue does not fit in one unsigned long? */

MAYBE_UNUSED
static inline unsigned long
modredc15ul_get_ul (const residueredc15ul_t s, 
		    const modulusredc15ul_t m MAYBE_UNUSED)
{
    return modredc2ul2_get_ul(s, (const struct __modulusredc2ul2_t *) m);
}


/* Returns the residue as a modintredc15ul_t */

MAYBE_UNUSED
static inline void
modredc15ul_get_int (modintredc15ul_t r, const residueredc15ul_t s, 
		     const modulusredc15ul_t m MAYBE_UNUSED)
{
    modredc2ul2_get_int(r, s, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline int
modredc15ul_equal (const residueredc15ul_t a, const residueredc15ul_t b, 
		   const modulusredc15ul_t m MAYBE_UNUSED)
{
    return modredc2ul2_equal(a, b, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline int
modredc15ul_is0 (const residueredc15ul_t a, const modulusredc15ul_t m MAYBE_UNUSED)
{
    return modredc2ul2_is0(a, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline int
modredc15ul_is1 (const residueredc15ul_t a, const modulusredc15ul_t m MAYBE_UNUSED)
{
    return modredc2ul2_is1(a, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_add (residueredc15ul_t r, const residueredc15ul_t a, 
		 const residueredc15ul_t b, const modulusredc15ul_t m)
{
    modredc2ul2_add(r, a, b, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_sub (residueredc15ul_t r, const residueredc15ul_t a, 
		 const residueredc15ul_t b, const modulusredc15ul_t m)
{
    modredc2ul2_sub(r, a, b, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_add1 (residueredc15ul_t r, const residueredc15ul_t a, 
		  const modulusredc15ul_t m)
{
    modredc2ul2_add1(r, a, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_add_ul (residueredc15ul_t r, const residueredc15ul_t a,
		    const unsigned long b, const modulusredc15ul_t m)
{
    modredc2ul2_add_ul(r, a, b, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_sub_ul (residueredc15ul_t r, const residueredc15ul_t a,
		    const unsigned long b, const modulusredc15ul_t m)
{
    modredc2ul2_sub_ul(r, a, b, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_neg (residueredc15ul_t r, const residueredc15ul_t a, 
		 const modulusredc15ul_t m)
{
    modredc2ul2_neg(r, a, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline void
modredc15ul_div2 (residueredc15ul_t r, const residueredc15ul_t a, 
		  const modulusredc15ul_t m)
{
    modredc2ul2_div2(r, a, (const struct __modulusredc2ul2_t *) m);
}


#ifdef WANT_ASSERT_EXPENSIVE
#if defined(__x86_64__)
#define ABORT_IF_CY "jnc 1f\n\tlea _GLOBAL_OFFSET_TABLE_(%%rip), %%rbx\n\tcall abort@plt\n1:\n\t"
#elif defined(__i386__)
#define ABORT_IF_CY "jnc 1f\n\tcall abort\n1:\n\t"
#endif
#else
#define ABORT_IF_CY
#endif

MAYBE_UNUSED
static inline void
modredc15ul_mul (residueredc15ul_t r, const residueredc15ul_t a, 
                 const residueredc15ul_t b, const modulusredc15ul_t m)
{
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM

  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_intlt (b, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu) * (%lu * 2^%d + %lu) / 2^%d) %% "
	  "(%lu * 2^%d + %lu)", 
          a[1], ULONG_BITS, a[0], b[1], ULONG_BITS, b[0], 2 * ULONG_BITS, 
	  m[0].m[1], ULONG_BITS, m[0].m[0]);
#endif

/* Since m1>0, m*u is maximal for m0=1 and u=2^64-1, so
   u*m is bounded by (2^96 - 2^64 + 1)*(2^64 - 1) = 
   2^160 - 2^128 - 2^96 - 1. Doesn't really save anything, tho */

  unsigned long dummy;
  __asm__ __VOLATILE (
    /* Product of low words */
// #define MODREDCUL15_VERBOSE_ASM 1
#ifdef MODREDCUL15_VERBOSE_ASM
    "# modredc15ul_mul(): asm output operands:\n\t"
    "# t0: %[t0], t1: %[t1], t2: %[t2]\n\t"
    "# modredc15ul_mul(): asm input operands:\n\t"
    "# a0: %[a0], a1: %[a1], b0: %[b0], b1: %[b1]\n\t"
    "# m0: %[m0], m1: %[m1], invm: %[invm]\n\t"
#endif
    "movq %[a0], %%rax\n\t"
    "mulq %[b0]\n\t"         /* rdx:rax = a0*b0 */
    "movq %%rdx, %[t0]\n\t"
    /* Compute u0*m, add to t0:rax */
    "imulq %[invm], %%rax\n\t"
    "movq %%rax, %[t2]\n\t"  /* t2 = u0 */
    "xorl %k[t1], %k[t1]\n\t"
    "mulq %[m0]\n\t"         /* rdx:rax = u0*m0 <= (2^64-1)^2 */
    "negq %%rax\n\t"         /* if low word != 0, carry to high word */
    "movq %[t2], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t0]\n\t"
    "setc %b[t1]\n\t"        /* t1:t0 = (a0*b0+u0*m0)/2^64 <= (2^64-1)^2/2^64 = 2*2^64-4 */
    "mulq %[m1]\n\t"         /* rdx:rax <= (2^64-1)*(2^32-1) = 2^96-2^64-2^32+1 */
    "addq %%rax, %[t0]\n\t"
    "movq %[a0], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0*b0+u0*m)/2^64 */
    ABORT_IF_CY              /* <= ((2^64-1)^2 + 2^160 - 2^128 - 2^96 - 1)/2^64
                                <= 2^96 - 2^32 - 2 */
    
    /* 2 products of low and high word */
    "xorl %k[t2], %k[t2]\n\t"
    "mulq %[b1]\n\t"         /* rdx:rax = a0*b1 <= (2^64-1)*(2^32-1) */
    "addq %%rax, %[t0]\n\t"
    "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0*b0+u0*m)/2^64 + a0*b1 */
    ABORT_IF_CY              /* <= 2^96 - 2^32 - 2 + 2*(2^64-1)*(2^32-1) 
                               = 2*2^96 - 2^64 - 2*2^32 - 1 */
    /* Free slot here */
    "mulq %[b0]\n\t"         /* rdx:rax = a1*b0 <= (2^64-1)*(2^32-1) */
    "addq %%rax, %[t0]\n\t"
    "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0*b0+u0*m)/2^64 + a0*b1 + a1*b0 */
    ABORT_IF_CY              /* <= 2*2^96 - 2^64 - 2*2^32 - 1 + (2^64-1)*(2^32-1)
                                =  3*2^96 - 2*2^64 - 3*2^32 */
    /* Free slot here */
    /* Product of high words */
    "imulq %[b1], %%rax\n\t" /* rax = a1*b1 <= (2^32-1)^2 = 2^64 - 2*2^32 + 1 */
    "addq %%rax, %[t1]\n\t"
    "setc %b[t2]\n\t"        /* t2:t1:rax = (a0*b0+u0*m)/2^64 + a0*b1 + a1*b0 + a1*b1*2^64
                                <= 3*2^96 - 2*2^64 - 3*2^32 + (2^32-1)^2*2^64
                                = 3*2^96 - 2*2^64 - 3*2^32 + 2^128 - 2*2^96 + 2^64
                                = 2^128 + 2^96 - 2^64 - 3*2^32 */
    "movq %[t0], %%rax\n\t"
    /* Compute u1*m, add to t2:t1:t0 */
    "imulq %[invm], %%rax\n\t"
    "movq %%rax, %[t0]\n\t" /* t0 = u1 */
    "mulq %[m0]\n\t"       /* rdx:rax = u1*m0 <= (2^64-1)^2 = 2^128 - 2*2^64 + 1 */
    "negq %%rax\n\t"       /* if low word != 0, carry to high word */
    "movq %[t0], %%rax\n\t"
    "adcq %%rdx, %[t1]\n\t"
    "adcq $0,%[t2]\n\t"    /* t2:t1:0 = (a0*b0+u0*m)/2^64 + a0*b1 + a1*b0 + a1*b1*2^64 + u1*m0 */
    ABORT_IF_CY            /* <= 2^128 + 2^96 - 2^64 - 3*2^32 + 2^128 - 2*2^64 + 1
                              = 2*2^128 + 2^96 - 3*2^64 - 3*2^32 + 1 */
                           /* t2:t1 = ((a0*b0+u0*m)/2^64 + a0*b1 + a1*b0  + u1*m0)/2^64 + a1*b1
                              <= 2*2^64 + 2^32 - 4 */

    "mulq %[m1]\n\t"       /* rdx:rax = u1*m1 <= (2^64-1)*(2^32-1) */
    "addq %%rax, %[t1]\n\t"
    "adcq %%rdx, %[t2]\n\t"/* t2:t1 = ((a0*b0+u0*m)/2^64 + a0*b1 + a1*b0 + u1*m)/2^64 + a1*b1 */
    ABORT_IF_CY            /* <= (2^128 + 2^96 - 2^64 - 3*2^32 + 2^160 - 2^128 - 2^96 - 1)/2^64
                              = (2^160 - 2^64 - 3*2^32 - 1)/2^64
                              <= 2^96 - 2 */

  /* t2:t1:0:0 = a*b + u0*m + u1*m*2^64
     t2:t1 <= (a*b + u0*m + u1*m*2^64) / 2^128
     <= (m^2 + 2^64*m + 2^64*(2^64-1)*m) / 2^128
     =  (m^2 + 2^64*m + (2^128-2^64)*m)/2^128
     =  m + (m^2)/2^128
     <= m + (2^96*m)/2^128
     <= m + m/2^32 */
    "movq %[t1], %%rax\n\t" /* See if result > m */
    "movq %[t2], %%rdx\n\t"
    "subq %[m0], %[t1]\n\t"
    "sbbq %[m1], %[t2]\n\t"
    "cmovc %%rax, %[t1]\n\t" /* No carry -> copy new result */
    "cmovc %%rdx, %[t2]\n\t"
    : [t0] "=&r" (dummy), [t1] "=&r" (r[0]), [t2] "=&r" (r[1])
    : [a0] "rme" (a[0]), [a1] "rme" (a[1]), [b0] "rm" (b[0]), [b1] "rm" (b[1]),
      [m0] "rm" (m[0].m[0]), [m1] "rm" (m[0].m[1]), [invm] "rm" (m[0].invm)
    : "%rax", "%rdx", "cc"
  );
#else
  unsigned long pl, ph, t[3];

  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_intlt (b, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu) * (%lu * 2^%d + %lu) / 2^%d) %% "
	  "(%lu * 2^%d + %lu)", 
          a[1], ULONG_BITS, a[0], b[1], ULONG_BITS, b[0], 2 * ULONG_BITS, 
	  m[0].m[1], ULONG_BITS, m[0].m[0]);
#endif
  
  /* Product of the two low words */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), a[0], b[0]); /* t1:t0 = a[0]*b[0] <= W^2 - 2W + 1 */

  /* One REDC step */
  modredc15ul_redc1 (t, t, m);  /* t1:t0 <= W^(3/2) + W - W^(1/2) - 3 */

  /* Products of one low and one high word  */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[1], b[0]);   /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), pl, ph); /* t1:t0 <= 2W^(3/2) - 2W^(1/2) - 2 */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[0], b[1]);   /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), pl, ph); /* t1:t0 <= 3W^(3/2) - 3W^(1/2) - W - 1 */

  t[2] = 0UL;
  pl = a[1] * b[1];                               /* pl <= (W^(1/2)-1)^2 = W - 2W^(1/2) + 1 */
  ularith_add_ul_2ul (&(t[1]), &(t[2]), pl);      /* t2:t1:t0 <= W^2 + W^(3/2) - 3W^(1/2) - 1 */

  modredc2ul2_redc1_wide_inplace(t, (const struct __modulusredc2ul2_t *) m);
  r[0] = t[1];
  r[1] = t[2];
  /* r1:r0 <= W^(3/2) + W - 2 */

  /* Result may be larger than m, but is < 2*m */
  ularith_sub_2ul_2ul_ge (&(r[0]), &(r[1]), m[0].m[0], m[0].m[1]);
#endif

#if defined(MODTRACE)
  printf (" == (%lu * 2^%d + %lu) /* PARI */ \n", r[1], ULONG_BITS, r[0]);
#endif
  ASSERT_EXPENSIVE (modredc15ul_intlt (r, m[0].m));
}


MAYBE_UNUSED
static inline void
modredc15ul_sqr (residueredc15ul_t r, const residueredc15ul_t a, 
                 const modulusredc15ul_t m)
{
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu)^2 / 2^%d) %% (%lu * 2^%d + %lu)", 
          a[1], ULONG_BITS, a[0], 2 * ULONG_BITS, m[0].m[1], ULONG_BITS, m[0].m[0]);
#endif
  
  unsigned long dummy;
  __asm__ __VOLATILE (
    /* Product of low words */
    "movq %[a0], %%rax\n\t"
    "mulq %%rax\n\t"         /* rdx:rax = a0*a0 */
    "movq %%rdx, %[t0]\n\t"
    /* Compute u0*m, add to t0:rax */
    "imulq %[invm], %%rax\n\t"
    "movq %%rax, %[t2]\n\t"  /* t2 = u0 */
    "xorl %k[t1], %k[t1]\n\t"
    "mulq %[m0]\n\t"         /* rdx:rax = u0*m0 <= (2^64-1)^2 */
    "negq %%rax\n\t"         /* if low word != 0, carry to high word */
    "movq %[t2], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t0]\n\t"
    "setc %b[t1]\n\t"        /* t1:t0 = (a0*a0+u0*m0)/2^64 <= (2^64-1)^2/2^64 = 2*2^64-4 */
    "mulq %[m1]\n\t"         /* rdx:rax <= (2^64-1)*(2^32-1) = 2^96-2^64-2^32+1 */
    "addq %%rax, %[t0]\n\t"
    "movq %[a0], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t1]\n\t"  /* t0:t1 = (a0*a0+u0*m)/2^64 */
    ABORT_IF_CY              /* <= ((2^64-1)^2 + 2^160 - 2^96 - 2^128 - 1)/2^64
                                <= 2^96 - 2^32 - 2 */
    
    /* Product of low and high word */
    "xorl %k[t2], %k[t2]\n\t"
    "mulq %[a1]\n\t"         /* rdx:rax = a0*a1 <= (2^64-1)*(2^32-1) */
    "shlq $1,%%rax\n\t"
    "rclq $1,%%rdx\n\t"
    ABORT_IF_CY
    "addq %%rax, %[t0]\n\t"
    "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0*a0+u0*m)/2^64 + 2*a0*a1 */
    ABORT_IF_CY              /* <= 2^96 - 2^32 - 2 + 2*(2^64-1)*(2^32-1)
                                =  3*2^96 - 2*2^64 - 3*2^32 */
    "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
    /* Free slot here */
    /* Product of high words */
    "imulq %%rax, %%rax\n\t" /* rax = a1*a1 <= (2^32-1)^2 = 2^64 - 2*2^32 + 1 */
    "addq %%rax, %[t1]\n\t"
    "setc %b[t2]\n\t"        /* t2:t1:rax = (a0*a0+u0*m)/2^64 + 2*a0*a1 + a1*a1*2^64
                                <= 3*2^96 - 2*2^64 - 3*2^32 + (2^32-1)^2*2^64
                                = 3*2^96 - 2*2^64 - 3*2^32 + 2^128 - 2*2^96 + 2^64
                                = 2^128 + 2^96 - 2^64 - 3*2^32 */
    "movq %[t0], %%rax\n\t"
    /* Compute u1*m, add to t2:t1:t0 */
    "imulq %[invm], %%rax\n\t"
    "movq %%rax, %[t0]\n\t" /* t0 = u1 */
    "mulq %[m0]\n\t"       /* rdx:rax = u1*m0 <= (2^64-1)^2 = 2^128 - 2*2^64 + 1 */
    "negq %%rax\n\t"       /* if low word != 0, carry to high word */
    "movq %[t0], %%rax\n\t"
    "adcq %%rdx, %[t1]\n\t"
    "adcq $0,%[t2]\n\t"    /* t2:t1:0 = (a0*a0+u0*m)/2^64 + 2*a0*a1 + a1*a1*2^64 + u1*m0 */
    ABORT_IF_CY            /* <= 2^128 + 2^96 - 2^64 - 3*2^32 + 2^128 - 2*2^64 + 1
                              = 2*2^128 + 2^96 - 3*2^64 - 3*2^32 + 1 */
                           /* t2:t1 = ((a0*a0+u*m)/2^64 + a0*a1 + a1*a0  + u*m0)/2^64 + a1*a1
                              <= 2*2^64 + 2^32 - 4 */

    "mulq %[m1]\n\t"       /* rdx:rax = u1*m1 <= (2^64-1)*(2^32-1) */
    "addq %%rax, %[t1]\n\t"
    "adcq %%rdx, %[t2]\n\t"/* t2:t1 = ((a0*a0+u0*m)/2^64 + 2*a0*a1 + u1*m)/2^64 + a1*a1 */
    ABORT_IF_CY            /* <= (2^128 + 2^96 - 2^64 - 3*2^32 + 2^160 - 2^128 - 2^96 - 1)/2^64
                              = (2^160 - 2^64 - 3*2^32 - 1)/2^64
                              <= 2^96 - 2 */

  /* t2:t1:0:0 = a*b + u0*m + u1*m*2^64
     t2:t1 <= (a*b + u0*m + u1*m*2^64) / 2^128
     <= (m^2 + 2^64*m + 2^64*(2^64-1)*m)/2^128
     =  (m^2 + 2^64*m + (2^128-2^64)*m)/2^128
     =  m + (m^2)/2^128
     <= m + (2^96*m)/2^128
     <= m + m/2^32 */
    "movq %[t1], %%rax\n\t" /* See if result > m */
    "movq %[t2], %%rdx\n\t"
    "subq %[m0], %[t1]\n\t"
    "sbbq %[m1], %[t2]\n\t"
    "cmovc %%rax, %[t1]\n\t" /* No carry -> copy new result */
    "cmovc %%rdx, %[t2]\n\t"
    : [t0] "=&r" (dummy), [t1] "=&r" (r[0]), [t2] "=&r" (r[1])
    : [a0] "rme" (a[0]), [a1] "rme" (a[1]), [m0] "rm" (m[0].m[0]), [m1] "rm" (m[0].m[1]), 
      [invm] "rm" (m[0].invm)
    : "%rax", "%rdx", "cc"
  );
#else
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu)^2 / 2^%d) %% (%lu * 2^%d + %lu)", 
          a[1], ULONG_BITS, a[0], 2 * ULONG_BITS, m[0].m[1], ULONG_BITS, m[0].m[0]);
#endif
  
  unsigned long pl, ph, t[3];
  
  /* Square of low word */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), a[0], a[0]); /* t1:t0 = a[0]*a[0] <= W^2 - 2W + 1 */

  /* One REDC step */
  modredc15ul_redc1 (t, t, m);  /* t1:t0 <= W^(3/2) + W - W^(1/2) - 3 */

  /* Product of low and high word  */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[1], a[0]);   /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), pl, ph); /* t1:t0 <= 2W^(3/2) - 2W^(1/2) - 2 */
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), pl, ph); /* t1:t0 <= 3W^(3/2) - 3W^(1/2) - W - 1 */

  /* Square of high word */
  t[2] = 0UL;
  pl = a[1] * a[1];                               /* pl <= (W^(1/2)-1)^2 = W - 2W^(1/2) + 1 */
  ularith_add_ul_2ul (&(t[1]), &(t[2]), pl);      /* t2:t1:t0 <= W^2 + W^(3/2) - 3W^(1/2) - 1 */

  modredc2ul2_redc1_wide_inplace(t, (const struct __modulusredc2ul2_t *) m);
  r[0] = t[1];
  r[1] = t[2];
  /* r1:r0 <= W^(3/2) + W - 2 */

  /* Result may be larger than m, but is < 2*m */
  ularith_sub_2ul_2ul_ge (&(r[0]), &(r[1]), m[0].m[0], m[0].m[1]);
#endif

#if defined(MODTRACE)
  printf (" == (%lu * 2^%d + %lu) /* PARI */ \n", r[1], ULONG_BITS, r[0]);
#endif
  ASSERT_EXPENSIVE (modredc15ul_intlt (r, m[0].m));
}


MAYBE_UNUSED
static inline int
modredc15ul_next (residueredc15ul_t r, const modulusredc15ul_t m)
{
    return modredc2ul2_next(r, (const struct __modulusredc2ul2_t *) m);
}


MAYBE_UNUSED
static inline int
modredc15ul_finished (const residueredc15ul_t r, const modulusredc15ul_t m)
{
    return modredc2ul2_finished(r, (const struct __modulusredc2ul2_t *) m);
}


/* Division by small integer n, where (n-1)*m may NOT overflow the most 
   significant word. Returns 1 if n is invertible modulo m, 0 if not. 
   
   w_mod_n is word base (e.g., 2^32 or  2^64) mod n
   inv_n contains -1/i (mod n) if i is coprime to n, or 0 if i is not coprime 
   to n, for 0 <= i < n
   c = n^(-1) (mod word base)
*/

static inline int
modredc15ul_divn (residueredc15ul_t r, const residueredc15ul_t a, 
		  const unsigned long n, const unsigned long w_mod_n, 
		  const unsigned long *inv_n, const unsigned long c,
		  const modulusredc15ul_t m)
{
    return modredc2ul2_divn(r, a, n, w_mod_n, inv_n, c, (const struct __modulusredc2ul2_t *) m);
}

static inline int
modredc15ul_inv (residueredc15ul_t r, const residueredc15ul_t A, 
		 const modulusredc15ul_t m) 
{
    return modredc2ul2_inv(r, A, (const struct __modulusredc2ul2_t *) m);
}


/* prototypes of non-inline functions */
#ifdef __cplusplus
extern "C" {
#endif

int modredc15ul_div3 (residueredc15ul_t, const residueredc15ul_t, 
		      const modulusredc15ul_t);
int modredc15ul_div5 (residueredc15ul_t, const residueredc15ul_t, 
		      const modulusredc15ul_t);
int modredc15ul_div7 (residueredc15ul_t, const residueredc15ul_t, 
		      const modulusredc15ul_t);
int modredc15ul_div11 (residueredc15ul_t, const residueredc15ul_t, 
		       const modulusredc15ul_t);
int modredc15ul_div13 (residueredc15ul_t, const residueredc15ul_t, 
		       const modulusredc15ul_t);
void modredc15ul_gcd (modintredc15ul_t, const residueredc15ul_t, 
		      const modulusredc15ul_t);
void modredc15ul_pow_ul (residueredc15ul_t, const residueredc15ul_t, 
			 unsigned long, const modulusredc15ul_t);
void modredc15ul_2pow_ul (residueredc15ul_t, unsigned long, 
                          const modulusredc15ul_t);
void modredc15ul_pow_mp (residueredc15ul_t, const residueredc15ul_t, 
			 const unsigned long *, int, 
			 const modulusredc15ul_t);
void modredc15ul_2pow_mp (residueredc15ul_t, const unsigned long *, int, 
			  const modulusredc15ul_t);
void modredc15ul_V_ul (residueredc15ul_t, const residueredc15ul_t, 
		       unsigned long, const modulusredc15ul_t);
void modredc15ul_V_mp (residueredc15ul_t, const residueredc15ul_t, 
		       const unsigned long *, int, 
		       const modulusredc15ul_t);
int modredc15ul_sprp (const residueredc15ul_t, const modulusredc15ul_t);
int modredc15ul_sprp2 (const modulusredc15ul_t);
int modredc15ul_isprime (const modulusredc15ul_t);
int modredc15ul_batchinv (residueredc15ul_t *, const residueredc15ul_t *,
                          size_t, const residueredc15ul_t,
                          const modulusredc15ul_t);
int modredc15ul_jacobi (const residueredc15ul_t, const modulusredc15ul_t);

#ifdef __cplusplus
}
#endif

#endif  /* CADO_MODREDC_15UL_H */
