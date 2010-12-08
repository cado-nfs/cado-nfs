/*****************************************************************
*                Functions for the factor base                  *
*****************************************************************/

#ifndef FB_H
#define FB_H

#include <stdio.h>
#include <gmp.h>

/* Data types */

typedef unsigned int fbprime_t; /* 32 bits should be enough for everyone */
#define FBPRIME_FORMAT "%u"
#define FBPRIME_MAX 4294967295U
#define FBPRIME_BITS 32
typedef fbprime_t fbroot_t;
#define FBROOT_FORMAT "%u"
typedef unsigned long largeprime_t; /* On IA32 they'll only get 32 bit 
                                       large primes */
#define LARGEPRIME_FORMAT "%lu"

/* Factor base entry with (possibly) several roots */
typedef struct {
  fbprime_t p;            /* A prime or a prime power */
  unsigned long invp;     /* -1/p (mod 2^wordsize) for REDC: although we need
			     only a 32-bit inverse in say redc_32, we need a
			     full-limb inverse on 64-bit machines for trial
			     division */
  unsigned char plog;     /* logarithm (to some suitable base) of this prime */
  unsigned char nr_roots; /* how many roots there are for this prime */
  unsigned char size;     /* The length of the struct in bytes */
  unsigned char dummy[1]; /* For dword aligning the roots. In fact uneeded, C99
                             guarantees proper alignment of roots[]. It's only a
                             precaution against -fpack-struct or other
                             ABI-breaking behaviours */
  fbroot_t roots[0];      /* the actual length of this array is determined
                             by nr_roots */
} factorbase_degn_t;

#define FB_END ((fbprime_t) 1)

void		fb_fprint_entry (FILE *, const factorbase_degn_t *);
void            fb_fprint (FILE *, const factorbase_degn_t *);
void            fb_sortprimes (fbprime_t *, const unsigned int);
unsigned char	fb_log (double, double, double);
factorbase_degn_t * 	fb_make_linear (const mpz_t *, const fbprime_t, 
					const fbprime_t, const double, 
					const int, const int, FILE *);
factorbase_degn_t *	fb_read (const char *, const double, const int);
factorbase_degn_t *     fb_read_addproj (const char *, const double, const int,
					 const fbprime_t *);
fbprime_t	*fb_extract_bycost (const factorbase_degn_t *, 
                                    const fbprime_t, const fbprime_t costlim);
size_t          fb_size (const factorbase_degn_t *);                   

/* Some inlined functions which need to be fast */
  
/* Hack to get around C's automatic multiplying constants to be added to 
   pointers by the pointer's base data type size */

__attribute__ ((unused))
static inline factorbase_degn_t *
fb_skip (const factorbase_degn_t *fb, const size_t s)
{
  return (factorbase_degn_t *)((char *)fb + s);
}
  
__attribute__ ((unused))
static inline factorbase_degn_t *
fb_next (const factorbase_degn_t *fb)
{
  return (factorbase_degn_t *)((char *)fb + fb->size);
}

__attribute__ ((unused))
static size_t
fb_entrysize (const factorbase_degn_t *fb)
{
  return (sizeof (factorbase_degn_t) + fb->nr_roots * sizeof (fbroot_t));
}

__attribute__ ((unused))
static size_t
fb_entrysize_uc (const unsigned char n)
{
  return (sizeof (factorbase_degn_t) + n * sizeof (fbroot_t));
}

#endif
