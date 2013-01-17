#ifndef POLYSELECT2L_STR_H
#define POLYSELECT2L_STR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <gmp.h>
#include <pthread.h>
#include "cado.h"
#include "utils.h"
#include "auxiliary.h"
#include "murphyE.h"

#define LEN_SPECIAL_Q 55
//#define DEBUG_HASH_TABLE

#ifndef INLINE
# if __GNUC__ && !__GNUC_STDC_INLINE__
#  define INLINE extern inline
# else
#  define INLINE inline
# endif
#endif

/* For the moment, this value is static. But it's highly critical for
   performance in the shash table cribble:
   * 10 (or 9) seems the best for an Intel nehalem (?).
   * 6 seems the best for Intel Core2 (?).
   * 7 seems the best for AMD (?).
   So, the next optimization will include real time tests to
   evaluate the best value.
   NOTA: The good range is between 6 and 10. Don't use values <= 4!
   Values >= 12 are not interesting.
*/
#define LN2SHASH_NBUCKETS 8

/* hash table slots */
typedef struct
{
  int64_t i;               /* contains the values of r such that p^2
                              divides N - (m0 + r)^2 */
  uint32_t p;              /* contains the primes */
} __slot_t;
typedef __slot_t slot_t;

/* hash table structure */
typedef struct
{
  slot_t *slot;
  unsigned int alloc;      /* total allocated size */
  unsigned int size;       /* number of entries in hash table */
#ifdef DEBUG_HASH_TABLE
  unsigned long coll;
  unsigned long coll_all;
#endif
} __hash_struct;
typedef __hash_struct hash_t[1];

#define SHASH_NBUCKETS (1<<LN2SHASH_NBUCKETS)

typedef struct
{
  uint64_t *current[SHASH_NBUCKETS+1]; /* +1 for guard */
  uint64_t *base[SHASH_NBUCKETS+1];    /* +1 for guard */
  uint64_t *mem;
  uint32_t alloc;      /* total allocated size */
  uint32_t balloc;     /* allocated size for each bucket */
} __shash_struct;
typedef __shash_struct shash_t[1];

/* thread structure */
typedef struct
{
  mpz_t N;
  unsigned int d;
  uint64_t ad;
  int thread;
} __tab_struct;
typedef __tab_struct tab_t[1];

/* structure to store P roots */
typedef struct
{
  unsigned long size;    /* used size */
  uint8_t *nr;     /* number of roots of x^d = N (mod p) */
  uint64_t **roots; /* roots of (m0+x)^d = N (mod p^2) */
} __proots_struct;
typedef __proots_struct proots_t[1];

/* structure to store q roots */
typedef struct
{
  unsigned int alloc;   /* allocated size */
  unsigned int size;    /* used size */
  unsigned int *q;
  unsigned int *nr;     /* number of roots of x^d = N (mod p) */
  uint64_t **roots; /* roots of (m0+x)^d = N (mod p^2) */
} __qroots_struct;
typedef __qroots_struct qroots_t[1];

/* structure to store information on N, d, ad, etc... */
typedef struct
{
  mpz_t N;
  unsigned long d;
  uint64_t ad;
  mpz_t Ntilde;
  mpz_t m0;
} _header_struct;
typedef _header_struct header_t[1];

/* inline functions */

#ifndef EMIT_ADDRESSABLE_shash_add
INLINE
#endif
void
shash_add (shash_t H, uint64_t i)
{
  /*
  unsigned int j;
  uint64_t **cur;
  */

  *(H->current[i & (SHASH_NBUCKETS - 1)])++ = i;
  /*
  if (UNLIKELY(H->current[i & (SHASH_NBUCKETS - 1)] >= H->base[i & (SHASH_NBUCKETS - 1)]))
    {
      fprintf (stderr, "A Shash bucket is full.\n");
      exit (1);
    }
  */
}

/* declarations */

extern const unsigned int SPECIAL_Q[];

unsigned long initPrimes (unsigned long, uint32_t**);
void printPrimes (uint32_t*, unsigned long);
void clearPrimes (uint32_t**);

void header_init (header_t, mpz_t, unsigned long, uint64_t);
void header_clear (header_t);

void proots_init (proots_t, unsigned long);
void proots_add (proots_t, unsigned long, uint64_t*, unsigned long);
void proots_print (proots_t, unsigned long);
void proots_clear (proots_t, unsigned long);

void qroots_init (qroots_t);
void qroots_realloc (qroots_t, unsigned long);
void qroots_add (qroots_t, unsigned int, unsigned int, uint64_t*);
void qroots_print (qroots_t);
void qroots_rearrange (qroots_t R);
void qroots_clear (qroots_t);

void hash_init (hash_t, unsigned int);
void shash_init (shash_t, unsigned int);
void hash_add (hash_t, unsigned long, int64_t, mpz_t, uint64_t,
               unsigned long, mpz_t, unsigned long, mpz_t);
int shash_find_collision (shash_t);
void gmp_hash_add (hash_t, uint32_t, int64_t, mpz_t, uint64_t,
		   unsigned long, mpz_t, uint64_t, mpz_t);
void hash_grow (hash_t);
void hash_clear (hash_t);
void shash_clear (shash_t);

void print_poly_info (mpz_t *, unsigned int d, mpz_t *, int);

#endif
