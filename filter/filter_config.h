#ifndef FILTER_CONFIG_H_
#define FILTER_CONFIG_H_


/* Duplicate removal */
/*********************/

// These 4 constants are used for defining the hash functions used by dup.
#define CA_DUP1 314159265358979323UL
#define CB_DUP1 271828182845904523UL

#define CA_DUP2 271828182845904523UL
#define CB_DUP2 577215664901532889UL

/* Purge (a.k.a singleton and clique removal) */
/**********************************************/

#define DEFAULT_PURGE_NTHREADS 1

// Default number of step of clique removal in purge.
#define DEFAULT_PURGE_NSTEPS 50

// How many more relations than prime ideals do we want?
// A value of 0.1 means 10% more relations than ideals after singleton removal.
// (Please update params.c90 if you change this value, and also
// scripts/cadofactor/README.)
#define DEFAULT_PURGE_REQUIRED_EXCESS 0.0

// Keep that many more relations (used by purge and merge).
#define DEFAULT_FILTER_EXCESS 160

/* Merge */
/*********/

/* Define DEBUG_STACK to replace the stack allocations of arrays of size
   MERGE_LEVEL_MAX x MERGE_LEVEL_MAX by calls to malloc. */
// #define DEBUG_STACK
/* Maximum level for a merge. */
#define MERGE_LEVEL_MAX 32
/* Maximum number of characters of a merge in the history file.
   For a k-merge, we have k-1 lines written in the same string.
   Each line consists of k numbers written in decimal:
   - if __SIZEOF_INDEX__=4, each number is a 32-bit value (<=10 digits)
   - if __SIZEOF_INDEX__=8, each number is a 64-bit value (<=20 digits)
   We also have k-1 spaces per line, an optional minus sign for the first
   number, a final \n, and (for DLP) an extra number of 10/20
   digits and two extra characters (" #"). Plus an extra \0 (end of string).
   This gives a total of:
   - (k-1)*(k*10*__SIZEOF_INDEX__/4 + k+1)+1 for the factorization
   - (k-1)*((k+1)*10*__SIZEOF_INDEX__/4 + k+3)+1 for DLP
   For MERGE_LEVEL_MAX=32 this gives:
   - __SIZEOF_INDEX__=4: 10944 for the factorization, 11316 for DLP
   - __SIZEOF_INDEX__=8: 20864 for the factorization, 21546 for DLP
*/
#ifndef FOR_DL
#define MERGE_CHAR_MAX ((MERGE_LEVEL_MAX - 1) * (MERGE_LEVEL_MAX * 10 * (__SIZEOF_INDEX__ / 4) + MERGE_LEVEL_MAX + 1) + 1)
#else
#define MERGE_CHAR_MAX ((MERGE_LEVEL_MAX - 1) * ((MERGE_LEVEL_MAX + 1) * 10 * (__SIZEOF_INDEX__ / 4) + MERGE_LEVEL_MAX + 3) + 1)
#endif

#ifndef FOR_DL
/* the default value 170 was determined experimentally on RSA-155
   with the git version 0a6a50c */
#define DEFAULT_MERGE_TARGET_DENSITY 170.0
#else /* for discrete log, a smaller density is better */
#define DEFAULT_MERGE_TARGET_DENSITY 100.0
#endif

#ifndef FOR_DL
#define DEFAULT_MERGE_SKIP 32 /* for factorization */
#else
#define DEFAULT_MERGE_SKIP 0  /* for discrete logarithm */
#endif


/* Some utilities */
/******************/
// (not sure they belong here...)

#define STR(s) XSTR(s)
#define XSTR(s) #s


#include "utils.h"  /* for ideal_merge_t */

static inline int
cmp_index (const void *p, const void *q)
{
  index_t x = *((index_t *)p);
  index_t y = *((index_t *)q);
  return (x <= y ? -1 : 1);
}

static inline int
cmp_ideal_merge (const void *p, const void *q)
{
  ideal_merge_t x = *((ideal_merge_t *)p);
  ideal_merge_t y = *((ideal_merge_t *)q);
  return (x.id <= y.id ? -1 : 1);
}

/* We also compare x[1] and y[1] to make the code deterministic
   since in case x[0] = y[0], qsort() may give different results on
   different machines */
static inline int
cmp_int2 (const void *p, const void *q)
{
  int *x = (int*) p;
  int *y = (int*) q;

  if (x[0] < y[0])
    return -1;
  else if (x[0] > y[0])
    return 1;
  else
    return (x[1] < y[1]) ? 1 : -1;
}

#endif /* FILTER_CONFIG_H_ */
