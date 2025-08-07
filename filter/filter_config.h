#ifndef CADO_FILTER_CONFIG_H
#define CADO_FILTER_CONFIG_H

#include "typedefs.h"

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

/* Maximum level for a merge. */
#define MERGE_LEVEL_MAX 32
/* Maximum number of characters of a merge in the history file.
   For a k-merge, we have up to k-1 lines written in the same string.
   Each line is of the following form:
   i1 i2 ... in
   which means that line i1 has to be added to i2, ..., in and then destroyed,
   or:
   -(i1+1) i2 ... in
   which means that line i1 has to be added to i2, ..., in, but not destroyed.
   Each line consists of numbers written in decimal:
   - if SIZEOF_INDEX=4, each number is a 32-bit value (<=10 digits)
   - if SIZEOF_INDEX=8, each number is a 64-bit value (<=20 digits)
   The worst case (in term of total length) is when we only have lines of
   the form "i1 i2" or "-(i1+1) i2". This is precisely what we get with the
   spanning tree algorithm. We then get exactly k-1 lines: first k-2 with a
   minus sign (row i1 not destroyed), and last one without minus sign (row
   destroyed).
/   We also a final \n per line, and (for DLP) an extra number of 10/20
   digits and two extra characters (" #"). Plus an extra \0 (end of string).
   This gives a maximum of:
   - (k-1)*(2*10*SIZEOF_INDEX/4 + 3) for the factorization
   - (k-1)*(3*10*SIZEOF_INDEX/4 + 5) for DLP
   For MERGE_LEVEL_MAX=32 this gives:
   - SIZEOF_INDEX=4: 713 for the factorization, 1085 for DLP
   - SIZEOF_INDEX=8: 1333 for the factorization, 2015 for DLP
*/
#ifndef FOR_DL
#define MERGE_CHAR_MAX ((MERGE_LEVEL_MAX - 1) * (2 * 10 * (SIZEOF_INDEX / 4) + 3))
#else
#define MERGE_CHAR_MAX ((MERGE_LEVEL_MAX - 1) * (3 * 10 * (SIZEOF_INDEX / 4) + 5))
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

#include "typedefs.h"  /* for ideal_merge_t */

#ifdef __cplusplus
extern "C" {
#endif

static inline int cmp_ideal_merge (const void *p, const void *q);
/* compare two index_t's */
static inline int cmp_index (const void *p, const void *q);
static inline int cmp_index2 (const void *p, const void *q);

#ifdef __cplusplus
}
#endif


static inline int
cmp_ideal_merge (const void *p, const void *q)
{
  ideal_merge_t x = *((ideal_merge_t *)p);
  ideal_merge_t y = *((ideal_merge_t *)q);
  return (x.id <= y.id ? -1 : 1);
}

/* compare two index_t's */
static inline int
cmp_index (const void *p, const void *q)
{
  index_t x = *((index_t *)p);
  index_t y = *((index_t *)q);
  return (x <= y ? -1 : 1);
}

/* Compare two pairs of index_t's.
   We also compare x[1] and y[1] to make the code deterministic
   since in case x[0] = y[0], qsort() may give different results on
   different machines. */
static inline int
cmp_index2 (const void *p, const void *q)
{
  index_t *x = (index_t*) p;
  index_t *y = (index_t*) q;

  if (x[0] < y[0])
    return -1;
  else if (x[0] > y[0])
    return 1;
  else
    return (x[1] < y[1]) ? 1 : -1;
}

static inline int cmp_typerow_t(const void * a, const void * b) {
#ifndef FOR_DL
    return cmp_index(a, b);
#else
    return cmp_ideal_merge(a, b);
#endif
}

#endif /* CADO_FILTER_CONFIG_H */
