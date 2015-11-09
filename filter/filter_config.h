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

#define DEFAULT_PURGE_NTHREADS 4

// Default number of step of clique removal in purge.
#define DEFAULT_PURGE_NSTEPS 50

// How many more relations than prime ideals do we want?
// A value of 0.1 means 10% more relations than ideals after singleton removal.
#define DEFAULT_PURGE_REQUIRED_EXCESS 0.0

// Keep that many more relations (used by purge and merge).
#define DEFAULT_FILTER_EXCESS 160

/* Merge */
/*********/

/* Maximum level for a merge. Such a large value is only useful when not using
 * BW. */
#define MERGE_LEVEL_MAX 256

#define DEFAULT_MERGE_MAXLEVEL 10
#define DEFAULT_MERGE_FORBW 0
#define DEFAULT_MERGE_RATIO 1.1
#define DEFAULT_MERGE_COVERNMAX 100.0
#define DEFAULT_MERGE_MKZTYPE 1 /* pure Markowitz */
#define DEFAULT_MERGE_WMSTMAX 7 /* relevant only if mkztype == 2 */

#ifndef FOR_DL
#define DEFAULT_MERGE_SKIP 32
#else
#define DEFAULT_MERGE_SKIP 0
#endif


/* Some utilities */
/******************/
// (not sure they belong to here...)

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
