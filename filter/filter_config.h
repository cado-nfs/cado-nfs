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

/* Maximum level for a merge. Such a large value is only useful when not using
 * BW. */
#define MERGE_LEVEL_MAX 64

#define DEFAULT_MERGE_MAXLEVEL 10
#ifndef FOR_DL
/* the default value 170 was determined experimentally on RSA-155
   with the git version 0a6a50c */
#define DEFAULT_MERGE_TARGET_DENSITY 170.0
#else /* for discrete log, a smaller density is better */
#define DEFAULT_MERGE_TARGET_DENSITY 100.0
#endif
#define DEFAULT_MERGE_MKZTYPE MKZTYPE_PURE /* pure Markowitz */
#define DEFAULT_MERGE_WMSTMAX 7 /* relevant only if mkztype == MKZTYPE_LIGHT */

#ifndef FOR_DL
#define DEFAULT_MERGE_SKIP 32
#else
#define DEFAULT_MERGE_SKIP 0
#endif


/* Some utilities */
/******************/
// (not sure they belong here...)

#define STR(s) XSTR(s)
#define XSTR(s) #s

#endif /* FILTER_CONFIG_H_ */
