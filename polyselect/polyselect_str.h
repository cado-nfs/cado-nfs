#ifndef POLYSELECT_STR_H
#define POLYSELECT_STR_H

#include <stdint.h>     // int64_t
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <gmp.h>
#include <limits.h> /* for ULONG_MAX */
#include "macros.h"     // LIKELY

#if ULONG_MAX == 4294967295UL
#define LEN_SPECIAL_Q 57
#else
#define LEN_SPECIAL_Q 59
#endif
//#define DEBUG_HASH_TABLE

/* hash table slots */
struct polyselect_hash_slot_s
{
  int64_t i;               /* contains the values of r such that p^2
                              divides N - (m0 + r)^2 */
  uint32_t p;              /* contains the primes */
};

/* hash table structure */
struct polyselect_hash_s
{
  struct polyselect_hash_slot_s *slot;
  unsigned int alloc;      /* total allocated size */
  unsigned int size;       /* number of entries in hash table */
#ifdef DEBUG_HASH_TABLE
  unsigned long coll;
  unsigned long coll_all;
#endif
};
typedef struct polyselect_hash_s polyselect_hash_t[1];
typedef struct polyselect_hash_s * polyselect_hash_ptr;
typedef const struct polyselect_hash_s * polyselect_hash_srcptr;

/* thread structure */
struct polyselect_thread_tab_s
{
  mpz_t N;
  unsigned int d;
  mpz_t ad;
  int thread;
};
typedef struct polyselect_thread_tab_s polyselect_thread_tab_t[1];
typedef struct polyselect_thread_tab_s * polyselect_thread_tab_ptr;
typedef const struct polyselect_thread_tab_s * polyselect_thread_tab_srcptr;

/* structure to store P roots */
struct polyselect_proots_s {
  unsigned long size;    /* used size */
  uint8_t *nr;     /* number of roots of x^d = N (mod p) */
  uint64_t **roots; /* roots of (m0+x)^d = N (mod p^2) */
};
typedef struct polyselect_proots_s polyselect_proots_t[1];
typedef struct polyselect_proots_s * polyselect_proots_ptr;
typedef const struct polyselect_proots_s * polyselect_proots_srcptr;

/* structure to store q roots */
struct polyselect_qroots_s {
  unsigned int alloc;   /* allocated size */
  unsigned int size;    /* used size */
  unsigned int *q;
  unsigned int *nr;     /* number of roots of x^d = N (mod p) */
  uint64_t **roots; /* roots of (m0+x)^d = N (mod p^2) */
};
typedef struct polyselect_qroots_s polyselect_qroots_t[1];
typedef struct polyselect_qroots_s * polyselect_qroots_ptr;
typedef const struct polyselect_qroots_s * polyselect_qroots_srcptr;

/* structure to store information on N, d, ad, etc... */
struct polyselect_poly_header_s {
  mpz_t N;
  unsigned long d;
  mpz_t ad;
  mpz_t Ntilde;
  mpz_t m0;
};
typedef struct polyselect_poly_header_s polyselect_poly_header_t[1];
typedef struct polyselect_poly_header_s * polyselect_poly_header_ptr;
typedef const struct polyselect_poly_header_s * polyselect_poly_header_srcptr;

/* structure for series of data */
struct polyselect_data_s {
  unsigned long size;  /* number of values */
  unsigned long alloc; /* allocated size */
  double *x;           /* values */
  double sum;          /* sum = x[0] + ... + x[size-1] */
  double var;          /* var = x[0]^2 + ... + x[size-1]^2 */
  double min, max;     /* minimum and maximum values */
};
typedef struct polyselect_data_s polyselect_data_t[1];
typedef struct polyselect_data_s * polyselect_data_ptr;
typedef const struct polyselect_data_s * polyselect_data_srcptr;

/* inline functions */

/* declarations */

extern const unsigned int SPECIAL_Q[];

#ifdef __cplusplus
extern "C" {
#endif

size_t expected_memory_usage_for_primes(unsigned long P);
unsigned long initPrimes (unsigned long, uint32_t**);
void printPrimes (uint32_t*, unsigned long);
void clearPrimes (uint32_t**);

void polyselect_poly_header_init (polyselect_poly_header_ptr, mpz_ptr, unsigned long, mpz_ptr);
void polyselect_poly_header_clear (polyselect_poly_header_ptr);
int polyselect_poly_header_skip (polyselect_poly_header_srcptr, unsigned long);

void polyselect_proots_init (polyselect_proots_ptr, unsigned long);
void polyselect_proots_add (polyselect_proots_ptr, unsigned long, uint64_t*, unsigned long);
void polyselect_proots_print (polyselect_proots_srcptr, unsigned long);
void polyselect_proots_clear (polyselect_proots_ptr, unsigned long);

void polyselect_qroots_init (polyselect_qroots_ptr);
void polyselect_qroots_realloc (polyselect_qroots_ptr, unsigned long);
void polyselect_qroots_add (polyselect_qroots_ptr, unsigned int, unsigned int, uint64_t*);
void polyselect_qroots_print (polyselect_qroots_srcptr);
void polyselect_qroots_rearrange (polyselect_qroots_ptr R);
void polyselect_qroots_clear (polyselect_qroots_ptr);

void polyselect_hash_init (polyselect_hash_ptr, unsigned int);

/* FIXME: what's the difference, really ? */
void polyselect_hash_add (polyselect_hash_ptr, unsigned long, int64_t, mpz_srcptr, mpz_srcptr, unsigned long, mpz_srcptr, unsigned long, mpz_srcptr);
void polyselect_hash_add_gmp (polyselect_hash_ptr, uint32_t, int64_t, mpz_srcptr, mpz_srcptr, unsigned long, mpz_srcptr, uint64_t, mpz_srcptr);

void polyselect_hash_grow (polyselect_hash_ptr);
void polyselect_hash_clear (polyselect_hash_ptr);

void polyselect_data_init (polyselect_data_ptr);
void polyselect_data_clear (polyselect_data_ptr);
void polyselect_data_add (polyselect_data_ptr, double);
double polyselect_data_mean (polyselect_data_srcptr);
double polyselect_data_var (polyselect_data_srcptr);


#ifdef __cplusplus
}
#endif

#endif
