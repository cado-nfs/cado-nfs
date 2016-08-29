#ifndef ROPT_PARAM_H
#define ROPT_PARAM_H


#include "ropt_str.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <pthread.h>
#include "portability.h"


/* -------------------------
   perhaps no need to change
   ------------------------- */


/* Stage 2 memory cutoff: at most, use 2^28 of int16_t
   which is about 256M memory. Large sieving region will
   be segmented into smaller chunks. */
#define SIZE_SIEVEARRAY 268435456

/* Max default sieve length of v */
#define SIZE_SIEVEARRAY_V_MAX 4194304

/* Min default sieve length of v */
#define SIZE_SIEVEARRAY_V_MIN 1000000

/* Max default sieve length of u */
#define SIZE_SIEVEARRAY_U 64

#define PI 3.14159265358979

#define SUP_ALPHA 4.843

#define MAX_LINE_LENGTH 4096

/* number of default parameters for various sublattices */
#define NUM_DEFAULT_SUBLATTICE 28

/* number of rows for size_total_sublattices */
#define NUM_DEFAULT_DIGITS 20

#define DEBUG 0

/* mu and sigma in function exp_alpha() */
#define MU 0.0

#define SIGMA 0.824

/* Don't change this. */
#define NUM_SUBLATTICE_PRIMES 9

/* Either rank stage 1 sublattices by E or by alpha, the former
   seems to be more accurate */
#define RANK_SUBLATTICE_BY_E 1

/* Toggle for tuning lognorm, it seems better on */
#define TUNE_LOGNORM_INCR 0

/* Top 32 (alpha) in each "SIZE_SIEVEARRAY": this does not
   affect the running-time since it is not dominating */
#define NUM_TOPALPHA_SIEVEARRAY 32

/* Top 32 (E) for each sublattice, may contain several
   SIZE_SIEVEARRAY */
#define NUM_TOPE_SUBLATTICE 32

/* Similar to above, but used in tune mode for speed */
#define TUNE_NUM_TOPALPHA_SIEVEARRAY 8
#define TUNE_NUM_TOPE_SUBLATTICE 8

/* Tuning parameter */
#define TUNE_EARLY_ABORT 1

#define TUNE_EARLY_ABORT_THR 2

#define TUNE_BOUND_ON_UV_TRIALS 64

#define TUNE_BOUND_ON_MOD_TRIALS 64

/* ratio of sublattices in Stage 1 ranking based on partial alpha values */
#define TUNE_RATIO_STAGE1_PART_ALPHA 4

/* cutoff from above by ranking based on full alpha values */
#define TUNE_RATIO_STAGE1_FULL_ALPHA 1

/* ratio between amount of fast tune and slow tune (a 'slow' tune is about 3-4
   times slower than a 'fast' tune). It could be set to be 4. However, setting
   a larger number will reduce the fineness of approximation. */
#define TUNE_RATIO_STAGE2_FAST_SLOW 2


/* -------------------------
   Parameters you may change
   ------------------------- */

/* maximum lognorm+exp_E increment for each rotation */
#define BOUND_LOGNORM_INCR_MAX 1.005

/* only used if TUNE_LOGNORM_INCR is 1 */
#define BOUND_LOGNORM_INCR_MAX_TUNESTEP 0.005


/* --- declarations --- */


extern unsigned int L1_cachesize;

/* set once in ropt_io.c and unchanged then */
extern unsigned int size_tune_sievearray;

extern const unsigned int primes[];

extern const unsigned char next_prime_idx[];

extern const unsigned int default_sublattice_pe[NUM_DEFAULT_SUBLATTICE][NUM_SUBLATTICE_PRIMES];

extern const unsigned long default_sublattice_prod[NUM_DEFAULT_SUBLATTICE];

extern const unsigned int s1_size_each_sublattice[NUM_SUBLATTICE_PRIMES][NUM_SUBLATTICE_PRIMES];

extern const unsigned int s1_size_each_sublattice_tune[NUM_SUBLATTICE_PRIMES];

extern const unsigned int size_total_sublattices[NUM_DEFAULT_DIGITS][4];

double exp_alpha (double logK);

extern unsigned int nthreads_sieve;

/* used as mutual exclusion lock for multithread of sieving */
extern pthread_mutex_t locksieve;

#endif
