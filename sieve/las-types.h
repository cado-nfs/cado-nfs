#ifndef LAS_TYPES_H_
#define LAS_TYPES_H_

#include <stdint.h>
#include "fb.h"
#include "trialdiv.h"
#include "las-config.h"
#include "cado_poly.h"
#include "ecm/facul.h"

struct sieve_info_s;
typedef struct sieve_info_s * sieve_info_ptr;
typedef const struct sieve_info_s * sieve_info_srcptr;

struct where_am_I_s;
typedef struct where_am_I_s * where_am_I_ptr;
typedef const struct where_am_I_s * where_am_I_srcptr;

#include "las-unsieve.h"
#include "las-smallsieve.h"

typedef struct {
    factorbase_degn_t * start;
    factorbase_degn_t * end;
} fb_interval;

struct sieve_side_info_s {
    unsigned char Bound[256]; /* zero for good lognorms, 127 otherwise */
    fbprime_t *trialdiv_primes;
    trialdiv_divisor_t *trialdiv_data;
    unsigned char lognorm_table[1 << NORM_BITS];
    factorbase_degn_t * fb;
    struct {
        factorbase_degn_t * pow2[2];
        factorbase_degn_t * pow3[2];
        factorbase_degn_t * td[2];
        factorbase_degn_t * rs[2];
        factorbase_degn_t * rest[2];
    } fb_parts[1];
    struct {
        int pow2[2];
        int pow3[2];
        int td[2];
        int rs[2];
        int rest[2];
    } fb_parts_x[1];


    /* These fields are used for the norm initialization essentially.
     * Only the scale is also relevant to part of the rest, since it
     * determines the logp contributions for factor base primes */
    double scale;      /* norm scale used on the algebraic side */
    double logmax;     /* norms on the alg-> side are < 2^alg->logmax */

    mpz_t *fij;       /* coefficients of F(a0*i+a1*j, b0*i+b1*j)  */
    double *fijd;     /* coefficients of F_q (divided by q on the special q side) */

    /* This updated by applying the special-q lattice transform to the
     * factor base. */
    small_sieve_data_t ssd[1];
    /* And this is just created as an extraction of the above */
    small_sieve_data_t rsd[1];
};

typedef struct sieve_side_info_s * sieve_side_info_ptr;
typedef const struct sieve_side_info_s * sieve_side_info_srcptr;
typedef struct sieve_side_info_s sieve_side_info[1];

// General information about the siever
struct sieve_info_s {
    cado_poly cpoly;

    // general operational flags
    int nb_threads;
    FILE *output;
    const char * outputname; /* keep track of whether it's gzipped or not */
    int verbose;
    int bench;
    int ratq;   // 0 means special q on alg side, otherwise on rat side

    // sieving area
    uint32_t I;
    uint32_t J;
    int logI; // such that I = 1<<logI

    // description of the q-lattice
    mpz_t q;
    mpz_t rho;
    int32_t a0, b0, a1, b1;

    // parameters for bucket sieving
    unsigned int td_thresh;
    int bucket_thresh;    // bucket sieve primes >= bucket_thresh
    int nb_buckets;
    double bucket_limit_multiplier;
    // unsigned int degree;   /* polynomial degree */
    sieve_side_info sides[2];
    double B;         /* bound for the norm computation */

    unsieve_aux_data us;

    facul_strategy_t *strategy;

    unsigned int degree;

    /* Only for the descent. Excerpts from the Euclidean algorithm
     */
    /* U[0]: middle remainders */
    /* U[1]: next remainders */
    mpz_t euclid_U[2][2];
};

typedef struct sieve_info_s sieve_info[1];

/* FIXME: This does not seem to work well */
#ifdef  __GNUC__
#define TYPE_MAYBE_UNUSED     __attribute__((unused));
#else
#define TYPE_MAYBE_UNUSED       /**/
#endif

struct where_am_I_s {
#ifdef TRACK_CODE_PATH
    fbprime_t p;        /* current prime or prime power, when applicable */
    fbprime_t r;        /* current root */
    int fb_idx;         /* index into the factor base si->sides[side]->fb
                           or into th->sides[side]->fb_bucket */
    unsigned int j;     /* row number in bucket */
    unsigned int x;     /* value in bucket */
    unsigned int N;     /* bucket number */
    int side;
    sieve_info_srcptr si;
#endif  /* TRACK_CODE_PATH */
} TYPE_MAYBE_UNUSED;
typedef struct where_am_I_s where_am_I[1];

#ifdef TRACK_CODE_PATH
#define WHERE_AM_I_UPDATE(w, field, value) (w)->field = (value)
#else
#define WHERE_AM_I_UPDATE(w, field, value) /**/
#endif

#endif	/* LAS_TYPES_H_ */
