#ifndef LAS_TYPES_H_
#define LAS_TYPES_H_

#include <stdint.h>
#include "fb.h"
#include "trialdiv.h"
#include "las-config.h"
#include "las-report-stats.h"
#include "bucket.h"
#include "cado_poly.h"
#include "ecm/facul.h"
#include "relation.h"
#include "las-forwardtypes.h"
#include "las-unsieve.h"
#include "las-qlattice.h"
#include "las-todo.h"
#include "las-smallsieve.h"
#include "las-dlog-base.h"
#include "ecm/batch.h"

/* {{{ siever_config */
/* The following structure lists the fields with an impact on the siever.
 * Different values for these fields will correspond to different siever
 * structures.
 */
struct siever_config_s {
    /* The bit size of the special-q. Counting in bits is no necessity,
     * we could imagine being more accurate */
    unsigned int bitsize;  /* bitsize == 0 indicates end of table */
    int side;
    int logI;
    unsigned long bucket_thresh;    // bucket sieve primes >= bucket_thresh
    unsigned long bucket_thresh1;   // primes above are 2-level bucket-sieved
    unsigned int td_thresh;
    unsigned int unsieve_thresh;
    struct {
        unsigned long lim; /* factor base bound */
        unsigned long powlim; /* bound on powers in the factor base */
        int lpb;           /* large prime bound is 2^lpbr */
        int mfb;           /* bound for residuals is 2^mfbr */
        int ncurves;       /* number of cofactorization curves */
        double lambda;     /* lambda sieve parameter */
    } sides[2][1];
};
typedef struct siever_config_s siever_config[1];
/* siever_config_*ptr defined in las-forwardtypes.h */

/* }}} */

/* This one wants to have siever_config defined */
#include "las-descent-trees.h"

/* {{{ descent_hint
 *
 * This is used for the descent. For each factor size, we provide a
 * reasonable siever_config value
 *
 * We also provide, based on experience, info relative to how long it
 * takes to finish the smoothing process for a prime factor of this size.
 */
struct descent_hint_s {
    siever_config conf;
    double expected_time;
    double expected_success;
};
typedef struct descent_hint_s descent_hint[1];
typedef struct descent_hint_s * descent_hint_ptr;
typedef const struct descent_hint_s * descent_hint_srcptr;

/* }}} */

struct root_s {
  unsigned char derivate; /* 0 = root of F; 1 = root of F'; and so on */
  double value;           /* root of F(i,1) or derivates of F. */
};
typedef struct root_s *root_ptr;
typedef struct root_s root_t[1];

/* These segments ((x, F(x)), (y, F(y))) are used in the smart normalization */
typedef struct sg_s {
  int begin, end;
  double f_begin, f_end;
} sg_t;

/* {{{ sieve_side_info */
struct sieve_side_info_s {
    unsigned char bound; /* A sieve array entry is a sieve survivor if it is
                            at most "bound" on each side */
    trialdiv_divisor_t *trialdiv_data;
    std::vector<fb_general_entry> *fb_smallsieved;
    struct {
        int pow2[2];
        int pow3[2];
        int td[2];
        int rs[2];
        int rest[2];
    } fb_parts_x[1];
    
    /* The reading, mapping or generating the factor base all create the
     * factor base in several pieces: small primes, and large primes split
     * into one piece for each thread.
     */
    fb_factorbase * fb;
    /* When threads pick up this sieve_info structure, they should check
     * their bucket allocation */
    double max_bucket_fill_ratio[FB_MAX_PARTS];


    /* These fields are used for the norm initialization essentially.
     * Only the scale is also relevant to part of the rest, since it
     * determines the logp contributions for factor base primes */
    double scale;      /* scale used for logarithms for fb and norm.
                          must be of form (int)x * 0.1 */
    char check_size[256]; /* non-zero for sizes which are ok */
    double cexp2[257]; /* for 2^X * scale + GUARD */
    double logmax;     /* norms are < 2^logmax */

    mpz_poly_t fij;   /* coefficients of F(a0*i+a1*j, b0*i+b1*j) (divided by 
                         q on the special-q side) */
    double *fijd;     /* coefficients of F_q (divided by q on the special q side) */
    unsigned int nroots; /* Number (+1) and values (+0.0) of the roots of */
    root_ptr roots;     /* F, F', F" and maybe F'" - cf las-norms.c,init_norms */

    /* This updated by applying the special-q lattice transform to the
     * factor base. */
    small_sieve_data_t ssd[1];
    /* And this is just created as an extraction of the above */
    small_sieve_data_t rsd[1];

    facul_strategy_t *strategy;
};
typedef struct sieve_side_info_s * sieve_side_info_ptr;
typedef const struct sieve_side_info_s * sieve_side_info_srcptr;
typedef struct sieve_side_info_s sieve_side_info[1];
/* }}} */

/* {{{ sieve_info
 *
 * General information about the siever, based on some input-dependent
 * configuration data with an impact on the output (as opposed to e.g.
 * file names, or verbosity flags, which do not affect the output).
 */
struct sieve_info_s {
    cado_poly_ptr cpoly; /* The polynomial pair */

    /* This conditions the validity of the sieve_info_side members
     * sides[0,1], as well as some other members */
    siever_config conf;

    las_todo_entry *doing;

    // sieving area
    uint32_t J;
    uint32_t I; // *conf has logI.

    // description of the q-lattice. The values here should remain
    // compatible with those in ->conf (this concerns notably the bit
    // size as well as the special-q side).
    qlattice_basis qbasis;

    // parameters for bucket sieving
    /* Actual number of buckets at toplevel used by current special-q */
    uint32_t nb_buckets[FB_MAX_PARTS];
    /* Max number of buckets if J=I/2, at each level */
    uint32_t nb_buckets_max[FB_MAX_PARTS];
    /* Largest level for which the corresponding fb_part is not empty */
    int toplevel;

    sieve_side_info sides[2];

    facul_strategies_t* strategies;
  
    /* Data for unsieving locations where gcd(i,j) > 1 */
    unsieve_aux_data_srcptr us;
    /* Data for divisibility tests p|i in lines where p|j */
    j_div_srcptr j_div;

    /* used in check_leftover_norm */
    mpz_t BB[2], BBB[2], BBBB[2];
};
typedef struct sieve_info_s sieve_info[1];

/* }}} */

/* {{{ las_info
 *
 * las_info holds general data, mostly unrelated to what is actually
 * computed within a sieve. las_info also contains outer data, which
 * lives outside the choice of one particular way to configure the siever
 * versus another.
 */
struct las_info_s {
    // general operational flags
    int nb_threads;
    FILE *output;
    const char * outputname; /* keep track of whether it's gzipped or not */
    const char * galois; /* a string to indicate which galois to use in las */
    int verbose;
    int suppress_duplicates;

    /* It's not ``general operational'', but global enough to be here */
    cado_poly cpoly;

    siever_config_srcptr default_config;

    /* This needs not be complete. The default_config field points here
     * if is. If not, the fields here are just used as a base for
     * initializing the other configurations */
    siever_config config_base;

    /* There may be several configured sievers. This is used mostly for
     * the descent.
     * TODO: For now these different sievers share nothing of their
     * factor bases, which is a shame. We should examine ways to get
     * around this limitation, but it is hard. One could imagine some,
     * though. Among them, given the fact that only _one_ siever is
     * active at a time, it might be possible to sort the factor base
     * again each time a new siever is used (but maybe it's too
     * expensive). Another way could be to work only on sharing
     * bucket-sieved primes.
     */
    sieve_info_ptr sievers;

    /* Descent lower steps */
    descent_hint * hint_table;
    unsigned int max_hint_bitsize[2];
    int * hint_lookups[2]; /* quick access indices into hint_table */

    /* This is an opaque pointer to C++ code. */
    void * descent_helper;

    las_todo_stack *todo;
    unsigned int nq_pushed;
    unsigned int nq_max;

    int random_sampling;
    gmp_randstate_t rstate;

    /* These are used for reading the todo list */
    mpz_t todo_q0;
    mpz_t todo_q1;
    FILE * todo_list_fd;

#ifdef  DLP_DESCENT
    las_dlog_base * dlog_base;
#endif
    descent_tree * tree;

    int batch; /* batch mode for cofactorization */
    cofac_list L; /* store (a,b) and corresponding cofactors in batch mode */
};

typedef struct las_info_s las_info[1];
/* }}} */

enum {
  OUTPUT_CHANNEL,
  ERROR_CHANNEL,
  STATS_CHANNEL,
  TRACE_CHANNEL,
  NR_CHANNELS /* This must be the last element of the enum */
};

#endif	/* LAS_TYPES_H_ */
