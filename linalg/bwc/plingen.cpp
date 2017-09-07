/* Copyright (C) 1999--2007 Emmanuel Thom'e --- see LICENSE file */
#include "cado.h"

#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/utsname.h>


#include <assert.h>

#include "portability.h"
#include "macros.h"
#include "utils.h"
#include "mpfq_layer.h"
#include "memusage.h"

/* lingen-matpoly is the default code. */
#include "lingen-matpoly.h"

#define ENABLE_MPI_LINGEN
#ifdef  HAVE_MPIR
#define ENABLE_CACHING_MPI_LINGEN
#endif

#ifdef  ENABLE_MPI_LINGEN
#include "lingen-bigmatpoly.h"
#endif

#ifdef  ENABLE_CACHING_MPI_LINGEN
#ifdef HAVE_MPIR
#include "lingen-bigmatpoly-ft.h"
#else
#error "ENABLE_CACHING_MPI_LINGEN requires MPIR"
#endif
#endif

#include "bw-common.h"		/* Handy. Allows Using global functions
                                 * for recovering parameters */
#include "plingen.h"
#include "plingen-tuning.h"
#include "logline.h"
#include "tree_stats.hpp"

/* Call tree for methods within this program:
 *
 * Two separate entry points.
 *
 * bw_biglingen_collective        when need to go collective
 *     |<--> bw_biglingen_recursive , loops to bw_biglingen_collective
 *     \->bw_lingen_single               when it makes sense to do this locally again
 * bw_lingen_single                      when computation can be done locally
 *    |<->bw_lingen_recursive
 *    |->bw_lingen_basecase
 *
 */

#define MPI_MY_SIZE_T   MPI_UNSIGNED_LONG

static unsigned int display_threshold = 10;
static int with_timings = 0;

/* This is an indication of the number of bytes we read at a time for A
 * (input) and F (output) */
static unsigned int io_block_size = 1 << 20;

/* If non-zero, then reading from A is actually replaced by reading from
 * a random generator */
static unsigned int random_input_length = 0;

static int split_input_file = 0;  /* unsupported ; do acollect by ourselves */
static int split_output_file = 0; /* do split by ourselves */

gmp_randstate_t rstate;

#ifdef HAVE_MPIR
static int caching = 1;
static unsigned int caching_threshold = 10;
#else
static int caching = 0;
#endif

static const char * checkpoint_directory;
static unsigned int checkpoint_threshold = 100;
static int save_gathered_checkpoints = 0;

static int allow_zero_on_rhs = 0;

int rank0_exit_code = EXIT_SUCCESS;


int global_flag_ascii = 0;
int global_flag_tune = 0;

struct bmstatus_s {
    dims d[1];
    unsigned int t;
    int * lucky;

    double t_basecase;
    double t_mp;
    double t_mul;
    double t_cp_io;

    unsigned int lingen_threshold;
    unsigned int lingen_mpi_threshold;
    int mpi_dims[2]; /* mpi_dims[0] = mpi[0] * thr[0] */
    MPI_Comm com[3]; /* [0]: MPI_COMM_WORLD, reordered.
                        [1]: row-wise
                        [2]: column-wise */
};
typedef struct bmstatus_s bmstatus[1];
typedef struct bmstatus_s *bmstatus_ptr;


tree_stats stats;

void plingen_decl_usage(param_list_ptr pl)
{
    param_list_decl_usage(pl, "ascii",
            "read and write data in ascii");
    param_list_decl_usage(pl, "timings",
            "provide timings on all output lines");
    param_list_decl_usage(pl, "tune",
            "activate tuning mode");
    param_list_decl_usage(pl, "allow_zero_on_rhs",
            "do not cry if the generator corresponds to a zero contribution on the RHS vectors");

    /* we must be square ! And thr is not supported. */
    param_list_decl_usage(pl, "mpi", "number of MPI nodes across which the execution will span, with mesh dimensions");
    param_list_decl_usage(pl, "thr", "number of threads (on each node) for the program, with mesh dimensions");

    param_list_decl_usage(pl, "nrhs",
            "number of columns to treat differently, as corresponding to rhs vectors");
    param_list_decl_usage(pl, "rhs",
            "file with rhs vectors (only the header is read)");

    param_list_decl_usage(pl, "afile",
            "input sequence file");
    param_list_decl_usage(pl, "random-input-with-length",
            "use surrogate for input");
    param_list_decl_usage(pl, "split-input-file",
            "work with split files on input");
    param_list_decl_usage(pl, "split-output-file",
            "work with split files on output");
    param_list_decl_usage(pl, "random_seed",
            "seed the random generator");
    param_list_decl_usage(pl, "ffile",
            "output generator file");

    param_list_decl_usage(pl, "caching",
            "whether we should use transform caching");
    param_list_decl_usage(pl, "caching-threshold",
            "threshold for transform caching");
    param_list_decl_usage(pl, "checkpoint-directory",
            "where to save checkpoints");
    param_list_decl_usage(pl, "checkpoint-threshold",
            "threshold for saving checkpoints");
    param_list_decl_usage(pl, "display-threshold",
            "threshold for outputting progress lines");
    param_list_decl_usage(pl, "io-block-size",
            "chunk size for reading the input or writing the output");

    param_list_decl_usage(pl, "lingen_mpi_threshold",
            "use MPI matrix operations above this size");
    param_list_decl_usage(pl, "lingen_threshold",
            "use recursive algorithm above this size");
    param_list_decl_usage(pl, "save_gathered_checkpoints",
            "save global checkpoints files, instead of per-job files");

    param_list_configure_switch(pl, "--tune", &global_flag_tune);
    param_list_configure_switch(pl, "--ascii", &global_flag_ascii);
    param_list_configure_switch(pl, "--timings", &with_timings);
    param_list_configure_alias(pl, "seed", "random_seed");

    plingen_tuning_decl_usage(pl);
}

/*{{{ basecase */

/* {{{ col sorting */
/* We sort only with respect to the global delta[] parameter. As it turns
 * out, we also access the column index in the same aray and sort with
 * respect to it, but this is only cosmetic.
 *
 * Note that unlike what is asserted in other coiped of the code, sorting
 * w.r.t. the local delta[] value is completely useless. Code which
 * relies on this should be fixed.
 */

typedef int (*sortfunc_t) (const void*, const void*);

static int lexcmp2(const int x[2], const int y[2])
{
    for(int i = 0 ; i < 2 ; i++) {
        int d = x[i] - y[i];
        if (d) return d;
    }
    return 0;
}

/* }}} */

static inline unsigned int expected_pi_length(dims * d, unsigned int len)/*{{{*/
{
    /* The idea is that we want something which may account for something
     * exceptional, bounded by probability 2^-64. This corresponds to a
     * column in e (matrix of size m*b) to be spontaneously equal to
     * zero. This happens with probability (#K)^-m.
     * The solution to
     * (#K)^(-m*x) > 2^-64
     * is m*x*log_2(#K) < 64
     *
     * We thus need to get an idea of the value of log_2(#K).
     *
     * (Note that we know that #K^abgroupsize(ab) < 2^64, but that bound
     * might be very gross).
     *
     * The same esitmate can be used to appreciate what we mean by
     * ``luck'' in the end. If a column happens to be zero more than
     * expected_pi_length(d,0) times in a row, then the cause must be
     * more than sheer luck, and we use it to detect generating rows.
     */
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    abdst_field ab MAYBE_UNUSED = d->ab;
    unsigned int res = iceildiv(len * m, b);
    mpz_t p;
    mpz_init(p);
    abfield_characteristic(ab, p);
    unsigned int l;
    if (mpz_cmp_ui(p, 1024) >= 0) {
        l = mpz_sizeinbase(p, 2);
        l *= abfield_degree(ab);    /* roughly log_2(#K) */
    } else {
        mpz_pow_ui(p, p, abfield_degree(ab));
        l = mpz_sizeinbase(p, 2);
    }
    mpz_clear(p);
    // unsigned int safety = iceildiv(abgroupsize(ab), m * sizeof(abelt));
    unsigned int safety = iceildiv(64, m * l);
    return res + safety;
}/*}}}*/

/* This destructively cancels the first len coefficients of E, and
 * computes the appropriate matrix pi which achieves this. The
 * elimination is done in accordance with the nominal degrees found in
 * delta.
 *
 * The result is expected to have degree ceil(len*m/b) coefficients, so
 * that E*pi is divisible by X^len.
 */

/* TODO: adapt for GF(2) */
static int bw_lingen_basecase(bmstatus_ptr bm, matpoly pi, matpoly E, unsigned int *delta) /*{{{*/
{
    stats.enter(__func__, E->size);
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    abdst_field ab = d->ab;
    ASSERT(E->m == m);
    ASSERT(E->n == b);

    ASSERT(pi->m == 0);
    ASSERT(pi->n == 0);
    ASSERT(pi->alloc == 0);

    /* Allocate something large enough for the result. This will be
     * soon freed anyway. Set it to identity. */
    unsigned int pi_room_base = expected_pi_length(d, E->size);

    matpoly_init(ab, pi, b, b, pi_room_base);
    pi->size = pi_room_base;

    /* Also keep track of the
     * number of coefficients for the columns of pi. Set pi to Id */

    unsigned int *pi_lengths = (unsigned int*) malloc(b * sizeof(unsigned int));
    for(unsigned int i = 0 ; i < b ; i++) {
        abset_ui(ab, matpoly_coeff(ab, pi, i, i, 0), 1);
        pi_lengths[i] = 1;
    }

    /* This is used below */
    int generator_found = 0;

    /* Keep a list of columns which have been used as pivots at the
     * previous iteration */
    unsigned int * pivots = (unsigned int*) malloc(m * sizeof(unsigned int));
    int * is_pivot = (int*) malloc(b * sizeof(int));
    memset(is_pivot, 0, b * sizeof(int));

    matpoly e;
    matpoly_init(ab, e, m, b, 1);
    e->size = 1;

    for (unsigned int t = 0; t < E->size ; t++, bm->t++) {

        /* {{{ Update the columns of e for degree t. Save computation
         * time by not recomputing those which can easily be derived from
         * previous iteration. Notice that the columns of e are exactly
         * at the physical positions of the corresponding columns of pi.
         */

        abelt_ur tmp_ur;
        abelt_ur_init(ab, &tmp_ur);

        abvec_ur e_ur;
        abvec_ur_init(ab, &e_ur, m*b);
        abvec_ur_set_zero(ab, e_ur, m*b);
        for(unsigned int j = 0 ; j < b ; j++) {
            if (is_pivot[j]) continue;
            /* We should never have to recompute from pi using discarded
             * columns. Discarded columns should always correspond to
             * pivots */
            ASSERT_ALWAYS(bm->lucky[j] >= 0);
            for(unsigned int s = 0 ; s < pi_lengths[j] ; s++) {
                for(unsigned int i = 0 ; i < m ; i++) {
                    for(unsigned int k = 0 ; k < b ; k++) {
                        abmul_ur(ab, tmp_ur,
                                matpoly_coeff(ab, E, i, k, t - s),
                                matpoly_coeff(ab, pi, k, j, s));
                        abelt_ur_add(ab,
                                abvec_ur_coeff_ptr(ab, e_ur, i*b+j),
                                abvec_ur_coeff_ptr_const(ab, e_ur, i*b+j),
                                tmp_ur);
                    }
                }
            }
        }
        unsigned int newluck = 0;
        for(unsigned int j = 0 ; j < b ; j++) {
            if (is_pivot[j]) continue;
            unsigned int nz = 0;
            for(unsigned int i = 0 ; i < m ; i++) {
                abreduce(ab, matpoly_coeff(ab, e, i, j, 0),
                        abvec_ur_coeff_ptr(ab, e_ur, i*b+j));
                nz += abcmp_ui(ab, matpoly_coeff(ab, e, i, j, 0), 0) == 0;
            }
            if (nz == m) {
                newluck++, bm->lucky[j]++;
            } else if (bm->lucky[j] > 0) {
                bm->lucky[j] = 0;
            }
        }
        abelt_ur_clear(ab, &tmp_ur);
        abvec_ur_clear(ab, &e_ur, m*b);
        if (newluck) {
            /* If newluck == n, then we probably have a generator. We add an
             * extra guarantee. newluck==n, for a total of k iterations in a
             * row, means m*n*k coefficients cancelling magically. We would
             * like this to be impossible by mere chance. Thus we want n*k >
             * luck_mini, which can easily be checked */

            int luck_mini = expected_pi_length(d, 0);
            unsigned int luck_sure = 0;

            printf("t=%d, canceled columns:", bm->t);
            for(unsigned int j = 0 ; j < b ; j++) {
                if (bm->lucky[j] > 0) {
                    printf(" %u", j);
                    luck_sure += bm->lucky[j] > luck_mini;
                }
            }

            if (newluck == n && luck_sure == n) {
                if (!generator_found) {
                    printf(", complete generator found, for sure");
                }
                generator_found = 1;
            }
            printf(".\n");
        }
        /* }}} */

        if (generator_found) break;

        int (*ctable)[2] = (int(*)[2]) malloc(b * 2 * sizeof(int));
        /* {{{ Now see in which order I may look at the columns of pi, so
         * as to keep the nominal degrees correct. In contrast with what
         * we used to do before, we no longer apply the permutation to
         * delta. So the delta[] array keeps referring to physical
         * indices, and we'll tune this in the end. */
        for(unsigned int j = 0; j < b; j++) {
            ctable[j][0] = delta[j];
            ctable[j][1] = j;
        }
        qsort(ctable, b, 2 * sizeof(int), (sortfunc_t) & lexcmp2);
        /* }}} */

        /* {{{ Now do Gaussian elimination */
        memset(is_pivot, 0, b * sizeof(int));
        unsigned int r = 0;
        /* Loop through logical indices */
        for(unsigned int jl = 0; jl < b; jl++) {
            unsigned int j = ctable[jl][1];
            unsigned int u = 0;
            /* {{{ Find the pivot */
            for( ; u < m ; u++) {
                if (abcmp_ui(ab, matpoly_coeff(ab, e, u, j, 0), 0) != 0)
                    break;
            }
            if (u == m) continue;
            assert(r < m);
            /* }}} */
            pivots[r++] = j;
            is_pivot[j] = 1;
            /* {{{ Cancel this coeff in all other columns. */
            abelt inv;
            abinit(ab, &inv);
            int rc = abinv(ab, inv, matpoly_coeff(ab, e, u, j, 0));
            if (!rc) {
                fprintf(stderr, "Error, found a factor of the modulus: ");
                abfprint(ab, stderr, inv);
                fprintf(stderr, "\n");
                exit(EXIT_FAILURE);
            }
            abneg(ab, inv, inv);
            for (unsigned int kl = jl + 1; kl < b ; kl++) {
                unsigned int k = ctable[kl][1];
                if (abcmp_ui(ab, matpoly_coeff(ab, e, u, k, 0), 0) == 0)
                    continue;
                // add lambda = e[u,k]*-e[u,j]^-1 times col j to col k.
                abelt lambda;
                abinit(ab, &lambda);
                abmul(ab, lambda, inv, matpoly_coeff(ab, e, u, k, 0));

                assert(delta[j] <= delta[k]);
                /* {{{ Apply on both e and pi */
                abelt tmp;
                abinit(ab, &tmp);
                for(unsigned int i = 0 ; i < m ; i++) {
                    /* TODO: Would be better if mpfq had an addmul */
                    abmul(ab, tmp, lambda, matpoly_coeff(ab, e, i, j, 0));
                    abadd(ab,
                            matpoly_coeff(ab, e, i, k, 0),
                            matpoly_coeff(ab, e, i, k, 0),
                            tmp);
                }
                if (bm->lucky[k] < 0) {
                    /* This column is already discarded, don't bother */
                    continue;
                }
                if (bm->lucky[j] < 0) {
                    /* This column is discarded. This is going to
                     * invalidate another column of pi. Not a problem,
                     * unless it's been marked as lucky previously ! */
                    ASSERT_ALWAYS(bm->lucky[k] <= 0);
                    printf("Column %u discarded from now on (through addition from column %u)\n", k, j);
                    bm->lucky[k] = -1;
                    continue;
                }
                for(unsigned int i = 0 ; i < b ; i++) {
                    /* Beware. One may be tempted to think that the code
                     * above is dubious. It is, in fact, not a problem.
                     * As long as the delta[] array undergoes no
                     * disturbing modification, everything is ok.
                     */
                    if (pi_lengths[k] < pi_lengths[j])
                        pi_lengths[k] = pi_lengths[j];
                    for(unsigned int s = 0 ; s < pi_lengths[j] ; s++) {
                        /* TODO: Would be better if mpfq had an addmul */
                        abmul(ab, tmp, lambda,
                                matpoly_coeff(ab, pi, i, j, s));
                        abadd(ab,
                                matpoly_coeff(ab, pi, i, k, s),
                                matpoly_coeff(ab, pi, i, k, s),
                                tmp);
                    }
                }
                abclear(ab, &tmp);
                /* }}} */
                abclear(ab, &lambda);
            }
            abclear(ab, &inv); /* }}} */
        }
        /* }}} */
        free(ctable);

        ASSERT_ALWAYS(r == m);

        /* {{{ Now for all pivots, multiply column in pi by x */
        for (unsigned int j = 0; j < b ; j++) {
            if (!is_pivot[j]) continue;
            if (pi_lengths[j] >= pi->alloc) {
                if (!generator_found) {
                    matpoly_realloc(ab, pi, pi->alloc + MAX(pi->alloc / (m+n), 1));
                    printf("t=%u, expanding allocation for pi (now %zu%%) ; lengths: ",
                            bm->t,
                            100 * pi->alloc / pi_room_base);
                    for(unsigned int j = 0; j < b; j++)
                        printf(" %u", pi_lengths[j]);
                    printf("\n");
                } else {
                    ASSERT_ALWAYS(bm->lucky[j] <= 0);
                    if (bm->lucky[j] == 0)
                        printf("t=%u, column %u discarded from now on\n",
                                bm->t, j);
                    bm->lucky[j] = -1;
                    pi_lengths[j]++;
                    delta[j]++;
                    continue;
                }
            }
            matpoly_multiply_column_by_x(ab, pi, j, pi_lengths[j]);
            pi_lengths[j]++;
            delta[j]++;
        }
        /* }}} */
    }
    pi->size = 0;
    for(unsigned int j = 0; j < b; j++) {
        if (pi_lengths[j] > pi->size)
            pi->size = pi_lengths[j];
    }
    pi->size = MIN(pi->size, pi->alloc);
    matpoly_clear(ab, e);
    free(is_pivot);
    free(pivots);
    free(pi_lengths);   /* What shall we do with this one ??? */

    return stats.leave(generator_found);
}/*}}}*/

/*}}}*/

/* TODO: adapt for GF(2) */
double avg_matsize(abdst_field ab, unsigned int m, unsigned int n, int ascii)/*{{{*/
{
    if (!ascii) {
        /* Easy case first. If we have binary input, then we know a priori
         * that the input data must have size a multiple of the element size.
         */
        size_t elemsize = abvec_elt_stride(ab, 1);
        size_t matsize = elemsize * m * n;
        return matsize;
    }

    /* Ascii is more complicated. We're necessarily fragile here.
     * However, assuming that each coefficient comes with only one space,
     * and each matrix with an extra space (this is how the GPU program
     * prints data -- not that this ends up having a considerable impact
     * anyway...), we can guess the number of bytes per matrix. */

    /* Formula for the average number of digits of an integer mod p,
     * written in base b:
     *
     * (k-((b^k-1)/(b-1)-1)/p)  with b = Ceiling(Log(p)/Log(b)).
     */
    double avg;
    mpz_t p, a;
    mpz_init(p);
    mpz_init(a);
    abfield_characteristic(ab, p);
    unsigned long k = ceil(log(mpz_get_d(p))/log(10));
    unsigned long b = 10;
    mpz_ui_pow_ui(a, b, k);
    mpz_sub_ui(a, a, 1);
    mpz_fdiv_q_ui(a, a, b-1);
    avg = k - mpz_get_d(a) / mpz_get_d(p);
    mpz_clear(p);
    mpz_clear(a);
    // printf("Expect roughly %.2f decimal digits for integers mod p.\n", avg);
    double matsize = (avg + 1) * m * n + 1;
    // printf("Expect roughly %.2f bytes for each sequence matrix.\n", matsize);
    return matsize;
}/*}}}*/

/* {{{ I/O helpers */
/* TODO: adapt for GF(2) */
/* {{{ matpoly_write
 * writes some of the matpoly data to f, either in ascii or binary
 * format. This can be used to write only part of the data (degrees
 * [k0..k1[). Returns the number of coefficients (i.e., matrices, so at
 * most k1-k0) successfully written, or
 * -1 on error (e.g. when some matrix was only partially written).
 */
int matpoly_write(abdst_field ab, FILE * f, matpoly_srcptr M, unsigned int k0, unsigned int k1, int ascii, int transpose)
{
    unsigned int m = transpose ? M->n : M->m;
    unsigned int n = transpose ? M->m : M->n;
    ASSERT_ALWAYS(k0 == k1 || (k0 < M->size && k1 <= M->size));
    abelt tmp;
    abinit(ab, &tmp);
    for(unsigned int k = k0 ; k < k1 ; k++) {
        int err = 0;
        int matnb = 0;
        for(unsigned int i = 0 ; !err && i < m ; i++) {
            for(unsigned int j = 0 ; !err && j < n ; j++) {
                absrc_elt x;
                x = transpose ? matpoly_coeff_const(ab, M, j, i, k)
                              : matpoly_coeff_const(ab, M, i, j, k);
                if (ascii) {
                    if (j) err = fprintf(f, " ") <= 0;
                    if (!err) err = abfprint(ab, f, x) <= 0;
                } else {
                    err = fwrite(x, abvec_elt_stride(ab, 1), 1, f) < 1;
                }
                if (!err) matnb++;
            }
            if (!err && ascii) err = fprintf(f, "\n") <= 0;
        }
        if (ascii) err = err || fprintf(f, "\n") <= 0;
        if (err) return (matnb == 0) ? (int) (k - k0) : -1;
    }
    abclear(ab, &tmp);
    return k1 - k0;
}

/* fw must be an array of FILE* pointers of exactly the same size as the
 * matrix to be written.
 */
int matpoly_write_split(abdst_field ab, FILE ** fw, matpoly_srcptr M, unsigned int k0, unsigned int k1, int ascii)
{
    ASSERT_ALWAYS(k0 == k1 || (k0 < M->size && k1 <= M->size));
    abelt tmp;
    abinit(ab, &tmp);
    for(unsigned int k = k0 ; k < k1 ; k++) {
        int err = 0;
        int matnb = 0;
        for(unsigned int i = 0 ; !err && i < M->m ; i++) {
            for(unsigned int j = 0 ; !err && j < M->n ; j++) {
                FILE * f = fw[i*M->n+j];
                absrc_elt x = matpoly_coeff_const(ab, M, i, j, k);
                if (ascii) {
                    err = abfprint(ab, f, x) <= 0;
                    if (!err) err = fprintf(f, "\n") <= 0;
                } else {
                    err = fwrite(x, abvec_elt_stride(ab, 1), 1, f) < 1;
                }
                if (!err) matnb++;
            }
        }
        if (err) return (matnb == 0) ? (int) (k - k0) : -1;
    }
    abclear(ab, &tmp);
    return k1 - k0;
}
/* }}} */

/* {{{ matpoly_read
 * reads some of the matpoly data from f, either in ascii or binary
 * format. This can be used to parse only part of the data (degrees
 * [k0..k1[, k1 being an upper bound). Returns the number of coefficients
 * (i.e., matrices, so at most k1-k0) successfully read, or
 * -1 on error (e.g. when some matrix was only partially read).
 *
 * Note that the matrix must *not* be in pre-init state. It must have
 * been already allocated.
 */

int matpoly_read(abdst_field ab, FILE * f, matpoly_ptr M, unsigned int k0, unsigned int k1, int ascii, int transpose)
{
    ASSERT_ALWAYS(!matpoly_check_pre_init(M));
    unsigned int m = transpose ? M->n : M->m;
    unsigned int n = transpose ? M->m : M->n;
    ASSERT_ALWAYS(k0 == k1 || (k0 < M->size && k1 <= M->size));
    abelt tmp;
    abinit(ab, &tmp);
    for(unsigned int k = k0 ; k < k1 ; k++) {
        int err = 0;
        int matnb = 0;
        for(unsigned int i = 0 ; !err && i < m ; i++) {
            for(unsigned int j = 0 ; !err && j < n ; j++) {
                abdst_elt x;
                x = transpose ? matpoly_coeff(ab, M, j, i, k)
                              : matpoly_coeff(ab, M, i, j, k);
                if (ascii) {
                    err = abfscan(ab, f, x) == 0;
                } else {
                    err = fread(x, abvec_elt_stride(ab, 1), 1, f) < 1;
                }
                if (!err) matnb++;
            }
        }
        if (err) return (matnb == 0) ? (int) (k - k0) : -1;
    }
    abclear(ab, &tmp);
    return k1 - k0;
}/* }}} */

/* }}} */

/*{{{ Checkpoints */

/* There's much copy-paste here */

struct cp_info {
    bmstatus_ptr bm;
    int level;
    unsigned int t0;
    unsigned int t1;
    int mpi;
    int rank;
    char * auxfile;
    char * sdatafile;
    char * gdatafile;
    const char * datafile;
    FILE * aux;
    FILE * data;
};

void cp_info_init_backend(struct cp_info * cp)
{
    int rc;
    cp->level = stats.depth;
    rc = asprintf(&cp->auxfile, "%s/pi.%d.%u.%u.aux",
            checkpoint_directory, cp->level, cp->t0, cp->t1);
    ASSERT_ALWAYS(rc >= 0);
    rc = asprintf(&cp->gdatafile, "%s/pi.%d.%u.%u.single.data",
            checkpoint_directory, cp->level, cp->t0, cp->t1);
    ASSERT_ALWAYS(rc >= 0);
    rc = asprintf(&cp->sdatafile, "%s/pi.%d.%u.%u.%d.data",
            checkpoint_directory, cp->level, cp->t0, cp->t1, cp->rank);
    ASSERT_ALWAYS(rc >= 0);
}

void cp_info_init(struct cp_info * cp, bmstatus_ptr bm, unsigned int t0, unsigned int t1)
{
    memset(cp, 0, sizeof(*cp));
    cp->bm = bm;
    cp->t0 = t0; cp->t1 = t1;
    cp_info_init_backend(cp);
    cp->datafile = cp->gdatafile;
}

void cp_info_init_mpi(struct cp_info * cp, bmstatus_ptr bm, unsigned int t0, unsigned int t1)
{
    memset(cp, 0, sizeof(*cp));
    cp->mpi = 1;
    cp->bm = bm;
    cp->t0 = t0; cp->t1 = t1;
    MPI_Comm_rank(bm->com[0], &(cp->rank));
    cp_info_init_backend(cp);
    /* a priori default. However load_mpi_checkpoint_file looks up both */
    cp->datafile = cp->sdatafile;
}

void cp_info_clear(struct cp_info * cp)
{
    free(cp->sdatafile);
    free(cp->gdatafile);
    free(cp->auxfile);
}

int cp_save_aux_file(struct cp_info * cp, size_t pi_size, unsigned int *delta, int done)
{
    dims * d = cp->bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    if (cp->rank) return 1;
    FILE * aux = fopen(cp->auxfile, "w");
    int rc;
    if (aux == NULL) {
        fprintf(stderr, "Warning: cannot open %s\n", cp->auxfile);
        return 0;
    }
    rc = fprintf(aux, "%zu\n", pi_size);
    if (rc <= 0) goto cp_save_aux_file_bailout;
    for(unsigned int i = 0 ; i < m + n ; i++) {
        rc = fprintf(aux, "%s%u", i?" ":"", delta[i]);
        if (rc <= 0) goto cp_save_aux_file_bailout;
    }
    rc = fprintf(aux, "\n");
    if (rc <= 0) goto cp_save_aux_file_bailout;
    for(unsigned int i = 0 ; i < m + n ; i++) {
        rc = fprintf(aux, "%s%d", i?" ":"", cp->bm->lucky[i]);
        if (rc <= 0) goto cp_save_aux_file_bailout;
    }
    rc = fprintf(aux, "\n");
    if (rc <= 0) goto cp_save_aux_file_bailout;
    rc = fprintf(aux, "%d\n", done);
    if (rc <= 0) goto cp_save_aux_file_bailout;
    rc = fclose(aux);
    if (rc == 0) return 1;
cp_save_aux_file_bailout:
    fclose(aux);
    unlink(cp->auxfile);
    return 0;
}

int cp_load_aux_file(struct cp_info * cp, size_t * p_pi_size, unsigned int *delta, int * p_done)
{
    dims * d = cp->bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    if (cp->rank) return 1;
    FILE * aux = fopen(cp->auxfile, "r");
    int rc;
    if (aux == NULL) {
        // fprintf(stderr, "Warning: cannot open %s\n", cp->auxfile);
        return 0;
    }
    rc = fscanf(aux, "%zu", p_pi_size);
    if (rc != 1) { fclose(aux); return 0; }
    for(unsigned int i = 0 ; i < m + n ; i++) {
        rc = fscanf(aux, "%u", delta + i);
        if (rc != 1) { fclose(aux); return 0; }
    }
    for(unsigned int i = 0 ; i < m + n ; i++) {
        rc = fscanf(aux, "%d", cp->bm->lucky + i);
        if (rc != 1) { fclose(aux); return 0; }
    }
    rc = fscanf(aux, "%d", p_done);
    if (rc != 1) { fclose(aux); return 0; }
    rc = fclose(aux);
    return rc == 0;
}

/* TODO: adapt for GF(2) */
int cp_load_data_file(struct cp_info * cp, matpoly pi, size_t pi_size)
{
    dims * d = cp->bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    FILE * data = fopen(cp->datafile, "rb");
    int rc;
    if (data == NULL) {
        fprintf(stderr, "Warning: cannot open %s\n", cp->datafile);
        return 0;
    }
    matpoly_init(ab, pi, m+n, m+n, pi_size);
    pi->size = pi_size;
    rc = matpoly_read(ab, data, pi, 0, pi->size, 0, 0);
    if (rc != (int) pi->size) { fclose(data); return 0; }
    rc = fclose(data);
    return rc == 0;
}

/* TODO: adapt for GF(2) */
/* I think we always have pi_size == pi->size, the only questionable
 * situation is when we're saving part of a big matrix */
int cp_save_data_file(struct cp_info * cp, matpoly pi, size_t pi_size)
{
    abdst_field ab = cp->bm->d->ab;
    FILE * data = fopen(cp->datafile, "wb");
    int rc;
    if (data == NULL) {
        fprintf(stderr, "Warning: cannot open %s\n", cp->datafile);
        unlink(cp->auxfile);
        return 0;
    }
    rc = matpoly_write(ab, data, pi, 0, pi_size, 0, 0);
    if (rc != (int) pi->size) goto cp_save_data_file_bailout;
    rc = fclose(data);
    if (rc == 0)
        return 1;
cp_save_data_file_bailout:
    fclose(data);
    unlink(cp->datafile);
    unlink(cp->auxfile);
    return 0;
}

int load_checkpoint_file(bmstatus_ptr bm, matpoly pi, unsigned int t0, unsigned int t1, unsigned int *delta, int * p_done)/*{{{*/
{
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    struct cp_info cp[1];
    cp_info_init(cp, bm, t0, t1);
    ASSERT_ALWAYS(matpoly_check_pre_init(pi));
    size_t pi_size;
    /* Don't output a message just now, since after all it's not
     * noteworthy if the checkpoint file does not exist. */
    int ok = cp_load_aux_file(cp, &pi_size, delta, p_done);
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Reading checkpoint %s", cp->datafile);
        ok = cp_load_data_file(cp, pi, pi_size);
        logline_end(&bm->t_cp_io,"");
        if (!ok)
            fprintf(stderr, "Warning: I/O error while reading %s\n", cp->datafile);
    }
    cp_info_clear(cp);
    if (ok) bm->t = t1;
    return ok;
}/*}}}*/

int save_checkpoint_file(bmstatus_ptr bm, matpoly pi, unsigned int t0, unsigned int t1, unsigned int *delta, int done)/*{{{*/
{
    /* corresponding t is bm->t - E->size ! */
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    struct cp_info cp[1];
    cp_info_init(cp, bm, t0, t1);
    logline_begin(stdout, SIZE_MAX, "Saving checkpoint %s%s",
            cp->datafile,
            cp->mpi ? " (MPI, scattered)" : "");
    int ok = cp_save_aux_file(cp, pi->size, delta, done);
    if (ok) ok = cp_save_data_file(cp, pi, pi->size);
    logline_end(&bm->t_cp_io,"");
    if (!ok && !cp->rank)
        fprintf(stderr, "Warning: I/O error while saving %s\n", cp->datafile);
    cp_info_clear(cp);
    return ok;
}/*}}}*/

#ifdef  ENABLE_MPI_LINGEN
int load_mpi_checkpoint_file_scattered(bmstatus_ptr bm, bigmatpoly xpi, unsigned int t0, unsigned int t1, unsigned int *delta, int * p_done)/*{{{*/
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm->com[0], &rank);
    dims * d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    struct cp_info cp[1];
    cp_info_init_mpi(cp, bm, t0, t1);
    ASSERT_ALWAYS(bigmatpoly_check_pre_init(xpi));
    size_t pi_size;
    int ok = cp_load_aux_file(cp, &pi_size, delta, p_done);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm->com[0]);
    MPI_Bcast(&pi_size, 1, MPI_MY_SIZE_T, 0, bm->com[0]);
    MPI_Bcast(delta, m + n, MPI_UNSIGNED, 0, bm->com[0]);
    MPI_Bcast(bm->lucky, m + n, MPI_INT, 0, bm->com[0]);
    MPI_Bcast(p_done, 1, MPI_INT, 0, bm->com[0]);
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Reading checkpoint %s (MPI, scattered)",
                cp->datafile);
        do {
            FILE * data = fopen(cp->datafile, "rb");
            int rc;
            ok = data != NULL;
            MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm->com[0]);
            if (!ok) {
                if (!rank)
                    fprintf(stderr, "Warning: cannot open %s\n", cp->datafile);
                if (data) free(data);
                break;
            }
            bigmatpoly_finish_init(ab, xpi, m+n, m+n, pi_size);
            bigmatpoly_set_size(xpi, pi_size);
            rc = matpoly_read(ab, data, bigmatpoly_my_cell(xpi), 0, xpi->size, 0, 0);
            ok = ok && rc == (int) xpi->size;
            rc = fclose(data);
            ok = ok && (rc == 0);
        } while (0);
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm->com[0]);
        logline_end(&bm->t_cp_io,"");
        if (!ok && !rank) {
            fprintf(stderr, "Warning: I/O error while reading %s\n",
                    cp->datafile);
        }
    }
    cp_info_clear(cp);
    if (ok) bm->t = t1;
    return ok;
}/*}}}*/

int save_mpi_checkpoint_file_scattered(bmstatus_ptr bm, bigmatpoly xpi, unsigned int t0, unsigned int t1, unsigned int *delta, int done)/*{{{*/
{
    /* corresponding t is bm->t - E->size ! */
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm->com[0], &rank);
    struct cp_info cp[1];
    cp_info_init_mpi(cp, bm, t0, t1);
    int ok = cp_save_aux_file(cp, xpi->size, delta, done);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm->com[0]);
    if (!ok && !rank) unlink(cp->auxfile);
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Saving checkpoint %s (MPI, scattered)",
                cp->datafile);
        ok = cp_save_data_file(cp, bigmatpoly_my_cell(xpi), xpi->size);
        logline_end(&bm->t_cp_io,"");
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm->com[0]);
        if (!ok) {
            if (cp->datafile) unlink(cp->datafile);
            if (!rank) unlink(cp->auxfile);
        }
    }
    if (!ok && !rank) {
        fprintf(stderr, "Warning: I/O error while saving %s\n", cp->datafile);
    }
    cp_info_clear(cp);
    return ok;
}/*}}}*/

int load_mpi_checkpoint_file_gathered(bmstatus_ptr bm, bigmatpoly xpi, unsigned int t0, unsigned int t1, unsigned int *delta, int * p_done)/*{{{*/
{
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm->com[0], &rank);
    dims * d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    struct cp_info cp[1];
    cp_info_init_mpi(cp, bm, t0, t1);
    cp->datafile = cp->gdatafile;
    size_t pi_size;
    int ok = cp_load_aux_file(cp, &pi_size, delta, p_done);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm->com[0]);
    MPI_Bcast(&pi_size, 1, MPI_MY_SIZE_T, 0, bm->com[0]);
    MPI_Bcast(delta, m + n, MPI_UNSIGNED, 0, bm->com[0]);
    MPI_Bcast(bm->lucky, m + n, MPI_INT, 0, bm->com[0]);
    MPI_Bcast(p_done, 1, MPI_INT, 0, bm->com[0]);
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Reading checkpoint %s (MPI, gathered)",
                cp->datafile);
        do {
            FILE * data = NULL;
            if (!rank) ok = (data = fopen(cp->datafile, "rb")) != NULL;
            MPI_Bcast(&ok, 1, MPI_INT, 0, bm->com[0]);
            if (!ok) {
                if (!rank)
                    fprintf(stderr, "Warning: cannot open %s\n", cp->datafile);
                if (data) free(data);
                break;
            }

            bigmatpoly_finish_init(ab, xpi, m+n, m+n, pi_size);
            bigmatpoly_set_size(xpi, pi_size);

            double avg = avg_matsize(ab, m + n, m + n, 0);
            unsigned int B = iceildiv(io_block_size, avg);

            /* This is only temp storage ! */
            matpoly pi;
            matpoly_init(ab, pi, m + n, m + n, B);
            pi->size = B;

            for(unsigned int k = 0 ; ok && k < xpi->size ; k += B) {
                unsigned int nc = MIN(B, xpi->size - k);
                if (!rank)
                    ok = matpoly_read(ab, data, pi, 0, nc, 0, 0) == (int) nc;
                MPI_Bcast(&ok, 1, MPI_INT, 0, bm->com[0]);
                bigmatpoly_scatter_mat_partial(ab, xpi, pi, k, nc);
            }
            matpoly_clear(ab, pi);

            if (!rank) {
                int rc = fclose(data);
                ok = ok && (rc == 0);
            }
        } while (0);
        MPI_Bcast(&ok, 1, MPI_INT, 0, bm->com[0]);
        logline_end(&bm->t_cp_io,"");
        if (!ok && !rank) {
            fprintf(stderr, "Warning: I/O error while reading %s\n",
                    cp->datafile);
        }
    }
    cp_info_clear(cp);
    if (ok) bm->t = t1;
    return ok;
}/*}}}*/

int save_mpi_checkpoint_file_gathered(bmstatus_ptr bm, bigmatpoly xpi, unsigned int t0, unsigned int t1, unsigned int *delta, int done)/*{{{*/
{
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm->com[0], &rank);
    dims * d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    struct cp_info cp[1];
    cp_info_init_mpi(cp, bm, t0, t1);
    cp->datafile = cp->gdatafile;
    logline_begin(stdout, SIZE_MAX, "Saving checkpoint %s (MPI, gathered)",
            cp->datafile);
    int ok = cp_save_aux_file(cp, xpi->size, delta, done);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm->com[0]);
    if (ok) {
        do {
            FILE * data = NULL;
            int rc;
            if (!rank) ok = (data = fopen(cp->datafile, "wb")) != NULL;
            MPI_Bcast(&ok, 1, MPI_INT, 0, bm->com[0]);
            if (!ok) {
                if (!rank)
                    fprintf(stderr, "Warning: cannot open %s\n", cp->datafile);
                if (data) free(data);
                break;
            }

            double avg = avg_matsize(ab, m + n, m + n, 0);
            unsigned int B = iceildiv(io_block_size, avg);

            /* This is only temp storage ! */
            matpoly pi;
            matpoly_init(ab, pi, m + n, m + n, B);
            pi->size = B;

            for(unsigned int k = 0 ; ok && k < xpi->size ; k += B) {
                unsigned int nc = MIN(B, xpi->size - k);
                bigmatpoly_gather_mat_partial(ab, pi, xpi, k, nc);
                if (!rank)
                    ok = matpoly_write(ab, data, pi, 0, nc, 0, 0) == (int) nc;
                MPI_Bcast(&ok, 1, MPI_INT, 0, bm->com[0]);
            }
            matpoly_clear(ab, pi);

            if (!rank) {
                rc = fclose(data);
                ok = ok && (rc == 0);
            }
        } while (0);
        MPI_Bcast(&ok, 1, MPI_INT, 0, bm->com[0]);
        if (!ok && !rank) {
            if (cp->datafile) unlink(cp->datafile);
            unlink(cp->auxfile);
        }
    }
    logline_end(&bm->t_cp_io,"");
    if (!ok && !rank) {
        fprintf(stderr, "Warning: I/O error while saving %s\n", cp->datafile);
    }
    cp_info_clear(cp);
    return ok;
}/*}}}*/

int load_mpi_checkpoint_file(bmstatus_ptr bm, bigmatpoly xpi, unsigned int t0, unsigned int t1, unsigned int *delta, int * p_done)/*{{{*/
{
    /* read scattered checkpoint with higher priority if available,
     * because we like distributed I/O. Otherwise, read gathered
     * checkpoint if we could find one.
     */
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm->com[0], &rank);
    struct cp_info cp[1];
    cp_info_init_mpi(cp, bm, t0, t1);
    int ok = 0;
    int aux_ok = rank || access(cp->auxfile, R_OK) == 0;
    int sdata_ok = access(cp->sdatafile, R_OK) == 0;
    int scattered_ok = aux_ok && sdata_ok;
    MPI_Allreduce(MPI_IN_PLACE, &scattered_ok, 1, MPI_INT, MPI_MIN, bm->com[0]);
    if (scattered_ok) {
        ok = load_mpi_checkpoint_file_scattered(bm, xpi, t0, t1, delta, p_done);
        if (ok) return ok;
    }
    int gdata_ok = rank || access(cp->gdatafile, R_OK) == 0;
    int gathered_ok = aux_ok && gdata_ok;
    MPI_Bcast(&gathered_ok, 1, MPI_INT, 0, bm->com[0]);
    if (gathered_ok) {
        ok = load_mpi_checkpoint_file_gathered(bm, xpi, t0, t1, delta, p_done);
    }
    cp_info_clear(cp);
    return ok;
}/*}}}*/

int save_mpi_checkpoint_file(bmstatus_ptr bm, bigmatpoly xpi, unsigned int t0, unsigned int t1, unsigned int *delta, int done)/*{{{*/
{
    if (save_gathered_checkpoints) {
        return save_mpi_checkpoint_file_gathered(bm, xpi, t0, t1, delta, done);
    } else {
        return save_mpi_checkpoint_file_scattered(bm, xpi, t0, t1, delta, done);
    }
}/*}}}*/
#endif  /* ENABLE_MPI_LINGEN */

/*}}}*/

/**********************************************************************/

/*{{{ Main entry points and recursive algorithm (with and without MPI) */

/* Forward declaration, it's used by the recursive version */
static int bw_lingen_single(bmstatus_ptr bm, matpoly pi, matpoly E, unsigned int *delta);

#ifdef  ENABLE_MPI_LINGEN
static int bw_biglingen_collective(bmstatus_ptr bm, bigmatpoly pi, bigmatpoly E, unsigned int *delta);
#endif

static int bw_lingen_recursive(bmstatus_ptr bm, matpoly pi, matpoly E, unsigned int *delta) /*{{{*/
{
    size_t z = E->size;
    stats.enter(__func__, E->size);
    dims * d = bm->d;
    abdst_field ab = d->ab;
    int done;

    /* XXX I think we have to start with something large enough to get
     * all coefficients of E_right correct */
    size_t half = E->size - (E->size / 2);

    // fprintf(stderr, "Enter %s\n", __func__);
    /* declare an lazy-alloc all matrices */
    matpoly E_left;
    matpoly pi_left;
    matpoly pi_right;
    matpoly E_right;
    matpoly_init(ab, E_left, 0, 0, 0);
    matpoly_init(ab, pi_left, 0, 0, 0);
    matpoly_init(ab, pi_right, 0, 0, 0);
    matpoly_init(ab, E_right, 0, 0, 0);

    matpoly_truncate(ab, E_left, E, half);
    done = bw_lingen_single(bm, pi_left, E_left, delta);
    ASSERT_ALWAYS(pi_left->size);
    matpoly_clear(ab, E_left);

    if (done) {
        matpoly_swap(pi_left, pi);
        matpoly_clear(ab, pi_left);
        // fprintf(stderr, "Leave %s\n", __func__);
        return stats.leave(1);
    }

    stats.begin_smallstep("MP");
    matpoly_rshift(ab, E, E, half - pi_left->size + 1);
    logline_begin(stdout, z, "t=%u MP(%zu, %zu) -> %zu",
            bm->t, E->size, pi_left->size, E->size - pi_left->size + 1);
    matpoly_mp(ab, E_right, E, pi_left);
    logline_end(&(bm->t_mp), "");
    stats.end_smallstep();

    done = bw_lingen_single(bm, pi_right, E_right, delta);
    matpoly_clear(ab, E_right);

    stats.begin_smallstep("MUL");
    logline_begin(stdout, z, "t=%u MUL(%zu, %zu) -> %zu",
            bm->t, pi_left->size, pi_right->size, pi_left->size + pi_right->size - 1);
    matpoly_mul(ab, pi, pi_left, pi_right);
    logline_end(&bm->t_mul, "");
    stats.end_smallstep();

    matpoly_clear(ab, pi_left);
    matpoly_clear(ab, pi_right);

    // fprintf(stderr, "Leave %s\n", __func__);
    return stats.leave(done);
}/*}}}*/

static int bw_lingen_single(bmstatus_ptr bm, matpoly pi, matpoly E, unsigned int *delta) /*{{{*/
{
    int rank;
    MPI_Comm_rank(bm->com[0], &rank);
    ASSERT_ALWAYS(!rank);
    unsigned int t0 = bm->t;
    unsigned int t1 = bm->t + E->size;

    int done;

    if (load_checkpoint_file(bm, pi, t0, t1, delta, &done))
        return done;

    // ASSERT_ALWAYS(E->size < bm->lingen_mpi_threshold);

    // fprintf(stderr, "Enter %s\n", __func__);
    if (E->size < bm->lingen_threshold) {
        bm->t_basecase -= seconds();
        done = bw_lingen_basecase(bm, pi, E, delta);
        bm->t_basecase += seconds();
    } else {
        done = bw_lingen_recursive(bm, pi, E, delta);
    }
    // fprintf(stderr, "Leave %s\n", __func__);

    save_checkpoint_file(bm, pi, t0, t1, delta, done);

    return done;
}/*}}}*/

#ifdef  ENABLE_MPI_LINGEN
static int bw_biglingen_recursive(bmstatus_ptr bm, bigmatpoly pi, bigmatpoly E, unsigned int *delta) /*{{{*/
{
    size_t z = E->size;
    stats.enter(__func__, E->size);

    dims * d = bm->d;
    abdst_field ab = d->ab;
    int done;

    int rank;
    MPI_Comm_rank(bm->com[0], &rank);

    // fprintf(stderr, "Enter %s\n", __func__);
    /* XXX I think we have to start with something large enough to get
     * all coefficients of E_right correct */
    size_t half = E->size - (E->size / 2);

    /* declare an lazy-alloc all matrices */
    bigmatpoly E_left;
    bigmatpoly E_right;
    bigmatpoly pi_left;
    bigmatpoly pi_right;
    bigmatpoly_init(ab, E_left, E, 0, 0, 0);
    bigmatpoly_init(ab, pi_left, pi, 0, 0, 0);
    bigmatpoly_init(ab, pi_right, pi, 0, 0, 0);
    bigmatpoly_init(ab, E_right, E, 0, 0, 0);

    bigmatpoly_truncate_loc(ab, E_left, E, half);
    done = bw_biglingen_collective(bm, pi_left, E_left, delta);
    bigmatpoly_clear(ab, E_left);

    if (done) {
        bigmatpoly_swap(pi_left, pi);
        bigmatpoly_clear(ab, pi_left);
        bigmatpoly_clear(ab, pi_right);
        bigmatpoly_clear(ab, E_right);
        // fprintf(stderr, "Leave %s\n", __func__);
        return stats.leave(1);
    }

    stats.begin_smallstep("MP");
    bigmatpoly_rshift(ab, E, E, half - pi_left->size + 1);
    logline_begin(stdout, z, "t=%u MPI-MP%s(%zu, %zu) -> %zu",
            bm->t, caching ? "-caching" : "",
            E->size, pi_left->size, E->size - pi_left->size + 1);
#ifdef  ENABLE_CACHING_MPI_LINGEN
    if (caching)
        bigmatpoly_mp_caching(ab, E_right, E, pi_left);
    else
#endif
        bigmatpoly_mp(ab, E_right, E, pi_left);
    logline_end(&bm->t_mp, "");
    stats.end_smallstep();

    done = bw_biglingen_collective(bm, pi_right, E_right, delta);
    bigmatpoly_clear(ab, E_right);

    stats.begin_smallstep("MUL");
    logline_begin(stdout, z, "t=%u MPI-MUL%s(%zu, %zu) -> %zu",
            bm->t, caching ? "-caching" : "",
            pi_left->size, pi_right->size, pi_left->size + pi_right->size - 1);
#ifdef  ENABLE_CACHING_MPI_LINGEN
    if (caching)
        bigmatpoly_mul_caching(ab, pi, pi_left, pi_right);
    else
#endif
        bigmatpoly_mul(ab, pi, pi_left, pi_right);
    logline_end(&bm->t_mul, "");
    stats.end_smallstep();

    bigmatpoly_clear(ab, pi_left);
    bigmatpoly_clear(ab, pi_right);

    // fprintf(stderr, "Leave %s\n", __func__);
    return stats.leave(done);
}/*}}}*/

static int bw_biglingen_collective(bmstatus_ptr bm, bigmatpoly pi, bigmatpoly E, unsigned int *delta)/*{{{*/
{
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    int done;
    int rank;
    int size;
    MPI_Comm_rank(bm->com[0], &rank);
    MPI_Comm_size(bm->com[0], &size);
    unsigned int t0 = bm->t;
    unsigned int t1 = bm->t + E->size;

    if (load_mpi_checkpoint_file(bm, pi, t0, t1, delta, &done))
        return done;

    // fprintf(stderr, "Enter %s\n", __func__);
    if (E->size >= bm->lingen_mpi_threshold)  {
        done = bw_biglingen_recursive(bm, pi, E, delta);
    } else {
        /* Fall back to local code */
        /* This entails gathering E locally, computing pi locally, and
         * dispathing it back. */

        matpoly sE, spi;
        matpoly_init(ab, sE, m, b, E->size);
        matpoly_init(ab, spi, 0, 0, 0);

        stats.begin_smallstep("gather(L+R)");
        bigmatpoly_gather_mat(ab, sE, E);
        stats.end_smallstep();

        /* Only the master node does the local computation */
        if (!rank)
            done = bw_lingen_single(bm, spi, sE, delta);

        stats.begin_smallstep("scatter(L+R)");
        bigmatpoly_scatter_mat(ab, pi, spi);
        MPI_Bcast(&done, 1, MPI_INT, 0, bm->com[0]);
        MPI_Bcast(delta, b, MPI_UNSIGNED, 0, bm->com[0]);
        MPI_Bcast(bm->lucky, b, MPI_UNSIGNED, 0, bm->com[0]);
        MPI_Bcast(&(bm->t), 1, MPI_UNSIGNED, 0, bm->com[0]);
        /* Don't forget to broadcast delta from root node to others ! */
        matpoly_clear(ab, spi);
        matpoly_clear(ab, sE);
        stats.end_smallstep();
    }
    // fprintf(stderr, "Leave %s\n", __func__);

    save_mpi_checkpoint_file(bm, pi, t0, t1, delta, done);

    return done;
}/*}}}*/
#endif  /* ENABLE_MPI_LINGEN */

/*}}}*/

/**********************************************************************/

void get_minmax_delta_on_solutions(bmstatus_ptr bm, unsigned int * res, unsigned int * delta)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int maxdelta = 0;
    unsigned int mindelta = UINT_MAX;
    for(unsigned int j = 0 ; j < m + n ; j++) {
        if (bm->lucky[j] <= 0)
            continue;
        if (delta[j] > maxdelta)
            maxdelta = delta[j];
        if (delta[j] < mindelta)
            mindelta = delta[j];
    }
    res[0] = mindelta;
    res[1] = maxdelta;
}/*}}}*/

unsigned int get_max_delta_on_solutions(bmstatus_ptr bm, unsigned int * delta)/*{{{*/
{
    unsigned int res[2];
    get_minmax_delta_on_solutions(bm, res, delta);
    return res[1];
}/*}}}*/

/**********************************************************************/

struct bm_io_s {/*{{{*/
    bmstatus_ptr bm;
    unsigned int t0;
    FILE ** fr; /* array of n files when split_input_file is supported
                   (which is not the case as of now),
                   or otherwise just a 1-element array */
    FILE ** fw; /* array of n^2 files (with split_output_file)
                   or otherwise just a 1-element array */
    char * iobuf;
    const char * input_file;
    const char * output_file;
    int ascii;
    /* This is only a rolling window ! */
    matpoly A, F;
    unsigned int (*fdesc)[2];
    /* This k is the coefficient in A(X) div X of the next coefficient to
     * be read. This is thus the total number of coefficients of A(X) div
     * X which have been read so far.
     * In writing mode, k is the number of coefficients of F which have
     * been written so far.
     */
    unsigned int k;

    unsigned int guessed_length;
};

typedef struct bm_io_s bm_io[1];
typedef struct bm_io_s * bm_io_ptr;/*}}}*/

/* The reading mode of bm_io is streaming, but with a look-back
 * functionality: we want to be able to access coefficients a few places
 * earlier, so we keep them in memory.
 *
 * The writing mode has a write-ahead feature. Coefficients of the result
 * are written at various times, some earlier than others. The time span
 * is the same as for the reading mode.
 */


/* We write the coefficients of the reversed polynomial \hat{F*\pi}, in
 * increasing degree order. Thus the first coefficients written
 * corresponds to high degree coefficients in \pi.
 * This is mostly for historical reasons, since in fact, we would prefer
 * having the coefficients in the reverse order.
 */

/* let mindelta and maxdelta be (as their name suggests) the minimum and
 * maximum over all deltas corresponding to solution columns.
 *
 * For 0<=i<n, coeff i,j,k of pi becomes coeff i,j,k' of the final f,
 * with k'=k-(maxdelta-delta[j]).
 *
 * For 0<=i<m, coeff n+i,j,k of pi becomes coeff c[i],j,k' of the final f,
 * with k'=k-(maxdelta-delta[j])-(t0-e[j]).
 *
 * Therefore the maximum write-behind distance is (maxdelta-mindelta)+t0.
 * We need one coeff more (because offset goes from 0 to maxoffset).
 */
unsigned int bm_io_set_write_behind_size(bm_io_ptr aa, unsigned int * delta)
{
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    abdst_field ab = d->ab;
    unsigned int res[2];
    get_minmax_delta_on_solutions(bm, res, delta);
    unsigned int mindelta = res[0];
    unsigned int maxdelta = res[1];
    unsigned int window = maxdelta - mindelta + aa->t0 + 1;
    if (d->nrhs) {
        /* in sm-outside-matrix mode for DLP, we form the matrix F
         * slightly differently, as some *rows* are shifted out before
         * writing */
        window++;
    }
    int rank;
    MPI_Comm_rank(aa->bm->com[0], &rank);
    if (rank) return window;
    matpoly_realloc(ab, aa->F, window);
    matpoly_zero(ab, aa->F);
    aa->F->size = window;
    return window;
}

/* Note that coefficients of F are written transposed */
void bm_io_write_one_F_coeff(bm_io_ptr aa, unsigned int deg)
{
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    abdst_field ab = d->ab;
    deg = deg % aa->F->alloc;
    // if (random_input_length) return;
    /* what the heck ?
    if (!aa->ascii || !random_input_length)
        matpoly_write(ab, aa->f, aa->F, deg, deg + 1, aa->ascii, 1);
        */
    if (random_input_length) {
        /* we used to do this only in the non-ascii case */
        matpoly_write(ab, aa->fw[0], aa->F, deg, deg + 1, aa->ascii, 1);
    } else if (split_output_file) {
        matpoly_write_split(ab, aa->fw, aa->F, deg, deg + 1, aa->ascii);
    } else {
        matpoly_write(ab, aa->fw[0], aa->F, deg, deg + 1, aa->ascii, 1);
    }

}

void bm_io_clear_one_F_coeff(bm_io_ptr aa, unsigned int deg)
{
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    abdst_field ab = d->ab;
    unsigned int n = d->n;
    deg = deg % aa->F->alloc;
    for(unsigned int j = 0 ; j < n ; j++) {
        for(unsigned int i = 0 ; i < n ; i++) {
            abset_zero(ab, matpoly_coeff(ab, aa->F, i, j, deg));
        }
    }
}

/* Transfers of matrix entries, from memory to memory ; this is possibly
 * almost a transparent layer, and the matpoly_* instances really _are_
 * transparent layers. However we want a unified interface for the case
 * where we read an MPI-scattered matrix by chunks.
 *
 * Given that we are talking memory-to-memory, the classes below are used
 * as follows. The *_write_task classes are at the beginning of the
 * computation, where the matrix A is read from disk, and the matrix E is
 * written to (possibly distributed) memory. The one which matters is E
 * here. Symmetrically, the *_read_task are for reading the computed
 * matrix pi from (possibly distributed) memory, and writing the final
 * data F to disk.
 *
 * Pay attention to the fact that the calls to these structures are
 * collective calls.
 */
#ifdef  ENABLE_MPI_LINGEN
class bigmatpoly_read_task { /* {{{ */
    /* This reads a bigmatpoly, by chunks, so that the memory footprint
     * remains reasonable. */
    bigmatpoly_srcptr xpi;
    matpoly pi; /* This is only temp storage ! */
    unsigned int B;
    unsigned int k0;
    abdst_field ab;
    int rank;

    public:
    bigmatpoly_read_task(bm_io_ptr aa, bigmatpoly_srcptr xpi) : xpi(xpi) {
        bmstatus_ptr bm = aa->bm;
        dims * d = bm->d;
        unsigned int m = d->m;
        unsigned int n = d->n;
        unsigned int b = m + n;
        MPI_Comm_rank(bm->com[0], &rank);
        ab = d->ab;

        /* Decide on the temp storage size */
        double avg = avg_matsize(ab, n, n, aa->ascii);
        B = iceildiv(io_block_size, avg);
        if (!rank && !random_input_length) {
            printf("Writing F to %s\n", aa->output_file);
            printf("Writing F by blocks of %u coefficients"
                    " (%.1f bytes each)\n", B, avg);
        }
        matpoly_init(ab, pi, b, b, B);

        k0 = UINT_MAX;
    }

    inline unsigned int chunk_size() const { return B; }

    inline unsigned int size() { return xpi->size; }

    inline absrc_elt coeff_const_locked(unsigned int i, unsigned int j, unsigned int k) {
        ASSERT(!rank);
        ASSERT(k0 != UINT_MAX && k - k0 < B);
        return matpoly_coeff_const(ab, pi, i, j, k - k0);
    }

    absrc_elt coeff_const(unsigned int i, unsigned int j, unsigned int k) {
        if (k0 == UINT_MAX || k - k0 >= B) {
            k0 = k - (k % B);
            matpoly_zero(ab, pi);
            pi->size = B;
            bigmatpoly_gather_mat_partial(ab, pi, xpi, k0, MIN(xpi->size - k0, B));
        }
        if (rank) return NULL;
        return coeff_const_locked(i, j, k);
    }

    ~bigmatpoly_read_task() {
        matpoly_clear(ab, pi);
    }
};      /* }}} */
class bigmatpoly_write_task { /* {{{ */
    /* This writes a bigmatpoly, by chunks, so that the memory footprint
     * remains reasonable.
     * Note that in any case, the coefficient indices mut be progressive
     * in the write.
     */
    bigmatpoly_ptr xE;
    matpoly E; /* This is only temp storage ! */
    unsigned int B;
    unsigned int k0;
    abdst_field ab;
    int rank;

    /* forbid copies */
    bigmatpoly_write_task(bigmatpoly_write_task&) {}

    public:
    bigmatpoly_write_task(bm_io_ptr aa, bigmatpoly_ptr xE) : xE(xE) {
        bmstatus_ptr bm = aa->bm;
        dims * d = bm->d;
        unsigned int m = d->m;
        unsigned int n = d->n;
        unsigned int b = m + n;
        ab = d->ab;
        MPI_Comm_rank(aa->bm->com[0], &rank);


        /* Decide on the MPI chunk size */
        double avg = avg_matsize(ab, m, n, 0);
        B = iceildiv(io_block_size, avg);

        if (!rank) {
            /* TODO: move out of here */
            if (aa->input_file)
                printf("Reading A from %s\n", aa->input_file);
            printf("Computing E by chunks of %u coefficients (%.1f bytes each)\n",
                    B, avg);
        }

        matpoly_init(ab, E, m, b, B);

        /* Setting E->size is rather artificial, since we use E essentially
         * as an area where we may write coefficients freely. The only aim is
         * to escape some safety checks involving ->size in matpoly_part */
        matpoly_zero(ab, E);

        k0 = UINT_MAX;
    }

    inline unsigned int chunk_size() const { return B; }

    // inline unsigned int size() { return xE->size; }
    inline void set_size(unsigned int s) {
        bigmatpoly_set_size(xE, s);
    }

    inline abdst_elt coeff_locked(unsigned int i, unsigned int j, unsigned int k) {
        ASSERT(k0 != UINT_MAX && k - k0 < B);
        ASSERT(!rank);

        E->size += (E->size == k - k0);
        ASSERT(E->size == k - k0 + 1);

        return matpoly_coeff(ab, E, i, j, k - k0);
    }

    abdst_elt coeff(unsigned int i, unsigned int j, unsigned int k) {
        if (k0 == UINT_MAX) {
            matpoly_zero(ab, E);
            E->size = 0;
            ASSERT(k == 0);
            k0 = 0;
        } else {
            if (k >= (k0 + B)) {
                /* We require progressive reads */
                ASSERT(k == k0 + B);
                bigmatpoly_scatter_mat_partial(ab, xE, E, k0, B);
                matpoly_zero(ab, E);
                k0 += B;
            }
        }
        ASSERT(k0 != UINT_MAX && k - k0 < B);
        if (rank) return NULL;
        return coeff_locked(i, j, k);
    }

    void finalize(unsigned int length) {
        if (length > k0) {
            /* Probably the last chunk hasn't been dispatched yet */
            ASSERT(length < k0 + B);
            bigmatpoly_scatter_mat_partial(ab, xE, E, k0, length - k0);
        }
        set_size(length);
    }

    ~bigmatpoly_write_task() {
        matpoly_clear(ab, E);
    }
    friend void matpoly_extract_column(abdst_field ab,
        bigmatpoly_write_task& dst, unsigned int jdst, unsigned int kdst,
        matpoly_ptr src, unsigned int jsrc, unsigned int ksrc);
};

void matpoly_extract_column(abdst_field ab,
        bigmatpoly_write_task& dst, unsigned int jdst, unsigned int kdst,
        matpoly_ptr src, unsigned int jsrc, unsigned int ksrc)
{
    dst.coeff_locked(0, jdst, kdst);
    matpoly_extract_column(ab, dst.E, jdst, kdst - dst.k0, src, jsrc, ksrc);
}

/* }}} */
#endif  /* ENABLE_MPI_LINGEN */
class matpoly_read_task {/*{{{*/
    /* This does the same, but for a simple matpoly. Of course this is
     * much simpler ! */
    matpoly_srcptr pi;
    abdst_field ab;

    public:
    matpoly_read_task(bm_io_ptr aa, matpoly_srcptr pi) :
            pi(pi),
            ab(aa->bm->d->ab)
    {
        if (!random_input_length) {
            printf("Writing F to %s\n", aa->output_file);
        }
    }

    inline unsigned int chunk_size() const { return 1; }

    inline unsigned int size() { return pi->size; }

    absrc_elt coeff_const_locked(unsigned int i, unsigned int j, unsigned int k) {
        return matpoly_coeff_const(ab, pi, i, j, k);
    }
    inline absrc_elt coeff_const(unsigned int i, unsigned int j, unsigned int k) {
        return coeff_const_locked(i, j, k);
    }

    ~matpoly_read_task() { }
};/*}}}*/
class matpoly_write_task { /* {{{ */
    matpoly_ptr E;
    abdst_field ab;

    /* forbid copies */
    matpoly_write_task(matpoly_write_task&) {}

    public:
    matpoly_write_task(bm_io_ptr aa, matpoly_ptr E) :
        E(E),
        ab(aa->bm->d->ab)
    {
            /* TODO: move out of here */
            if (aa->input_file)
                printf("Reading A from %s\n", aa->input_file);
    }

    inline unsigned int chunk_size() const { return 1; }
    inline void set_size(unsigned int s) { E->size = s; }
    // inline unsigned int size() { return E->size; }

    abdst_elt coeff_locked(unsigned int i, unsigned int j, unsigned int k) {
        E->size += (E->size == k);
        ASSERT(E->size == k + 1);
        return matpoly_coeff(ab, E, i, j, k);
    }
    inline abdst_elt coeff(unsigned int i, unsigned int j, unsigned int k) {
        return coeff_locked(i, j, k);
    }
    void finalize(unsigned int length) { set_size(length); }
    ~matpoly_write_task() { }
    friend void matpoly_extract_column(abdst_field ab,
        matpoly_write_task& dst, unsigned int jdst, unsigned int kdst,
        matpoly_ptr src, unsigned int jsrc, unsigned int ksrc);
};
void matpoly_extract_column(abdst_field ab,
        matpoly_write_task& dst, unsigned int jdst, unsigned int kdst,
        matpoly_ptr src, unsigned int jsrc, unsigned int ksrc)
{
    matpoly_extract_column(ab, dst.E, jdst, kdst, src, jsrc, ksrc);
}

/* }}} */


template<class Reader>
void bm_io_compute_final_F(bm_io_ptr aa, Reader& pi, unsigned int * delta)/*{{{*/
{
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    abdst_field ab = d->ab;
    int rank;
    MPI_Comm_rank(bm->com[0], &rank);
    int leader = rank == 0;


    /* We are not interested by pi->size, but really by the number of
     * coefficients for the columns which give solutions. */
    unsigned int maxdelta = get_max_delta_on_solutions(bm, delta);

    if (leader) printf("Final f(X)=f0(X)pi(X) has degree %u\n", maxdelta);

    unsigned int window = aa->F->alloc;

    /* Which columns of F*pi will make the final generator ? */
    unsigned int * sols = (unsigned int *) malloc(n * sizeof(unsigned int));
    for(unsigned int j = 0, jj=0 ; j < m + n ; j++) {
        if (bm->lucky[j] <= 0)
            continue;
        sols[jj++]=j;
    }

    double tt0 = wct_seconds();
    double next_report_t = tt0 + 10;
    unsigned next_report_s = pi.size() / 100;

    /*
     * first compute the rhscontribs. Use that to decide on a renumbering
     * of the columns, because we'd like to get the "primary" solutions
     * first, in a sense. Those are the ones with fewer zeroes in the
     * rhscontrib part. So we would want that matrix to have its columns
     * sorted in decreasing weight order
     *
     * An alternative, possibly easier, is to have a function which
     * decides the solution ordering precisely based on the inspection of
     * this rhscoeffs matrix (?). But how should we spell that info when
     * we give it to mksol ??
     */

    /* This **modifies** the "sols" array */
    if (d->nrhs) {
        /* The data is only gathered at rank 0. So the stuff we compute
         * can only be meaningful there. On the other hand, it is
         * necessary that we call *collectively* the pi.coeff() routine,
         * because we need synchronisation of all ranks.
         */

        matpoly rhs;
        if (leader) {
            matpoly_init(ab, rhs, d->nrhs, n, 1);
            rhs->size = 1;
        }

        /* Now redo the exact same loop as above, this time
         * adding the contributions to the rhs matrix. */
        for(unsigned int ipi = 0 ; ipi < m + n ; ipi++) {
            for(unsigned int jF = 0 ; jF < n ; jF++) {
                unsigned int jpi = sols[jF];
                unsigned int iF, offset;
                if (ipi < n) {
                    iF = ipi;
                    offset = 0;
                } else {
                    iF = aa->fdesc[ipi-n][1];
                    offset = aa->t0 - aa->fdesc[ipi-n][0];
                }
                if (iF >= d->nrhs) continue;
                ASSERT_ALWAYS(delta[jpi] >= offset);
                unsigned kpi = delta[jpi] - offset;

                ASSERT_ALWAYS(d->nrhs);
                ASSERT_ALWAYS(iF < d->nrhs);
                ASSERT_ALWAYS(jF < n);
                absrc_elt src = pi.coeff_const(ipi, jpi, kpi);

                if (leader) {
                    abdst_elt dst = matpoly_coeff(ab, rhs, iF, jF, 0);
                    abadd(ab, dst, dst, src);
                }
            }
        }

        if (leader) {
            /* Now comes the time to prioritize the different solutions. Our
             * goal is to get the unessential solutions last ! */
            int (*sol_score)[2];
            sol_score = (int (*)[2]) malloc(n * 2 * sizeof(int));
            memset(sol_score, 0, n * 2 * sizeof(int));
            /* score per solution is the number of non-zero coefficients,
             * that's it. Since we have access to lexcmp2, we want to use it.
             * Therefore, desiring the highest scoring solutions first, we
             * negate the hamming weight.
             */
            for(unsigned int jF = 0 ; jF < n ; jF++) {
                sol_score[jF][1] = jF;
                for(unsigned int iF = 0 ; iF < d->nrhs ; iF++) {
                    int z = !abis_zero(ab, matpoly_coeff(ab, rhs, iF, jF, 0));
                    sol_score[jF][0] -= z;
                }
            }
            qsort(sol_score, n, 2 * sizeof(int), (sortfunc_t) & lexcmp2);

            if (leader) {
                printf("Reordered solutions:\n");
                for(unsigned int i = 0 ; i < n ; i++) {
                    printf(" %d (col %d in pi, weight %d on rhs vectors)\n", sol_score[i][1], sols[sol_score[i][1]], -sol_score[i][0]);
                }
            }

            /* We'll now modify the sols[] array, so that we get a reordered
             * F, too (and mksol/gather don't have to care about our little
             * tricks */
            {
                matpoly rhs2;
                matpoly_init(ab, rhs2, d->nrhs, n, 1);
                rhs2->size = 1;
                for(unsigned int i = 0 ; i < n ; i++) {
                    matpoly_extract_column(ab, rhs2, i, 0, rhs, sol_score[i][1], 0);
                }
                matpoly_swap(rhs2, rhs);
                matpoly_clear(ab, rhs2);
                if (sol_score[0][0] == 0) {
                    if (allow_zero_on_rhs) {
                        printf("Note: all solutions have zero contribution on the RHS vectors -- we will just output right kernel vectors (ok because of allow_zero_on_rhs=1)\n");
                    } else {
                        fprintf(stderr, "ERROR: all solutions have zero contribution on the RHS vectors -- we will just output right kernel vectors (maybe use allow_zero_on_rhs ?)\n");
                        rank0_exit_code = EXIT_FAILURE;
                    }
                }
                /* ugly: use sol_score[i][0] now to provide the future
                 * "sols" array. We'll get rid of sol_score right afterwards
                 * anyway.
                 */
                for(unsigned int i = 0 ; i < n ; i++) {
                    sol_score[i][0] = sols[sol_score[i][1]];
                }
                for(unsigned int i = 0 ; i < n ; i++) {
                    sols[i] = sol_score[i][0];
                }
            }

            free(sol_score);

            /* Because we've done the swap above, the rhs matrix is now a
             * matrix with dimension nrhs * n, as required.
             *
             * XXX hmm. there are nrhs contributions, but n solutions. So
             * we'd better get that transposed, right ?
             */
            matpoly_transpose_dumb(ab, rhs, rhs);

            if (split_output_file) {
                FILE ** fw = (FILE**) malloc(d->n * d->nrhs * sizeof(FILE*));
                for(unsigned int i = 0 ; i < d->n ; i++) {
                    for(unsigned int j = 0 ; j < d->nrhs ; j++) {
                        char * tmp;
                        int rc = asprintf(&tmp, "%s.sols%u-%u.%u-%u.rhs",
                                aa->output_file, i, i+1, j, j+1);
                        ASSERT_ALWAYS(rc >= 0);
                        fw[i*d->nrhs+j] = fopen(tmp, aa->ascii ? "w" : "wb");
                        ASSERT_ALWAYS(fw[i*d->nrhs+j]);
                        free(tmp);
                    }
                }
                matpoly_write_split(ab, fw, rhs, 0, 1, aa->ascii);
                for(unsigned int i = 0 ; i < d->n * d->nrhs ; i++) {
                    fclose(fw[i]);
                }
            } else {
                /* we used to write F transposed, but we no longer do
                 * that.
                 */
                char * tmp;
                int rc = asprintf(&tmp, "%s.rhs", aa->output_file);
                ASSERT_ALWAYS(rc >= 0);
                FILE * f = fopen(tmp, aa->ascii ? "w" : "wb");
                ASSERT_ALWAYS(f);
                matpoly_write(ab, f, rhs, 0, 1, aa->ascii, 0);
                fclose(f);
                printf("Note: contributions to RHS coefficients saved to the rhs file %s\n", tmp);
                free(tmp);
            }
            matpoly_clear(ab, rhs);
        }
    }

    /* we need to read pi backwards. The number of coefficients in pi is
     * pilen = maxdelta + 1 - t0. Hence the first interesting index is
     * maxdelta - t0. However, for notational ease, we'll access
     * coefficients from index pi->size downwards. The latter is always
     * large enough.
     */

    ASSERT_ALWAYS(pi.size() >= maxdelta + 1 - aa->t0);

    for(unsigned int s = 0 ; s < pi.size() ; s++) {
        unsigned int kpi = pi.size() - 1 - s;

        /* as above, we can't just have ranks > 0 do nothing. We need
         * synchronization of the calls to pi.coeff()
         */

        /* This call is here only to trigger the gather call. This one
         * must therefore be a collective call. Afterwards we'll use a
         * _locked call which is okay to call only at rank 0.
         */
        pi.coeff_const(0, 0, kpi);

        if (rank) continue;

        /* Coefficient kpi + window of F has been totally computed,
         * because of previous runs of this loop (which reads the
         * coefficients of pi).
         */
        if (kpi + window == maxdelta) {
            /* Avoid writing zero coefficients. This can occur !
             * Example:
             * tt=(2*4*1200) mod 1009, a = (tt cat tt cat * tt[0..10])
             */
            for(unsigned int j = 0 ; j < n ; j++) {
                int z = 1;
                for(unsigned int i = 0 ; z && i < n ; i++) {
                    absrc_elt src = matpoly_coeff_const(ab, aa->F, i, j, 0);
                    z = abis_zero(ab, src);
                }

                if (z) {
                    /* This is a bit ugly. Given that we're going to
                     * shift one column of F, we'll have a
                     * potentially deeper write-back buffer. Columns
                     * which seemed to be ready still are, but they
                     * will now be said so only at the next step.
                     */
                    printf("Reduced solution column #%u from"
                            " delta=%u to delta=%u\n",
                            sols[j], delta[sols[j]], delta[sols[j]]-1);
                    window++;
                    matpoly_realloc(ab, aa->F, window);
                    aa->F->size = window;
                    delta[sols[j]]--;
                    /* shift this column */
                    for(unsigned int k = 1 ; k < window ; k++) {
                        matpoly_extract_column(ab,
                                aa->F, j, k-1, aa->F, j, k);
                    }
                    matpoly_zero_column(ab, aa->F, j, window - 1);
                    break;
                }
            }
        }
        /* coefficient of degree maxdelta-kpi-window is now complete */
        if (kpi + window <= maxdelta) {
            bm_io_write_one_F_coeff(aa, (maxdelta-kpi) - window);
            bm_io_clear_one_F_coeff(aa, (maxdelta-kpi) - window);
        }
        /* Now use our coefficient. This might tinker with
         * coefficients up to degree kpi-(window-1) in the file F */

        if (kpi > maxdelta + aa->t0 ) {
            /* this implies that we always have kF > delta[jpi],
             * whence we expect a zero contribution */
            for(unsigned int i = 0 ; i < m + n ; i++) {
                for(unsigned int jF = 0 ; jF < n ; jF++) {
                    unsigned int jpi = sols[jF];
                    absrc_elt src = pi.coeff_const_locked(i, jpi, kpi);
                    ASSERT_ALWAYS(abis_zero(ab, src));
                }
            }
            continue;
        }

        for(unsigned int ipi = 0 ; ipi < m + n ; ipi++) {
            for(unsigned int jF = 0 ; jF < n ; jF++) {
                unsigned int jpi = sols[jF];
                unsigned int iF, offset;
                if (ipi < n) {
                    /* Left part of the initial F is x^0 times
                     * identity. Therefore, the first n rows of pi
                     * get multiplied by this identity matrix, this
                     * is pretty simple.
                     */
                    iF = ipi;
                    offset = 0;
                } else {
                    /* next m rows of the initial F are of the form
                     * x^(some value) times some canonical basis
                     * vector. Therefore, the corresponding row in pi
                     * ends up contributing to some precise row in F,
                     * and with an offset which is dictated by the
                     * exponent of x.
                     */
                    iF = aa->fdesc[ipi-n][1];
                    offset = aa->t0 - aa->fdesc[ipi-n][0];
                }
                unsigned int subtract = maxdelta - delta[jpi] + offset;
                ASSERT(subtract < window);
                if (maxdelta < kpi + subtract) continue;
                unsigned int kF = (maxdelta - kpi) - subtract;
                unsigned int kF1 = kF - (iF < d->nrhs);
                abdst_elt dst;
                if (kF1 == UINT_MAX) {
                    /* this has been addressed in the first pass,
                     * earlier.
                     */
                    continue;
                } else {
                    dst = matpoly_coeff(ab, aa->F, iF, jF, kF1 % window);
                }
                absrc_elt src = pi.coeff_const_locked(ipi, jpi, kpi);
                ASSERT_ALWAYS(kF <= delta[jpi] || abis_zero(ab, src));
                abadd(ab, dst, dst, src);
            }
        }

        if (leader && s > next_report_s) {
            double tt = wct_seconds();
            if (tt > next_report_t) {
                printf(
                        "Written %u coefficients (%.1f%%) in %.1f s\n",
                        s, 100.0 * s / pi.size(), tt-tt0);
                next_report_t = tt + 10;
                next_report_s = s + pi.size() / 100;
            }
        }
    }
    /* flush the pipe */
    if (leader && window <= maxdelta) {
        for(unsigned int s = window ; s-- > 0 ; )
            bm_io_write_one_F_coeff(aa, maxdelta - s);
    }

    free(sols);
}/*}}}*/

/* read 1 coefficient into the sliding window of input coefficients of
 * the input series A. The io_window parameter controls the size of the
 * sliding window. There are in fact two behaviours:
 *  - io_window == 0: there is no sliding window, really, and the new
 *    coefficient is appended as the last coefficient of aa->A.
 *  - io_window > 0: there, we really have a sliding window. Coeffs
 *    occupy places in a circular fashion within the buffer.
 */
int bm_io_read1(bm_io_ptr aa, unsigned int io_window)/*{{{*/
{
    bmstatus_ptr bm = aa->bm;
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    matpoly_ptr A = aa->A;
    static unsigned int generated_random_coefficients = 0;

    unsigned int pos = aa->k;
    if (io_window) {
        pos = pos % io_window;
        ASSERT_ALWAYS(A->size == io_window);
    } else {
        ASSERT_ALWAYS(A->size == aa->k);
        if (aa->k >= A->alloc) {
            matpoly_realloc(aa->bm->d->ab, A, A->alloc + 1);
        }
        A->size++;
    }
    if (random_input_length) {
        if (generated_random_coefficients >= random_input_length)
            return 0;

        for (unsigned int i = 0; i < m ; i++) {
            for (unsigned int j = 0; j < n ; j++) {
                abdst_elt x = matpoly_coeff(ab, A, i, j, pos);
                abrandom(ab, x, rstate);
            }
        }
        generated_random_coefficients++;
        aa->k++;
        return 1;
    }

    for (unsigned int i = 0; i < m ; i++) {
        for (unsigned int j = 0; j < n ; j++) {
            abdst_elt x = matpoly_coeff(ab, A, i, j, pos);
            int rc;
            if (aa->ascii) {
                rc = abfscan(ab, aa->fr[0], x);
                /* rc is the number of bytes read -- non-zero on success */
            } else {
                size_t elemsize = abvec_elt_stride(ab, 1);
                rc = fread(x, elemsize, 1, aa->fr[0]);
                rc = rc == 1;
                abnormalize(ab, x);
            }
            if (!rc) {
                if (i == 0 && j == 0) {
                    return 0;
                }
                fprintf(stderr,
                        "Parse error while reading coefficient (%d,%d,%d)%s\n",
                        i, j, 1 + aa->k,
                        aa->ascii ? "" : " (forgot --ascii?)");
                exit(EXIT_FAILURE);
            }
        }
    }
    aa->k++;
    return 1;
}/*}}}*/

void bm_io_init(bm_io_ptr aa, bmstatus_ptr bm, const char * input_file, const char * output_file, int ascii)/*{{{*/
{
    memset(aa, 0, sizeof(*aa));
    aa->bm = bm;
    aa->ascii = ascii;
    aa->input_file = input_file;
    aa->output_file = output_file;
}/*}}}*/

void bm_io_begin_read(bm_io_ptr aa)/*{{{*/
{
    bmstatus_ptr bm = aa->bm;
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    int rank;
    MPI_Comm_rank(aa->bm->com[0], &rank);
    if (rank) return;

    matpoly_init(ab, aa->A, m, n, 1);

    if (random_input_length) {
        /* see below. I think it would be a bug to not do that */
        bm_io_read1(aa, 0);
        return;
    }

    if (split_input_file) {
        fprintf(stderr, "--split-input-file not supported yet\n");
        exit(EXIT_FAILURE);
    }
    aa->fr = (FILE**) malloc(sizeof(FILE*));
    aa->fr[0] = fopen(aa->input_file, aa->ascii ? "r" : "rb");

    DIE_ERRNO_DIAG(aa->fr[0] == NULL, "fopen", aa->input_file);
    aa->iobuf = (char*) malloc(2 * io_block_size);
    setbuffer(aa->fr[0], aa->iobuf, 2 * io_block_size);

    /* read the first coefficient ahead of time. This is because in most
     * cases, we'll discard it. Only in the DL case, we will consider the
     * first coefficient as being part of the series. Which means that
     * the coefficient reads in the I/O loop will sometimes correspond to
     * the coefficient needed at that point in time, while we will also
     * (in the DL case) need data from the previous read.
     */
    if (!bm_io_read1(aa, 0)) {
        fprintf(stderr, "Read error from %s\n", aa->input_file);
        exit(EXIT_FAILURE);
    }
}/*}}}*/

void bm_io_end_read(bm_io_ptr aa)/*{{{*/
{
    int rank;
    MPI_Comm_rank(aa->bm->com[0], &rank);
    if (rank) return;
    matpoly_clear(aa->bm->d->ab, aa->A);
    if (random_input_length) return;
    fclose(aa->fr[0]);
    free(aa->fr);
    aa->fr = NULL;
    free(aa->iobuf);
    aa->iobuf = 0;
}/*}}}*/

void bm_io_begin_write(bm_io_ptr aa)/*{{{*/
{
    bmstatus_ptr bm = aa->bm;
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int n = d->n;
    int rank;
    MPI_Comm_rank(aa->bm->com[0], &rank);
    if (rank) return;

    matpoly_init(ab, aa->F, n, n, aa->t0 + 1);

    if (random_input_length) {
        /* It's horrible. I really want "checksum" to be printed at the
         * exact right moment, hence the following hack...  */
        aa->fw = (FILE**) malloc(sizeof(FILE*));
        aa->fw[0] = popen("(echo \"checksum:\" ; sha1sum) | xargs echo", "w");
        DIE_ERRNO_DIAG(aa->fw[0] == NULL, "popen", "random sink");
    } else if (split_output_file) {
        aa->fw = (FILE**) malloc(n*n*sizeof(FILE*));
        for(unsigned int i = 0 ; i < n ; i++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                FILE * f;
                char * str;
                int rc = asprintf(&str, "%s.sols%d-%d.%d-%d",
                        aa->output_file,
                        j, j + 1,
                        i, i + 1);
                ASSERT_ALWAYS(rc >= 0);
                f = fopen(str, aa->ascii ? "w" : "wb");
                DIE_ERRNO_DIAG(f == NULL, "fopen", str);
                /* I doubt the iobuf makes sense in this case */
                /*
                   aa->iobuf = (char*) malloc(2 * io_block_size);
                   setbuffer(f, aa->iobuf, 2 * io_block_size);
                   */
                aa->fw[i*n+j] = f;
            }
        }
    } else {
        aa->fw = (FILE**) malloc(sizeof(FILE*));
        aa->fw[0] = fopen(aa->output_file, aa->ascii ? "w" : "wb");
        DIE_ERRNO_DIAG(aa->fw[0] == NULL, "fopen", aa->output_file);
        aa->iobuf = (char*) malloc(2 * io_block_size);
        setbuffer(aa->fw[0], aa->iobuf, 2 * io_block_size);
    }
}/*}}}*/

void bm_io_end_write(bm_io_ptr aa)/*{{{*/
{
    int rank;
    MPI_Comm_rank(aa->bm->com[0], &rank);
    if (rank) return;
    matpoly_clear(aa->bm->d->ab, aa->F);
    if (random_input_length) {
        pclose(aa->fw[0]);
    } else if (split_output_file) {
        for(unsigned int i = 0 ; i < aa->bm->d->n ; i++) {
            for(unsigned int j = 0 ; j < aa->bm->d->n ; j++) {
                fclose(aa->fw[i*aa->bm->d->n+j]);
            }
        }
    } else {
        fclose(aa->fw[0]);
    }
    free(aa->fw);
    aa->fw = NULL;
    if (aa->iobuf) free(aa->iobuf);
    aa->iobuf = 0;
}/*}}}*/

void bm_io_clear(bm_io_ptr aa)/*{{{*/
{
    if (aa->fdesc) free(aa->fdesc);
    memset(aa, 0, sizeof(*aa));
}/*}}}*/

void bm_io_guess_length(bm_io_ptr aa)/*{{{*/
{
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    abdst_field ab = d->ab;
    int rank;
    MPI_Comm_rank(aa->bm->com[0], &rank);

    if (random_input_length) {
        aa->guessed_length = random_input_length;
        return;
    }

    if (!rank) {
        struct stat sbuf[1];
        int rc = fstat(fileno(aa->fr[0]), sbuf);
        if (rc < 0) {
            perror(aa->input_file);
            exit(EXIT_FAILURE);
        }

        size_t filesize = sbuf->st_size;

        if (!aa->ascii) {
            size_t avg = avg_matsize(ab, m, n, aa->ascii);
            if (filesize % avg) {
                fprintf(stderr, "File %s has %zu bytes, while its size should be amultiple of %zu bytes (assuming binary input; perhaps --ascii is missing ?).\n", aa->input_file, filesize, avg);
                exit(EXIT_FAILURE);
            }
            aa->guessed_length = filesize / avg;
        } else {
            double avg = avg_matsize(ab, m, n, aa->ascii);
            double expected_length = filesize / avg;
            if (!rank)
                printf("Expect roughly %.2f items in the sequence.\n", expected_length);

            /* First coefficient is always lighter, so we add a +1. */
            aa->guessed_length = 1 + expected_length;
        }
#if 0
        /* we don't have the struct bw at hand here... */
        if (bw->end || bw->start) {
            printf("(Note: from bw parameters, we expect %u).\n",
                    bw->end - bw->start);
        }
        printf(".\n");
#endif
    }
    MPI_Bcast(&(aa->guessed_length), 1, MPI_UNSIGNED, 0, aa->bm->com[0]);
}/*}}}*/

void bm_io_compute_initial_F(bm_io_ptr aa) /*{{{ */
{
    bmstatus_ptr bm = aa->bm;
    dims *d = bm->d;
    abdst_field ab = d->ab;
    unsigned int m = d->m;
    unsigned int n = d->n;
    matpoly_ptr A = aa->A;

    unsigned int (*fdesc)[2] = (unsigned int(*)[2])malloc(2 * m * sizeof(unsigned int));

    int rank;
    MPI_Comm_rank(aa->bm->com[0], &rank);
    if (!rank) {
        /* read the first few coefficients. Expand A accordingly as we are
         * doing the read */

        ASSERT(A->m == m);
        ASSERT(A->n == n);

        abelt tmp;
        abinit(ab, &tmp);

        /* First try to create the initial F matrix */
        printf("Computing t0\n");

        /* We want to create a full rank m*m matrix M, by extracting columns
         * from the first coefficients of A */

        matpoly M;
        matpoly_init(ab, M, m, m, 1);
        M->size = 1;

        /* For each integer i between 0 and m-1, we have a column, picked
         * from column cnum[i] of coeff exponent[i] of A which, once reduced modulo
         * the other ones, has coefficient at row pivots[i] unequal to zero.
         */
        unsigned int *pivots = (unsigned int*) malloc(m * sizeof(unsigned int));
        unsigned int *exponents = (unsigned int*) malloc(m * sizeof(unsigned int));
        unsigned int *cnum = (unsigned int*) malloc(m * sizeof(unsigned int));
        unsigned int r = 0;

        for (unsigned int k = 0; r < m ; k++) {
            /* read a new coefficient */
            bm_io_read1(aa, 0);

            for (unsigned int j = 0; r < m && j < n; j++) {
                /* Extract a full column into M (column j, degree k in A) */
                /* adjust the coefficient degree to take into account the
                 * fact that for SM columns, we might in fact be
                 * interested by the _previous_ coefficient */
                matpoly_extract_column(ab, M, r, 0, A, j, k + (j >= bm->d->nrhs));

                /* Now reduce it modulo all other columns */
                for (unsigned int v = 0; v < r; v++) {
                    unsigned int u = pivots[v];
                    /* the v-th column in the M is known to
                     * kill coefficient u (more exactly, to have a -1 as u-th
                     * coefficient, and zeroes for the other coefficients
                     * referenced in the pivots[0] to pivots[v-1] indices).
                     */
                    /* add M[u,r]*column v of M to column r of M */
                    for(unsigned int i = 0 ; i < m ; i++) {
                        if (i == u) continue;
                        abmul(ab, tmp,
                                  matpoly_coeff(ab, M, i, v, 0),
                                  matpoly_coeff(ab, M, u, r, 0));
                        abadd(ab, matpoly_coeff(ab, M, i, r, 0),
                                  matpoly_coeff(ab, M, i, r, 0),
                                  tmp);
                    }
                    abset_zero(ab,
                            matpoly_coeff(ab, M, u, r, 0));
                }
                unsigned int u = 0;
                for( ; u < m ; u++) {
                    if (abcmp_ui(ab, matpoly_coeff(ab, M, u, r, 0), 0) != 0)
                        break;
                }
                if (u == m) {
                    printf("[X^%d] A, col %d does not increase rank (still %d)\n",
                           k + (j >= bm->d->nrhs), j, r);

                    /* we need at least m columns to get as starting matrix
                     * with full rank. Given that we have n columns per
                     * coefficient, this means at least m/n matrices.
                     */

                    if (k * n > m + 40) {
                        printf("The choice of starting vectors was bad. "
                               "Cannot find %u independent cols within A\n", m);
                        exit(EXIT_FAILURE);
                    }
                    continue;
                }

                /* Bingo, it's a new independent col. */
                pivots[r] = u;
                cnum[r] = j;
                exponents[r] = k;

                /* Multiply the column so that the pivot becomes -1 */
                int rc = abinv(ab, tmp, matpoly_coeff(ab, M, u, r, 0));
                if (!rc) {
                    fprintf(stderr, "Error, found a factor of the modulus: ");
                    abfprint(ab, stderr, tmp);
                    fprintf(stderr, "\n");
                    exit(EXIT_FAILURE);
                }
                abneg(ab, tmp, tmp);
                for(unsigned int i = 0 ; i < m ; i++) {
                    abmul(ab, matpoly_coeff(ab, M, i, r, 0),
                              matpoly_coeff(ab, M, i, r, 0),
                              tmp);
                }

                r++;

                // if (r == m)
                    printf
                        ("[X^%d] A, col %d increases rank to %d (head row %d)%s\n",
                         k + (j >= bm->d->nrhs), j, r, u,
                         (j < bm->d->nrhs) ? " (column not shifted because of the RHS)":"");
            }
        }

        if (r != m) {
            printf("This amount of data is insufficient. "
                   "Cannot find %u independent cols within A\n", m);
            exit(EXIT_FAILURE);
        }

        aa->t0 = exponents[r - 1] + 1;
        /* We always have one extra coefficient of back log */
        ASSERT_ALWAYS(aa->t0 == aa->k - 1);

        printf("Found satisfying init data for t0=%d\n", aa->t0);

        /* We've also got some adjustments to make: room for one extra
         * coefficient is needed in A. Reading of further coefficients will
         * pickup where they stopped, and will always leave the last t0+2
         * coefficients readable. */
        matpoly_realloc(ab, A, aa->t0 + 2);
        A->size++;


        for(unsigned int j = 0 ; j < m ; j++) {
            fdesc[j][0] = exponents[j];
            fdesc[j][1] = cnum[j];
            ASSERT_ALWAYS(exponents[j] < aa->t0);
        }
        free(pivots);
        free(exponents);
        free(cnum);
        matpoly_clear(ab, M);
        abclear(ab, &tmp);
    }
    aa->fdesc = fdesc;
    MPI_Bcast(aa->fdesc, 2*m, MPI_UNSIGNED, 0, aa->bm->com[0]);
    MPI_Bcast(&(aa->t0), 1, MPI_UNSIGNED, 0, aa->bm->com[0]);
    bm->t = aa->t0;
}				/*}}} */

template<class Writer>
void bm_io_compute_E(bm_io_ptr aa, Writer& E, unsigned int expected, unsigned int allocated)/*{{{*/
{
    // F0 is exactly the n x n identity matrix, plus the
    // X^(s-exponent)e_{cnum} vectors. fdesc has the (exponent, cnum)
    // pairs
    bmstatus_ptr bm = aa->bm;
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    abdst_field ab = d->ab;
    int rank;
    MPI_Comm_rank(aa->bm->com[0], &rank);

    unsigned int window = aa->t0 + 2;

    double tt0 = wct_seconds();
    double next_report_t = tt0 + 10;
    unsigned next_report_k = expected / 100;

    int over = 0;
    unsigned int final_length = UINT_MAX;

    for(unsigned int kE = 0 ; kE + aa->t0 < allocated ; kE ++) {

        if (kE % E.chunk_size() || kE + aa->t0 >= expected) {
            MPI_Bcast(&over, 1, MPI_INT, 0, bm->com[0]);
            if (over) break;
        }
        if (rank == 0 && !over) {
            over = !bm_io_read1(aa, window);
            if (over) {
                printf("EOF met after reading %u coefficients\n", aa->k);
                final_length = kE;
            }
        }

        E.coeff(0, 0, kE);

        if (over || rank)
            continue;

        if (kE + aa->t0 > allocated) {
            fprintf(stderr, "Going way past guessed length%s ???\n", aa->ascii ? " (more than 5%%)" : "");
        }

        for(unsigned int j = 0 ; j < n ; j++) {
            /* If the first columns of F are the identity matrix, then
             * in E we get data from coefficient kE+t0 in A. More
             * generally, if it's x^q*identity, we read
             * coeficient of index kE + t0 - q.
             *
             * Note that we prefer to take q=0 anyway, since a
             * choice like q=t0 would create duplicate rows in E,
             * and that would be bad.
             */
            unsigned int kA = kE + aa->t0 + (j >= bm->d->nrhs);
            matpoly_extract_column(ab, E, j, kE, aa->A, j, kA % window);
        }

        for(unsigned int jE = n ; jE < m + n ; jE++) {
            unsigned int jA = aa->fdesc[jE-n][1];
            unsigned int offset = aa->fdesc[jE-n][0];
            unsigned int kA = kE + offset + (jA >= bm->d->nrhs);
            matpoly_extract_column(ab, E, jE, kE, aa->A, jA, kA % window);
        }
        /* This is because aa->k integrates some backlog because of
         * the SM / non-SM distinction (for DL) */
        ASSERT_ALWAYS(kE + 1 + aa->t0 + 1 == aa->k);
        if (aa->k > next_report_k) {
            double tt = wct_seconds();
            if (tt > next_report_t) {
                printf(
                        "Read %u coefficients (%.1f%%)"
                        " in %.1f s (%.1f MB/s)\n",
                        aa->k, 100.0 * aa->k / expected,
                        tt-tt0, aa->k * avg_matsize(ab, m, n, aa->ascii) / (tt-tt0)/1.0e6);
                next_report_t = tt + 10;
                next_report_k = aa->k + expected / 100;
            }
        }
    }
    MPI_Bcast(&final_length, 1, MPI_UNSIGNED, 0, bm->com[0]);
    ASSERT_ALWAYS(final_length != UINT_MAX);
    E.finalize(final_length);
}/*}}}*/

void usage()
{
    fprintf(stderr, "Usage: ./plingen [options, to be documented]\n");
    fprintf(stderr,
            "General information about bwc options is in the README file\n");
    exit(EXIT_FAILURE);
}

void bmstatus_init(bmstatus_ptr bm, unsigned int m, unsigned int n)/*{{{*/
{
    /* Fields will not all be significant right now. The initialization
     * is a long process */
    memset(bm, 0, sizeof(bmstatus));
    bm->d->m = m;
    bm->d->n = n;
    bm->lucky = (int*) malloc((m + n) * sizeof(int));
    memset(bm->lucky, 0, (m + n) * sizeof(int));
}/*}}}*/

void bmstatus_clear(bmstatus_ptr bm)/*{{{*/
{
    free(bm->lucky);
    memset(bm, 0, sizeof(bmstatus));
}/*}}}*/


unsigned int count_lucky_columns(bmstatus_ptr bm)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int b = m + n;
    int luck_mini = expected_pi_length(d, 0);
    MPI_Bcast(bm->lucky, b, MPI_UNSIGNED, 0, bm->com[0]);
    unsigned int nlucky = 0;
    for(unsigned int j = 0 ; j < b ; nlucky += bm->lucky[j++] >= luck_mini) ;
    return nlucky;
}/*}}}*/

int check_luck_condition(bmstatus_ptr bm)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;
    unsigned int nlucky = count_lucky_columns(bm);

    int rank;
    MPI_Comm_rank(bm->com[0], &rank);

    if (!rank) {
        printf("Number of lucky columns: %u (%u wanted)\n", nlucky, n);
    }

    if (nlucky == n)
        return 1;

    if (!rank) {
        fprintf(stderr, "Could not find the required set of solutions (nlucky=%u)\n", nlucky);
    }
    if (random_input_length) {
        static int once=0;
        if (once++) {
            if (!rank) {
                fprintf(stderr, "Solution-faking loop crashed\n");
            }
            MPI_Abort(bm->com[0], EXIT_FAILURE);
        }
        if (!rank) {
            printf("Random input: faking successful computation\n");
        }
        for(unsigned int j = 0 ; j < n ; j++) {
            bm->lucky[(j * 1009) % (m+n)] = expected_pi_length(d, 0);
        }
        return check_luck_condition(bm);
    }

    return 0;
}/*}}}*/

void display_deltas(bmstatus_ptr bm, unsigned int * delta)/*{{{*/
{
    dims * d = bm->d;
    unsigned int m = d->m;
    unsigned int n = d->n;

    int rank;
    MPI_Comm_rank(bm->com[0], &rank);

    if (!rank) {
        printf("Final, t=%u: delta =", bm->t);
        for(unsigned int j = 0; j < m + n; j++) {
            printf(" %u", delta[j]);
            if (bm->lucky[j] < 0) {
                printf("(*)");
            }
        }
        printf("\n");
    }
}/*}}}*/

void print_node_assignment(MPI_Comm comm)
{
    int rank;
    int size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    struct utsname me[1];
    int rc = uname(me);
    if (rc < 0) { perror("uname"); MPI_Abort(comm, 1); }
    size_t sz = 1 + sizeof(me->nodename);
    char * global = (char*) malloc(size * sz);
    memset(global, 0, size * sz);
    memcpy(global + rank * sz, me->nodename, sizeof(me->nodename));

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            global, sz, MPI_BYTE, comm);
    if (rank == 0) {
        char name[80];
        int len=80;
        MPI_Comm_get_name(comm, name, &len);
        name[79]=0;
        for(int i = 0 ; i < size ; i++) {
            printf("%s rank %d: %s\n", name, i, global + i * sz);
        }
    }
    free(global);
}


int main(int argc, char *argv[])
{
    bmstatus bm;
    dims * d = bm->d;

    param_list pl;

    bw_common_init(bw, &argc, &argv);
    param_list_init(pl);

    bw_common_decl_usage(pl);
    plingen_decl_usage(pl);
    logline_decl_usage(pl);

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);
    /* {{{ interpret our parameters */
    gmp_randinit_default(rstate);

    int rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    logline_init_timer();

    param_list_parse_int(pl, "allow_zero_on_rhs", &allow_zero_on_rhs);
    param_list_parse_uint(pl, "random-input-with-length", &random_input_length);
    param_list_parse_int(pl, "split-output-file", &split_output_file);
    param_list_parse_int(pl, "split-input-file", &split_input_file);
    param_list_parse_int(pl, "caching", &caching);

    const char * afile = param_list_lookup_string(pl, "afile");

    if (bw->m == -1) {
	fprintf(stderr, "no m value set\n");
	exit(EXIT_FAILURE);
    }
    if (bw->n == -1) {
	fprintf(stderr, "no n value set\n");
	exit(EXIT_FAILURE);
    }
    if (!global_flag_tune && !(afile || random_input_length)) {
        fprintf(stderr, "No afile provided\n");
        exit(EXIT_FAILURE);
    }

    /* we allow ffile and ffile to be both NULL */
    const char * tmp = param_list_lookup_string(pl, "ffile");
    char * ffile = NULL;
    if (tmp) {
        ffile = strdup(tmp);
    } else if (afile) {
        int rc = asprintf(&ffile, "%s.gen", afile);
        ASSERT_ALWAYS(rc >= 0);
    }
    ASSERT_ALWAYS((afile==NULL) == (ffile == NULL));

#ifndef ENABLE_CACHING_MPI_LINGEN
    if (caching) {
        fprintf(stderr, "--caching=1 only supported with ENABLE_CACHING_MPI_LINGEN\n");
        exit(EXIT_FAILURE);
    }
#endif

    bmstatus_init(bm, bw->m, bw->n);

    const char * rhs_name = param_list_lookup_string(pl, "rhs");
    if (!random_input_length) {
        if (!rhs_name) {
            fprintf(stderr, "When using plingen, you must either supply --random-input-with-length, or provide a rhs, or possibly provide rhs=none\n");
        } else if (strcmp(rhs_name, "none") == 0) {
            rhs_name = NULL;
        }
    }
    if ((rhs_name != NULL) && param_list_parse_uint(pl, "nrhs", &(bm->d->nrhs))) {
        fprintf(stderr, "the command line arguments rhs= and nrhs= are incompatible\n");
        exit(EXIT_FAILURE);
    }
    if (rhs_name) {
        if (!rank)
            get_rhs_file_header(rhs_name, NULL, &(bm->d->nrhs), NULL);
        MPI_Bcast(&bm->d->nrhs, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    }

    unsigned int m = d->m;
    unsigned int n = d->n;
    abdst_field ab = d->ab;

    abfield_init(ab);
    abfield_specify(ab, MPFQ_PRIME_MPZ, bw->p);

    bm->lingen_threshold = 10;
    bm->lingen_mpi_threshold = 1000;
    param_list_parse_uint(pl, "lingen_threshold", &(bm->lingen_threshold));
    param_list_parse_uint(pl, "display-threshold", &(display_threshold));
#ifdef HAVE_MPIR
    param_list_parse_uint(pl, "caching-threshold", &(caching_threshold));
#endif
    param_list_parse_uint(pl, "lingen_mpi_threshold", &(bm->lingen_mpi_threshold));
    param_list_parse_uint(pl, "io-block-size", &(io_block_size));
    gmp_randseed_ui(rstate, bw->seed);
    if (bm->lingen_mpi_threshold < bm->lingen_threshold) {
        bm->lingen_mpi_threshold = bm->lingen_threshold;
        fprintf(stderr, "Argument fixing: setting lingen_mpi_threshold=%u (because lingen_threshold=%u)\n",
                bm->lingen_mpi_threshold, bm->lingen_threshold);
    }
    checkpoint_directory = param_list_lookup_string(pl, "checkpoint-directory");
    param_list_parse_uint(pl, "checkpoint-threshold", &checkpoint_threshold);
    param_list_parse_int(pl, "save_gathered_checkpoints", &save_gathered_checkpoints);



#if defined(FAKEMPI_H_)
    bm->lingen_mpi_threshold = UINT_MAX;
#endif

    /* }}} */

    /* {{{ Parse MPI args. Make bm->com[0] a better mpi communicator */
    bm->mpi_dims[0] = 1;
    bm->mpi_dims[1] = 1;
    param_list_parse_intxint(pl, "mpi", bm->mpi_dims);
    {
        /* Display node index wrt MPI_COMM_WORLD */
        print_node_assignment(MPI_COMM_WORLD);

        /* Reorder all mpi nodes so that each node gets the given number
         * of jobs, but close together.
         */
        int mpi[2] = { bm->mpi_dims[0], bm->mpi_dims[1], };
        int thr[2] = {1,1};
        param_list_parse_intxint(pl, "thr", thr);

#ifdef  FAKEMPI_H_
        if (mpi[0]*mpi[1] > 1 || thr[0]*thr[1] > 1) {
            fprintf(stderr, "non-trivial option mpi= or thr= can't be used with fakempi. Please do an MPI-enabled build (MPI=1)\n");
            exit(EXIT_FAILURE);
        }
#endif
        /* Asssume numbering in MPI_COMM_WORLD is all mpi jobs from the
         * same node together. So we pick them by bunches of size
         * thr[0]*thr[1].
         */
        if (!rank)
            printf("size=%d mpi=%dx%d thr=%dx%d\n", size, mpi[0], mpi[1], thr[0], thr[1]);
        ASSERT_ALWAYS(size == mpi[0] * mpi[1] * thr[0] * thr[1]);
        /* The mpi and thr command line argument lead to the same number
         * of working processes than with krylov. Except that here, we're
         * not really doing threads, but running real mpi jobs with their
         * own address space (hence mpirun will need a bigger -n
         * argument).  In a sense, here we strive to follow the semantics
         * used in bwc otherwise.  In fact, we should probably use the
         * pi_comm structures as well, so that we can say that we
         * consistently use the same infrastructure.  For simplicity of
         * the code, we don't.
         */
        bm->mpi_dims[0] *= thr[0];
        bm->mpi_dims[1] *= thr[1];
        if (bm->mpi_dims[0] != bm->mpi_dims[1]) {
            if (!rank)
                fprintf(stderr, "The current plingen code is limited to square splits ; here, we received a %d x %d split, which will not work\n",
                    bm->mpi_dims[0], bm->mpi_dims[1]);
            abort();
        } else if ((m % bm->mpi_dims[0] != 0) || (n % bm->mpi_dims[0] != 0)) {
            if (!rank)
                fprintf(stderr, "The process grid dimensions must divide gcd(m,n)\n");
            abort();
        }
        int tl_grank = rank % (thr[0] * thr[1]); // thread-level global rank
        int tl_irank = tl_grank / thr[1];
        int tl_jrank = tl_grank % thr[1];
        int ml_grank = rank / (thr[0] * thr[1]); // thread-level global rank
        int ml_irank = ml_grank / mpi[1];
        int ml_jrank = ml_grank % mpi[1];
        int irank = ml_irank * thr[0] + tl_irank;
        int jrank = ml_jrank * thr[1] + tl_jrank;
        int newrank = irank * mpi[1] * thr[1] + jrank;

        MPI_Comm_split(MPI_COMM_WORLD, 0, newrank, &(bm->com[0]));
        /* MPI Api has some very deprecated prototypes */
        MPI_Comm_set_name(bm->com[0], (char*) "world");

        char commname[12];
        snprintf(commname, 12, "row%d\n", irank);
        MPI_Comm_split(MPI_COMM_WORLD, irank, jrank, &(bm->com[1]));
        MPI_Comm_set_name(bm->com[1], commname);

        snprintf(commname, 12, "col%d\n", jrank);
        MPI_Comm_split(MPI_COMM_WORLD, jrank, irank, &(bm->com[2]));
        MPI_Comm_set_name(bm->com[2], commname);

        print_node_assignment(bm->com[0]);
    }
    /* }}} */


    /* plingen tuning accepts some arguments. We look them up so as to
     * avoid failures down the line */
    plingen_tuning_lookup_parameters(pl);
    
    logline_interpret_parameters(pl);

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    if (global_flag_tune) {
        plingen_tuning(bm->d->ab, bm->d->m, bm->d->n, bm->com[0], pl);
        MPI_Finalize();
        return 0;
    }


    /* We now have a protected structure for a_reading task which does
     * the right thing concerning parallelism among MPI nodes (meaning
     * that non-root nodes essentially do nothing while the master job
     * does the I/O stuff) */
    bm_io aa;
    bm_io_init(aa, bm, afile, ffile, global_flag_ascii);
    bm_io_begin_read(aa);
    bm_io_guess_length(aa);

    bm_io_compute_initial_F(aa);

    /* this is somewhat ugly, too */
    unsigned int t0 = bm->t;
    unsigned int * delta = (unsigned int*) malloc((m + n) * sizeof(unsigned int));
    for(unsigned int j = 0 ; j < m + n ; delta[j++]=t0);

    int go_mpi = aa->guessed_length >= bm->lingen_mpi_threshold;

    if (go_mpi && !rank) {
        if (size) {
            printf("Expected length %u exceeds MPI threshold %u, going MPI now.\n", aa->guessed_length, bm->lingen_mpi_threshold);
        } else {
            printf("Expected length %u exceeds MPI threshold %u, but the process is not running  in an MPI context.\n", aa->guessed_length, bm->lingen_mpi_threshold);
        }
    }

    if (size > 1) {
#ifdef  ENABLE_MPI_LINGEN
        unsigned int m = aa->bm->d->m;
        unsigned int n = aa->bm->d->n;
        unsigned int b = m + n;
        unsigned int guess = aa->guessed_length;
        size_t safe_guess = aa->ascii ? ceil(1.05 * guess) : guess;
        bigmatpoly model;
        bigmatpoly xpi, xE;

        bigmatpoly_init_model(model, bm->com, bm->mpi_dims[0], bm->mpi_dims[1]);
        bigmatpoly_init(ab, xE, model, 0, 0, 0);
        bigmatpoly_init(ab, xpi, model, 0, 0, 0);   /* pre-init for now */
        bigmatpoly_finish_init(ab, xE, m, b, safe_guess);

        bigmatpoly_write_task writer(aa, xE);

        bm_io_compute_E(aa, writer, guess, safe_guess);

        bm_io_end_read(aa);

        bw_biglingen_collective(bm, xpi, xE, delta);

        display_deltas(bm, delta);
        if (!rank) printf("(pi->alloc = %zu)\n", bigmatpoly_my_cell(xpi)->alloc);

        if (check_luck_condition(bm)) {
            bm_io_begin_write(aa);
            bm_io_set_write_behind_size(aa, delta);
            bigmatpoly_read_task xpi_reader(aa, xpi);
            bm_io_compute_final_F(aa, xpi_reader, delta);
            bm_io_end_write(aa);
        }

        /* clear everything */
        bigmatpoly_clear(ab, xE);
        bigmatpoly_clear(ab, xpi);
        bigmatpoly_clear_model(model);
#else
        /* The ENABLE_MPI_LINGEN flag should be turned on for a proper
         * MPI run.
         */
        ASSERT_ALWAYS(0);
#endif
    } else {
        unsigned int m = aa->bm->d->m;
        unsigned int n = aa->bm->d->n;
        unsigned int b = m + n;
        unsigned int guess = aa->guessed_length;
        size_t safe_guess = aa->ascii ? ceil(1.05 * guess) : guess;

        matpoly pi, E;


        matpoly_init(ab, E, m, b, safe_guess);
        matpoly_init(ab, pi, 0, 0, 0);  /* pre-init for now */


        matpoly_write_task E_writer(aa, E);

        bm_io_compute_E(aa, E_writer, guess, safe_guess);

        bm_io_end_read(aa);

        bw_lingen_single(bm, pi, E, delta);

        display_deltas(bm, delta);
        if (!rank) printf("(pi->alloc = %zu)\n", pi->alloc);

        if (check_luck_condition(bm)) {
            bm_io_begin_write(aa);
            bm_io_set_write_behind_size(aa, delta);
            matpoly_read_task pi_reader(aa, pi);
            bm_io_compute_final_F(aa, pi_reader, delta);
            bm_io_end_write(aa);
        }

        /* clear everything */
        matpoly_clear(ab, E);
        matpoly_clear(ab, pi);

    }



    if (!rank) {
        printf("t_basecase = %.2f\n", bm->t_basecase);
        printf("t_mp = %.2f\n", bm->t_mp);
        printf("t_mul = %.2f\n", bm->t_mul);
        printf("t_cp_io = %.2f\n", bm->t_cp_io);
    }

    bm_io_clear(aa);
    free(delta);
    abfield_clear(ab);
    bmstatus_clear(bm);
    if (ffile) free(ffile);

    gmp_randclear(rstate);

    param_list_clear(pl);
    bw_common_clear(bw);

    return rank0_exit_code;
}

/* vim:set sw=4 sta et: */
