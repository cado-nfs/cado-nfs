#ifndef BW_COMMON_H_
#define BW_COMMON_H_

/* This file gathers the parameters that are relevant to all BW
 * programs -- most programs parse the corresponding values -- for some,
 * the values are unused, which is not a lot to worry about.
 */

#include <gmp.h>
#include <stdio.h>
#include <stdint.h>  // for uint32_t
#include "params.h"

/* relevant to secure.c as well ; really, this is almost a hard-coded
 * constant.  Don't imagine it's easily tunable. If one really insists,
 * perhaps a b128/secure would work, but this has not been checked. */
#define NCHECKS_CHECK_VECTOR_GF2    64
#define NCHECKS_CHECK_VECTOR_GFp    1

#define MAX_NUMBER_OF_CHECK_STOPS       32

struct bw_params {

    /* m,n blocking factors as in the textbook BW description */
    int m,n;

    /* We have 0<=nrhs<=n. This corresponds to solving inhomogeneous
     * systems, but in reality we consider this parameter only when
     * related to SM handling in the DL computation, when the SM blocks
     * are part of the block Wiedemann starting vectors. This influences
     * the number of columns which are shifted in the generator
     * computation ("considering A(X) div X"), so that we force the
     * computed relation to have non-zero coefficients for the SM columns
     * (hence, not shifted in A(X)), while those corresponding to random
     * vectors _are_ shifted.
     *
     * mksol and gather are also affected by this parameter.
     *
     * By default we have nrhs=0
     */
    // we make this a phantom parameter. It's deduced from reading the
    // rhs file, and we have no compelling need for allotting room for it
    // in the global struct, given its rare use.
    // int nrhs;

    /* modulus */
    mpz_t p;

    /* _at the moment_ this corresponds to the checking & checkpointing
     * interval. Quite clearly, the two could be separated.
     */
    int interval;

    /* defined, but unused */
    int verbose;

    /* Whether the current job/thread may print to stdout */
    int can_print;

    /* This indicates the starting iteration -- only for krylov and mksol */
    int start;
    int end;

    /*
     * relevant for krylov mostly, but not only.
     *
     * for krylov: indicates the coordinate range in [0..n[ that is
     * relevant for us. ys[1]-ys[0] defines the ``local'' blocking factor
     * n'.
     *
     * for dispatch: only the difference ys[1]-ys[0] is used, and fed as
     * a parameter to the cache building routines (those have an interest
     * in knowing the width of the vectors which will be worked on.
     */
    int ys[2];

    /* for mksol and gather, indicate which sub-range of the solution
     * space we intend to compute. This must be aligned on multiples of
     * the quantity ys[1]-ys[0] derived from the supplied values ys.
     */
    unsigned int solutions[2];

    /* dir is a boolean flag equal to 1 if we are looking for the right
     * nullspace of the matrix. In matmul_top speak, it indicates where the
     * source vector is.
     */
    int dir;

    /* secure, prep, and lingen are not deterministic. */
    int seed;

    /* If enabled, do not check the intermediate results while running
     * (the possibility of doing the checks offline still exists,
     * though). This makes it possible to run krylov without having run
     * the secure program first.
     */
    int skip_online_checks;

    /* If keep_rolling_checkpoints is defined to a positive integer X,
     * keep only the last X vector checkpoints, and remove the others.
     * Otherwise (if the parameter is zero), all checkpoints are kept.
     * See also next parameter.
     *
     * Note that having skip_online_checks and keep_rolling_checkpoints
     * at the same time is dangerous.
     */
    int keep_rolling_checkpoints;

    /* This moderates the previous behaviour. If a checkpoint to be
     * discarded as per the previous flag turns out to have been
     * modified less than X seconds ago, keep it anyway.
     */
    int keep_checkpoints_younger_than;

    /* In case the previous flag is enabled, still keep all the
     * checkpoints which are multiple of the length given here.
     */
    int checkpoint_precious;

    /* The check_stops[] array contains the number of known check
     * vectors, against which data is checked while running.
     * Note that krylov and mksol store only one vector for the check
     * stops, and do not afford storing many. This is of course to save
     * memory. The secure program computes all requested check vectors.
     * The yet-to-be-written offline check program checks against
     * everything it can.
     *
     * The interval specified by the interval parameter is automatically
     * counted as a check stop, unless skip_online_checks is true.
     */
    int number_of_check_stops;
    int check_stops[MAX_NUMBER_OF_CHECK_STOPS];

    int original_argc;
    char ** original_argv;

    double wct_base;
};

extern struct bw_params bw[1];
extern const char * bw_dirtext[];

#ifdef __cplusplus
extern "C" {
#endif

/* Typical use pattern:
     
    bw_common_init(bw, &argc, &argv);
    param_list_init(pl);

    bw_common_decl_usage(pl);
    // more decl_usage functions.

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);
    // bw_common_cheat_parameter_lookup(pl);
    // more interpret_parameters functions.

    param_list_warn_unused(pl);
    param_list_clear(pl);

    // program !

    bw_common_clear(bw);

*/

void bw_common_decl_usage(param_list);
void bw_common_parse_cmdline(struct bw_params * bw, param_list pl, int * p_argc, char *** p_argv);
void bw_common_interpret_parameters(struct bw_params * bw, param_list pl);
int bw_common_init(struct bw_params * bw, int * p_argc, char *** p_argv);
int bw_common_clear(struct bw_params * bw);


/* utility function. */
extern int get_rhs_file_header(const char * filename, uint32_t * p_nrows, unsigned int * p_nrhs, mpz_ptr p_p);
extern int get_rhs_file_header_stream(FILE * f, uint32_t * p_nrows, unsigned int * p_nrhs, mpz_ptr p_p);

extern unsigned int bw_set_length_and_interval_krylov(struct bw_params * bw, unsigned int dims[2]);
extern unsigned int bw_set_length_and_interval_mksol(struct bw_params * bw, unsigned int dims[2]);
#ifdef __cplusplus
}
#endif

#endif	/* BW_COMMON_H_ */
