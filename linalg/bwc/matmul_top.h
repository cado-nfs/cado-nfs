#ifndef MATMUL_TOP_H_
#define MATMUL_TOP_H_

#include <stddef.h>
#include <stdint.h>

#include "select_mpi.h"
#include "intersections.h"
#include "parallelizing_info.h"
#include "matmul.h"
#include "params.h"
#include "balancing.h"
#include "misc.h"
#include "mpfq/mpfq_vbase.h"
#include "async.h"

// the ``all_v'' field collects all the pointers to the per-thread vector
// values. There are exactly pi->wr[0]->ncores such pointers in
// mmt->wr[0]. In some cases, these pointers may be equal (for source
// vectors, never used as destination), and in some cases not.

struct mmt_vec_s;
typedef struct mmt_vec_s * mmt_vec_ptr;

struct mmt_vec_s {
    mpfq_vbase_ptr abase;
    parallelizing_info_ptr pi;
    int d;
    pi_datatype_ptr pitype;
    void * v;
    mmt_vec_ptr * siblings;     /* pi->wr[d]->ncores siblings ;
                                   only in case all cores in the communicator
                                   have their own data area v */
    mmt_vec_ptr * mpals;        /* pi->m->ncores siblings, always */
    mmt_vec_ptr * wrpals[2];    /* pi->wr[0]->ncores and pi->wr[1]->ncores
                                   siblings, always. */
    unsigned int n;             // total size in items
    unsigned int i0;
    unsigned int i1;
    size_t rsbuf_size;          // auto-expanded on demand.
    void * rsbuf[2];            // only for RS_CHOICE == RS_CHOICE_MINE
    int consistency;    /* 0 == inconsistent ; 1 == partial ; 2 == full */
};
typedef struct mmt_vec_s mmt_vec[1];

/* some handy macros to access portions of mmt_vec's */
#define SUBVEC(v,w,offset) (v)->abase->vec_subvec(v->abase, (v)->w, offset)
#define SUBVEC_const(v,w,offset) (v)->abase->vec_subvec_const(v->abase, (v)->w, offset)

/* yikes. yet another list structure. Wish we had the STL. */
struct permutation_data_s {
    size_t n;
    size_t alloc;
    unsigned int (*x)[2];
};
typedef struct permutation_data_s permutation_data[1];
typedef struct permutation_data_s * permutation_data_ptr;
/* all methods are private */

struct matmul_top_matrix_s {
    // global stuff.
    //
    // n[0] is the number of rows, includding padding.
    // n[1] is the number of columns, includding padding.
    //
    // n0[] is without padding.
    //
    // Note that a communicator in direction 0 handles in total a dataset
    // of size n[1] (the number of items in a matrix row is the number of
    // columns.
    unsigned int n[2];
    unsigned int n0[2]; // n0: unpadded.

    // this really ends up within the mm field. It's not a complete file
    // name though. We lack the implementation extension, the possible
    // transposition tag, as well as the .bin extension.
    char * locfile;

    /* These two are global to all threads and jobs (well, each thread
     * has its own pointer though, it's not a shared_malloc. It could be,
     * but it isn't).
     *
     * For random matrices, both strings below are NULL.
     */
    char * mname;
    char * bname;

    balancing bal;

    /* These are local excerpts of the balancing permutation: arrays of
     * pairs (row index in the (sub)-matrix) ==> (row index in the
     * original matrix), but only when both coordinates of this pair are
     * in the current row and column range. This can be viewed as the set
     * of non-zero positions in the permutation matrix if it were split
     * just like the current matrix is. */
    permutation_data_ptr perm[2];       /* rowperm, colperm */

    matmul_ptr mm;
};

typedef struct matmul_top_matrix_s matmul_top_matrix[1];
typedef struct matmul_top_matrix_s * matmul_top_matrix_ptr;
typedef struct matmul_top_matrix_s const * matmul_top_matrix_srcptr;

struct matmul_top_data_s {
    parallelizing_info_ptr pi;
    mpfq_vbase_ptr abase;
    pi_datatype_ptr pitype;
    /* These n[] and n0[] correspond to the dimensions of the product
     *
     * n[0] is matrices[0]->n[0]
     * n[1] is matrices[nmatrices-1]->n[1]
     */
    unsigned int n[2];
    unsigned int n0[2]; // n0: unpadded.
    int nmatrices;
    matmul_top_matrix * matrices;
};

typedef struct matmul_top_data_s matmul_top_data[1];
typedef struct matmul_top_data_s * matmul_top_data_ptr;
typedef struct matmul_top_data_s const * matmul_top_data_srcptr;

/* These are flags for the distributed vectors. For the moment we have
 * only one flag */
#define THREAD_SHARED_VECTOR    1

#ifdef __cplusplus
extern "C" {
#endif

extern void matmul_top_init(matmul_top_data_ptr mmt,
        mpfq_vbase_ptr abase,
        parallelizing_info_ptr pi,
        param_list pl,
        int optimized_direction);


extern void matmul_top_decl_usage(param_list_ptr pl);
extern void matmul_top_lookup_parameters(param_list_ptr pl);
extern void matmul_top_report(matmul_top_data_ptr mmt, double scale);
extern void matmul_top_clear(matmul_top_data_ptr mmt);
extern unsigned int matmul_top_rank_upper_bound(matmul_top_data_ptr mmt);
#if 0
extern void matmul_top_fill_random_source(matmul_top_data_ptr mmt, int d);
#endif
extern size_t mmt_my_own_size_in_items(mmt_vec_ptr v);
extern size_t mmt_my_own_size_in_bytes(mmt_vec_ptr v);
extern size_t mmt_my_own_offset_in_items(mmt_vec_ptr v);
extern size_t mmt_my_own_offset_in_bytes(mmt_vec_ptr v);
extern void * mmt_my_own_subvec(mmt_vec_ptr v);

extern void mmt_vec_set_random_through_file(mmt_vec_ptr v, const char * name, unsigned int itemsondisk, gmp_randstate_t rstate);
/* do not use this function if you want consistency when the splitting
 * changes ! */
extern void mmt_vec_set_random_inconsistent(mmt_vec_ptr v, gmp_randstate_t rstate);
extern void mmt_vec_truncate(matmul_top_data_ptr mmt, mmt_vec_ptr v);
extern void mmt_vec_set_x_indices(mmt_vec_ptr y, uint32_t * gxvecs, int m, unsigned int nx);

extern void matmul_top_mul_cpu(matmul_top_data_ptr mmt, int midx, int d, mmt_vec_ptr w, mmt_vec_ptr v);
extern void matmul_top_comm_bench(matmul_top_data_ptr mmt, int d);
extern void matmul_top_mul_comm(mmt_vec_ptr w, mmt_vec_ptr v);

/* v is both input and output. w is temporary */
extern void matmul_top_mul(matmul_top_data_ptr mmt, mmt_vec *w, struct timing_data * tt);

extern void mmt_vec_init(matmul_top_data_ptr mmt, mpfq_vbase_ptr abase, pi_datatype_ptr pitype, mmt_vec_ptr v, int d, int flags, unsigned int n);
extern void mmt_vec_clear(matmul_top_data_ptr mmt, mmt_vec_ptr v);
extern void mmt_own_vec_set(mmt_vec_ptr w, mmt_vec_ptr v);
extern void mmt_own_vec_set2(mmt_vec_ptr z, mmt_vec_ptr w, mmt_vec_ptr v);
extern void mmt_full_vec_set(mmt_vec_ptr w, mmt_vec_ptr v);
extern void mmt_full_vec_set_zero(mmt_vec_ptr v);
#if 0
extern void matmul_top_fill_random_source_generic(matmul_top_data_ptr mmt, size_t stride, mmt_vec_ptr v, int d);
#endif
extern int mmt_vec_load_stream(pi_file_handle f, mmt_vec_ptr v, unsigned int itemsondisk);
extern int mmt_vec_save_stream(pi_file_handle f, mmt_vec_ptr v, unsigned int itemsondisk);
extern int mmt_vec_load(mmt_vec_ptr v, const char * name, unsigned int itemsondisk);
extern int mmt_vec_save(mmt_vec_ptr v, const char * name, unsigned int itemsondisk);
extern void mmt_vec_reduce_mod_p(mmt_vec_ptr v);
extern void mmt_vec_clear_padding(mmt_vec_ptr v, size_t unpadded, size_t padded);

extern void mmt_vec_broadcast(mmt_vec_ptr v);
extern void mmt_vec_reduce(mmt_vec_ptr w, mmt_vec_ptr v);
extern void mmt_vec_reduce_sameside(mmt_vec_ptr v);
extern void mmt_vec_allreduce(mmt_vec_ptr v);
extern void mmt_vec_twist(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_vec_untwist(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_vec_apply_T(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_vec_unapply_T(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_vec_apply_S(matmul_top_data_ptr mmt, int midx, mmt_vec_ptr y);
extern void mmt_vec_unapply_S(matmul_top_data_ptr mmt, int midx, mmt_vec_ptr y);
extern void mmt_vec_apply_P(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_vec_unapply_P(matmul_top_data_ptr mmt, mmt_vec_ptr y);
extern void mmt_apply_identity(mmt_vec_ptr w, mmt_vec_ptr v);
extern void indices_twist(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_TOP_H_ */
