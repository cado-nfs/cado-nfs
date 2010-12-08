#ifndef MATMUL_TOP_H_
#define MATMUL_TOP_H_

#include <stdint.h>

#include "abase.h"
#include "select_mpi.h"
#include "intersections.h"
#include "parallelizing_info.h"
#include "matmul.h"
#include "params.h"
#include "balancing.h"

/* Don't touch this. */
#define CONJUGATED_PERMUTATIONS

/* A wiring is the information concerning our matrix in one of its two
 * dimensions. Everything is made so that inverting the direction flag
 * should make it possible to compute with the transposed matrix easily.
 * 
 * Wiring number 0 is ``horizontal''. In detail, what this means is:
 *
 * n is the total number of rows ; the total number of instances of this
 * ``wiring'', so to say.
 * v is the LEFT vector (result of matrix-times-vector, or input of
 *   vector-times-matrix).
 * fences denote the row indices within the matrix where the different
 *   jobs are split.
 * x denotes the intersection data for communication across rows.
 * i0 denotes the row indices of our local matrix chunk.
 *
 * Horizontal indices i0 to i1 of the LEFT column vector are stored locally.
 * --> v is a chunk of i1-i0 data elements.
 *
 * Communication with the left vector. The intent is to obtain a right
 * vector.
 * --> reduce is done across rows. The destination for row index i is
 *   column index i on conjugation mode(*)
 * XXX TBC.
 *
 * (*) it's somewhat different in non-conjugation mode, since a
 * permutation has to be taken (in which case allreduce is probably
 * preferrable).
 *
 *
 *
 * Note also that we have the pi_wiring to consider as well. Those
 * have nothing to do with the matrix, and are relevant with regard to
 * the parallelization of the work.
 */

// the ``all_v'' field collects all the pointers to the per-thread vector
// values. There are exactly pi->wr[0]->ncores such pointers in
// mmt->wr[0]. In some cases, these pointers may be equal (for source
// vectors, never used as destination), and in some cases not.

struct mmt_vec_s {
    abt * v;
    abt * * all_v;
    // internal use
    int flags;
};
typedef struct mmt_vec_s mmt_vec[1];
typedef struct mmt_vec_s * mmt_vec_ptr;


struct mmt_wiring_s {
    mmt_vec v;
    unsigned int i0;
    unsigned int i1;
#ifdef  CONJUGATED_PERMUTATIONS
    /* This is likely to be relevant only for conjugated matrices */
    struct isect_info * x;
    unsigned int xlen;

    // this second intersection list is used for reduction, where the row
    // threads collectively build the result vector from their respective
    // parts -- deciding where the computed sum eventually goes is
    // determined by the intersection of i0+k/n*(i1-i0)..i0+(k+1)/n*(i1-i0)
    // with the vertical fences.
    //
    // note that the way reduction is performed is likely to cause some
    // cache invalidates unless some fences end up correctly aligned.
    struct isect_info * y;
    unsigned int ylen;
#endif
};
typedef struct mmt_wiring_s mmt_wiring[1];
typedef struct mmt_wiring_s * mmt_wiring_ptr;
typedef struct mmt_wiring_s const * mmt_wiring_srcptr;

/* TODO: It's unhandy to have abase here, since it forces info-file to ba
 * abase-dependent, which it isn't in reality... Perhaps split the interface
 * in two ?
 */
struct matmul_top_data_s {
    abobj_ptr abase;
    parallelizing_info_ptr pi;
    matmul_ptr mm;

    // w <- M v, or v^T <- w^T M would make sense. Of course v^T and w^T
    // are represented in exactly the same way.

    // global stuff. n[] is going to be misleading in any situation. It's
    // the size of the matrix in the given direction. Meaning that n[0]
    // is the horizontal size -- the size of the rows --. This means the
    // number of columns.
    unsigned int n[2];
    // unsigned int ncoeffs_total;

    // local stuff
    // // unsigned int ncoeffs;
    
    // this really ends up within the mm field.
    char * locfile;

    mmt_wiring wr[2];

#ifdef  CONJUGATED_PERMUTATIONS
    // fences[0]: horizontal fences == row indices (horizontal separating line)
    // fences[1]: vertical fences == col indices (vertical separating line)
    // fences[0] contains 1 + pi->wr[1]->totalsize values
    // fences[1] contains 1 + pi->wr[0]->totalsize values
    //
    // it turns out later on that this numbering is slightly orthogonal
    // to the one used in most of the rest of the code, so maybe it's a
    // bad decision. Fortunately, fences[] is rarely used.
    unsigned int * fences[2];
#endif
    balancing bal;
};

/* THREAD_MULTIPLE_VECTOR is when several threads will be writing to the
 * vector simultaneously. This implies having separate areas.
 *
 * TODO: Think about waiving this restriction in the spirit of bucket
 * sieving, maybe */
#define THREAD_MULTIPLE_VECTOR    0
#define THREAD_SHARED_VECTOR    1

typedef struct matmul_top_data_s matmul_top_data[1];
typedef struct matmul_top_data_s * matmul_top_data_ptr;
typedef struct matmul_top_data_s const * matmul_top_data_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

extern void matmul_top_init(matmul_top_data_ptr mmt,
        abobj_ptr abase,
        parallelizing_info_ptr pi,
        int const * flags,
        param_list pl,
        int optimized_direction);


extern void matmul_top_clear(matmul_top_data_ptr mmt, abobj_ptr abase);
#if 0
extern void matmul_top_fill_random_source(matmul_top_data_ptr mmt, int d);
#endif
extern void matmul_top_load_vector(matmul_top_data_ptr mmt, const char * name, int d, unsigned int iter);
extern void matmul_top_save_vector(matmul_top_data_ptr mmt, const char * name, int d, unsigned int iter);
extern void matmul_top_mul_cpu(matmul_top_data_ptr mmt, int d);
extern void matmul_top_mul_comm(matmul_top_data_ptr mmt, int d);
static inline void matmul_top_mul(matmul_top_data_ptr mmt, int d)
{
    matmul_top_mul_cpu(mmt, d);
    matmul_top_mul_comm(mmt, d);
}


/* Now some of the generic interface calls. By design, not everything is
 * possible with these calls. In particular, nothing critical is doable.
 * As a general convention, the only thing these calls need to know about
 * the abase is the stride value, which corresponds to abbytes(abase, 1).
 * Besides that, the generic calls take the vector pointer as argument.
 * Specifying NULL as vector argument is equivalent to taking the default
 * vector defined in the mmt data for the given direction flag. */

/* Same as mmt_vec, but not bound to the given abase. Used for generic
 * calls */
struct mmt_generic_vec_s {
    abase_generic_ptr v;
    abase_generic_ptr * all_v;
    // internal use
    int flags;
};
typedef struct mmt_generic_vec_s mmt_generic_vec[1];
typedef struct mmt_generic_vec_s * mmt_generic_vec_ptr;

extern void matmul_top_vec_init_generic(matmul_top_data_ptr mmt, size_t stride, mmt_generic_vec_ptr v, int d, int flags);
extern void matmul_top_vec_clear_generic(matmul_top_data_ptr mmt, size_t stride MAYBE_UNUSED, mmt_generic_vec_ptr v, int d);
extern void matmul_top_fill_random_source_generic(matmul_top_data_ptr mmt, size_t stride, mmt_generic_vec_ptr v, int d);
extern void matmul_top_load_vector_generic(matmul_top_data_ptr mmt, size_t stride, mmt_generic_vec_ptr v, const char * name, int d, unsigned int iter);
extern void matmul_top_save_vector_generic(matmul_top_data_ptr mmt, size_t stride, mmt_generic_vec_ptr v, const char * name, int d, unsigned int iter);

/* These two do not really belong here, but comes as a useful complement
 */
extern void vec_init_generic(pi_wiring_ptr, size_t, mmt_generic_vec_ptr, int, unsigned int);
extern void vec_clear_generic(pi_wiring_ptr, size_t, mmt_generic_vec_ptr, unsigned int);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_TOP_H_ */
