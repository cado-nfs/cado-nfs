#ifndef RAW_MATRIX_U32_H_
#define RAW_MATRIX_U32_H_

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct matrix_u32_s {
    // input arguments.
    const char * mfile;    // matrix file name
    const char * bfile;    // balancing file name ; NULL will mean auto-detect
    /* The 'transpose' flag just tells that the in-memory representation
     * of the sparse matrix below is given by columns and not by rows.
     * This is sometimes the format that is preferred by some MM layers,
     * and cado-nfs decides to present the local matrices to the MM
     * layers precisely in their preferred format.
     */
    int transpose;
    int withcoeffs;
    // output arguments.
    uint32_t * p;
    size_t size;
};
typedef struct matrix_u32_s matrix_u32[1];
typedef struct matrix_u32_s * matrix_u32_ptr;

/* The interface around matrix_u32 is pretty thin, as this header defines
 * no function...
 */
/* Constructors:
 *
 * It's the entire responsibility of the caller. Structures of type
 * matrix_u32 are to be created empty and zeroed out by the caller. The
 * fields mfile and bfile are also expected from the caller.
 */
/*
 * Initialization:
 *
 * The main initializer function is balancing_get_matrix_u32 defined in
 * balancing_workhorse.c
 * Alternatively, there is also the random_matrix_get_u32 function,
 * defined in random_matrix.c.
 */
/*
 * Destructors:
 *
 * There is no destructor function for this type. The global assumption
 * is that the mm layer which is fed with this structure has the right to
 * claim ownership of the data referenced from this structure. If it does
 * not, then it is this layer's responsibility to release the memory. In
 * a sense, matmul_build_cache *is* the destructor...
 *
 * To illustrate this, the layer mm_basic reuses (and in fact, also
 * modifies) the m->p data, and thus it remains in memory.  In contrast,
 * the bucket layer discards this data, as the post-processed matrix form
 * is more efficient eventually.
 */
 
#ifdef __cplusplus
}
#endif

#endif	/* RAW_MATRIX_U32_H_ */
