#ifndef CADO_MERGE_COMPUTE_WEIGHTS_H
#define CADO_MERGE_COMPUTE_WEIGHTS_H

#include "merge_replay_matrix.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void compute_weights_backend (filter_matrix_t *mat, index_t j0);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_MERGE_COMPUTE_WEIGHTS_H */
