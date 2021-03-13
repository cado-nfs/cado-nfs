#ifndef MERGE_COMPUTE_WEIGHTS_H_
#define MERGE_COMPUTE_WEIGHTS_H_

#include "merge_replay_matrix.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void compute_weights_backend (filter_matrix_t *mat, index_t j0);

#ifdef __cplusplus
}
#endif

#endif	/* MERGE_COMPUTE_WEIGHTS_H_ */
