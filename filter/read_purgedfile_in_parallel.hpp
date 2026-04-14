#ifndef CADO_READ_PURGEDFILE_IN_PARALLEL_H
#define CADO_READ_PURGEDFILE_IN_PARALLEL_H

#include <stdint.h>
#include "merge_replay_matrix.h"

/* maximal number of threads when reading purged file
 *
 * note that this is mostly a limitation of the filesystem more than
 * anything else.
 */
#define MAX_IO_THREADS 16

#ifdef __cplusplus
extern "C" {
#endif

uint64_t read_purgedfile_in_parallel(filter_matrix_t * mat, const char * filename);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_READ_PURGEDFILE_IN_PARALLEL_H */
