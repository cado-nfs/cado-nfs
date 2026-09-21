#ifndef CADO_READ_PURGEDFILE_IN_PARALLEL_HPP
#define CADO_READ_PURGEDFILE_IN_PARALLEL_HPP

#include <cstdint>

#include <string>

#include "merge_replay_matrix.h"

/* maximal number of threads when reading purged file
 *
 * note that this is mostly a limitation of the filesystem more than
 * anything else.
 */
inline constexpr unsigned int MAX_IO_THREADS = 16;

/* Fill mat->rows with the rows of the purged file, and return how many
 * were read. The file must be seekable: it is read by several threads at
 * once, each with its own reading head. */
uint64_t read_purgedfile_in_parallel(filter_matrix_t * mat,
                                     std::string const & filename);

#endif	/* CADO_READ_PURGEDFILE_IN_PARALLEL_HPP */
