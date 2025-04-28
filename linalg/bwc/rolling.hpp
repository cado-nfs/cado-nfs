#ifndef CADO_ROLLING_HPP
#define CADO_ROLLING_HPP

#include <string>

/* Removes all vector files in bwc's wdir (current directory, in fact,
 * since bwc programs always run in the current directory) according to
 * the rules specified by the parameters keep_rolling_checkpoints and
 * checkpoint_precious.
 *
 * If v == 0, consider all checkpoints. Otherwise, restrict to those
 * whose index is <= v
 *
 * The bal argument is used only to compose the filename according to the
 * checksum of the current balancing permutation.
 */

void keep_rolling_checkpoints(std::string const & stem, unsigned int v);

#endif	/* ROLLING_HPP_ */
