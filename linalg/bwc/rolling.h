#ifndef ROLLING_H_
#define ROLLING_H_

#ifdef __cplusplus
#include <string>
#endif

#ifdef __cplusplus
extern "C" {
#endif

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
extern void keep_rolling_checkpoints(const char * stem, unsigned int v);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
static inline void keep_rolling_checkpoints(std::string const & stem, unsigned int v)
{
    keep_rolling_checkpoints(stem.c_str(), v);
}
#endif

#endif	/* ROLLING_H_ */
