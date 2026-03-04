#ifndef UTILS_REPLACEMENT_LATCH_HPP_
#define UTILS_REPLACEMENT_LATCH_HPP_

/* This is meant to provide an equivalent for std::latch for
 * libstdc++10, which doesn't have it (only version 11.1 has
 * std::latch)
 *
 */

#ifndef __cpp_lib_latch
#include <barrier>

namespace std {
    using latch = cado_std_replacement::barrier;
}
#endif

#endif	/* UTILS_REPLACEMENT_LATCH_HPP_ */
