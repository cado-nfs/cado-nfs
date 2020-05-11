#ifndef LAS_WHERE_AM_I_PROD_HPP_
#define LAS_WHERE_AM_I_PROD_HPP_

// IWYU pragma: private, include "las-where-am-i.hpp"

#include <cstdint>
#include "las-where-am-i-proxy.hpp"        // for where_am_I
#include "macros.h"     // IWYU pragma: keep

struct where_am_I::impl {
};

#define WHERE_AM_I_UPDATE(w, field, value) /**/

static inline int test_divisible(where_am_I const & w MAYBE_UNUSED) { return 1; }
static inline int trace_on_spot_N(unsigned int N MAYBE_UNUSED) { return 0; }
static inline int trace_on_spot_Nx(unsigned int N MAYBE_UNUSED, unsigned int x MAYBE_UNUSED) { return 0; }
static inline int trace_on_range_Nx(unsigned int N MAYBE_UNUSED, unsigned int x0 MAYBE_UNUSED, unsigned int x1 MAYBE_UNUSED) { return 0; }
static inline int trace_on_spot_x(unsigned int x MAYBE_UNUSED) { return 0; }
static inline int trace_on_spot_ab(int64_t a MAYBE_UNUSED, uint64_t b MAYBE_UNUSED) { return 0; }
static inline int trace_on_spot_ij(int i MAYBE_UNUSED, unsigned int j MAYBE_UNUSED) { return 0; }

static inline void sieve_increase(unsigned char *S, const unsigned char logp, where_am_I const & w MAYBE_UNUSED)
{
    *S += logp;
}


#endif	/* LAS_WHERE_AM_I_PROD_HPP_ */
