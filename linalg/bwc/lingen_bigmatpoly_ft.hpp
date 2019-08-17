#ifndef BIGMATPOLY_FT_HPP_
#define BIGMATPOLY_FT_HPP_

#include "mpfq_layer.h"
#include "lingen_matpoly_ft.hpp"
#include "lingen_bigmatpoly.hpp"
#include "flint-fft/fft.h"

#include "select_mpi.h"
#include "lingen_substep_schedule.hpp"
#include "tree_stats.hpp"

void bigmatpoly_mul_caching_adj(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, unsigned int adj, const struct lingen_substep_schedule * S);

void bigmatpoly_mp_caching_adj(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, unsigned int adj, const struct lingen_substep_schedule * S);

static inline void bigmatpoly_mul_caching(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, const struct lingen_substep_schedule * S) { return bigmatpoly_mul_caching_adj(t, c, a, b, UINT_MAX, S); }

static inline void bigmatpoly_mp_caching(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, const struct lingen_substep_schedule * S) { return bigmatpoly_mp_caching_adj(t, c, a, b, UINT_MAX, S); }

#endif	/* BIGMATPOLY_FT_HPP_ */
