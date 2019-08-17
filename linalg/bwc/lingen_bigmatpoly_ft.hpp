#ifndef BIGMATPOLY_FT_HPP_
#define BIGMATPOLY_FT_HPP_

#include "mpfq_layer.h"
#include "lingen_matpoly_ft.hpp"
#include "lingen_bigmatpoly.hpp"
#include "flint-fft/fft.h"

#include "select_mpi.h"
#include "lingen_substep_schedule.hpp"
#include "tree_stats.hpp"

/* This defines an MPI-shared polynomial matrix type */

void bigmatpoly_mul_caching_adj(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, unsigned int adj, lingen_call_companion::mul_or_mp_times * M);

void bigmatpoly_mp_caching_adj(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, unsigned int adj, lingen_call_companion::mul_or_mp_times * M);

static inline void bigmatpoly_mul_caching(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, lingen_call_companion::mul_or_mp_times * M) { return bigmatpoly_mul_caching_adj(t, c, a, b, UINT_MAX, M); }

static inline void bigmatpoly_mp_caching(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, lingen_call_companion::mul_or_mp_times * M) { return bigmatpoly_mp_caching_adj(t, c, a, b, UINT_MAX, M); }

#endif	/* BIGMATPOLY_FT_HPP_ */
