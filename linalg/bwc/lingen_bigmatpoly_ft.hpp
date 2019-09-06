#ifndef BIGMATPOLY_FT_HPP_
#define BIGMATPOLY_FT_HPP_

#include "lingen_matpoly_ft.hpp"
#include "lingen_bigmatpoly.hpp"

#include "select_mpi.h"
#include "lingen_substep_schedule.hpp"
#include "tree_stats.hpp"

/* The class is only used as a namespace, really */
template<typename fft_type>
class bigmatpoly_ft : public matpoly_ft<fft_type> {
public:
    static void mp_caching(tree_stats & stats, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, lingen_call_companion::mul_or_mp_times * M);
    static void mul_caching(tree_stats & stats, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, lingen_call_companion::mul_or_mp_times * M);
};

#ifdef SELECT_MPFQ_LAYER_u64k1
extern template class bigmatpoly_ft<gf2x_fake_fft>;
extern template class bigmatpoly_ft<gf2x_cantor_fft>;
extern template class bigmatpoly_ft<gf2x_ternary_fft>;
#else
extern template class bigmatpoly_ft<fft_transform_info>;
#endif

#endif	/* BIGMATPOLY_FT_HPP_ */
