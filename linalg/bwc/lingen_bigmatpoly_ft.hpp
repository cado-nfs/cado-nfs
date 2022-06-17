#ifndef BIGMATPOLY_FT_HPP_
#define BIGMATPOLY_FT_HPP_

#include "lingen_matpoly_ft.hpp"
#include "lingen_bigmatpoly.hpp"
#include "lingen_fft_select.hpp" // IWYU pragma: keep
#include "lingen_call_companion.hpp"
/*
#ifdef LINGEN_BINARY
struct gf2x_fake_fft_info;
struct gf2x_cantor_fft_info;
struct gf2x_ternary_fft_info;
#else
struct fft_transform_info;
#endif
*/

class tree_stats;

/* The class is only used as a namespace, really */
template<typename fft_type>
class bigmatpoly_ft : public matpoly_ft<fft_type> {
public:
    static bigmatpoly mp_caching(tree_stats & stats, bigmatpoly const & a, bigmatpoly const & b, lingen_call_companion::mul_or_mp_times * M);
    static bigmatpoly mul_caching(tree_stats & stats, bigmatpoly const & a, bigmatpoly const & b, lingen_call_companion::mul_or_mp_times * M);
};

#ifdef LINGEN_BINARY
extern template class bigmatpoly_ft<gf2x_fake_fft_info>;
extern template class bigmatpoly_ft<gf2x_cantor_fft_info>;
extern template class bigmatpoly_ft<gf2x_ternary_fft_info>;
#else
extern template class bigmatpoly_ft<fft_transform_info>;
#endif

#endif	/* BIGMATPOLY_FT_HPP_ */
