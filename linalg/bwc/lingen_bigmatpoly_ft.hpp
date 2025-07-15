#ifndef CADO_LINGEN_BIGMATPOLY_FT_HPP
#define CADO_LINGEN_BIGMATPOLY_FT_HPP

#include "lingen_matpoly_ft.hpp"
#include "lingen_bigmatpoly.hpp"
#include "lingen_fft_select.hpp"
#include "lingen_call_companion.hpp"

class tree_stats;

/* The class is only used as a namespace, really */
template<typename fft_type>
class bigmatpoly_ft : public matpoly_ft<fft_type> {
public:
    static constexpr bool is_binary = is_binary_fft<fft_type>::value;
    static 
        auto
        mp_caching(tree_stats & stats, bigmatpoly<is_binary> const & a, bigmatpoly<is_binary> const & b, lingen_call_companion::mul_or_mp_times * M)
        -> bigmatpoly<is_binary>;
    static
    auto
    mul_caching(tree_stats & stats, bigmatpoly<is_binary> const & a, bigmatpoly<is_binary> const & b, lingen_call_companion::mul_or_mp_times * M)
    ->
    bigmatpoly<is_binary>;
};

#ifdef LINGEN_BINARY
extern template class bigmatpoly_ft<gf2x_fake_fft_info>;
extern template class bigmatpoly_ft<gf2x_cantor_fft_info>;
extern template class bigmatpoly_ft<gf2x_ternary_fft_info>;
#else
extern template class bigmatpoly_ft<fft_transform_info>;
#endif

#endif	/* LINGEN_BIGMATPOLY_FT_HPP_ */
