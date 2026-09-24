#ifndef CADO_UTILS_ARITHXX_BATCH_Q_TO_FP_IMPL_HPP
#define CADO_UTILS_ARITHXX_BATCH_Q_TO_FP_IMPL_HPP

#include <cstddef>
#include <cstdint>

#include <span>
#include <type_traits>
#include <vector>

#include "arithxx_batch_Q_to_Fp.hpp"
#include "macros.h"
#include "modint.hpp"
#include "u64arith.h"

/* note that the modredc64 impl has an explicit specialization of this
 */
template<typename layer>
arithxx_details::batch_Q_to_Fp_context<layer>::batch_Q_to_Fp_context(Integer const & num, Integer const & den)
        : remainder(num % den)
        , quotient((num - remainder).divexact(den))
        , D(den)
    { }

template<typename layer>
bool arithxx_details::batch_Q_to_Fp_context<layer>::operator()(
        std::span<uint64_t> r, std::span<uint64_t const> p, int const k) const
{
    ASSERT_ALWAYS(r.size() == p.size());

    /* We use -rem (mod den) here. batchinv_redc() does not
       mandate its c parameter to be fully reduced, which occurs here in the
       case of rem == 0. */
    Integer const c(D.m - remainder);

    if constexpr (std::is_same_v<Integer, Integer64>) {
        /* den fits in one word: the inverses go straight to r */
        if (!D.batchinv_redc(r, p, c))
            return false;
    } else {
        auto const ri = D.batchinv_redc(p, c);
        if (ri.empty() && !p.empty())
            return false;
        for (size_t i = 0; i < p.size(); i++)
            r[i] = ri[i][0];
    }

    for (size_t i = 0; i < p.size(); i++)
        r[i] = u64arith_post_process_inverse(r[i], p[i],
                remainder[0], -D.invm, quotient[0], k);

    return true;
}

template<typename layer>
std::vector<uint64_t> arithxx_details::batch_Q_to_Fp_context<layer>::operator()(std::vector<uint64_t> const & p, int const k) const
{
    std::vector<uint64_t> r(p.size());
    if (!(*this)(std::span<uint64_t>(r), std::span<uint64_t const>(p), k))
        return {};
    return r;
}


#endif	/* UTILS_ARITHXX_BATCH_Q_TO_FP_IMPL_HPP_ */
