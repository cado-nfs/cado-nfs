#ifndef CADO_UTILS_ARITHXX_BATCH_Q_TO_FP_IMPL_HPP
#define CADO_UTILS_ARITHXX_BATCH_Q_TO_FP_IMPL_HPP

#include <cstddef>
#include <cstdint>

#include <vector>

#include "arithxx_batch_Q_to_Fp.hpp"

/* note that the modredc64 impl has an explicit specialization of this
 */
template<typename layer>
arithxx_details::batch_Q_to_Fp_context<layer>::batch_Q_to_Fp_context(Integer const & num, Integer const & den)
        : remainder(num % den)
        , quotient((num - remainder).divexact(den))
        , D(den)
    { }

template<typename layer>
std::vector<uint64_t> arithxx_details::batch_Q_to_Fp_context<layer>::operator()(std::vector<uint64_t> const & p, int const k) const
{
    /* We use -rem (mod den) here. batchinv_ul() does not
       mandate its c parameter to be fully reduced, which occurs here in the
       case of rem == 0. */
    auto r = D.batchinv_redc(p, Integer(D.m - remainder));
    if (r.empty())
        return {};

    std::vector<uint64_t> ri(p.size());

    for (size_t i = 0; i < p.size(); i++)
        ri[i] = u64arith_post_process_inverse(r[i][0], p[i],
                remainder[0], -D.invm, quotient[0], k);

    return ri;
}


#endif	/* UTILS_ARITHXX_BATCH_Q_TO_FP_IMPL_HPP_ */
