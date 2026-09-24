#ifndef CADO_UTILS_ARITHXX_REDC_IMPL_HPP
#define CADO_UTILS_ARITHXX_REDC_IMPL_HPP

#include <cstdint>

#include <span>
#include <vector>

#include "arithxx_redc.hpp"
#include "arithxx_batch_Q_to_Fp.hpp"

template<typename layer>
    std::vector<uint64_t>
arithxx_details::redc<layer>::batch_Q_to_Fp(Integer const & num,
        Integer const & den, int k,
            std::vector<uint64_t> const & p)
{
    return batch_Q_to_Fp_context<layer>(num, den)(p, k);
}

template<typename layer>
bool
arithxx_details::redc<layer>::batch_Q_to_Fp(std::span<uint64_t> r,
        Integer const & num, Integer const & den, int k,
        std::span<uint64_t const> p)
{
    return batch_Q_to_Fp_context<layer>(num, den)(r, p, k);
}


#endif	/* UTILS_ARITHXX_REDC_IMPL_HPP_ */
