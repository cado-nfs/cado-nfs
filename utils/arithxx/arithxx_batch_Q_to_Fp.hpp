#ifndef CADO_UTILS_ARITHXX_BATCH_Q_TO_FP_HPP
#define CADO_UTILS_ARITHXX_BATCH_Q_TO_FP_HPP

#include <cstdint>

#include <span>
#include <vector>

namespace arithxx_details {

template<typename layer>
struct batch_Q_to_Fp_context {
    typedef typename layer::Integer Integer;
    typedef typename layer::Modulus Modulus;
    Integer remainder, quotient;
    Modulus D;
    batch_Q_to_Fp_context(Integer const & num, Integer const & den);
    /* put num/den/2^k mod p[i] in r[i], which must have the size of p.
     * Return false if den is not invertible modulo one of the p[i]. */
    bool operator()(std::span<uint64_t> r, std::span<uint64_t const> p, int k=0) const;
    /* same, or an empty vector if den is not invertible */
    std::vector<uint64_t> operator()(std::vector<uint64_t> const & p, int k=0) const;
};

}


#endif	/* UTILS_ARITHXX_BATCH_Q_TO_FP_HPP_ */
