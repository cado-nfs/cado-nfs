#ifndef CADO_LINGEN_BW_DIMENSIONS_HPP
#define CADO_LINGEN_BW_DIMENSIONS_HPP

#include "cxx_mpz.hpp"
#include "lingen_matpoly_select.hpp"

template<bool is_binary>
struct bw_dimensions {
    unsigned int m, n, nrhs = 0;
    typename matpoly<is_binary>::arith_hard ab;
    bw_dimensions(unsigned int m, unsigned int n, cxx_mpz const & p)
        : m(m)
        , n(n)
        , ab { p, 1U }
    {}
};


#endif	/* LINGEN_BW_DIMENSIONS_HPP_ */
