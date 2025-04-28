#ifndef CADO_LINGEN_BW_DIMENSIONS_HPP
#define CADO_LINGEN_BW_DIMENSIONS_HPP

#include "arith-hard.hpp"       // IWYU pragma: keep
#include "cxx_mpz.hpp"
#include "lingen_matpoly_select.hpp"   // for matpoly

struct bw_dimensions {
    unsigned int m, n, nrhs = 0;
    matpoly::arith_hard ab;
    bw_dimensions(unsigned int m, unsigned int n, cxx_mpz const & p)
        : m(m)
        , n(n)
        , ab { p,1U }
    {}
};


#endif	/* LINGEN_BW_DIMENSIONS_HPP_ */
