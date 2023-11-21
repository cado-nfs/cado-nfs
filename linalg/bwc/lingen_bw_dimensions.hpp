#ifndef LINGEN_BW_DIMENSIONS_HPP_
#define LINGEN_BW_DIMENSIONS_HPP_

#include "arith-hard.hpp"
#include "cxx_mpz.hpp"
#include "lingen_matpoly_select.hpp"   // for matpoly

struct bw_dimensions {
    unsigned int m,n,nrhs;
    matpoly::arith_hard ab;
    bw_dimensions(unsigned int m, unsigned int n, cxx_mpz const & p)
        : m(m)
        , n(n)
        , nrhs(0)
        , ab(p,1)
    {}
};


#endif	/* LINGEN_BW_DIMENSIONS_HPP_ */
