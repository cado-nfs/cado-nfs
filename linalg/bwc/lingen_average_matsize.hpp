#ifndef CADO_LINGEN_AVERAGE_MATSIZE_HPP
#define CADO_LINGEN_AVERAGE_MATSIZE_HPP

#include "lingen_matpoly_select.hpp"

template<bool is_binary>
double average_matsize(typename matpoly<is_binary>::arith_hard * ab, unsigned int m, unsigned int n, int ascii);

#ifdef LINGEN_BINARY
template<> double average_matsize<true>(typename matpoly<true>::arith_hard * ab, unsigned int m, unsigned int n, int ascii);
#else
template<> double average_matsize<false>(typename matpoly<false>::arith_hard * ab, unsigned int m, unsigned int n, int ascii);
#endif

#endif	/* LINGEN_AVERAGE_MATSIZE_HPP_ */
