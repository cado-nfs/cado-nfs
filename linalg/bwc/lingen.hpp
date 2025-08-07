#ifndef CADO_LINALG_BWC_LINGEN_HPP
#define CADO_LINALG_BWC_LINGEN_HPP

#include "lingen_bigmatpoly.hpp"
#include "lingen_bmstatus.hpp"
#include "lingen_matpoly_select.hpp"

template<bool is_binary>
matpoly<is_binary> bw_lingen_single(bmstatus<is_binary> & bm, matpoly<is_binary> & E);

template<bool is_binary>
bigmatpoly<is_binary> bw_biglingen_collective(bmstatus<is_binary> & bm, bigmatpoly<is_binary> & E);

#endif	/* LINALG_BWC_LINGEN_HPP_ */
