#ifndef LINALG_BWC_LINGEN_HPP_
#define LINALG_BWC_LINGEN_HPP_

#include "lingen_bigmatpoly.hpp"          // for bigmatpoly, bigmatpoly_model
#include "lingen_bmstatus.hpp"            // for bmstatus
#include "lingen_matpoly_select.hpp"      // for matpoly, matpoly::memory_guard

matpoly bw_lingen_single(bmstatus & bm, matpoly & E);
bigmatpoly bw_biglingen_collective(bmstatus & bm, bigmatpoly & E);


#endif	/* LINALG_BWC_LINGEN_HPP_ */
