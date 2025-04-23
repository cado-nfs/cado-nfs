#ifndef CADO_MATMUL_TOP_COMM_HPP
#define CADO_MATMUL_TOP_COMM_HPP

#include "matmul_top_vec.hpp"

extern void mmt_vec_allreduce(mmt_vec & v);
extern void mmt_vec_broadcast(mmt_vec & v);
extern void mmt_vec_reduce(mmt_vec & w, mmt_vec & v);
extern void mmt_vec_reduce_mod_p(mmt_vec & v);
extern void mmt_vec_reduce_sameside(mmt_vec & v);


#endif	/* MATMUL_TOP_COMM_HPP_ */
