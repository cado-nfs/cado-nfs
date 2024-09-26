#ifndef MATMUL_MF_HPP_
#define MATMUL_MF_HPP_

#include "matmul.hpp"
#include "raw_matrix_u32.h"

void mf_prepare_matrix_u32(matmul_ptr mm, matrix_u32_ptr m, const char * file, int withcoeffs);

#endif	/* MATMUL_MF_HPP_ */
