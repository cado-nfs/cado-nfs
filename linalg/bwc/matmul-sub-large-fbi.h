#ifndef MATMUL_SUB_LARGE_FBI_H_
#define MATMUL_SUB_LARGE_FBI_H_

#include "mpfq_layer.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void matmul_sub_large_fbi(abdst_field x, abelt ** sb, const abelt * z, const uint8_t * q, unsigned int n);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_SUB_LARGE_FBI_H_ */
