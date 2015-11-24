#ifndef MATMUL_SUB_SMALL1_H_
#define MATMUL_SUB_SMALL1_H_

#include "mpfq_layer.h"

#ifdef __cplusplus
extern "C" {
#endif

const uint16_t * matmul_sub_small1(abdst_field x, abelt * where, const abelt * from, const uint16_t * q, unsigned long n);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_SUB_SMALL1_H_ */

