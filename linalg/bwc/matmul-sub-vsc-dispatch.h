#ifndef MATMUL_SUB_VSC_DISPATCH_H_
#define MATMUL_SUB_VSC_DISPATCH_H_

#include "mpfq_layer.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void matmul_sub_vsc_dispatch(abdst_field x, abelt * dst, abelt const * src, const uint16_t * q, unsigned long count);


#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_SUB_VSC_DISPATCH_H_ */
