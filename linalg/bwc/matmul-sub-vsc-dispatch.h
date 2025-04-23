#ifndef CADO_MATMUL_SUB_VSC_DISPATCH_H
#define CADO_MATMUL_SUB_VSC_DISPATCH_H

/* This is only compiled when the underlying layer is b64 */

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* the arith_hard underlying type must be uint64_t, of course... */
extern void matmul_sub_vsc_dispatch_asm(uint64_t * dst, uint64_t const * src, const uint16_t * q, unsigned long count);


#ifdef __cplusplus
}
#endif

#endif	/* CADO_MATMUL_SUB_VSC_DISPATCH_H */
