#ifndef CADO_MATMUL_SUB_VSC_COMBINE_H
#define CADO_MATMUL_SUB_VSC_COMBINE_H

/* This is only compiled when the underlying layer is b64 */

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* the arith_hard underlying type must be uint64_t, of course... */
extern void matmul_sub_vsc_combine_asm(uint64_t * dst, const uint64_t ** mptrs, const uint8_t * q, unsigned long count, unsigned int defer);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_MATMUL_SUB_VSC_COMBINE_H */
