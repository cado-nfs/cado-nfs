#ifndef CADO_MATMUL_SUB_LARGE_FBI_H
#define CADO_MATMUL_SUB_LARGE_FBI_H

/* This is only compiled when the underlying layer is b64 */

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* the arith_hard underlying type must be uint64_t, of
 * course... */
extern void matmul_sub_large_fbi_asm(uint64_t ** sb, const uint64_t * z, const uint8_t * q, unsigned int n);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_MATMUL_SUB_LARGE_FBI_H */
