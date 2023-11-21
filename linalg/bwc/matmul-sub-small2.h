#ifndef MATMUL_SUB_SMALL2_H_
#define MATMUL_SUB_SMALL2_H_

/* This is only compiled when the underlying layer is b64 */

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* the arith_hard underlying type must be uint64_t, of
 * course... */
const uint16_t * matmul_sub_small2_asm(uint64_t * where, uint64_t const * from, const uint16_t * q, unsigned long n);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_SUB_SMALL2_H_ */

