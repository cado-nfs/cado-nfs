#ifndef XDOTPROD_H_
#define XDOTPROD_H_

#include <stdint.h>              // for uint32_t
#include "matmul_top.h"

/* This interface is relevant to both krylov and mksol, since it's used
 * both for checking and computing A files (obviously only krylov is
 * concerned by the latter aspect). */

#ifdef __cplusplus
extern "C" {
#endif

void x_dotprod(void * dst, uint32_t * xv, unsigned int m, unsigned int nx, mmt_vec_ptr v, int sign);

#ifdef __cplusplus
}
#endif

#endif	/* XDOTPROD_H_ */
