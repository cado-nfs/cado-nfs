#ifndef UTILSINT64_T_H
#define UTILSINT64_T_H

#include <stdint.h>
#include "utils.h"
#include "mod_ul_default.h"
#include "mod_ul.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Compute the modular inverse of xx mod mm.
 *
 * xx: the number for which we want the inverse.
 * mm: the modulo.
 */
uint64_t invmod_uint64(uint64_t xx, uint64_t mm);

void int64_fdiv_qr(int64_t * q, int64_t * r, int64_t n, int64_t d);

static inline void swap_int64(int64_t * a, int64_t * b)
{
  int64_t tmp = * a;
  * a = * b;
  * b = tmp;
}

static inline int64_t int64_sgn(int64_t a)
{
  if (a == 0) {
    return 0;
  } else if (a < 0) {
    return -1;
  } else {
    return 1;
  }
}

void int64_xgcd(int64_t * g, int64_t * u, int64_t * v, int64_t a, int64_t b);

void int64_gcdext(int64_t * e, int64_t * s, int64_t * t, int64_t a, int64_t b);

#ifdef __cplusplus
}
#endif
#endif /* UTILSINT64_T_H */
