#ifndef CADO_GCD_H
#define CADO_GCD_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

uint64_t gcd_int64 (int64_t a, int64_t b);
uint64_t gcd_uint64 (uint64_t a, uint64_t b);
unsigned long gcd_ul (unsigned long a, unsigned long b);
unsigned long xgcd_ul (unsigned long * xa, unsigned long a, unsigned long b);
unsigned long invert_ul (unsigned long a, unsigned long b);
uint64_t bin_gcd_int64 (int64_t a, int64_t b);
uint64_t bin_gcd_int64_safe (int64_t a, int64_t b);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_GCD_H */
