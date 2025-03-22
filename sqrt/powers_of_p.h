#ifndef CADO_POWERS_OF_P_H
#define CADO_POWERS_OF_P_H

#include <gmp.h>

// this is a trampoline to C++. Entry points are in C.
#ifdef __cplusplus
extern "C" {
#endif

void * power_lookup_table_init(unsigned long p);
void power_lookup_table_clear(void * t);
mpz_srcptr power_lookup(void * t, int i);
mpz_srcptr power_lookup_const(const void * t, int i);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_POWERS_OF_P_H */
