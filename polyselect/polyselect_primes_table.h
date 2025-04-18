#ifndef CADO_POLYSELECT_PRIMES_TABLE_H
#define CADO_POLYSELECT_PRIMES_TABLE_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct polyselect_primes_table_s {
    uint32_t *Primes;
    size_t lenPrimes;
};

typedef struct polyselect_primes_table_s polyselect_primes_table[1];
typedef struct polyselect_primes_table_s * polyselect_primes_table_ptr;
typedef const struct polyselect_primes_table_s * polyselect_primes_table_srcptr;

extern void polyselect_primes_table_init(polyselect_primes_table_ptr pt, unsigned long P);
extern void polyselect_primes_table_clear(polyselect_primes_table_ptr pt);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_POLYSELECT_PRIMES_TABLE_H */
