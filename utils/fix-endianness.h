#ifndef FIX_ENDIANNESS_H_
#define FIX_ENDIANNESS_H_

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "cado-endian.h"
#ifdef HAVE_SYS_ENDIAN_H
#include <sys/endian.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

size_t fread32_little(uint32_t * ptr, size_t nmemb, FILE * stream);
size_t fwrite32_little(const uint32_t * ptr, size_t nmemb, FILE * stream);
size_t fread64_little(uint64_t * ptr, size_t nmemb, FILE * stream);
size_t fwrite64_little(const uint64_t * ptr, size_t nmemb, FILE * stream);

#ifdef __cplusplus
}

#ifdef HAVE_BSWAP32
/* Freebsd already has bswap32.  */
static inline uint32_t bswap32(uint32_t x)
{
    // x is 3210
    uint32_t lohi = (x >> 16) | (x << 16);
    // lohi is 1032
    uint32_t m = UINT32_C(0xff00ff00);
    // retrieve x as .1.3  OR  0.2
    return ((lohi & m) >> 8) | ((lohi << 8) & m);
}
#endif

#ifdef CADO_LITTLE_ENDIAN
static inline uint32_t bfix32(uint32_t x) { return x; }
#elif defined(CADO_BIG_ENDIAN)
static inline uint32_t bfix32(uint32_t x) { return bswap32(x); }
#else
#error "implement me"
#endif
#endif

#endif	/* FIX_ENDIANNESS_H_ */
