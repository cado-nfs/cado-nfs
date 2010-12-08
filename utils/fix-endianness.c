#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "macros.h"

#include "fix-endianness.h"
#include "cado-endian.h"

#ifdef CADO_LITTLE_ENDIAN
size_t fread32_little(uint32_t * ptr, size_t nmemb, FILE * stream)
{
    return fread(ptr, sizeof(uint32_t), nmemb, stream);
}
size_t fwrite32_little(const uint32_t * ptr, size_t nmemb, FILE * stream)
{
    return fwrite(ptr, sizeof(uint32_t), nmemb, stream);
}
size_t fread64_little(uint64_t * ptr, size_t nmemb, FILE * stream)
{
    return fread(ptr, sizeof(uint64_t), nmemb, stream);
}
size_t fwrite64_little(const uint64_t * ptr, size_t nmemb, FILE * stream)
{
    return fwrite(ptr, sizeof(uint64_t), nmemb, stream);
}
#elif defined(CADO_BIG_ENDIAN)

/* This is slow and stupid. */
size_t fread32_little(uint32_t * ptr, size_t nmemb, FILE * stream)
{
    uint32_t * p = malloc(nmemb * sizeof(uint32_t));
    ASSERT_ALWAYS(sizeof(uint32_t) == 4);
    size_t n = fread(p, sizeof(uint32_t), nmemb, stream);
    int twist[sizeof(uint32_t)];
    for(unsigned int j = 0 ; j < sizeof(uint32_t) ; j++) {
        const uint32_t c = 0x03020100;
        twist[j] = ((uint8_t *)&c)[j];
        // byte j, in machine order, contains value twist[j], which is (by
        // construction of c) to be found at offset twist[j] in a
        // little-endian 32-bit integer.
    }
    for(size_t i = 0 ; i < n ; i++) {
        uint8_t * dst = (uint8_t *) (ptr + i);
        const uint8_t * src = (const uint8_t *) (p + i);
        for(unsigned int j = 0 ; j < sizeof(uint32_t) ; j++) {
            dst[j] = src[twist[j]];
        }
    }
    free(p);
    return n;
}

size_t fwrite32_little(const uint32_t * ptr, size_t nmemb, FILE * stream)
{
    uint32_t * p = malloc(nmemb * sizeof(uint32_t));
    ASSERT_ALWAYS(sizeof(uint32_t) == 4);
    int twist[sizeof(uint32_t)];
    for(unsigned int j = 0 ; j < sizeof(uint32_t) ; j++) {
        const uint32_t c = 0x03020100;
        twist[j] = ((uint8_t *)&c)[j];
        // byte j, in machine order, contains value twist[j], which is (by
        // construction of c) to be found at offset twist[j] in a
        // little-endian 32-bit integer.
    }
    for(size_t i = 0 ; i < nmemb ; i++) {
        uint8_t * dst = (uint8_t *) (p + i);
        const uint8_t * src = (const uint8_t *) (ptr + i);
        for(unsigned int j = 0 ; j < sizeof(uint32_t) ; j++) {
            dst[twist[j]] = src[j];
        }
    }
    size_t n = fwrite(p, sizeof(uint32_t), nmemb, stream);
    free(p);
    return n;
}


size_t fread64_little(uint64_t * ptr, size_t nmemb, FILE * stream)
{
    uint64_t * p = malloc(nmemb * sizeof(uint64_t));
    ASSERT_ALWAYS(sizeof(uint64_t) == 8);
    size_t n = fread(p, sizeof(uint64_t), nmemb, stream);
    int twist[sizeof(uint64_t)];
    for(unsigned int j = 0 ; j < sizeof(uint64_t) ; j++) {
        const uint64_t c = UINT64_C(0x0706050403020100);
        twist[j] = ((uint8_t *)&c)[j];
        // byte j, in machine order, contains value twist[j], which is (by
        // construction of c) to be found at offset twist[j] in a
        // little-endian 64-bit integer.
    }
    for(size_t i = 0 ; i < n ; i++) {
        uint8_t * dst = (uint8_t *) (ptr + i);
        const uint8_t * src = (const uint8_t *) (p + i);
        for(unsigned int j = 0 ; j < sizeof(uint64_t) ; j++) {
            dst[j] = src[twist[j]];
        }
    }
    free(p);
    return n;
}

size_t fwrite64_little(const uint64_t * ptr, size_t nmemb, FILE * stream)
{
    uint64_t * p = malloc(nmemb * sizeof(uint64_t));
    ASSERT_ALWAYS(sizeof(uint64_t) == 8);
    int twist[sizeof(uint64_t)];
    for(unsigned int j = 0 ; j < sizeof(uint64_t) ; j++) {
        const uint64_t c = UINT64_C(0x0706050403020100);
        twist[j] = ((uint8_t *)&c)[j];
        // byte j, in machine order, contains value twist[j], which is (by
        // construction of c) to be found at offset twist[j] in a
        // little-endian 64-bit integer.
    }
    for(size_t i = 0 ; i < nmemb ; i++) {
        uint8_t * dst = (uint8_t *) (p + i);
        const uint8_t * src = (const uint8_t *) (ptr + i);
        for(unsigned int j = 0 ; j < sizeof(uint64_t) ; j++) {
            dst[twist[j]] = src[j];
        }
    }
    size_t n = fwrite(p, sizeof(uint64_t), nmemb, stream);
    free(p);
    return n;
}
#else
#error "implement me"
#endif
