#include "cado.h" // IWYU pragma: keep

#include <stddef.h>
#include <stdint.h>

#include "macros.h"        // cado_crc_lfsr
#include "crc.h"        // cado_crc_lfsr

/* This computes checksum or arbitrary data ranges. The data is piped
 * through an LFSR over GF(2^32) with a suitable defining polynomial.
 */

uint32_t cado_crc_lfsr_turn1(cado_crc_lfsr_ptr l, uint32_t c)
{
    static uint32_t const twist = 3486328325U;

    uint32_t w = 0;

    // FIXME. We used to have 31U^l->r, which very much seems to be a typo
    // for 31U&l->r. Unfortunately, ``fixing'' this would mean changing
    // behaviour on our main platform, so for the time being we change
    // this in a compatible way so as to reach the same output on non-x86
    // hardware (which wraps around shift counts, not a guaranteed
    // behaviour everywhere).
    w = (c >> (31U&(31U^l->r))) ^ (c << (31U&-l->r));
    l->i--;
    l->r+=11;
    l->i &= 31U;

    w ^= l->t[ l->i             ];
    w ^= l->t[(l->i + 22U) & 31U];
    w ^= l->t[(l->i +  2U) & 31U];
    w ^= l->t[(l->i +  1U) & 31U];

    w = w >> 1U ^ (twist & -(w&1U));
    l->t[l->i] = w;

    return w;
}

void cado_crc_lfsr_init(cado_crc_lfsr_ptr l)
{
    l->i = 0;
    l->r = 0;
    for(int k = 0 ; k < 32 ; k++) l->t[k] = k;
    for(int k = 0 ; k < 96 ; k++) {
        cado_crc_lfsr_turn1(l, 0xdeadbeef / (k+1));
    }
}

void cado_crc_lfsr_clear(cado_crc_lfsr_ptr l MAYBE_UNUSED) { }

uint32_t cado_crc_lfsr_turn(cado_crc_lfsr_ptr l, const void * data, size_t count)
{
    const uint8_t * ptr = (const uint8_t *) data;
    uint32_t w = 0;

    for( ; count-- ; ) {
        w = cado_crc_lfsr_turn1(l, *ptr++);
    }
    return w;
}

/* This version yields the same checksum regardless of the endianness */
uint32_t cado_crc_lfsr_turn32_little(cado_crc_lfsr_ptr l, const uint32_t * data, size_t count)
{
    ASSERT_ALWAYS(sizeof(uint32_t) == 4);
    ASSERT_ALWAYS(count % 4 == 0);
    int twist[sizeof(uint32_t)];
    for(int j = 0 ; j < (int) sizeof(uint32_t) ; j++) {
        const uint32_t c = 0x03020100;
        twist[((uint8_t *)&c)[j]]=j;
    }
    uint32_t w = 0;

    for(unsigned int i = 0 ; i * sizeof(uint32_t) < count ; i++) {
        const uint8_t * ptr = (const uint8_t *) (data + i);
        for(size_t j = 0 ; j < sizeof(uint32_t) ; j++) {
            w = cado_crc_lfsr_turn1(l, ptr[twist[j]]);
        }
    }
    return w;
}

uint32_t crc32(const void * data, size_t count)
{
    cado_crc_lfsr l;
    cado_crc_lfsr_init(l);
    uint32_t w = cado_crc_lfsr_turn(l, data, count);
    cado_crc_lfsr_clear(l);
    return w;
}
