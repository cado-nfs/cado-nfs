#include "utils.h"

/* This computes checksum or arbitrary data ranges. The data is piped
 * through an LFSR over GF(2^32) with a suitable defining polynomial.
 */

uint32_t cado_crc_lfsr_turn1(cado_crc_lfsr l, uint32_t c)
{
    static uint32_t const twist = 3486328325;

    uint32_t w = 0;

    // FIXME. We used to have 31^l->r, which very much seems to be a typo
    // for 31&l->r. Unfortunately, ``fixing'' this would mean changing
    // behaviour on our main platform, so for the time being we change
    // this in a compatible way so as to reach the same output on non-x86
    // hardware (which wraps around shift counts, not a guaranteed
    // behaviour everywhere).
    w = (c >> (31&(31^l->r))) ^ (c << (31&-l->r));
    l->i--;
    l->r+=11;
    l->i &= 31;

    w ^= l->t[ l->i           ];
    w ^= l->t[(l->i + 22) & 31];
    w ^= l->t[(l->i +  2) & 31];
    w ^= l->t[(l->i +  1) & 31];

    w = w >> 1 ^ (twist & -(w&1));
    l->t[l->i] = w;

    return w;
}

void cado_crc_lfsr_init(cado_crc_lfsr l)
{
    l->i = 0;
    l->r = 0;
    for(int k = 0 ; k < 32 ; k++) l->t[k] = k;
    for(int k = 0 ; k < 96 ; k++) {
        cado_crc_lfsr_turn1(l, 0xdeadbeef / (k+1));
    }
}

void cado_crc_lfsr_clear(cado_crc_lfsr l MAYBE_UNUSED) { }

uint32_t cado_crc_lfsr_turn(cado_crc_lfsr l, const void * data, unsigned int count)
{
    const uint8_t * ptr = (const uint8_t *) data;
    uint32_t w = 0;

    for( ; count-- ; ) {
        w = cado_crc_lfsr_turn1(l, *ptr++);
    }
    return w;
}

/* This version yields the same checksum regardless of the endianness */
uint32_t cado_crc_lfsr_turn32_little(cado_crc_lfsr l, const uint32_t * data, unsigned int count)
{
    ASSERT_ALWAYS(sizeof(uint32_t) == 4);
    int twist[sizeof(uint32_t)];
    for(unsigned int j = 0 ; j < sizeof(uint32_t) ; j++) {
        const uint32_t c = 0x03020100;
        twist[((uint8_t *)&c)[j]]=j;
    }
    uint32_t w = 0;

    for( ; count-- ; data++) {
        const uint8_t * ptr = (const uint8_t *) data;
        for(unsigned int j = 0 ; j < sizeof(uint32_t) ; j++) {
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
