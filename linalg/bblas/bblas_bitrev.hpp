#ifndef BBLAS_BITREV_HPP_
#define BBLAS_BITREV_HPP_

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"

#include <cstdint>
#include "macros.h"

static inline uint64_t MAYBE_UNUSED bitrev(uint64_t a)
{
    a = (a >> 32) ^ (a << 32);
    uint64_t m;
    m = UINT64_C(0x0000ffff0000ffff);
    a = ((a >> 16) & m) ^ ((a << 16) & ~m);
    m = UINT64_C(0x00ff00ff00ff00ff);
    a = ((a >> 8) & m) ^ ((a << 8) & ~m);
    m = UINT64_C(0x0f0f0f0f0f0f0f0f);
    a = ((a >> 4) & m) ^ ((a << 4) & ~m);
    m = UINT64_C(0x3333333333333333);
    a = ((a >> 2) & m) ^ ((a << 2) & ~m);
    m = UINT64_C(0x5555555555555555);
    a = ((a >> 1) & m) ^ ((a << 1) & ~m);
    return a;
}
/* like bitrev, but keep nibbles intact */
static inline uint64_t MAYBE_UNUSED nibrev(uint64_t a)
{
    a = (a >> 32) ^ (a << 32);
    uint64_t m;
    m = UINT64_C(0x0000ffff0000ffff);
    a = ((a >> 16) & m) ^ ((a << 16) & ~m);
    m = UINT64_C(0x00ff00ff00ff00ff);
    a = ((a >> 8) & m) ^ ((a << 8) & ~m);
    m = UINT64_C(0x0f0f0f0f0f0f0f0f);
    a = ((a >> 4) & m) ^ ((a << 4) & ~m);
    return a;
}

#endif	/* BBLAS_BITREV_HPP_ */
