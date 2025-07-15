#include "cado.h" // IWYU pragma: keep

#include <cstddef>

#include "lingen_round_operand_size.hpp"

size_t lingen_round_operand_size(size_t x, int bits) {
    /* round x up to the next size that has all but its six most significant
     * bits set to 0.
     */
    if (x == 0) return x;
    x -= 1;
    size_t y = x >> bits;
    for(int const i = 1 ; y ; y >>= i) x |= y;
    x += 1;
    return x;
}
