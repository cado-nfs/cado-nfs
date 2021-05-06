#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>

int main() {
    volatile __mmask8 k = _cvtu32_mask8(14);
    volatile __mmask8 l = _cvtu32_mask8(125);
    return _cvtmask8_u32(_kxor_mask8(k, l)) == (14 ^ 125) ? EXIT_SUCCESS : EXIT_FAILURE;
}
