/* This source file is our test case for ssse3 support. */
#include <stdint.h>
#include <string.h>
/* we could use tmmintrin.h with pre-4.5.0 gcc. See
 * https://stackoverflow.com/questions/11228855/header-files-for-x86-simd-intrinsics
 * Note however that we tend to use x86intrin anyway in the code.
 */
#include <x86intrin.h>

int main()
{
    volatile uint32_t a0 = 0x03020100;
    volatile uint32_t a1 = 0x1F1E1D1C;
    __m128i x = _mm_setr_epi32(a0, 0x07060504, 0x0B0A0908, 0x0F0E0D0C);
    __m128i y = _mm_setr_epi32(0x13121110, 0x17161514, 0x1B1A1918, a1);
    uint64_t a[2], b[2] = { 0x0C0B0A0908070605, 0x14131211100F0E0D };
    y = _mm_alignr_epi8(y, x, 0x5);
    memcpy (a, &y, 16);
    return(a[0] != b[0] || a[1] != b[1]);
}
