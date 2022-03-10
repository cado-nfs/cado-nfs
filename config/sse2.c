/* This source file is our test case for sse-2 support. If this can be
 * compiled, then we're happy because sse-2 will be usable in cantor/ ;
 */
#include <stdint.h>
#include <emmintrin.h>

int main(int argc, char *argv[])
{
    volatile int a0 = 17;
    volatile int a1 = 42;
    __m128i foo = _mm_setr_epi32(argc, argc + 1, argc + 2, argc + 3);
    __m128i bar = _mm_setr_epi32(argc + 3, argc + 2, argc + 1, argc);
    __m128i x = _mm_setr_epi32(a1, 0, a0, 0);
    __m128d g = _mm_set_pd((double) a1, (double) a0);

    x = _mm_srl_epi64(x, _mm_setr_epi32(2,0,0,0));
    foo = _mm_mullo_epi16(foo, bar);
    foo = _mm_slli_epi64(foo, 1);
    foo = _mm_xor_si128(bar, _mm_unpacklo_epi32(foo, bar));
    foo = _mm_srli_epi64(foo, 1);
    foo = _mm_mullo_epi16(foo, bar);
    foo = _mm_shuffle_epi32(foo, 78);
    foo = _mm_xor_si128(bar, _mm_unpacklo_epi32(foo, bar));
    foo = _mm_srli_si128(foo, 1);
    foo = _mm_xor_si128(foo, x);
    foo = _mm_xor_si128(foo, _mm_castpd_si128 (g));

    return _mm_extract_epi16(foo, 0) & (argc - 1);
}
