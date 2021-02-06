/* This source file is our test case for mmx support. If this can be
 */
#include <stdint.h>
#include <stdlib.h>
#include <mmintrin.h>

int main(int argc, char *argv[])
{
    volatile uint64_t a1 = 42;
    __m64 mmask = _mm_set1_pi8(a1);
    __m64 foo = (__m64) a1;

    return _m_to_int(_mm_cmpeq_pi8(foo, mmask)) != 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
