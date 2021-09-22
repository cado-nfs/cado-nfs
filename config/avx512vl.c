#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>

__mmask8 foo(__m256i a, __m256i b)
{
    return _mm256_cmpneq_epu32_mask(a, b);
}

int main() {
    __m512i * x = _mm_malloc(120 * sizeof(__m512i), sizeof(__m512i));
    memset(x, 0, 120 * sizeof(__m512i));
    for(int i = 0 ; i < 100 ; i += 10) {
        for(int j = 0 ; j < 10 ; j++) {
            x[i+j]=_mm512_add_epi64(x[i+j], x[i+j+1]);
            _mm512_storeu_si512 (x + i + j + 10, x[i + j]);
        }
        _mm_sfence();
    }
    _mm_empty();
    _mm_free(x);
    return 0;
}
