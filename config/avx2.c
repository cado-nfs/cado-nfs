#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>

int main() {
    __m256i * x = _mm_malloc(120 * sizeof(__m256i), sizeof(__m256i));
    memset(x, 0, 120 * sizeof(__m256i));
    for(int i = 0 ; i < 100 ; i += 10) {
        for(int j = 0 ; j < 10 ; j++) {
            x[i+j]=_mm256_add_epi64(x[i+j], x[i+j+1]);
            _mm256_stream_si256 (x + i + j + 10, x[i + j]);
        }
        _mm_sfence();
    }
    _mm_empty();
    _mm_free(x);
    return 0;
}
