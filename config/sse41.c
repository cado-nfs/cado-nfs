#include <stdint.h>
#include <smmintrin.h>

int main() {
    volatile __v2di x = { (uint64_t) 42, (uint64_t) 17 };
    volatile __v2di y = { (uint64_t) 41, (uint64_t) 17 };
    x = _mm_cmpeq_epi64(x, y);
    /* the following test is for emulated 32-bit on physical 64-bit */
    if (sizeof(unsigned long) != 8)
      abort ();
    return 0;
}

