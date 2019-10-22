#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>

int main() {
    unsigned __int128 r;
    const uint64_t a = UINT64_C(0x123456789ABCDEF), b = UINT64_C(0x4242424242424242);
    uint64_t r0, r1;
    
    r = (unsigned __int128) a * b;
    r0 = r;
    r1 = r >> 64;
    
    if (r0 == UINT64_C(0x782D15307F00B59E) && r1 == UINT64_C(0x4B6347F977C2DA)) {
        exit(EXIT_SUCCESS);
    } else { 
        fprintf(stderr, "r0 = %" PRIu64 ", r1 = %" PRIu64 "\n", r0, r1);
        exit(EXIT_FAILURE);
    }
}
