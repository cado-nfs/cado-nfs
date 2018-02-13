#include <arm_neon.h>
#include <stdlib.h>
#include <stdio.h>

int do_test()
{
    const uint32x4_t a = {1, 17, 42, 4294967000U}, b = {5, 4, 3, 1000};
    uint32x4_t c;
    int ok;
    c = vaddq_u32(a, b);
    ok = (c[0] == 6 && c[1] == 21 && c[2] == 45 && c[3] == 704);
    if (!ok) {
        fprintf(stderr, "Error, c = {%d, %d, %d, %d}\n",
            c[0], c[1], c[2], c[3]);
    }
    return ok;
}

int main()
{
    if (do_test())
        exit(EXIT_SUCCESS);
    else
        exit(EXIT_FAILURE);
}
