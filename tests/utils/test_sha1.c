#include "cado.h" // IWYU pragma: keep
#include <stdint.h>     // uint8_t
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "sha1.h"

struct {
    const char * input;
    const char * digest;
} test_vectors[] = {
    { "abc", "a9993e364706816aba3e25717850c26c9cd0d89d" },
};

void stdio_util(FILE * in, FILE * out) {
    SHA1_CTX ctx;
    uint8_t mac[20];
    SHA1Init(&ctx);
    for (int x ; (x = getc(in)) != EOF ; ) {
        unsigned char c = x;
        SHA1Update(&ctx, &c, 1);
    }
    SHA1Final(mac, &ctx);
    for(int i = 0 ; i < 20 ; i++) {
        printf("%02x", (unsigned int) mac[i]);
    }
    fputs("\n", out);
}

int main(int argc, char * argv[])
{
    if (argc == 2 && strcmp(argv[1], "--test") == 0) {
        unsigned int ntests = sizeof(test_vectors)/sizeof(test_vectors[0]);
        for(unsigned int i = 0 ; i < ntests ; ++i) {
            SHA1_CTX ctx;
            uint8_t mac[20];
            const char * in = test_vectors[i].input;
            const char * xout = test_vectors[i].digest;
            SHA1Init(&ctx);
            SHA1Update(&ctx, (const unsigned char *) in, strlen(in));
            SHA1Final(mac, &ctx);
            char out[41];
            for(int i = 0 ; i < 20 ; i++) {
                snprintf(out + 2*i, 3, "%02x", (unsigned int) mac[i]);
            }
            if (strcmp(xout, out) != 0) {
                fprintf(stderr, "test failure, computed sha1(%s)=%s, expected %s\n",
                        in, out, xout);
                exit(EXIT_FAILURE);
            }
        }
    } else {
        stdio_util(stdin, stdout);
    }

    return 0;
}
