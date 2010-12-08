#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <string.h>

/* This is a stupid program for checking a matrix-times-vector product:
 *
 * ./a.out <matrix file> <vector file> : compute M * v
 * ./a.out -t <matrix file> <vector file> : compute v * M
 *
 * dimensions are guessed from the vector size. data files are expected
 * in binary format, and printed to stdout in ascii hex.
 */

int main(int argc, char * argv[])
{
    struct stat sbuf[1];
    int mul_left = 0;
    for( ; argc > 3 ; ) {
        if (strcmp(argv[1], "-t") == 0) {
            mul_left = 1;
            argv++, argc--;
        }
    }
    if (argc != 3) {
        fprintf(stderr, "Usage: short-matmul <matrix file> <vector file>\n");
        exit(1);
    }
    FILE * m = fopen(argv[1], "r");
    if (m == NULL) {
        perror(argv[1]);
        exit(1);
    }
    if (stat(argv[2], sbuf) < 0) {
        perror(argv[2]);
        exit(1);
    }
    unsigned long ncols = sbuf->st_size;
    unsigned long nrows;
    assert(ncols % sizeof(uint64_t) == 0);
    ncols /= sizeof(uint64_t);
    FILE * v = fopen(argv[2], "r");
    if (v == NULL) {
        perror(argv[1]);
        exit(1);
    }
    if (mul_left == 0) {
        fprintf(stderr, "Found %lu columns in %s\n", ncols, argv[2]);
        nrows = 0;
    } else {
        fprintf(stderr, "Found %lu rows in %s\n", ncols, argv[2]);
        nrows = ncols;
        ncols += ncols / 10;
    }
    uint64_t * vec;
    vec = malloc(ncols * sizeof(uint64_t));
    memset(vec, 0, ncols * sizeof(uint64_t));
    assert(vec);
    
    if (mul_left == 0) {
        unsigned long l = fread(vec, sizeof(uint64_t), ncols, v);
        if (l < ncols) {
            fprintf(stderr, "%s: short read\n", argv[2]);
            exit(1);
        }
        fclose(v);
        for( ; ; ) {
            uint32_t rowlen;
            int k = fread(&rowlen, sizeof(uint32_t), 1, m);
            if (k != 1) {
                if (!feof(m)) {
                    fprintf(stderr, "%s: short read\n", argv[1]);
                    exit(1);
                }
                break;
            }
            uint64_t w = 0;
            for( ; rowlen-- ; ) {
                uint32_t c;
                int k = fread(&c, sizeof(uint32_t), 1, m);
                if (k != 1) {
                    fprintf(stderr, "%s: short read\n", argv[1]);
                    exit(1);
                }
                assert(c < ncols);
                w ^= vec[c];
            }
            printf("%016" PRIx64 "\n", w);
        }
    } else {
        unsigned long rrows;
        for(rrows = 0 ; ; rrows++) {
            uint32_t rowlen;
            int k = fread(&rowlen, sizeof(uint32_t), 1, m);
            if (k != 1) {
                if (!feof(m)) {
                    fprintf(stderr, "%s: short read\n", argv[1]);
                    exit(1);
                }
                break;
            }
            uint64_t w = 0;
            k = fread(&w, sizeof(uint64_t), 1, v);
            if (k != 1) {
                fprintf(stderr, "%s: short read\n", argv[1]);
                exit(1);
            }
            for( ; rowlen-- ; ) {
                uint32_t c;
                int k = fread(&c, sizeof(uint32_t), 1, m);
                if (k != 1) {
                    fprintf(stderr, "%s: short read\n", argv[1]);
                    exit(1);
                }
                if (c >= ncols) {
                    unsigned long old_ncols = ncols;
                    ncols += ncols / 10;
                    fprintf(stderr, "re-allocating vec array, now %lu cols\n", ncols);
                    vec = realloc(vec, ncols * sizeof(uint64_t));
                    memset(vec + old_ncols, 0, (ncols - old_ncols) * sizeof(uint64_t));
                }
                vec[c] ^= w;
            }
        }
        assert(nrows == rrows);
        unsigned long k;
        for(k = 0 ; k < nrows ; k++) {
            printf("%016" PRIx64 "\n", vec[k]);
        }
    }

    fclose(m);
    free(vec);
}

