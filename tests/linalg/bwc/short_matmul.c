#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <gmp.h>

#include "gmp-hacks.h"
#include "gmp_aux.h"
#include "macros.h"

/* This is a stupid program for checking a matrix-times-vector product:
 *
 * ./a.out <matrix file> <vector file> : compute M * v
 * ./a.out -t <matrix file> <vector file> : compute v * M
 *
 * dimensions are guessed from the vector size. data files are expected
 * in binary format, and printed in the same binary format to stdout.
 * Pipe through xxd -ps -c8 to get binary hex (or xxd -ps -c8 | uniq -c
 * to verify that you reach 0 as expected for the final kernel vector).
 *
 * For non-binary matrices, do:
 *
 * ./a.out -p <prime> <matrix file> <vector file> : compute M * v
 * ./a.out -p <prime> -t <matrix file> <vector file> : compute v * M
 *
 * setting prime>2 will cause the matrix to be read with coefficients (as
 * signed 32-bit integers). The result will then be reduced modulo the
 * prime.
 *
 * (setting prime==2 has the same effect as not setting it at all).
 */

void short_matmul_binary(FILE * out, FILE * v,  const char * uri, int mul_left)
{
    struct stat sbuf[1];
    if (fstat(fileno(v), sbuf) < 0) {
        perror("input vector");
        exit(1);
    }

    unsigned long ncols = sbuf->st_size;
    unsigned long nrows;
    ASSERT(ncols % sizeof(uint64_t) == 0);
    ncols /= sizeof(uint64_t);
    if (mul_left == 0) {
        fprintf(stderr, "Found %lu columns in input vector\n", ncols);
        nrows = 0;
    } else {
        fprintf(stderr, "Found %lu rows in input vector\n", ncols);
        nrows = ncols;
        ncols += ncols / 10;
    }
    uint64_t * vec;
    vec = malloc(ncols * sizeof(uint64_t));
    memset(vec, 0, ncols * sizeof(uint64_t));
    ASSERT(vec);

    FILE * m = fopen(uri, "rb");

    if (mul_left == 0) {
        unsigned long l = fread(vec, sizeof(uint64_t), ncols, v);
        if (l < ncols) {
            fprintf(stderr, "input vector: short read\n");
            exit(1);
        }
        fclose(v);
        for( ; ; ) {
            uint32_t x;
            int k = fread(&x, sizeof(uint32_t), 1, m);
            uint32_t rowlen = x;
            if (k != 1) {
                if (k) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                break;
            }
            uint64_t w = 0;
            for(uint32_t i = 0 ; i < rowlen ; i++) {
                int k = fread(&x, sizeof(uint32_t), 1, m);
                if (k != 1) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                ASSERT(x < ncols);
                w ^= vec[x];
            }
            // printf("%016" PRIx64 "\n", w);
            fwrite(&w, sizeof(uint64_t), 1, out);
        }
    } else {
        unsigned long rrows;
        for(rrows = 0 ; ; rrows++) {
            uint32_t x;
            int k = fread(&x, sizeof(uint32_t), 1, m);
            uint32_t rowlen = x;
            if (k != 1) {
                if (k) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                break;
            }
            uint64_t w = 0;
            k = fread(&w, sizeof(uint64_t), 1, v);
            if (k != 1) {
                fprintf(stderr, "matrix: short read\n");
                exit(1);
            }
            for(uint32_t i = 0 ; i < rowlen ; i++) {
                // coverity[tainted_argument]
                k = fread(&x, sizeof(uint32_t), 1, m);
                if (k != 1) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                if (x >= ncols) {
                    unsigned long old_ncols = ncols;
                    for( ; x >= ncols ; ncols += ncols / 10);
                    fprintf(stderr, "re-allocating vec array, now %lu cols\n", ncols);
                    vec = realloc(vec, ncols * sizeof(uint64_t));
                    memset(vec + old_ncols, 0, (ncols - old_ncols) * sizeof(uint64_t));
                }
                vec[x] ^= w;
            }
        }
        fclose(v);
        ASSERT_ALWAYS(nrows == rrows);
        unsigned long k;
        for(k = 0 ; k < nrows ; k++) {
            uint64_t w = vec[k];
            // printf("%016" PRIx64 "\n", w);
            fwrite(&w, sizeof(uint64_t), 1, out);
        }
    }

    fclose(m);
    free(vec);
}

void short_matmul_prime(FILE * out, FILE * v,  const char * uri, int mul_left, mpz_srcptr p)
{
    struct stat sbuf[1];
    if (fstat(fileno(v), sbuf) < 0) {
        perror("input vector");
        exit(1);
    }

    unsigned long ncols = sbuf->st_size;
    unsigned long nrows;
    size_t coeffsize = sizeof(mp_limb_t) * mpz_size(p);
    ASSERT(ncols % coeffsize == 0);
    ncols /= coeffsize;
    if (mul_left == 0) {
        fprintf(stderr, "Found %lu columns in input vector\n", ncols);
        nrows = 0;
    } else {
        fprintf(stderr, "Found %lu rows in input vector\n", ncols);
        nrows = ncols;
        ncols += ncols / 10;
    }
    mp_limb_t * vec = malloc(ncols * coeffsize);
    ASSERT(vec);
    memset(vec, 0, ncols * coeffsize);

    mp_limb_t * varea = malloc(coeffsize);
    ASSERT(varea);
    memset(varea, 0, coeffsize);

    FILE * m = fopen(uri, "rb");

    mpz_t tmp, w;
    mpz_init(tmp);
    mpz_init(w);
    if (mul_left == 0) {
        unsigned long l = fread(vec, coeffsize, ncols, v);
        if (l < ncols) {
            fprintf(stderr, "input vector: short read\n");
            exit(1);
        }
        fclose(v);
        for( ; ; ) {
            uint32_t x;
            int32_t c;
            // coverity[tainted_argument]
            int k = fread(&x, sizeof(uint32_t), 1, m);
            uint32_t rowlen = x;
            if (k != 1) {
                if (k) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                break;
            }
            mpz_set_ui(w, 0);
            for(uint32_t i = 0 ; i < rowlen ; i++) {
                k = fread(&x, sizeof(uint32_t), 1, m);
                if (k != 1) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                k = fread(&c, sizeof(int32_t), 1, m);
                if (k != 1) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                ASSERT(x < ncols);
                MPZ_SET_MPN(tmp, vec + x * mpz_size(p), mpz_size(p));
                mpz_addmul_si(w, tmp, c);
            }
            mpz_mod(w, w, p);
            MPN_SET_MPZ(varea, mpz_size(p), w);
            // printf("%016" PRIx64 "\n", w);
            fwrite(varea, coeffsize, 1, out);
        }
    } else {
        unsigned long rrows;
        for(rrows = 0 ; ; rrows++) {
            uint32_t x;
            int32_t c;
            int k = fread(&x, sizeof(uint32_t), 1, m);
            uint32_t rowlen = x;
            if (k != 1) {
                if (k) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                break;
            }
            k = fread(varea, coeffsize, 1, v);
            if (k != 1) {
                fprintf(stderr, "matrix: short read\n");
                exit(1);
            }
            MPZ_SET_MPN(w, varea, mpz_size(p));
            for(uint32_t i = 0 ; i < rowlen ; i++) {
                k = fread(&x, sizeof(uint32_t), 1, m);
                if (k != 1) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                k = fread(&c, sizeof(int32_t), 1, m);
                if (k != 1) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                if (x >= ncols) {
                    unsigned long old_ncols = ncols;
                    for( ; x >= ncols ; ncols += ncols / 10);
                    fprintf(stderr, "re-allocating vec array, now %lu cols\n", ncols);
                    vec = realloc(vec, ncols * coeffsize);
                    memset(vec + old_ncols * mpz_size(p),
                            0, (ncols - old_ncols) * coeffsize);
                }
                MPZ_SET_MPN(tmp, vec + x * mpz_size(p), mpz_size(p));
                mpz_addmul_si(tmp, w, c);
                mpz_mod(tmp, tmp, p);
                MPN_SET_MPZ(vec + x * mpz_size(p), mpz_size(p), tmp);
            }
        }
        fclose(v);
        ASSERT_ALWAYS(nrows == rrows);
        unsigned long k;
        for(k = 0 ; k < nrows ; k++) {
            // printf("%016" PRIx64 "\n", w);
            fwrite(vec + k * mpz_size(p),coeffsize, 1, out);
        }
    }
    mpz_clear(tmp);
    mpz_clear(w);

    fclose(m);
    free(varea);
    free(vec);
}

void usage()
{
    fprintf(stderr, "Usage: short-matmul [-t] <matrix file> <vector file> <vector output file>\n");
    exit(1);
}
int main(int argc, char * argv[])
{
    int mul_left = 0;
    mpz_t p;
    mpz_init_set_ui(p, 2);
    for( ; argc > 4 ; ) {
        if (strcmp(argv[1], "-t") == 0) {
            mul_left = 1;
            argv++, argc--;
        } else if (strcmp(argv[1], "-p") == 0) {
            argv++, argc--;
            mpz_set_str(p, argv[1], 0);
            argv++, argc--;
        } else {
            usage();
        }
    }
    if (argc != 4) {
        usage();
    }
    FILE * v = fopen(argv[2], "rb");
    if (v == NULL) {
        perror(argv[2]);
        exit(1);
    }
    

    FILE * w = fopen(argv[3], "wb");
    if (w == NULL) {
        perror(argv[2]);
        exit(1);
    }

    /* fclose(v) is done within short_matmul */
    if (mpz_cmp_ui(p, 2) == 0) {
        short_matmul_binary(w, v, argv[1], mul_left);
    } else {
        short_matmul_prime(w, v, argv[1], mul_left, p);
    }
    fclose(w);
    mpz_clear(p);
}

