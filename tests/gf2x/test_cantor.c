#include "cado.h" // IWYU pragma: keep
#include <inttypes.h>
#include <stdlib.h>
#include <gmp.h>
#include "macros.h"
#include "gf2x.h"
#include "gf2x-fft.h"
#include "tests_common.h"

/* using the selected FFT engine, multiply two n by n matrices of
 * polynomials with Fl and Gl words in each entry, respectively */
long gf2x_cantor_fft_matmul(unsigned long ** H, unsigned long ** F, size_t Fl, unsigned long ** G, size_t Gl, int n)
{
    gf2x_cantor_fft_info_t order;

    gf2x_cantor_fft_info_init(order, Fl * GF2X_WORDSIZE, Gl * GF2X_WORDSIZE);

    gf2x_cantor_fft_ptr f = gf2x_cantor_fft_alloc(order, n*n);

    size_t sizes[3];
    gf2x_cantor_fft_info_get_alloc_sizes(order, sizes);
    gf2x_cantor_fft_ptr temp1 = malloc(sizes[1] * sizeof(gf2x_cantor_fft_elt));
    gf2x_cantor_fft_ptr temp2 = malloc(sizes[2] * sizeof(gf2x_cantor_fft_elt));

    for(int i = 0 ; i < n ; i++)
        for(int j = 0 ; j < n ; j++)
            gf2x_cantor_fft_dft(order, gf2x_cantor_fft_get(order, f, i*n+j), F[i*n+j], Fl * GF2X_WORDSIZE, temp1);

    gf2x_cantor_fft_ptr g = gf2x_cantor_fft_alloc(order, n*n);

    for(int i = 0 ; i < n ; i++)
        for(int j = 0 ; j < n ; j++)
            gf2x_cantor_fft_dft(order, gf2x_cantor_fft_get(order, g, i*n+j), G[i*n+j], Gl * GF2X_WORDSIZE, temp1);

    gf2x_cantor_fft_ptr h = gf2x_cantor_fft_alloc(order, n*n);

    for(int i = 0 ; i < n ; i++)
        for(int j = 0 ; j < n ; j++) {
            gf2x_cantor_fft_compose(order,
                    gf2x_cantor_fft_get(order, h, i*n+j),
                    gf2x_cantor_fft_get_const(order, (gf2x_cantor_fft_srcptr) f, i*n/*+k*/),
                    gf2x_cantor_fft_get_const(order, (gf2x_cantor_fft_srcptr) g, /*k*n+*/j),
                    temp2);
            for(int k = 1 ; k < n ; k++)
                gf2x_cantor_fft_addcompose(order,
                        gf2x_cantor_fft_get(order, h, i*n+j),
                        gf2x_cantor_fft_get_const(order, (gf2x_cantor_fft_srcptr) f, i*n+k),
                        gf2x_cantor_fft_get_const(order, (gf2x_cantor_fft_srcptr) g, k*n+j), temp2, temp1);
        }

    gf2x_cantor_fft_free(order, g, n*n);
    gf2x_cantor_fft_free(order, f, n*n);

    for(int i = 0 ; i < n ; i++)
        for(int j = 0 ; j < n ; j++)
            gf2x_cantor_fft_ift(order, H[i*n+j], (Fl+Gl) * GF2X_WORDSIZE - 1,  gf2x_cantor_fft_get(order, h, i*n+j), temp1);

    long res = gf2x_cantor_fft_info_order(order);

    free(temp1);
    free(temp2);

    gf2x_cantor_fft_free(order, h, n*n);
    gf2x_cantor_fft_info_clear(order);

    return res;
}



int main(int argc, const char **argv)
{
    const char *argv0 = argv[0];
    tests_common_cmdline(&argc, &argv, PARSE_SEED);
    if (argc != 2) {
        fprintf(stderr, "usage: %s <option> N\n", argv0);
        fprintf(stderr,
                "  where N is the number of limbs of operands\n");
        exit(1);
    }
    int N = atoi(argv[1]);

    if (N <= 0)
        return EXIT_FAILURE;

#ifdef __COVERITY__
    __coverity_mark_pointee_as_sanitized__(p, LOOP_BOUND);
#endif
    
    unsigned long *f, *g, *h1, *h2;
    f = (unsigned long *) malloc(N * sizeof(unsigned long));
    g = (unsigned long *) malloc(N * sizeof(unsigned long));
    h1 = (unsigned long *) malloc(2 * N * sizeof(unsigned long));
    h2 = (unsigned long *) malloc(2 * N * sizeof(unsigned long));
    ASSERT_ALWAYS(f != NULL);
    ASSERT_ALWAYS(g != NULL);
    ASSERT_ALWAYS(h1 != NULL);
    ASSERT_ALWAYS(h2 != NULL);

    for (int i = 0; i < N; ++i) {
        f[i] = gmp_urandomb_ui(state, 64);
        g[i] = gmp_urandomb_ui(state, 64);
    }

    gf2x_cantor_fft_matmul(&h1,&f,N,&g,N,1);
    gf2x_mul(h2,f,N,g,N);
    for (int i = 0; i < 2*N; ++i) {
        ASSERT_ALWAYS(h1[i] == h2[i]);
    }

    free(f);
    free(g);
    free(h1);
    free(h2);

    return EXIT_SUCCESS;
}
