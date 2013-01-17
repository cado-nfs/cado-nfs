#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "bwc_config.h"
#include "macros.h"
#include "utils.h"

int verbose = 0;

/* Returns something centered on 1. */
static inline double dist_func()
{
    long int y = random();
    y += y == 0;
    double x = (double) y / RAND_MAX;
    /* With this distribution, bizarrely we tend to get a density which
     * is larger than expected... would be accounted for by our
     */
    return (x + 2 / sqrt(x) - 3) / 1.5;
    // return 2 * x;
}

int gen_row(double lambda, int n, int * ptr)
{
    int c = 0;
    for(int i = 0 ; ; ) {
        int e = 1 + (int) (lambda * dist_func() * (i+1));
        i += e;
        c += 1;
        if (i < 0 || i > n)
            break;
        if (ptr) *ptr++ = i - 1;
    }
    return c-1;
}

/* The process above sets i_0 = 0, and:
 * i_k, knowing i_{k-1} is on average i_{k-1} + 1 + (i_{k-1}+1) * lambda.
 * And E(i_k) = E(E(i_k|i_{k-1})). So if u_k is E(i_k) we have:
 * u_k = 1 + lambda + (1+lambda) * u_{k-1}.
 * From which the general expression for u_k follows.
 *
 * We use it boldly to reverse the sense of the expected value, so as to
 * say that the expected density it the first integer k such that the
 * expected value for the k-th pick exceeds n (which looks very wrong).
 */
double expected_density_from_lambda(double lambda, int n)
{
#if 0
    double u = log(1+lambda*n / (1 + lambda));
    double v = log(1+lambda);
    return u/v;
#else
    double r = 0;
    for(int i = 0 ; i < 10 ; i++) {
        r += gen_row(lambda, n, NULL);
    }
    return r / 10.0;
#endif
}

double get_corresponding_lambda(int dens, int n)
{
    double l0 = 1.0e-2;
    double l1 = n;

    double e0, e1;

    ASSERT_ALWAYS((e0 = expected_density_from_lambda(l0,n)) > dens);
    ASSERT_ALWAYS((e1 = expected_density_from_lambda(l1,n)) < dens);
    for( ; l1/l0 > 1.01 ; ) {
        double lc = sqrt(l0 * l1);
        double ec = expected_density_from_lambda(lc,n);
        // printf("l=%.2f e=%.2f\n", lc, ec);
        if (ec > dens) {
            l0 = lc;
            e0 = ec;
        } else {
            l1 = lc;
            e1 = ec;
        }
    }
    return l0;
}

void printrows(int nrows, int ncols, double lambda, int maxc)
{
    int * ptr = malloc((ncols + 1) * sizeof(int));
    int total_coeffs = 0;
    double tot_sq = 0;
    for(int i = 0 ; i < nrows ; i++) {
        int c = gen_row(lambda, ncols, ptr);
        printf("%d", c);
        for(int j = 0 ; j < c ; j++) {
            printf(" %d", ptr[j]);
            if (maxc) {
                int co = rand() % (2 * maxc + 1) - maxc;
                printf(" %d", co);
            }
        }
        printf("\n");
        total_coeffs += c;
        tot_sq += (double) c * (double) c;
    }
    free(ptr);
    if (verbose) {
        double e = (double) total_coeffs / nrows;
        double s = (double) tot_sq / nrows;
        double sdev = sqrt(s - e*e);
        fprintf(stderr, "Actual density per row avg %.2f sdev %.2f\n",
                e, sdev);
    }
}

void usage()
{
    fprintf(stderr, "Usage: ./random <nrows> [<ncols>] [<density>] [options]\n"
            "Options:\n"
            "\t-d <density> : another way for specifying the density per row\n"
            "\t-s <seed> : seed\n"
            "\t-c <maxc> : add coefficients\n"
            "\t-v : turn verbosity on\n"
            "\t--kleft <d>: ensure at least a left kernel of dimension d\n"
            "\t--kright <d>: ditto for right kernel\n");
    exit(1);
}


int main(int argc, char * argv[])
{
    param_list pl;

    param_list_init(pl);

    argv++, argc--;
    int wild = 0;

    int maxcoeff = 0;

    int wild_args[3] = { -1, -1, -1 }; // nrows ncols coeffs_per_row

    param_list_configure_alias(pl, "density", "-d");
    param_list_configure_alias(pl, "seed", "-s");
    param_list_configure_switch(pl, "-v", &verbose);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

        if (argv[0][0] != '-' && wild < 3) {
            char * tmp;
            wild_args[wild] = strtoul(argv[0], &tmp, 0);
            if (*tmp != '\0') {
                fprintf(stderr, "Parse error for parameter %s\n", argv[0]);
                exit(1);
            }
            argv++, argc--;
            wild++;
            continue;
        }

        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage();
    }

    int seed = 0;
    param_list_parse_int(pl, "seed", &seed);

    srand(seed ? seed : time(NULL));

    param_list_parse_int(pl, "c", &maxcoeff);

    int nrows = wild_args[0];
    if (nrows < 0) {
        fprintf(stderr, "Please specify nrows\n");
        exit(1);
    }
    int ncols = wild_args[1];
    if (ncols < 0) {
        ncols = nrows;
    }
    int dens = wild_args[2];      // average coeffs per row
    if (dens < 0) {
        param_list_parse_int(pl, "density", &dens);
    }

    ASSERT_ALWAYS(ncols > 10 && nrows > 10);

    int kernel_left = 0;
    int kernel_right = 1;

    param_list_parse_int(pl, "kleft", &kernel_left);
    param_list_parse_int(pl, "kright", &kernel_right);
    
    // default is lambda = 1
    double lambda = 1;

    if (dens > 0) {
        lambda = get_corresponding_lambda(dens, ncols);
        if (verbose) {
            fprintf(stderr, "Selected lambda=%.2f\n", lambda);
        }
    }

    printf("%d %d\n", nrows, ncols);

    if (ncols > nrows)
        kernel_right -= ncols - nrows;
    if (kernel_right < 0)
        kernel_right = 0;

    if (nrows > ncols)
        kernel_left -= nrows - ncols;
    if (kernel_left < 0)
        kernel_left = 0;

    if (kernel_right > nrows / 4) {
        fprintf(stderr, "Warning, right kernel is large."
                " Could trigger misbehaviours\n");
    }

    if (kernel_left > ncols / 4) {
        fprintf(stderr, "Warning, left kernel is large."
                " Could trigger misbehaviours\n");
    }

    printrows(nrows - kernel_right, ncols - kernel_left, lambda, maxcoeff);
    for(int i = 0 ; i < kernel_right ; i++) {
        printf("0\n");
    }

    param_list_clear(pl);

    return 0;
}
