#include "cado.h" // IWYU pragma: keep

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <complex.h>

#include "polyroots.h"

int main(int argc, char const * argv[])
{
    int degree = argc - 2;
    printf("degree %d\n", degree);
    double coeffs[MAX_POLYROOTS_ROOTFINDER_DEGREE] = {0,};
    for(int i = 0 ; i <= degree ; i++) {
        coeffs[i] = atof(argv[1+i]);
    }
    double _Complex roots[MAX_POLYROOTS_ROOTFINDER_DEGREE];
    poly_roots_double(coeffs, degree, roots);
    for(int i = 0 ; i < degree ; i++) {
        printf("%f+i*%f\n", creal(roots[i]), cimag(roots[i]));
    }
}
