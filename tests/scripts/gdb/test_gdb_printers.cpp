#include "cado.h" // IWYU pragma: keep

#include <gmp.h>
#include "cxx_mpz.hpp"
#include "mpz_mat.h"
#include "mpz_poly.h"
// #include "mpz_poly_bivariate.hpp"
#include "macros.h"

static void foo()
{
}

int main(int argc, char const * argv[])
{
    cxx_mpz A = 1728;
    if (argc >= 2)
        mpz_set_str(A, argv[1], 0);

    cxx_mpz & rA = A;
    cxx_mpz const & rcA = A;
    ASSERT_ALWAYS(rA == A);
    ASSERT_ALWAYS(rcA == A);


    cxx_mpq C(17, 42);
    ASSERT_ALWAYS(mpq_cmp_ui(C, 0, 1) > 0);

    cxx_mpz_mat M0(3, 3);
    for(int i = 0 ; i < 3 ; i++)
        for(int j = 0 ; j < 3 ; j++)
            mpz_set_ui(mpz_mat_entry(M0, i, j), 10 * i + j);

    cxx_mpq_mat M1;
    mpq_mat_set_mpz_mat_denom(M1, M0, cxx_mpz(17));

    cxx_mpz_poly f("1+2*x-(x-1)^2");

    // cxx_mpz_poly_bivariate ff("1+2*x*y-(x-y)^2");

    foo();
}
