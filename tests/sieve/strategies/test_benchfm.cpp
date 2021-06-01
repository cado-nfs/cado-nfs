#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <gmp.h>
#include "generate_factoring_method.hpp"
#include "fm.h"                           // for fm_create, fm_free, fm_set_...
#include "tab_fm.h"                       // for tabular_fm_add_fm, tabular_...

// coverity[root_function]
int main ()
{
    gmp_randstate_t state;
    gmp_randinit_default(state);
    mpz_t seed;
    mpz_init_set_ui(seed, 42);
    gmp_randseed(state, seed);
    mpz_clear(seed);

    tabular_fm_t* tab = tabular_fm_create ();

    // We exercise the code, without really checking that it works
    fm_t* fm = fm_create();
    unsigned long elem[4] = {1, 0, 20, 100};
    fm_set_method (fm, elem, 4);
    tabular_fm_add_fm (tab, fm);
    bench_proba(state, tab, 20, 20, 5);
    bench_time(state, tab, 100);

    fm_free (fm);
    tabular_fm_free (tab);
    gmp_randclear(state);
    return EXIT_SUCCESS;
}
