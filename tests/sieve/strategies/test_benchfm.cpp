#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdio>

#include <gmp.h>

#include "generate_factoring_method.hpp"
#include "fm.hpp"
#include "tab_fm.hpp"

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
    for(unsigned int i = 0 ; i < tab->size ; i++) {
        printf("method %u:", i);
        for(unsigned int j = 0 ; j < tab->tab[i]->len_time ; j++) {
            printf(" %3.2f", tab->tab[i]->time[j]);
        }
        printf("\n");
    }

    fm_free (fm);
    tabular_fm_free (tab);
    gmp_randclear(state);
    return EXIT_SUCCESS;
}
