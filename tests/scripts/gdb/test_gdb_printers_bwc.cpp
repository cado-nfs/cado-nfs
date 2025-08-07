#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include "bw-common.h"
#include "cxx_mpz.hpp"
#include "lingen_matpoly_select.hpp"
#include "macros.h"
#include "params.h"

static void foo()
{
}

int main(int argc, char const * argv[])
{
#ifdef LINGEN_BINARY
    constexpr bool is_binary = true;
#else
    constexpr bool is_binary = false;
#endif
    cxx_mpz p;
    cxx_param_list pl;

    bw_common_init(bw, &argc, &argv);

    const char * argv0 = argv[0];
    argv++,argc--;
    /* read all command-line parameters */
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf(stderr, "Unexpected argument %s\n", argv[0]);
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }
    if constexpr (!is_binary) {
        if (!param_list_parse(pl, "prime", p)) {
            fprintf(stderr, "--prime is mandatory\n");
            param_list_print_command_line (stdout, pl);
            exit(EXIT_FAILURE);
        }
    } else {
        p = 2;
    }
    unsigned int m = is_binary ? 64 : 3;
    unsigned int n = is_binary ? 128 : 5;

    param_list_parse(pl, "m", m);
    param_list_parse(pl, "n", n);
    if (param_list_warn_unused(pl))
        exit(EXIT_FAILURE);

    matpoly<is_binary>::arith_hard K(p, 1U);
    matpoly<is_binary> M(&K, m, n, 0);

    M.realloc(3);
    M.zero();
    M.set_size(3);

#if defined(ARITH_MODP)
    for(unsigned int i = 0 ; i < m ; i++)
        for(unsigned int j = 0 ; j < n ; j++)
            K.set(*M.part(i, j), i * 10 + j);
#elif defined(ARITH_MOD2)
    for(unsigned int i = 0 ; i < m ; i++)
        for(unsigned int j = 0 ; j < n ; j++)
            *M.part(i, j) = i * 10 + j;
#endif
    M.multiply_column_by_x(0, 2);
    M.multiply_column_by_x(1, 2);
    M.multiply_column_by_x(1, 2);

    ASSERT_ALWAYS(M.nrows() == m);
    foo();

    bw_common_clear(bw);
}
