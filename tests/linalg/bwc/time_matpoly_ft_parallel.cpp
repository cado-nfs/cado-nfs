#include "cado.h"
#include <gmp.h>
#include "lingen_matpoly_ft.hpp"
#include "gmp_aux.h"
#include "lingen_substep_characteristics.hpp"
#include "timing.h" // seconds

template<typename fft_type>
struct matpoly_checker_ft {
    abfield ab;
    unsigned int m;
    unsigned int n;
    unsigned int L;

    gmp_randstate_t rstate;
    unsigned long seed;

    matpoly::memory_guard dummy;
    typename matpoly_ft<fft_type>::memory_guard dummy_ft;

    matpoly_checker_ft(cxx_mpz const & p, unsigned int m, unsigned int n, unsigned int L, gmp_randstate_t rstate0)
        : m(m)
        , n(n)
        , L(L)
        , seed(gmp_urandomm_ui(rstate0, ULONG_MAX))
        , dummy(SIZE_MAX)
        , dummy_ft(SIZE_MAX)
    {
        abfield_init(ab);
        abfield_specify(ab, MPFQ_PRIME_MPZ, (mpz_srcptr) p);
        gmp_randinit_default(rstate);
        gmp_randseed_ui(rstate, seed);
    }
    ~matpoly_checker_ft() {
        gmp_randclear(rstate);
        abfield_clear(ab);
    }
    private:
    static inline int max_threads() {
#ifdef HAVE_OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }
    public:

    void doit_mul(lingen_platform const & P, std::ostream& os) {
        typedef lingen_substep_characteristics<op_mul<fft_type>> X_t;
        size_t LE = L / 2;
        size_t Lpi = m * LE / (m + n);
        X_t X(ab, rstate,
                0,      /* not used -- this field should probably go away
                           anyway */
                m+n, m+n, m+n,
                Lpi, Lpi, Lpi+Lpi-1);
        X.report_size_stats_human(os);
        X.fill_tvec(os, typename X_t::microbench_dft(P, X));
        X.fill_tvec(os, typename X_t::microbench_ift(P, X));
        X.fill_tvec(os, typename X_t::microbench_conv(P, X));
    }

    void doit_mp(lingen_platform const & P, std::ostream& os) {
        typedef lingen_substep_characteristics<op_mp<fft_type>> X_t;
        size_t LE = L / 2;
        size_t Lpi = m * LE / (m + n);
        X_t X(ab, rstate,
                0,      /* not used -- this field should probably go away
                           anyway */
                m, m+n, m+n,
                LE + Lpi, Lpi, LE + 1);
        X.report_size_stats_human(os);
        X.fill_tvec(os, typename X_t::microbench_dft(P, X));
        X.fill_tvec(os, typename X_t::microbench_ift(P, X));
        X.fill_tvec(os, typename X_t::microbench_conv(P, X));
    }

    void doit(lingen_platform const & P, std::ostream& os) {
        doit_mp(P, os);
        doit_mul(P, os);
    }
};

void declare_usage(cxx_param_list & pl)
{
#ifndef SELECT_MPFQ_LAYER_u64k1
    param_list_decl_usage(pl, "prime", "(mandatory) prime defining the base field");
#else
    param_list_decl_usage(pl, "prime", "(unused) prime defining the base field -- we only use 2");
#endif
    param_list_decl_usage(pl, "m", "dimension m");
    param_list_decl_usage(pl, "n", "dimension n");
    param_list_decl_usage(pl, "L", "length (step length in recursive algorithm)");
    param_list_decl_usage(pl, "seed", "random seed");
    lingen_platform::declare_usage(pl);
}

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);

    cxx_mpz p;
    gmp_randstate_t rstate;

    unsigned int m = 4;
    unsigned int n = 2;
    unsigned int L = 10000;
    unsigned long seed = 0;

    cxx_param_list pl;

    declare_usage(pl);
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
#ifndef SELECT_MPFQ_LAYER_u64k1
    if (!param_list_parse_mpz(pl, "prime", (mpz_ptr) p)) {
        fprintf(stderr, "--prime is mandatory\n");
        param_list_print_command_line (stdout, pl);
        exit(EXIT_FAILURE);
    }
#else
    mpz_set_ui(p, 2);   /* unused anyway */
    param_list_parse_mpz(pl, "prime", (mpz_ptr) p);
#endif
    param_list_parse_uint(pl, "m", &m);
    param_list_parse_uint(pl, "n", &n);
    param_list_parse_uint(pl, "L", &L);
    param_list_parse_ulong(pl, "seed", &seed);

    lingen_platform P(MPI_COMM_WORLD, pl);

    if (param_list_warn_unused(pl))
        exit(EXIT_FAILURE);
#ifdef SELECT_MPFQ_LAYER_u64k1
    if (m & 63) {
        unsigned int nm = 64 * iceildiv(m, 64);
        printf("Round m=%u to m=%u\n", m, nm);
        m = nm;
    }
    if (n & 63) {
        unsigned int nn = 64 * iceildiv(n, 64);
        printf("Round n=%u to n=%u\n", n, nn);
        n = nn;
    }
#endif

    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, seed);

#ifdef SELECT_MPFQ_LAYER_u64k1
    {
        matpoly_checker_ft<gf2x_fake_fft_info> checker_ft(p, m, n, L, rstate);
        checker_ft.doit(P, std::cout);
    }
    {
        matpoly_checker_ft<gf2x_cantor_fft_info> checker_ft(p, m, n, L, rstate);
        checker_ft.doit(P, std::cout);
    }
    {
        matpoly_checker_ft<gf2x_ternary_fft_info> checker_ft(p, m, n, L, rstate);
        checker_ft.doit(P, std::cout);
    }
#else
    {
        matpoly_checker_ft<fft_transform_info> checker_ft(p, m, n, L, rstate);
        checker_ft.doit(P, std::cout);
    }
#endif

    gmp_randclear(rstate);

    MPI_Finalize();
}
