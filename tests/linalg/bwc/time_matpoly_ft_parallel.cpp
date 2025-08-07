#include "cado.h" // IWYU pragma: keep

#include <gmp.h>

#include "lingen_matpoly_ft.hpp"
#include "gmp_aux.h"
#include "lingen_substep_characteristics.hpp"
#include "lingen_mul_substeps.hpp"
#include "lingen_substep_schedule.hpp"
#include "lingen_fft_select.hpp"
#include "timing.h" // seconds
#include "params.h"

template<typename fft_type>
struct matpoly_checker_ft {
    static constexpr bool is_binary = is_binary_fft<fft_type>::value;
    using matpoly_type = matpoly<is_binary>;
    matpoly_type::arith_hard ab;
    unsigned int m;
    unsigned int n;
    unsigned int L;

    cxx_gmp_randstate rstate;
    unsigned long seed;

    matpoly_type::memory_guard dummy;
    typename matpoly_ft<fft_type>::memory_guard dummy_ft;

    matpoly_checker_ft(cxx_mpz const & p, unsigned int m, unsigned int n, unsigned int L, cxx_gmp_randstate & rstate0)
        : ab(p, 1U)
        , m(m)
        , n(n)
        , L(L)
        , seed(gmp_urandomm_ui(rstate0, ULONG_MAX))
        , dummy(SIZE_MAX)
        , dummy_ft(SIZE_MAX)
    {
        gmp_randseed_ui(rstate, seed);
    }
    private:
    static int max_threads() {
#ifdef HAVE_OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }
    public:

    void doit_mul(lingen_platform const & P, std::ostream& os) {
        using X_t = lingen_substep_characteristics<is_binary>;
        using op_t = op_mul<is_binary, fft_type>;
        size_t LE = L / 2;
        size_t Lpi = m * LE / (m + n);
        X_t X(&ab, rstate,
                0,      /* not used -- this field should probably go away
                           anyway */
                op_mul_or_mp_base::OP_MUL,
                m+n, m+n, m+n,
                Lpi, Lpi, Lpi+Lpi-1);
        auto op_generic = X.instantiate(encode_fft_type<fft_type>);
        auto const & op = dynamic_cast<op_t const &>(*op_generic);
        X.report_size_stats_human(os, op);
        constexpr unsigned int mesh = 1;
        X.fill_tvec(os, microbench_dft<op_t>(op, P, mesh, X));
        X.fill_tvec(os, microbench_ift<op_t>(op, P, mesh, X));
        X.fill_tvec(os, microbench_conv<op_t>(op, P, mesh, X));
    }

    void doit_mp(lingen_platform const & P, std::ostream& os) {
        using X_t = lingen_substep_characteristics<is_binary>;
        using op_t = op_mul<is_binary, fft_type>;
        size_t LE = L / 2;
        size_t Lpi = m * LE / (m + n);
        X_t X(&ab, rstate,
                0,      /* not used -- this field should probably go away
                           anyway */
                op_mul_or_mp_base::OP_MP,
                m, m+n, m+n,
                LE + Lpi, Lpi, LE + 1);
        auto op_generic = X.instantiate(encode_fft_type<fft_type>);
        auto const & op = dynamic_cast<op_t const &>(*op_generic);
        X.report_size_stats_human(os, op);
        constexpr unsigned int mesh = 1;
        X.fill_tvec(os, microbench_dft<op_t>(op, P, mesh, X));
        X.fill_tvec(os, microbench_ift<op_t>(op, P, mesh, X));
        X.fill_tvec(os, microbench_conv<op_t>(op, P, mesh, X));
    }

    void doit(lingen_platform const & P, std::ostream& os) {
        doit_mp(P, os);
        doit_mul(P, os);
    }
};

template<bool is_binary>
static void declare_usage(cxx_param_list & pl)
{
    if constexpr (!is_binary)
        param_list_decl_usage(pl, "prime", "(mandatory) prime defining the base field");
    else
        param_list_decl_usage(pl, "prime", "(unused) prime defining the base field -- we only use 2");
    param_list_decl_usage(pl, "m", "dimension m");
    param_list_decl_usage(pl, "n", "dimension n");
    param_list_decl_usage(pl, "L", "length (step length in recursive algorithm)");
    param_list_decl_usage(pl, "seed", "random seed");
    lingen_platform::declare_usage(pl);
}

int main(int argc, char const * argv[])
{
#ifdef LINGEN_BINARY
    constexpr bool is_binary = true;
#else
    constexpr bool is_binary = false;
#endif
    MPI_Init(&argc, (char ***) &argv);

    cxx_mpz p;
    cxx_gmp_randstate rstate;

    unsigned int m = 4;
    unsigned int n = 2;
    unsigned int L = 10000;
    unsigned long seed = 0;

    cxx_param_list pl;

    declare_usage<is_binary>(pl);
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
        if (!param_list_parse_mpz(pl, "prime", (mpz_ptr) p)) {
            fprintf(stderr, "--prime is mandatory\n");
            param_list_print_command_line (stdout, pl);
            exit(EXIT_FAILURE);
        }
    } else {
        mpz_set_ui(p, 2);   /* unused anyway */
        param_list_parse_mpz(pl, "prime", (mpz_ptr) p);
    }
    param_list_parse_uint(pl, "m", &m);
    param_list_parse_uint(pl, "n", &n);
    param_list_parse_uint(pl, "L", &L);
    param_list_parse_ulong(pl, "seed", &seed);

    lingen_platform const P(MPI_COMM_WORLD, pl);

    if (param_list_warn_unused(pl))
        exit(EXIT_FAILURE);
    if constexpr (is_binary) {
        if (m & 63) {
            unsigned int const nm = 64 * iceildiv(m, 64);
            printf("Round m=%u to m=%u\n", m, nm);
            m = nm;
        }
        if (n & 63) {
            unsigned int const nn = 64 * iceildiv(n, 64);
            printf("Round n=%u to n=%u\n", n, nn);
            n = nn;
        }
    }

    gmp_randseed_ui(rstate, seed);

#ifdef LINGEN_BINARY
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

    MPI_Finalize();
}
