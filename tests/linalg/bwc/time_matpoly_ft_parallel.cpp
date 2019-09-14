#include "cado.h"
#include <gmp.h>
#include "lingen_matpoly_ft.hpp"
#include "utils.h"
#include "gmp_aux.h"
#include "lingen_substep_characteristics.hpp"

struct matpoly_checker_base {
    abfield ab;
    unsigned int m;
    unsigned int n;
    unsigned int len1;
    unsigned int len2;
    matpoly::memory_guard dummy;

    gmp_randstate_t rstate;
    /* tests are free to seed and re-seed the checker's private rstate, a
     * priori with this random seed which was taken from the initial
     * random state */
    unsigned long seed;

    matpoly_checker_base(cxx_mpz const & p, unsigned int m, unsigned int n, unsigned int len1, unsigned int len2, gmp_randstate_t rstate0)
        : m(m)
        , n(n)
        , len1(len1)
        , len2(len2)
        , dummy(SIZE_MAX)
        , seed(gmp_urandomm_ui(rstate0, ULONG_MAX))
    {
        abfield_init(ab);
        abfield_specify(ab, MPFQ_PRIME_MPZ, (mpz_srcptr) p);
        gmp_randinit_default(rstate);
        gmp_randseed_ui(rstate, seed);
    }
    mpz_srcptr prime() const {
        return abfield_characteristic_srcptr(ab);
    }

    matpoly_checker_base(matpoly_checker_base const & o)
        : m(o.m)
        , n(o.n)
        , len1(o.len1)
        , len2(o.len2)
        , dummy(SIZE_MAX)
        , seed(o.seed)
    {
        cxx_mpz p;
        abfield_init(ab);
        abfield_specify(ab, MPFQ_PRIME_MPZ, o.prime());
        gmp_randinit_default(rstate);
        gmp_randseed_ui(rstate, seed);
    }
    ~matpoly_checker_base() {
        gmp_randclear(rstate);
        abfield_clear(ab);
    }
};

template<typename fft_type>
struct matpoly_checker_ft : public matpoly_checker_base {
    typename matpoly_ft<fft_type>::memory_guard dummy_ft;

    matpoly_checker_ft(matpoly_checker_base const & base)
        : matpoly_checker_base(base)
        , dummy_ft(SIZE_MAX)
    {}
    matpoly_checker_ft(cxx_mpz const & p, unsigned int m, unsigned int n, unsigned int len1, unsigned int len2, gmp_randstate_t rstate0)
        : matpoly_checker_base(p, m, n, len1, len2, rstate0)
        , dummy_ft(SIZE_MAX)
    {}
    private:
    static inline int max_threads() {
#ifdef HAVE_OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }
    public:
    double time_dft(unsigned int nparallel, unsigned int cap = 25) {
        ASSERT_ALWAYS(nparallel <= (unsigned int) max_threads());
        fft_type fti = fft_type::polynomial_mul_info(prime(), len1, len2, n);
        double tt = 0;
        unsigned int n = 1;
        unsigned int v = 2;
        // 64M is bigger than any L3 cache
        size_t z;
        char buf[20];
        for( ; (z = nparallel * v * fti.size0_bytes()) < (1UL << 26) ; v <<= 1) ;
        printf("ft area (%u * %u), total size %s\n", nparallel, v, size_disp(z, buf));
            fflush(stdout);

        matpoly a(ab, nparallel, v, len1);
        a.fill_random(len1, rstate);
        matpoly_ft<fft_type> ta(a.m, a.n, fti);
        unsigned int ci = 0;
        for(unsigned int k = 0 ; k < cap ; k++) {
            tt = 0;
            n = 1 << k;
            tt = -wct_seconds();
            for(unsigned int i = 0 ; i < n ; i++, ci++) {
                submatrix_range R(0, ci % a.n,a.m,1);
                matpoly_ft<fft_type>::dft(ta.view(R), a.view(R));
            }
            tt = (tt + wct_seconds());
            printf("# dft time for %u threads: %.3g / %u = %.3g\n",
                    nparallel, tt, n, tt / n);
            fflush(stdout);
        }
        return tt / n;
    }

    void doit_mul(size_t L) {
        typedef lingen_substep_characteristics<op_mul<fft_type>> X_t;
        size_t LE = L / 2;
        size_t Lpi = m * LE / (m + n);
        X_t X(ab, rstate,
                0,      /* not used -- this field should probably go away
                           anyway */
                m+n, m+n, m+n,
                Lpi, Lpi, Lpi+Lpi-1);
        X.fill_tvec(&X_t::measure_dft_raw , "dft");
        X.fill_tvec(&X_t::measure_ift_raw , "ift");
        X.fill_tvec(&X_t::measure_conv_raw, "conv");
    }

    void doit_mp(size_t L) {
        typedef lingen_substep_characteristics<op_mp<fft_type>> X_t;
        size_t LE = L / 2;
        size_t Lpi = m * LE / (m + n);
        X_t X(ab, rstate,
                0,      /* not used -- this field should probably go away
                           anyway */
                m, m+n, m+n,
                LE + Lpi, Lpi, LE + 1);
        X.fill_tvec(&X_t::measure_dft_raw , "dft");
        X.fill_tvec(&X_t::measure_ift_raw , "ift");
        X.fill_tvec(&X_t::measure_conv_raw, "conv");
    }

    void doit(size_t L) {
        doit_mp(L);
        doit_mul(L);
    }

    /*
    int mp_and_mp_caching_are_consistent() {
        matpoly P(ab, m,   n, len1);
        matpoly Q(ab, n, n, len2);
        matpoly M0;
        matpoly M1;

        P.fill_random(len1, rstate);
        Q.fill_random(len2, rstate);

        M0.mp(P, Q);
        matpoly_ft<fft_type>::mp_caching(M1, P, Q, NULL);
        return M0.cmp(M1) == 0;
    }
    */
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
    param_list_decl_usage(pl, "len1", "length 1");
    param_list_decl_usage(pl, "len2", "length 2");
    param_list_decl_usage(pl, "seed", "random seed");
    param_list_decl_usage(pl, "thr", "threads (integer)");
}

int main(int argc, char * argv[])
{
    cxx_mpz p;
    gmp_randstate_t rstate;

    unsigned int m = 4;
    unsigned int n = 2;
    unsigned int len1 = 1000;
    unsigned int len2 = 600;
    unsigned int thr = 1;
    unsigned long seed = 0;

    cxx_param_list pl;

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
    param_list_parse_uint(pl, "len1", &len1);
    param_list_parse_uint(pl, "len2", &len2);
    param_list_parse_uint(pl, "thr", &thr);
    param_list_parse_ulong(pl, "seed", &seed);
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

    matpoly_checker_base checker(p, m, n, len1, len2, rstate);

    size_t L = len1;
#ifdef SELECT_MPFQ_LAYER_u64k1
    {
        matpoly_checker_ft<gf2x_fake_fft_info> checker_ft(checker);
        checker_ft.time_dft(thr);
        checker_ft.doit(L);
    }
    {
        matpoly_checker_ft<gf2x_cantor_fft_info> checker_ft(checker);
        checker_ft.time_dft(thr);
        checker_ft.doit(L);
    }
    {
        matpoly_checker_ft<gf2x_ternary_fft_info> checker_ft(checker);
        checker_ft.time_dft(thr);
        checker_ft.doit(L);
    }
#else
    {
        matpoly_checker_ft<fft_transform_info> checker_ft(checker);
        checker_ft.time_dft(thr);
        checker_ft.doit(L);
    }
#endif

    gmp_randclear(rstate);
}
