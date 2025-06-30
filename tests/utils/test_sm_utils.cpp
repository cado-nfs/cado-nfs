#include "cado.h" // IWYU pragma: keep

#include <cinttypes>
#include <cstdint> // for uint64_t, int64_t
#include <cstdio>
#include <cstdlib> // for EXIT_FAILURE, EXIT_SUCCESS

#include <istream>
#include <fstream>
#include <vector>

#include <gmp.h>
#include "fmt/base.h"

#include "cxx_mpz.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "sm_utils.hpp"
#include "utils_cxx.hpp"

#define FREQ 2 // when possible one time out of FREQ we try sm_single_rel

/* Return number of errors */
static int test_sm(std::string const & datafile)
{
    std::ifstream is(datafile);

    int err = 0;
    do {
        int degF, degN, degD, nb_ab, nbSM;
        unsigned int nb_test_single_rel = 0;
        cxx_mpz_poly F, N, Nc, D, Dc, SM, SMc;
        cxx_mpz tmp, ell;
        int64_t a, e[MAX_LEN_RELSET];
        uint64_t b, len_relset, r[MAX_LEN_RELSET];
        std::vector<pair_and_sides> ab_polys;
        sm_relset relset;
        is >> expect("in") >> degF;
        if (is.eof())
            break;

        for (int i = 0; i <= degF; i++) {
            is >> tmp;
            mpz_poly_setcoeff(F, i, tmp);
        }
        is >> ell;

        sm_side_info sm_info(F, ell, false);

        is >> nb_ab;

        for (int i = 0; i < nb_ab; i++) {
            is >> a >> b;
            ab_polys.emplace_back(a, b, 0, 1);
        }

        is >> len_relset;
        ASSERT_ALWAYS(len_relset <= MAX_LEN_RELSET);

        for (uint64_t i = 0; i < len_relset; i++) {
            is >> r[i] >> e[i];
        }

        /* make sure that all input was parsed correctly */
        is >> std::ws;
        ASSERT_ALWAYS(is);

        is >> expect("out") >> degN;
        ASSERT_ALWAYS(degN >= -1);
#ifdef __COVERITY__
        __coverity_mark_pointee_as_sanitized(&degN, GENERIC);
#endif

        for (int i = 0; i <= degN; i++) {
            is >> tmp;
            mpz_poly_setcoeff(N, i, tmp);
        }

        is >> degD;
        ASSERT_ALWAYS(0 <= degD);

        for (int i = 0; i <= degD; i++) {
            is >> tmp;
            mpz_poly_setcoeff(D, i, tmp);
        }

        is >> nbSM;
        ASSERT_ALWAYS(0 <= nbSM && nbSM <= degF);

        for (int i = 0; i < nbSM; i++) {
            is >> tmp;
            mpz_poly_setcoeff(SM, i, tmp);
        }

        is >> std::ws;

        //     /* Real tests begin here */
        //     /* artificially duplicate data, to test both sides */
        //     mpz_poly_ptr FF[2];
        //     FF[0] = &F[0]; FF[1] = &F[0];
        //     cxx_mpz_poly SMc2;
        //     mpz_poly_ptr SSMc[2];
        //     SSMc[0] = &SMc[0]; SSMc[1] = &SMc2[0];
        if (len_relset == 1 && e[0] == 1 && nb_test_single_rel % FREQ == 0) {
            nb_test_single_rel++;
            sm_info.compute_piecewise(SMc, ab_polys[r[0]].ab);
        } else {
            std::vector<mpz_poly_srcptr> const FF {F, F};
            relset = sm_build_one_relset(r, e, len_relset, ab_polys, FF,
                                         sm_info.ell2);
            Nc = relset.num[0];
            Dc = relset.denom[0];
            mpz_poly_reduce_frac_mod_f_mod_mpz(relset.num[0], relset.denom[0],
                                               F, sm_info.ell2);
            sm_info.compute_piecewise(SMc, relset.num[0]);
        }
        // mpz_poly_clear(SMc2);

        /* In case of error, print all relevant information */
        if (SM != SMc) {
            err++;
            fmt::print(stderr,
                    "### ERROR: computation of SM is wrong with:\nF = ");
            mpz_poly_fprintf(stderr, F);
            fmt::print(stderr, "ell = {}\nell2 = {}\n\n", ell, sm_info.ell2);
            sm_info.print(stderr);
            fprintf(stderr, "# Relation-set is:\n%" PRIu64 "", len_relset);
            for (uint64_t i = 0; i < len_relset; i++)
                fmt::print(stderr, " {}:{}", r[i], e[i]);
            fprintf(stderr, "\n# (a,b) pairs are:\n");
            for (int i = 0; i < nb_ab; i++) {
                cxx_mpz tmp;
                mpz_neg(tmp, mpz_poly_coeff_const(ab_polys[i].ab, 1));
                fmt::print(stderr, "{} {},{}\n",
                        i, cxx_mpz(mpz_poly_coeff_const(ab_polys[i].ab, 0)),
                        tmp);
            }
            if (N != Nc) {
                fmt::print(stderr,
                        "# Expected numerator in fraction corresponding to "
                        "the relation-set:\n");
                fmt::print(stderr, "{}\n", N);
                fprintf(stderr, "# Instead computed numerator is:\n");
                fmt::print(stderr, "{}\n", Nc);
            }
            if (D != Dc) {
                fmt::print(stderr,
                        "# Expected denominator in fraction corresponding to "
                        "the relation-set:\n");
                fmt::print(stderr, "{}\n", D);
                fmt::print(stderr, "# Instead computed denominator is:\n");
                fmt::print(stderr, "{}\n", Dc);
            }
            fprintf(stderr, "# Values of SM should be:\n");
            for (int i = 0; i < nbSM; i++) {
                fmt::print(stderr, "{} ", cxx_mpz(mpz_poly_coeff_const(SM, i)));
            }
            fprintf(stderr, "\n# but computed values of SM are:\n");
            for (int i = 0; i < nbSM; i++) {
                fmt::print(stderr, "{} ", cxx_mpz(mpz_poly_coeff_const(SMc, i)));
            }
            fprintf(stderr, "\n#######################\n");
        }
    } while (true);

    return err;
}

int main(int argc, char const * argv[])
{
    ASSERT_ALWAYS(argc == 2);
    char const * datafilename = argv[1];

    int const err = test_sm(datafilename);

    if (err)
        fprintf(stderr, "# %d erro%s found\n", err, (err == 1) ? "r" : "rs");
    return (err) ? EXIT_FAILURE : EXIT_SUCCESS;
}
