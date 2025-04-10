#include "cado.h" // IWYU pragma: keep

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <vector>

#include <gmp.h>

#include "arith/modredc_15ul.h"
#include "arith/modredc_2ul2.h"
#include "arith/modredc_ul.h"
#include "cxx_mpz.hpp"
#include "decomp.hpp"
#include "facul.hpp"
#include "facul_ecm.h"
#include "facul_method.hpp"
#include "facul_strategies.hpp"
#include "fm.hpp"
#include "generate_factoring_method.hpp"
#include "generate_strategies.hpp"
#include "macros.h"
#include "random_distributions.hpp"
#include "strategy.hpp"
#include "tab_decomp.hpp"
#include "tab_fm.hpp"
#include "tab_strategy.hpp"
#include "timing.h"

static double EPSILON_DBL = LDBL_EPSILON;

// convert type: from strategy_t to facul_strategy_oneside.
facul_strategy_oneside convert_strategy_to_facul_strategy(strategy_t const * t)
{
    tabular_fm_t * tab_fm = strategy_get_tab_fm(t);
    unsigned int nb_methods = tab_fm->size;

    facul_strategy_oneside strategy;

    strategy.lpb = UINT_MAX;
    strategy.BB = 0;
    strategy.BBB = 0;

    for (unsigned int i = 0; i < nb_methods; i++) {
        fm_t const * fm = tab_fm->tab[i];
        auto const method = facul_method_code(fm->method[0]);
        auto const curve = ec_parameterization_t(fm->method[1]);
        unsigned long const B1 = fm->method[2];
        unsigned long const B2 = fm->method[3];

        int const verbose = 0;
        if (method == PM1_METHOD || method == PP1_27_METHOD ||
            method == PP1_65_METHOD) {
            strategy.methods.emplace_back(
                facul_method::parameters(method, B1, B2), verbose);
        } else if (method == EC_METHOD) {
            /* XXX see with generate_fm which uses BOUND_SIGMA */
            unsigned long const sigma = curve == MONTY16 ? 1 : (2 + rand());
            int const extra_primes = 0;

            strategy.methods.emplace_back(
                facul_method::parameters(method, B1, B2, curve, sigma,
                                         extra_primes),
                verbose);
        } else {
            exit(EXIT_FAILURE);
        }
        ASSERT_ALWAYS(strategy.methods.back().plan != nullptr);
    }

    return strategy;
}

/* Given an list of decompositions, together with their population count,
 * pick one at random.
 */
static decomp pick_decomp_uniform(tabular_decomp const & t,
                                  gmp_randstate_ptr rstate)
{
    double total = 0;
    for (auto const & d: t)
        total += d.nb_elem;
    double r = total * random_uniform(rstate);
    for (auto const & d: t) {
        if (r <= d.nb_elem)
            return d;
        r -= d.nb_elem;
    }
    ASSERT_ALWAYS(0);
}

/*
  This function compute directly the probabity and the time to find a
  factor in a good decomposition in a cofactor of length r.

  duplicated in sieve/strategies/test_strategy_cado.cpp ?? XXX
 */
weighted_success bench_proba_time_st(gmp_randstate_t state,
                                     facul_strategy_oneside const & strategy,
                                     tabular_decomp const & init_tab, int r,
                                     int lpb)
{
    unsigned int const nb_success_max = 10000;
    unsigned int const nb_test_max = 10 * nb_success_max; // 20 for five percent

    unsigned nb_success = 0;
    double time = 0;

    for (unsigned int nb_test = 0; nb_test < nb_test_max; nb_test++) {
        auto const d = pick_decomp_uniform(init_tab, state);

        /* It's a bit weird. We're special-casing d[0] here, even
         * though it's just one factor among many. It seems that we
         * should rather form a number that is exactly along the
         * given distribution, right? XXX
         */
        cxx_mpz const N = generate_composite_integer(state, d[0], r);

        auto const starttime = microseconds();
        facul_result const f = facul(N, strategy);
        auto const endtime = microseconds();

        /* time += (end_test.tv_sec - st_test.tv_sec)*1000000L */
        /* 	+ end_test.tv_usec - st_test.tv_usec; */
        time += double(endtime - starttime);

        /* This seems really, really bogus! The test
         * (len_factor <= lpb && (r-len_factor) <= lpb)
         * can't be correct: we have decompositions that we're
         * interested in, with three factors, and our restriction to
         * d[0] + the rest above completely clashes with the
         * simplification that is done here.
         */

        // We should probably only be doing:
        // nb_success += f.status == FACUL_SMOOTH;
        //
        // And instead, we do:
        if (f.primes.size() != 2)
            continue;
        if (mpz_sizeinbase(f.primes[0], 2) > (size_t) lpb)
            continue;
        if (mpz_sizeinbase(f.primes[1], 2) > (size_t) lpb)
            continue;
        nb_success++;
        if (nb_success >= nb_success_max)
            return { nb_success, time, nb_test };
    }

    return { nb_success, time, nb_test_max };
}

int get_nb_word(int r)
{
    /*
     * We add 0.5 to the length of one word, because for our times the
     * length is inclusive. For example, if MODREDCUL_MAXBITS = 64
     * bits, a cofactor is in one word if is length is less OR equal to
     * 64 bits. So, if you don't add 0.5 to MODREDCUL_MAXBITS, you
     * miss the equality case and thus insert an error in your math.
     */
    double half_word = (MODREDCUL_MAXBITS + 0.5) / 2.0;
    int number_half_wd = floor(r / half_word);
    int ind = (number_half_wd < 2) ? 0 : number_half_wd - 1;
    return ind;
}

/*
 * This function is like bench_time() but only compute the time for one
 * length. In fact, i remove the unnecessary computation to reduce the
 * cost of this binary.
 */
void bench_time_mini(gmp_randstate_t state, tabular_fm_t * fm, int r)
{
    unsigned int len = fm->size; // number of methods!
    // precompute 4 arrays for our bench_time!
    //{{
    int len_n;
    if (r <= MODREDCUL_MAXBITS)
        len_n = MODREDCUL_MAXBITS;
    else if (r <= MODREDC15UL_MAXBITS)
        len_n = MODREDC15UL_MAXBITS;
    else if (r <= MODREDC2UL2_MAXBITS)
        len_n = MODREDC2UL2_MAXBITS;
    else
        len_n = MODREDC2UL2_MAXBITS + 30;
    int nb_test = 100000;
    std::vector<cxx_mpz> N;
    N.reserve(nb_test);
    for (int i = 0; i < nb_test; i++) {
        N.push_back(generate_prime_factor(state, len_n));
    }
    //}}
    int ind = get_nb_word(r);
    for (unsigned int i = 0; i < len; i++) {
        fm_t * elem = tabular_fm_get_fm_rw(fm, i);
        unsigned long const * param = fm_get_method(elem);
        auto const method = facul_method_code(param[0]);
        auto const curve = ec_parameterization_t(param[1]);
        unsigned long const B1 = param[2];
        unsigned long const B2 = param[3];
        if (B1 != 0 || B2 != 0) {
            facul_strategy_oneside st = generate_fm(method, B1, B2, curve);

            double res[4];
            ASSERT_ALWAYS(ind < 4);
            res[ind] = bench_time_fm_onelength(st, N, nb_test);
            fm_set_time(elem, res, 4);
        } else {
            double time[4] = {0, 0, 0, 0};
            fm_set_time(elem, time, 4);
        }
    }
}
/*
  This function is like bench_proba() but only compute the necessary
  probabilities. In fact, i remove the unnecessary computation to
  reduce the cost of this binary.
*/

static void bench_proba_mini(gmp_randstate_t state, tabular_fm_t * fm,
                             unsigned int const * val_p, unsigned int len_val_p,
                             unsigned int len_p_min)
{
    unsigned int len = fm->size; // number of methods!
    unsigned int p_max = 100;
    double * proba = (double *)calloc(p_max, sizeof(double));
    ASSERT(proba != NULL);

    //{{Will contain the precomputes of composite integer!
    std::vector<std::vector<cxx_mpz>> N(len_val_p);
    unsigned int nb_test_max = 10000;

    //}}

    for (unsigned int i = 0; i < len; i++) {
        fm_t * elem = tabular_fm_get_fm_rw(fm, i);
        unsigned long const * param = fm_get_method(elem);
        auto const method = facul_method_code(param[0]);
        auto const curve = ec_parameterization_t(param[1]);
        unsigned long const B1 = param[2];
        unsigned long const B2 = param[3];

        facul_strategy_oneside st = generate_fm(method, B1, B2, curve);

        unsigned int max_index = 0;
        for (unsigned int j = 0; j < len_val_p; j++) {
            unsigned int ind_proba = val_p[j] - len_p_min;
            if (B1 == 0 && B2 == 0) {
                proba[ind_proba] = 0;
            } else {
                unsigned int len_p = len_p_min + ind_proba;
                unsigned int len_n = 60 + len_p;
                proba[ind_proba] =
                    bench_proba_fm(st, state, len_p, len_n, N[j], nb_test_max);
            }
            if (ind_proba > max_index)
                max_index = ind_proba;
        }
        fm_set_proba(elem, proba, max_index + 1, len_p_min);
    }
    free(proba);
}

/************************************************************************/
/*     MAIN                                                             */
/************************************************************************/

// coverity[root_function]
int main()
{
    gmp_randstate_t state;
    gmp_randinit_default(state);

    unsigned int const fbb = 15; // 25;
    unsigned int const lpb = 19; // 29;
    unsigned int const r = 35;   // 55;

    // input decomp!
    tabular_decomp init_tab;
    init_tab.emplace_back(decomp(100, {17U, 18U}));
    init_tab.emplace_back(decomp(200, {15U, 20U}));

    unsigned int const len_val_factor = 4;
    unsigned int const val_factor[4] = {15, 17, 18, 20};

    // fm
    fm_t * pm1 = fm_create();
    unsigned long elem1[4] = {PM1_METHOD, 0, 50, 500};
    fm_set_method(pm1, elem1, 4);
    fm_t * pp1 = fm_create();
    unsigned long elem2[4] = {PP1_65_METHOD, 0, 70, 700};
    fm_set_method(pp1, elem2, 4);
    fm_t * ecm = fm_create();
    unsigned long elem3[4] = {EC_METHOD, BRENT12, 80, 1000};
    fm_set_method(ecm, elem3, 4);
    tabular_fm_t * tab = tabular_fm_create();
    tabular_fm_add_fm(tab, pm1);
    tabular_fm_add_fm(tab, pp1);
    tabular_fm_add_fm(tab, ecm);
    // bench our method
    bench_proba_mini(state, tab, val_factor, len_val_factor, fbb);
    // bench_time_mini (state, tab, r);
    /* printf ("mini\n"); */
    /* tabular_fm_print (tab); */

    /* bench_proba (state, tab, fbb); */
    /* bench_time(state, tab); */
    /* printf ("\n all \n"); */
    /* tabular_fm_print (tab); */
    // generate some examples of strategies!
    strategy_t * strat1 = strategy_create();
    strategy_add_fm(strat1, tab->tab[0]); // pm1
    strategy_add_fm(strat1, tab->tab[1]); // pp1
    strategy_add_fm(strat1, tab->tab[2]); // ecm
    const double prob1 = compute_proba_strategy(init_tab, strat1, fbb, lpb);
    // double time1 = compute_time_strategy(init_tab, strat1, r);

    // Create our strategy to use facul().
    const auto st = convert_strategy_to_facul_strategy(strat1);
    // bench our strategies!
    const auto res = bench_proba_time_st(state, st, init_tab, r, lpb);
    const double prob2 = res.prob;
    // double time2 = res[1];

    printf("prob1 = %lf, prob2 = %lf\n", prob1, prob2);
    // printf ("time1 = %lf, time2 = %lf\n", time1, time2);

    const double precision_p = 0.05;
    if ((prob1 - prob2) > precision_p || (prob2 - prob1) > precision_p) {
        fprintf(stderr, "error with the test(1)\n");
        return EXIT_FAILURE;
    }

    //{{test the function: bench_proba_time_pset()
    const auto c = int(tab->tab[0]->method[3] / tab->tab[0]->method[2]);
    const int param[6] = {
        (int)tab->tab[0]->method[2], (int)tab->tab[0]->method[2], 1, c, c, 1};

    tabular_fm_t * tmp = bench_proba_time_pset(
        facul_method_code(tab->tab[0]->method[0]),
        ec_parameterization_t(tab->tab[0]->method[1]),
        state, 17, 20, 18 * 3, param);

    /*check this probability: the probability to find a prime number
    of length 17 or 18 bits must be more than this to find 18 bits and
    less than this for 18 bits.*/
    fm_t * pm1_bis = tmp->tab[1];
    /* fm_print (tab->tab[0]); */
    /* fm_print (tmp->tab[1]); */
    if (pm1_bis->proba[0] < tab->tab[0]->proba[20 - fbb] ||
        pm1_bis->proba[0] > tab->tab[0]->proba[17 - fbb]) {
        fprintf(stderr, "error with the test(2)\n");
        return EXIT_FAILURE;
    }
    tabular_fm_free(tmp);

    //}}

    // test the generation of our strategy for one pair of cofactor!!!
    //{{
    strategy_set_proba(strat1, 0.4);
    strategy_set_time(strat1, 40);
    tabular_strategy_t * strat_r0 = tabular_strategy_create();
    tabular_strategy_t * strat_r1 = tabular_strategy_create();

    // firstly, create the zero strategy!
    strategy_t * zero_st = strategy_create();
    unsigned long tab0[4] = {PM1_METHOD, 0, 0, 0};
    fm_t * zero_fm = fm_create();
    fm_set_method(zero_fm, tab0, 4);
    strategy_add_fm(zero_st, zero_fm);
    strategy_set_proba(zero_st, 0.0);

    tabular_strategy_add_strategy(strat_r0, zero_st);
    tabular_strategy_add_strategy(strat_r0, strat1);
    tabular_strategy_add_strategy(strat_r1, zero_st);
    tabular_strategy_t * res2 = generate_strategy_r0_r1(strat_r0, strat_r1);

    // check proba + time!
    if (!(res2->size == 1 && res2->tab[0]->proba < EPSILON_DBL &&
          res2->tab[0]->time < EPSILON_DBL)) {
        fprintf(stderr, "error with the test(3)\n");
        return EXIT_FAILURE;
    }
    strategy_set_proba(strat_r1->tab[0], 1.0);
    tabular_strategy_t * res3 = generate_strategy_r0_r1(strat_r0, strat_r1);

    // check proba + time!
    if (!(res3->size == 2 && res3->tab[0]->proba < EPSILON_DBL &&
          res3->tab[0]->time < EPSILON_DBL &&
          (res3->tab[1]->proba - strat1->proba) < EPSILON_DBL &&
          (res3->tab[1]->time - strat1->time) < EPSILON_DBL)) {
        fprintf(stderr, "error with the test(3)\n");
        return EXIT_FAILURE;
    }

    //}}
    // free
    strategy_free(strat1);
    tabular_fm_free(tab);
    fm_free(pm1);
    fm_free(pp1);
    fm_free(ecm);
    gmp_randclear(state);

    fm_free(zero_fm);
    strategy_free(zero_st);
    tabular_strategy_free(res2);
    tabular_strategy_free(res3);
    tabular_strategy_free(strat_r0);
    tabular_strategy_free(strat_r1);

    return EXIT_SUCCESS;
}
