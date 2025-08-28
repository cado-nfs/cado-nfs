#include "cado.h" // IWYU pragma: keep

#include <ctime>
#include <ostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>
#include <stdexcept>

#include <gmp.h>

#include "fmt/ranges.h" // used to print std::vector<> // IWYU pragma: keep

#include "getprime.h"
#include "imaginary_quadratic_class_groups.hpp"
#include "params.h"     // param_list
#include "roots_mod.h"
#include "cado_poly.h"
#include "verbose.h"
#include "prime_power_factorization.hpp"
#include "cxx_mpz.hpp"
#include "mpz_poly.h"

#include <functional>

/* boost::hash_combine
 * Copyright 2005-2014 Daniel James.
 * Distributed under the Boost Software License, Version 1.0. (See accompanying
 * file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
template <class T>
static inline void hash_combine(std::size_t & seed, const T & v)
{
  const std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template<>
struct std::hash<cxx_mpz>
{
    size_t operator()(cxx_mpz const & v) const
    {
        const size_t n = mpz_size(v);
        size_t seed = n;
        hash_combine(seed, mpz_sgn(v));
        for (size_t i = 0; i < n; ++i) {
            hash_combine(seed, mpz_getlimbn(v, i));
        }
        return seed;
    }
};

template<>
struct std::hash<imaginary_quadratic_form>
{
    size_t operator()(imaginary_quadratic_form const & f) const
    {
        const std::hash<cxx_mpz> hasher;
        size_t seed = hasher(f.a);
        hash_combine(seed, f.b);
        hash_combine(seed, f.c);
        return seed;
    }
};

/* it's almost an interface that we'd like to port to
 * cado::prime_power_factorization, but the tricks about zero-exponents
 * are specific to this code here.
 */
class cxx_mpz_factored
{
    public:

    using valuations = std::vector<unsigned int>;

    private:

    std::vector<cxx_mpz> P;
    cxx_mpz v;
    valuations vals;

    public:

    cxx_mpz_factored(std::vector<cxx_mpz> const & primes,
                     valuations const & vals)
        : P(primes), v(1U), vals(vals)
    {
        for (size_t i = 0; i < P.size(); ++i) {
            auto e = vals[i];
            cxx_mpz pe;
            mpz_pow_ui(pe, primes[i], e);
            mpz_mul(v, v, pe);
        }
    }


    explicit cxx_mpz_factored(cado::prime_power_factorization const & F)
        : v(1U)
    {
        for(auto const & [ p, e ] : F) {
            cxx_mpz pe;
            P.emplace_back(p);
            vals.emplace_back(e);
            mpz_pow_ui(pe, p, e);
            mpz_mul(v, v, pe);
        }
    }

    /* This one looks very special */
    explicit cxx_mpz_factored(std::vector<cxx_mpz> const & primes)
        : P(primes), v(1), vals(primes.size(), 0)
    {
    }


    bool operator==(cxx_mpz_factored const & o) const
    {
        return v == o.v;
    }

    void lcm(valuations const & other_vals)
    {
        ASSERT_ALWAYS(vals.size() == other_vals.size());
        for (size_t i = 0; i < vals.size(); ++i) {
            for ( ; vals[i] < other_vals[i]; ++vals[i]) {
                mpz_mul(v, v, P[i]);
            }
        }
    }

    cxx_mpz const & value() const
    {
        return v;
    }

    std::vector<cxx_mpz> const & primes() const
    {
        return P;
    }

    cxx_mpz const & prime(size_t i) const
    {
        return P[i];
    }

    unsigned int valuation(size_t i) const
    {
        return vals[i];
    }

    unsigned int valuation_in(cxx_mpz const & p) const
    {
        for (size_t i = 0; i < nfactors(); ++i) {
            if (P[i] == p) {
                return vals[i];
            }
        }
        return 0U;
    }

    size_t nfactors() const
    {
        return P.size();
    }

    friend std::ostream & operator<<(std::ostream & o,
                                     cxx_mpz_factored const & n)
    {
        for (size_t i = 0; i < n.nfactors(); ++i) {
            o << (i ? " * " : "") << n.P[i] << "^" << n.vals[i];
        }
        return o << " = " << n.v;
    }
};

namespace fmt {
    template <> struct formatter<cxx_mpz_factored>: ostream_formatter {};
}


class imaginary_quadratic_cl_structure
{
    imaginary_quadratic_class_group cl;
    cxx_mpz_factored group_order;
    mutable cxx_gmp_randstate randgen;
    unsigned int const Dmod8;
    unsigned int const Dmod2;
    uint64_t B; /* Bound on norm of prime forms that generate the whole group */

    static cxx_mpz disc_from_poly(cxx_mpz_poly const & pol)
    {
        cxx_mpz D;
        mpz_neg(D, mpz_poly_coeff_const(pol, 0));
        if (mpz_fdiv_ui(D, 4) == 2 || mpz_fdiv_ui(D, 4) == 3) {
            mpz_mul_2exp(D, D, 2);
        }
        return D;
    }

    /* Return the smallest b in [0,p] such that (p,b,...) is a form or ULONG_MAX
     * if it does not exist.
     * Assumes p is prime.
     */
    unsigned long primeform_b_coeff(unsigned long p) const
    {
        if (p == 2U) {
            if (Dmod8 == 1U) {
                return 1U;
            } else if (Dmod8 == 0U) {
                return 0U;
            } else if (Dmod8 == 4U) {
                return 2U;
            } else { /* As D = 0 or 1 mod 4, here D = 5 mod 8 */
                return ULONG_MAX; /* no roots mod 8 */
            }
        } else {
            uint64_t rr[2];
            unsigned long Dmodp = mpz_fdiv_ui(cl.discriminant(), p);
            if (!Dmodp) { /* D is 0 mod p */
                return Dmod2 ? p : 0U;
            }
            unsigned int nr = roots_mod_uint64(rr, Dmodp, 2, p, randgen);
            /* nr is 0 or 2 as the case of 1 square root was already handled */
            if (nr == 0U) {
                return ULONG_MAX; /* no roots mod p */
            } else {
                /* rr[1] = p-rr[0], so their parity is different */
                return (rr[0] % 2U == Dmod2) ? rr[0] : rr[1];
            }
        }
    }

    bool compute_order_inner_rec(cxx_mpz_factored::valuations & v,
                                 unsigned int i,
                                 imaginary_quadratic_form f) const
    {
        cxx_mpz const & pi = group_order.prime(i);
        unsigned int ei = group_order.valuation(i);
        for (unsigned int j = 0; j <= ei; ++j) {
            if (f.is_one()) {
                v[i] = j;
                return true;
            }
            if (i+1 < group_order.nfactors()) {
                if (compute_order_inner_rec(v, i+1, f)) {
                    v[i] = j;
                    return true;
                }
            }
            f = f^pi;
        }
        return false;
    }

    void compute_order(cxx_mpz_factored::valuations & v,
                       imaginary_quadratic_form const & f) const
    {
        v.assign(group_order.nfactors(), 0U);
        if (!compute_order_inner_rec(v, 0, f)) {
            throw std::runtime_error("could not compute the order of element");
        }
    }

    struct Exponent
    {
        /* The exponent as a factored integer and a form whose order is the
         * exponent
         */
        cxx_mpz_factored exponent;
        imaginary_quadratic_form g;

        explicit Exponent(imaginary_quadratic_cl_structure const & G)
            : exponent(G.group_order.primes())
            , g(G.cl.one())
        {
        }

        void update(imaginary_quadratic_form const & f,
                    cxx_mpz_factored::valuations const & v)
        {
            cxx_mpz_factored::valuations exp_f, exp_g;
            for (size_t i = 0; i < exponent.nfactors(); ++i) {
                if (exponent.valuation(i) > v[i]) {
                    exp_g.push_back(0U);
                    exp_f.push_back(v[i]);
                } else {
                    exp_g.push_back(exponent.valuation(i));
                    exp_f.push_back(0U);
                }
            }
            exponent.lcm(v);
            auto r = g^(cxx_mpz_factored(exponent.primes(), exp_g).value());
            g = r * f^(cxx_mpz_factored(exponent.primes(), exp_f).value());
        }
    };

    struct pSylow
    {
        cxx_mpz pmpz;
        std::vector<std::pair<unsigned int, imaginary_quadratic_form>> groups;

        /* Assumes the p part of the group order fits in a unsigned int. */
        void compute_groups(std::unordered_set<imaginary_quadratic_form> const & G,
                            unsigned int exponent_pval,
                            imaginary_quadratic_form const & ge)
        {
            ASSERT_ALWAYS(mpz_fits_uint_p(pmpz));
            unsigned int p = mpz_get_ui(pmpz);
            std::vector<std::set<imaginary_quadratic_form>> by_order;
            by_order.resize(exponent_pval+1);
            std::vector<unsigned int> power_of_p {1U};
            for (unsigned int i = 0; i < exponent_pval; i++) {
                power_of_p.push_back(power_of_p.back() * p);
            }

            for (auto const & g: G) {
                /* compute order */
                auto f = g;
                unsigned int i = 0;
                for (; !f.is_one() && i < exponent_pval; ++i) {
                    f = f^pmpz;
                }
                ASSERT_ALWAYS(f.is_one());
                by_order[i].emplace(g);
            }
            ASSERT_ALWAYS(by_order[0].size() == 1);
            auto one = *by_order[0].begin();

            /* Compute the size of the different part of the p-Sylow.
             * The size of the ith group is p^groups[i].first.
             */
            unsigned int target_size = 1U;
            unsigned int curr_size = 1U;
            groups.clear();
            for (unsigned int i = 1; i <= exponent_pval; ++i) {
                target_size += by_order[i].size();
                for (unsigned int j = 0; curr_size < target_size;
                                                        curr_size *= p, ++j) {
                    if (i == 1) {
                        /* Set the generator to one for now, will be computed
                         * later.
                         */
                        groups.emplace_back(1U, one);
                    } else {
                        ASSERT_ALWAYS(groups[j].first+1 == i);
                        groups[j].first = i;
                    }
                }
            }

            /* Compute a generator for each group. */
            std::unordered_set<imaginary_quadratic_form> already_gen = {one};
            auto is_in_already_gen =
                [&already_gen] (imaginary_quadratic_form const & h) {
                    return already_gen.contains(h);
                };
            bool first = true;
            for (auto & [e, g]: groups) {
                if (first) {
                    /* To be coherent, for the first group (one of the largest),
                     * use for generator the power of the one stored in the
                     * corresponding Exponent struct
                     */
                    g = ge;
                } else {
                    auto it = find_if_not(by_order[e].begin(),
                                          by_order[e].end(), is_in_already_gen);
                    ASSERT_ALWAYS(it != by_order[e].end());
                    g = *it;
                }
                auto f = g;
                std::unordered_set<imaginary_quadratic_form> F;
                for (unsigned int i = 0; i < power_of_p[e]; ++i) {
                    for (auto const & h: already_gen) {
                        F.emplace(f*h);
                    }
                    f = f*g;
                }
                already_gen.merge(F);
                first = false;
            }
            ASSERT_ALWAYS(already_gen.size() == G.size());
        }

        friend std::ostream & operator<<(std::ostream & o, pSylow const & S)
        {
            for (size_t i = 0; i < S.groups.size(); ++i) {
                auto const & [sg_ord_pval, g] = S.groups[i];
                o << (i ? " x " : "")
                  << "Z/" << S.pmpz << "^" << sg_ord_pval << "Z";
            }
            for (size_t i = 0; i < S.groups.size(); ++i) {
                auto const & [sg_ord_pval, g] = S.groups[i];
                o << "\ng" << S.pmpz << "_" << i << " = " << g;
            }
            return o;
        }
    };

    public:

    imaginary_quadratic_cl_structure(cxx_cado_poly const & cpoly,
                                     cxx_mpz_factored group_order,
                                     unsigned int seed,
                                     uint64_t generators_bound = 0U)
        : cl(disc_from_poly(cpoly->pols[0]))
        , group_order(std::move(group_order))
        , Dmod8(mpz_fdiv_ui(cl.discriminant(), 8))
        , Dmod2(Dmod8 % 2)
        , B(generators_bound)
    {
        if (B == 0U) {
            B = cl.bound_generating_set_grh();
        }
        verbose_fmt_print(0, 2, "# Will use all prime forms up to {}\n", B);
        gmp_randseed_ui(randgen, seed);
    }

    imaginary_quadratic_form one() const
    {
        return cl.one();
    }

    Exponent exponent() const
    {
        verbose_fmt_print(0, 1, "\n# Computing exponent of {}\n", cl);
        verbose_fmt_print(0, 1, "# group_order = {}\n", group_order);

        Exponent E(*this);

        cxx_mpz_factored::valuations v;
        prime_info pi;
        prime_info_init(pi);

        for (unsigned long p = 2; p < B; p = getprime_mt(pi)) {
            if (E.exponent == group_order) {
                break; /* early abort if exponent is the group order */
            }

            unsigned long r = primeform_b_coeff(p);
            if (r == ULONG_MAX) {
                continue; /* skip this prime: no prime form exists for p */
            }
            verbose_fmt_print(0, 2, "# prime ideal (p, r) = ({}, {})\n", p, r);
            try {
                imaginary_quadratic_form fp(cl, p, r);
                verbose_fmt_print(0, 3, "# corresponding form: {}\n", fp);

                /* check if order of fp divides the exponent => skip it. */
                auto g = fp^E.exponent.value();
                if (g.is_one()) {
                    verbose_fmt_print(0, 2, "# skipping it: order divides the "
                                            "exponent\n");
                    continue;
                }

                compute_order(v, fp);
                verbose_fmt_print(0, 2, "# element has order {}\n", v);
                E.update(fp, v);
            } catch (imaginary_quadratic_form::not_primitive const & e) {
              /* We skipped primes that divided the conductor: cado-nfs.py
               * forbids them below the lpb so it should be safe to skip them.
               */
            }
        }

        prime_info_clear(pi);
        return E;
    }

     /* Assumes the p part of the group order fits in a unsigned int. */
    pSylow compute_p_sylow_naive(unsigned int p, Exponent const & E) const
    {
        verbose_fmt_print(0, 1, "# Computing {}-sylow of {}\n", p, cl);

        cxx_mpz pmpz (p);
        const unsigned int group_order_pval = group_order.valuation_in(p);
        const unsigned int exp_pval = E.exponent.valuation_in(p);
        unsigned long p_group_order = 1U;
        for (unsigned int i = 0; i < group_order_pval; ++i, p_group_order *= p);

        /* Compute cofac such that exponent = p^exp_pval * cofac */
        cxx_mpz cofac = E.exponent.value();
        mpz_remove(cofac, cofac, pmpz);
        verbose_fmt_print(0, 2, "# cofac = {}\n", cofac);

        /* G will contains all elements of the p-Sylow group (all elements whose
         * order is a power of p).
         */
        std::unordered_set<imaginary_quadratic_form> G;

        /* Use E.g (an element whose order is the exponent) to start to fill G:
         * E.g^cofact is an element of order p_group_order.
         */
        auto const ge = E.g^cofac;
        verbose_fmt_print(0, 3, "# using form {} of order {}^{}\n", ge, p,
                                                                    exp_pval);
        auto f = ge;
        for (unsigned long i = 0U; i < p_group_order; ++i) {
            G.emplace(f);
            f = f*ge;
        }
        verbose_fmt_print(0, 2, "# G contains {} elements\n", G.size());

        prime_info pi;
        prime_info_init(pi);

        for (unsigned long q = 2; q < B; q = getprime_mt(pi)) {
            if (G.size() == p_group_order) {
                break; /* early abort if we have all forms of order p^i */
            }

            unsigned long r = primeform_b_coeff(q);
            if (r == ULONG_MAX) {
                continue; /* skip this prime: no prime form exists for q */
            }
            verbose_fmt_print(0, 2, "# prime ideal (p, r) = ({}, {})\n", q, r);
            try {
                imaginary_quadratic_form fq(cl, q, r);
                verbose_fmt_print(0, 3, "# corresponding form: {}\n", fq);

                auto g = fq^cofac;

                std::unordered_set<imaginary_quadratic_form> H;
                unsigned long ord = 1U;
                auto f = g;
                for (; ord <= p_group_order; ++ord) {
                    if (f.is_one()) {
                        break;
                    } else if (!G.contains(f)) {
                        for (auto const & h: G) {
                            H.emplace(f*h);
                        }
                    }
                    f = f*g;
                }
                verbose_fmt_print(0, 2, "# order={}\n", ord);
                G.merge(H);
                verbose_fmt_print(0, 2, "# G contains {} elements\n", G.size());
            } catch (imaginary_quadratic_form::not_primitive const & e) {
              /* We skipped primes that divided the conductor: cado-nfs.py
               * forbids them below the lpb so it should be safe to skip them.
               */
            }
        }

        pSylow S(p, {});
        S.compute_groups(G, exp_pval, ge);

        prime_info_clear(pi);

        return S;
    }

    std::vector<pSylow> pSylow_groups(Exponent const & E,
                                      unsigned int naive_pSylow_bound) const;
};

namespace fmt {
    template <> struct formatter<imaginary_quadratic_cl_structure::pSylow>: ostream_formatter {};
}

std::vector<imaginary_quadratic_cl_structure::pSylow>
imaginary_quadratic_cl_structure::pSylow_groups(Exponent const & E,
                                                unsigned int naive_pSylow_bound)
                                                const
{
    std::vector<imaginary_quadratic_cl_structure::pSylow> S;
    for (size_t i = 0; i < group_order.nfactors(); ++i) {
        cxx_mpz const & pmpz = group_order.prime(i);
        const unsigned int group_order_pval = group_order.valuation(i);
        const unsigned int exp_pval = E.exponent.valuation(i);

        verbose_fmt_print(0, 2, "# group order has valuation {} in {}\n"
                                "# exponent has valuation {} in {}\n",
                                group_order_pval, pmpz, exp_pval, pmpz);

        /* Compute cofac such that exponent = p^exp_pval * cofac */
        cxx_mpz cofac = E.exponent.value();
        mpz_remove(cofac, cofac, pmpz);
        verbose_fmt_print(0, 2, "# cofac = {}\n", cofac);

        if (0 < exp_pval && exp_pval < group_order_pval) {
            cxx_mpz pe;
            mpz_pow_ui(pe, pmpz, group_order.valuation(i));

            if (pe < naive_pSylow_bound) {
                unsigned int p = mpz_get_ui(pmpz);
                auto psylow = compute_p_sylow_naive(p, E);
                verbose_fmt_print(0, 1, "{}-Sylow: {}\n", p, psylow);
                S.emplace_back(std::move(psylow));
            } else {
                verbose_fmt_print(0, 1, "{}-Sylow: too large, skipped\n",
                                        pmpz);
                continue;
            }
        } else { /* exp_pval == group_order or exp_pval = 0 */
            pSylow psylow {pmpz, {{exp_pval, E.g^cofac}}};
            verbose_fmt_print(0, 1, "{}-Sylow: {}\n", pmpz, psylow);
            S.emplace_back(std::move(psylow));
        }

        unsigned int actual_pval = 0;
        for (auto const & sg: S.back().groups) {
            actual_pval += sg.first;
        }
        if (actual_pval != group_order_pval) {
            verbose_fmt_print(0, 1, "# Warning: for p={}, from parameters "
                                    "expected valuation of p in the group "
                                    "order to be {} , got {} instead\n", pmpz,
                                    group_order_pval, actual_pval);
        }

    }
    return S;
}


static void
usage (cxx_param_list & pl, const char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}


struct command_line
{
    std::string polyfilename;
    cado::prime_power_factorization factored_order;

    uint64_t bound = 0;
    unsigned int seed;
    unsigned int pSylow_bound = 40000U; /* allow naive pSylow computations up
                                         * to 2^15, 3^10, 5^6, 7^5, 11^4, 13^4,
                                         * 17^3, ..., 31^3, 37^2, ..., 199^2
                                         */
    int verbosity_level = 1; /* each -v on command line increases it by 1 */

    static void declare_usage(cxx_param_list & pl) {
        param_list_usage_header(pl, "Compute the exponent of the group given "
                                    "the factorization of its order by "
                                    "computing the order of all elements "
                                    "corresponding to columns of the matrix\n");
        param_list_decl_usage(pl, "poly", "input polynomial file");
        param_list_decl_usage(pl, "order", "factorization of the group order");
        param_list_decl_usage(pl, "B", "bound to build the generating set "
                                       "(default is 6*log(|D|)^2)");
        param_list_decl_usage(pl, "seed", "seed for random generator");
        param_list_decl_usage(pl, "pSylow-bound", "Use naive pSylow computation"
                                                  " if the p-part of the group "
                                                  "order is below this bound "
                                                  "(default: 40000).");
        param_list_decl_usage(pl, "v", "enable verbose output");
        verbose_decl_usage(pl);
    }

    void configure_switches(cxx_param_list & pl)
    {
        param_list_configure_switch(pl, "-v", &verbosity_level);
    }

    void lookup_parameters(cxx_param_list & pl)
    {
        param_list_parse(pl, "poly", polyfilename);
        param_list_parse(pl, "order", factored_order);
        param_list_parse(pl, "pSylow-bound", pSylow_bound);
        param_list_parse(pl, "B", bound);
        if (!param_list_parse(pl, "seed", seed)) {
            seed = time(NULL);
        }
    }

    void check_inconsistencies(const char * argv0, cxx_param_list & pl) const
    {
        if (polyfilename.empty()) {
            fmt::print(stderr, "Error, missing -poly\n");
            usage(pl, argv0);
        }
        if (factored_order.empty()) {
            fmt::print(stderr, "Error, -primes must not be empty\n");
            usage(pl, argv0);
        }
        for(auto const & [ p, e ] : factored_order) {
            if (e <= 0) {
                fmt::print(stderr, "Error, exponents in the factored order must be positive\n");
                usage(pl, argv0);
            }
        }
    }
};

int main(int argc, char const * argv[])
{
    const char *argv0 = argv[0];

    cxx_param_list pl;
    command_line cmdline;

    command_line::declare_usage(pl);
    cmdline.configure_switches(pl);

    argv++, argc--;
    if (argc == 0)
      usage(pl, argv0);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv))
            continue;
        fmt::print(stderr, "Unhandled parameter {}\n", argv[0]);
        usage(pl, argv0);
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    cmdline.lookup_parameters(pl);
    if (param_list_warn_unused(pl)) {
        usage(pl, argv0);
    }
    cmdline.check_inconsistencies(argv0, pl);

    verbose_output_init(2);
    verbose_output_add(0, stdout, cmdline.verbosity_level);
    verbose_output_add(1, stderr, 1);

    cxx_cado_poly cpoly;
    if (!cado_poly_read (cpoly, cmdline.polyfilename.c_str()))
    {
        fmt::print(stderr, "Error reading polynomial file\n");
        exit(EXIT_FAILURE);
    }

    /* check that poly is x^2-d with d < 0 */
    if (!cado_poly_is_imaginary_quadratic(cpoly)) {
        fmt::print(stderr, "This binary is only working for 1 poly of the form "
                           "x^2-d, with negative d\n");
        exit(EXIT_FAILURE);
    }

    /* To help the parsing from Python script. */
    const std::string prefix("structure: ");

    /* Set the group order from the parameters */
    cxx_mpz_factored group_order(cmdline.factored_order);

    /* Set the class group */
    imaginary_quadratic_cl_structure cl(cpoly, group_order, cmdline.seed,
                                                            cmdline.bound);

    /* Compute the exponent */
    auto E = cl.exponent();
    verbose_fmt_print(0, 1, "{}exponent = {}\n{}gen = {}\n",
                            prefix, E.exponent, prefix, E.g);

    /* Compute the p-Sylow groups if possible */
    auto S = cl.pSylow_groups(E, cmdline.pSylow_bound);

    /* If possible print the structure (with generators) of the group */
    if (S.size() == E.exponent.nfactors()) {
        std::vector<imaginary_quadratic_form> gens { E.g };
        verbose_fmt_print(0, 1, "{}group = Z/{}Z",
                                prefix, E.exponent.value());
        cxx_mpz ord, pe;
        imaginary_quadratic_form gi = cl.one();
        size_t i = 1;
        do {
            ord = 1U;
            gi = cl.one();
            for (auto const & s: S) {
                if (i < s.groups.size()) {
                    mpz_pow_ui(pe, s.pmpz, s.groups[i].first);
                    mpz_mul(ord, ord, pe);
                    gi = gi*s.groups[i].second;
                }
            }
            if (ord > 1U) {
                verbose_fmt_print(0, 1, " x Z/{}Z", ord);
                gens.emplace_back(gi);
            }
            ++i;
        } while (ord > 1U);
        verbose_fmt_print(0, 1, "\n");
        for (size_t i = 0; i < gens.size(); ++i) {
            verbose_fmt_print(0, 1, "{}gen_{} = {}\n", prefix, i, gens[i]);
        }
    }

    verbose_output_clear();

    return S.size() == E.exponent.nfactors() ? EXIT_SUCCESS : EXIT_FAILURE;
}
