#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>

#include <memory>
#include <type_traits>
#include <vector>

#include <gmp.h>

#include "arithxx/mod_mpz_new.hpp"
#include "arithxx/modredc126.hpp"
#include "arithxx/modredc64.hpp"
#include "arithxx/modredc96.hpp"
#include "cxx_mpz.hpp"
#include "ecm.hpp"
#include "facul.hpp"
#include "facul_doit.hpp"
#include "facul_method.hpp"
#include "macros.h"
#include "modset.hpp"
#include "mpqs_doit.h"
#include "pm1.hpp"
#include "pp1.hpp"

#ifdef ENABLE_UNSAFE_FACUL_STATS
extern unsigned long stats_called[];
extern unsigned long stats_found_n[];
extern int stats_current_index;
#endif

/* {{{ integer operations, for the Integer types of the arithxx layers */
namespace {

double int_get_double(cxx_mpz const & x) { return mpz_get_d(x); }
template<typename I> double int_get_double(I const & x) { return double(x); }

size_t int_bits(cxx_mpz const & x) { return mpz_sizeinbase(x, 2); }
template<typename I> size_t int_bits(I const & x) { return x.bits(); }

cxx_mpz int_divexact(cxx_mpz const & n, cxx_mpz const & f)
{
    cxx_mpz q;
    mpz_divexact(q, n, f);
    return q;
}
template<typename I> I int_divexact(I const & n, I const & f) { return n.divexact(f); }

/* The factor x of a number over some layer, over the smallest layer that
 * fits. */
std::unique_ptr<FaculModulusBase> facul_modulus(cxx_mpz const & x)
{
    return std::unique_ptr<FaculModulusBase>(FaculModulusBase::init_mpz(x));
}
std::unique_ptr<FaculModulusBase> facul_modulus(Integer64 const & x)
{
    return std::make_unique<FaculModulus<arithxx_modredc64>>(x);
}
std::unique_ptr<FaculModulusBase> facul_modulus(Integer128 const & x)
{
    size_t const bits = x.bits();
    if (bits <= 64)
        return std::make_unique<FaculModulus<arithxx_modredc64>>(Integer64(x[0]));
    if (bits <= 96)
        return std::make_unique<FaculModulus<arithxx_modredc96>>(x);
    return std::make_unique<FaculModulus<arithxx_modredc126>>(x);
}

} /* namespace */
/* }}} */

template<typename layer>
facul_status facul_doit_onefm(std::vector<cxx_mpz> & factors,
        typename layer::Modulus const & m,
        facul_method const & method,
        std::vector<std::unique_ptr<FaculModulusBase>> & composites,
        unsigned long const lpb, double const BB, double const BBB)
{
    using Integer = typename layer::Integer;

    ASSERT_ALWAYS(composites.empty());

    Integer n = m.getmod();
    Integer f;
    int bt = 0;

    switch(method.method) {
        case PM1_METHOD:
            bt = pm1<layer>(f, m, *(pm1_plan_t const *) (method.plan));
            break;
        case PP1_27_METHOD:
            bt = pp1_27<layer>(f, m, *(pp1_plan_t const *) (method.plan));
            break;
        case PP1_65_METHOD:
            bt = pp1_65<layer>(f, m, *(pp1_plan_t const *) (method.plan));
            break;
        case EC_METHOD:
            bt = ecm<layer>(f, m, *(ecm_plan_t const *) (method.plan));
            break;
        case MPQS_METHOD:
            {
                cxx_mpz fz;
                cxx_mpz const nz(n);
                mpqs_doit(fz, nz, 0);
                f = Integer(fz);
            }
            break;
        case NO_METHOD:
            ASSERT_ALWAYS(0);
    }

    if (f == uint64_t(1)) {
        if (bt == 0)
            /* No factor found, no backtracking... this was a simple miss. */
            return FACUL_MAYBE;
        else
            /* Backtracking was used, so there was a point where all
             * factors had been found simultaneously, but backing up to
             * the previous checkpoint resulted in no factors being
             * found. We could try to do some more clever backtracking to
             * discover the factors yet. TODO. For now, just continue to
             * the next method. */
            return FACUL_MAYBE;
    } else if (f == n) {
#ifdef ENABLE_UNSAFE_FACUL_STATS
        if (stats_current_index < STATS_LEN)
            stats_found_n[stats_current_index]++;
#endif  /* ENABLE_UNSAFE_FACUL_STATS */
        if (bt == 0)
            /* Input number was found without any backtracking happening?
             * Find out when this can occur and how to get a chance of
             * finding the factors yet. TODO. */
            return FACUL_MAYBE;
        else
            /* see above */
            return FACUL_MAYBE;
    }

    /* So we found a non-trivial factor. See if it is prime, if the
       cofactor is prime, and if one of them is, whether they are too
       large for our smoothness bounds */

    /* A quick test if the factor is <= fbb^2 and >2^lpb */
    double const f_dbl = int_get_double(f);
    bool fprime = f_dbl < BB;
    if (fprime && int_bits(f) > lpb)
        /* A prime > 2^lpb, not smooth */
        return FACUL_NOT_SMOOTH;
    else if (2 * lpb < int_bits(f) && f_dbl < BBB)
        /* if L^2 < f < B^3, it cannot be smooth */
        return FACUL_NOT_SMOOTH;

    /* Compute the cofactor */
    n = int_divexact(n, f);

    /* Do the same tests, and see if the cofactor is something non smooth
     */
    double const n_dbl = int_get_double(n);
    bool cfprime = n_dbl < BB;
    if (cfprime && int_bits(n) > lpb)
        return FACUL_NOT_SMOOTH;
    else if (2 * lpb < int_bits(n) && n_dbl < BBB)
        return FACUL_NOT_SMOOTH;

    /* At this point, if fprime or cfprime are false, it means that the
     * might be composite, but we need to check that.
     *
     * We still have a chance to abort if we find out that one of them is
     * indeed prime, and out of range.
     *
     * Determine now for certain if the factor is prime
     */

    std::unique_ptr<FaculModulusBase> fm;
    if (!fprime) {
        fm = facul_modulus(f);
        fprime = fm->isprime();
        if (fprime && int_bits(f) > lpb)
            return FACUL_NOT_SMOOTH;
    }

    std::unique_ptr<FaculModulusBase> cfm;

    if (!cfprime) {
        cfm = facul_modulus(n);
        cfprime = cfm->isprime();
        if (cfprime && int_bits(n) > lpb)
            return FACUL_NOT_SMOOTH;
    }

    /* So each of factor and cofactor is either a prime < 2^lpb,
       or is composite */

    if (fprime) {
        factors.emplace_back(cxx_mpz(f));
    } else {
        composites.emplace_back(std::move(fm));
    }

    if (cfprime) {
        factors.emplace_back(cxx_mpz(n));
    } else {
        composites.emplace_back(std::move(cfm));
    }

    return composites.empty() ? FACUL_SMOOTH : FACUL_MAYBE;
}

facul_status facul_doit_onefm(std::vector<cxx_mpz> & factors,
        cxx_mpz const & n,
        facul_method const & method,
        std::vector<std::unique_ptr<FaculModulusBase>> & composites,
        unsigned long const lpb, double const BB, double const BBB)
{
    std::unique_ptr<FaculModulusBase> const fm(FaculModulusBase::init_mpz(n));
    return fm->facul_doit_onefm(factors, method, composites, lpb, BB, BBB);
}

#define FACUL_DOIT_INSTANTIATE(layer)                                   \
template facul_status facul_doit_onefm<layer>(std::vector<cxx_mpz> &,   \
        layer::Modulus const &, facul_method const &,                   \
        std::vector<std::unique_ptr<FaculModulusBase>> &,               \
        unsigned long, double, double);

FACUL_DOIT_INSTANTIATE(arithxx_modredc64)
FACUL_DOIT_INSTANTIATE(arithxx_modredc96)
FACUL_DOIT_INSTANTIATE(arithxx_modredc126)
FACUL_DOIT_INSTANTIATE(arithxx_mod_mpz_new)
