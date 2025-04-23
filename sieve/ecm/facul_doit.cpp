/* This file is _NOT_ a standalone compilation unit. It is included by
 * facul_doit_blah.cpp ; therefore it's _normal_ if this file does not
 * seem to be self-contained.
 */

/* We don't even need to #include "cado.h", because the
 * facul_doit_blah.cpp has to do so already. However, our policy
 * detection script will insist on seeing this include, and it doesn't do
 * any harm
 */
#include "cado.h" // IWYU pragma: keep

#ifndef FACUL_DOIT_READY_TO_INCLUDE_IMPL_CODE
#error "This file must not be used as a standalone compilation unit"
#endif

#include <cstdlib>
#include <memory>

#include "facul_method.hpp"
#include "facul_doit.hpp"
#include "modset.hpp"
#include "pm1.h"
#include "pp1.h"
#include "facul_ecm.h"
#include "facul.hpp"
#include "mpqs.h"
#include "cxx_mpz.hpp"
#include "arith/mod_ul.h"
#include "macros.h"
#include "utils_cxx.hpp"

#ifdef ENABLE_UNSAFE_FACUL_STATS
extern unsigned long stats_called[];
extern unsigned long stats_found_n[];
extern int stats_current_index;
#endif

#define mod_intget_mpz MOD_RENAME(int_get_mpz)
void
mod_intget_mpz(mpz_t z, const modint_t x) {
#ifdef MOD_SIZE
    mpz_import(z, MOD_SIZE, -1, sizeof(unsigned long), 0, 0, x);
#else
    mpz_set(z, x);
#endif
}

#define mod_intget_cxx_mpz MOD_RENAME(int_get_cxx_mpz)
cxx_mpz mod_intget_cxx_mpz(const modint_t x) {
    cxx_mpz c;
    mod_intget_mpz(c, x);
    return c;
}

static inline FaculModulusBase *
modset_init (modint_t m)
{
  return FaculModulusBase::MOD_APPEND_TYPE(init)(m);
}

/*****************************************************************************/
/*                       STRATEGY BOOK                                       */
/*****************************************************************************/


/* This new version of facul_doit_onefm applies one factoring method
 * 'method' on an integer 'm'.
 *
 * It appends the prime factors it finds to the [factors] argument
 * The [composites] argument is first cleared, and then all composites
 * that are found are stored there.
 *
 * Note that it is possible to find composites even if no prime factor is
 * found!
 *
 * This returns FACUL_NOT_SMOOTH, FACUL_MAYBE, or FACUL_SMOOTH
 */
facul_status facul_doit_onefm(std::vector<cxx_mpz> & factors,
        const modulus_t m,
        facul_method const & method,
        std::vector<std::unique_ptr<FaculModulusBase>> & composites,
        unsigned long lpb, double BB, double BBB)
{
    ASSERT_ALWAYS(composites.empty());

    modint_t n;
    modint_t f;
    mod_intinit(n);
    mod_intinit(f);
    auto dt_n = call_dtor([&](){mod_intclear(n);});
    auto dt_f = call_dtor([&](){mod_intclear(f);});
    int bt;

    mod_getmod_int (n, m);
    mod_intset_ul (f, 1UL);

    switch(method.method) {
        case PM1_METHOD:
            bt = pm1 (f, m, (pm1_plan_t *) (method.plan));
            break;
        case PP1_27_METHOD:
            bt = pp1_27 (f, m, (pp1_plan_t *) (method.plan));
            break;
        case PP1_65_METHOD:
            bt = pp1_65 (f, m, (pp1_plan_t *) (method.plan));
            break;
        case EC_METHOD:
            bt = ecm (f, m, (ecm_plan_t *) (method.plan));
            break;
        case MPQS_METHOD:
            bt = mpqs (f, m);
            break;
        case NO_METHOD:
            ASSERT_ALWAYS(0);
    }


    if (mod_intequal_ul (f, 1UL)) {
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
    } else if (mod_intequal (f, n)) {
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
    double f_dbl = mod_intget_double (f);
    bool fprime = f_dbl < BB;
    if (fprime && mod_intbits (f) > lpb)
        /* A prime > 2^lpb, not smooth */
        return FACUL_NOT_SMOOTH;
    else if (2 * lpb < mod_intbits (f) && f_dbl < BBB)
        /* if L^2 < f < B^3, it cannot be smooth */
        return FACUL_NOT_SMOOTH;

    /* Compute the cofactor */
    mod_intdivexact (n, n, f);

    /* Do the same tests, and see if the cofactor is something non smooth
     */
    double n_dbl = mod_intget_double (n);
    bool cfprime = n_dbl < BB;
    if (cfprime && mod_intbits (n) > lpb)
        return FACUL_NOT_SMOOTH;
    else if (2 * lpb < mod_intbits (n) && n_dbl < BBB)
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
        fm.reset(modset_init(f));
        fprime = fm->isprime ();
        if (fprime && mod_intbits (f) > lpb)
            return FACUL_NOT_SMOOTH;
    }

    std::unique_ptr<FaculModulusBase> cfm;

    if (!cfprime) {
        cfm.reset(modset_init(n));
        cfprime = cfm->isprime();
        if (cfprime && mod_intbits (n) > lpb)
            return FACUL_NOT_SMOOTH;
    }

    /* So each of factor and cofactor is either a prime < 2^lpb, 
       or is composite */

    if (fprime) {
        factors.emplace_back(mod_intget_cxx_mpz(f));
    } else {
        composites.emplace_back(std::move(fm));
    }

    if (cfprime) {
        factors.emplace_back(mod_intget_cxx_mpz(n));
    } else {
        composites.emplace_back(std::move(cfm));
    }

    return composites.empty() ? FACUL_SMOOTH : FACUL_MAYBE;
}
