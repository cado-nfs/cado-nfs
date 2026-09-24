#include "cado.h" // IWYU pragma: keep

#include <cstdint>

#include "facul_arithxx.hpp"
#include "pm1.hpp"
#include "pp1_stage2.hpp"

/* #define PARI */

#ifdef PARI
#include "fmt/base.h"
#endif

/* Do we want backtracking when processing factors of 2 in E? */
#ifndef PM1_BACKTRACKING
/* Default is "yes." Set to 0 for "no." */
#define PM1_BACKTRACKING 1
#endif

/* Looks for a factor of the modulus m, using the P-1 algorithm.
 * The parameters of P-1 are given in plan. The factor found (or 1) is
 * stored in f. Returns 1 if backtracking was used, 0 otherwise.
 */
template<typename layer>
static int pm1_run(typename layer::Integer & f,
        typename layer::Modulus const & m,
        pm1_plan_t const & plan)
{
    using Residue = typename layer::Residue;

    Residue x(m), t(m), X(m), one(m), two(m);
    int bt = 0;

    m.set1(one);
    m.add(two, one, one);

    /* Stage 1, a simple exponentiation ... */
    m.pow2(x, plan.E);
    /* ... except for the backtracking part for the 2's in the exponent */
    m.set(t, x);
    for (unsigned int i = 0; i < plan.exp2; i++) {
        m.sqr(x, x);
#if PM1_BACKTRACKING
        if (m.is1(x)) {
            m.set(x, t);
            bt = 1;
            break;
        }
        m.set(t, x);
#endif
    }

#ifdef PARI
    fmt::print("E = B1_exponent ({}); x = Mod(2, {})^E; x == {} /* PARI */\n",
            plan.B1, m.getmod(), m.get(x));
#endif

    m.sub(t, x, one);
    m.gcd(f, t);

    if (f > uint64_t(1) || plan.B1 >= plan.stage2.B2)
        return 0;

    /* Compute X = x + 1/x. The modulus is odd, so that x, a power of 2,
     * is invertible. */
    // coverity[check_return]
    m.inv(X, x);
    m.add(X, X, x);

#ifdef PARI
    fmt::print("X = x+1/x; X == {} /* PARI */\n", m.get(X));
#endif

    bt = pp1_stage2<layer>(t, X, plan.stage2, two, m);
    m.gcd(f, t);

    return bt;
}

int pm1(modintredcul_t f, const modulusredcul_t m, const pm1_plan_t * plan)
{
    return facul_arithxx_run(f, m, [&]<typename layer>(auto & g, auto const & mm) {
            return pm1_run<layer>(g, mm, *plan); });
}

int pm1(modintredc15ul_t f, const modulusredc15ul_t m, const pm1_plan_t * plan)
{
    return facul_arithxx_run(f, m, [&]<typename layer>(auto & g, auto const & mm) {
            return pm1_run<layer>(g, mm, *plan); });
}

int pm1(modintredc2ul2_t f, const modulusredc2ul2_t m, const pm1_plan_t * plan)
{
    return facul_arithxx_run(f, m, [&]<typename layer>(auto & g, auto const & mm) {
            return pm1_run<layer>(g, mm, *plan); });
}

int pm1(modintmpz_t f, const modulusmpz_t m, const pm1_plan_t * plan)
{
    return facul_arithxx_run(f, m, [&]<typename layer>(auto & g, auto const & mm) {
            return pm1_run<layer>(g, mm, *plan); });
}
