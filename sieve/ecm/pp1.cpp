#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <utility>

#include "bytecode.h"
#include "arithxx/modredc64.hpp"
#include "arithxx/modredc96.hpp"
#include "arithxx/modredc126.hpp"
#include "arithxx/mod_mpz_new.hpp"
#include "pp1.hpp"
#include "pp1_stage2.hpp"

/* Interpret the PRAC bytecode for P+1: X <- V_k(X), where k is the
 * product that the bytecode encodes. */
template<typename layer>
static void pp1_stage1(typename layer::Residue & X, bytecode_const bc,
        typename layer::Residue const & two,
        typename layer::Modulus const & m)
{
    using Residue = typename layer::Residue;

    /* we need 5 points: 3 for PRAC + 2 temporary points */
    Residue R[5] = { Residue(m), Residue(m), Residue(m), Residue(m), Residue(m) };

    /* current point (here starting point) go into R[0] at init */
    m.set(R[0], X);

    for (bool finished = false; !finished; bc++) {
        switch (*bc) {
            case PRAC_SWAP: /* [ = 's' ] Swap R[0], R[1] */
                std::swap(R[0], R[1]);
                break;
            case PRAC_SUBBLOCK_INIT: /* [ = 'i' ] Start of a sub-block */
                m.set(R[1], R[0]);
                m.set(R[2], R[0]);
                m.V_dbl(R[0], R[0], two);
                break;
            case PRAC_SUBBLOCK_FINAL: /* [ = 'f' ] End of a sub-block */
                m.V_dadd(R[0], R[0], R[1], R[2]);
                break;
            case PRAC_BLOCK_FINAL: /* [ = 'F' ] End of the block */
                m.V_dadd(R[1], R[0], R[1], R[2]);
                finished = true;
                break;
            case 1:
                m.V_dadd(R[3], R[0], R[1], R[2]);
                m.V_dadd(R[4], R[3], R[0], R[1]);
                m.V_dadd(R[1], R[1], R[3], R[0]);
                m.set(R[0], R[4]);
                break;
            case 2:
                m.V_dadd(R[1], R[0], R[1], R[2]);
                m.V_dbl(R[0], R[0], two);
                break;
            case 3:
                m.V_dadd(R[3], R[1], R[0], R[2]);
                m.set(R[2], R[1]);
                m.set(R[1], R[3]);
                break;
            case 4:
                m.V_dadd(R[1], R[1], R[0], R[2]);
                m.V_dbl(R[0], R[0], two);
                break;
            case 5:
                m.V_dadd(R[2], R[2], R[0], R[1]);
                m.V_dbl(R[0], R[0], two);
                break;
            case 6:
                m.V_dbl(R[3], R[0], two);
                m.V_dadd(R[4], R[0], R[1], R[2]);
                m.V_dadd(R[4], R[3], R[4], R[2]);
                m.set(R[2], R[4]);
                m.V_dadd(R[4], R[3], R[0], R[0]);
                m.set(R[0], R[4]);
                std::swap(R[1], R[2]);
                break;
            case 7:
                m.V_dadd(R[3], R[0], R[1], R[2]);
                m.V_dadd(R[4], R[3], R[0], R[1]);
                m.set(R[1], R[4]);
                m.V_dbl(R[3], R[0], two);
                m.set(R[4], R[0]);
                m.V_dadd(R[0], R[0], R[3], R[4]);
                break;
            case 8:
                m.V_dadd(R[3], R[0], R[1], R[2]);
                m.V_dadd(R[2], R[2], R[0], R[1]);
                std::swap(R[1], R[3]);
                m.V_dbl(R[3], R[0], two);
                m.V_dadd(R[4], R[0], R[3], R[0]);
                m.set(R[0], R[4]);
                break;
            case 9:
                m.V_dadd(R[2], R[2], R[1], R[0]);
                m.V_dbl(R[1], R[1], two);
                break;
            case 10:
                /* Combined final add of old subchain and init of new
                 * subchain [=fi] */
                m.V_dadd(R[1], R[0], R[1], R[2]);
                m.set(R[2], R[1]);
                m.V_dbl(R[0], R[1], two);
                break;
            case 11:
                /* Combined rule 3 and rule 0 [=\x3s] */
                m.set(R[3], R[0]);
                m.V_dadd(R[0], R[1], R[0], R[2]);
                m.set(R[2], R[1]);
                m.set(R[1], R[3]);
                break;
            case 12:
                /* Combined rule 3, then subchain end/start [=\x3fi] */
                m.V_dadd(R[3], R[1], R[0], R[2]);
                m.V_dadd(R[2], R[0], R[3], R[1]);
                m.set(R[1], R[2]);
                m.V_dbl(R[0], R[2], two);
                break;
            case 13:
                /* Combined rule 3, swap, rule 3 and swap, merged a bit
                 * [=\x3s\x3s] */
                m.set(R[3], R[1]);
                m.V_dadd(R[1], R[1], R[0], R[2]);
                m.set(R[2], R[0]);
                m.V_dadd(R[0], R[0], R[1], R[3]);
                break;
            default:
                printf("Fatal error in %s at %s:%d -- unknown bytecode 0x%02x\n",
                        __func__, __FILE__, __LINE__, *bc);
                abort();
        }
    }

    m.set(X, R[1]);
}

/* P+1 from the starting value b. The factor found (or 1) is stored in f.
 * Returns 1 if backtracking was used, 0 otherwise.
 */
template<typename layer>
static int pp1_run(typename layer::Integer & f,
        typename layer::Residue & b,
        typename layer::Residue const & two,
        typename layer::Modulus const & m,
        pp1_plan_t const & plan)
{
    typename layer::Residue t(m);
    int bt = 0;

    /* stage 1 */
    pp1_stage1<layer>(b, plan.bc, two, m);

    /* Backtracking for the 2's in the exponent */
    m.set(t, b);
    for (unsigned int i = 0; i < plan.exp2; i++) {
        m.V_dbl(b, b, two);
#if PP1_BACKTRACKING
        if (m.equal(b, two)) {
            m.set(b, t);
            bt = 1;
            break;
        }
        m.set(t, b);
#endif
    }
    m.sub(t, b, two);
    m.gcd(f, t);

    /* stage 2 */
    if (bt == 0 && f == uint64_t(1) && plan.B1 < plan.stage2.B2) {
        bt = pp1_stage2<layer>(t, b, plan.stage2, two, m);
        m.gcd(f, t);
    }

    return bt;
}

/* P+1 starting from 2/7 */
template<typename layer>
int pp1_27(typename layer::Integer & f,
        typename layer::Modulus const & m,
        pp1_plan_t const & plan)
{
    typename layer::Residue b(m), two(m);
    m.set1(two);
    m.add(two, two, two);
    m.set(b, two);
    m.div7(b, b);
    return pp1_run<layer>(f, b, two, m, plan);
}

/* P+1 starting from 6/5 */
template<typename layer>
int pp1_65(typename layer::Integer & f,
        typename layer::Modulus const & m,
        pp1_plan_t const & plan)
{
    typename layer::Residue b(m), two(m);
    m.set1(two);
    m.add(two, two, two);
    m.set(b, two);
    m.add(b, b, two);
    m.add(b, b, two);
    m.div5(b, b);
    return pp1_run<layer>(f, b, two, m, plan);
}

template int pp1_27<arithxx_modredc64>(arithxx_modredc64::Integer &, arithxx_modredc64::Modulus const &, pp1_plan_t const &);
template int pp1_27<arithxx_modredc96>(arithxx_modredc96::Integer &, arithxx_modredc96::Modulus const &, pp1_plan_t const &);
template int pp1_27<arithxx_modredc126>(arithxx_modredc126::Integer &, arithxx_modredc126::Modulus const &, pp1_plan_t const &);
template int pp1_27<arithxx_mod_mpz_new>(arithxx_mod_mpz_new::Integer &, arithxx_mod_mpz_new::Modulus const &, pp1_plan_t const &);
template int pp1_65<arithxx_modredc64>(arithxx_modredc64::Integer &, arithxx_modredc64::Modulus const &, pp1_plan_t const &);
template int pp1_65<arithxx_modredc96>(arithxx_modredc96::Integer &, arithxx_modredc96::Modulus const &, pp1_plan_t const &);
template int pp1_65<arithxx_modredc126>(arithxx_modredc126::Integer &, arithxx_modredc126::Modulus const &, pp1_plan_t const &);
template int pp1_65<arithxx_mod_mpz_new>(arithxx_mod_mpz_new::Integer &, arithxx_mod_mpz_new::Modulus const &, pp1_plan_t const &);
