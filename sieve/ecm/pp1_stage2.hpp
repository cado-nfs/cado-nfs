#ifndef CADO_SIEVE_ECM_PP1_STAGE2_HPP
#define CADO_SIEVE_ECM_PP1_STAGE2_HPP

#include <vector>

#include "macros.h"
#include "stage2.h"

// #define PP1_STAGE2_DEBUG /* define to print debug information for stage 2 */

#ifdef PP1_STAGE2_DEBUG
#include "fmt/base.h"
#endif

/* Do we want backtracking in stage 2? */
#ifndef PP1_BACKTRACKING
/* Default is "yes." Set to 0 for "no." */
#define PP1_BACKTRACKING 1
#endif

/* Stage 2 of P+1, which is also that of P-1 (with X = x + 1/x).
 *
 * Computes r = prod_{(u,v) in plan->pairs} (V_{v*w}(X) - V_u(X)), where
 * V_n is the Lucas sequence V_n(X, 1), V_0 = 2, V_1 = X,
 * V_{n+2} = X * V_{n+1} - V_n. See stage2_make_plan in stage2.c.
 *
 * Returns 1 if backtracking was used (a zero product was replaced by the
 * product before the last increase of v), 0 otherwise.
 */
template<typename layer>
int pp1_stage2(typename layer::Residue & r,
        typename layer::Residue const & X,
        stage2_plan_t const & plan,
        typename layer::Residue const & two,
        typename layer::Modulus const & m)
{
    using Residue = typename layer::Residue;

    ASSERT(plan.w % 6 == 0); /* see stage2_make_plan */
    ASSERT(plan.vmin <= plan.vmax);

    Residue Xw(m), t(m); /* V_w (X) and a temp */
    /* V_{v*w} (X) for vmin <= v <= vmax, and V_u (X) for u in U */
    std::vector<Residue> Xvw, Xu;
    Xvw.reserve(plan.vmax - plan.vmin + 1);
    for (unsigned int i = 0; i < plan.vmax - plan.vmin + 1; i++)
        Xvw.emplace_back(m);
    Xu.reserve(plan.U_len);
    for (unsigned int i = 0; i < plan.U_len; i++)
        Xu.emplace_back(m);
    int bt = 0;

#ifdef PP1_STAGE2_DEBUG
    fmt::print("# {}: w={}\n", __func__, plan.w);
    fmt::print("# {}: input X={}\n", __func__, m.get(X));
#endif

    { /***************************** baby step **************************/
        /* Compute V_u (X) for u in U. Compute all the u, 1 <= u < d/2,
         * gcd(u,w)=1 with two arithmetic progressions 1+6k and 5+6k (this
         * assumes 6|w). We need two values of each progression (1, 7 and
         * 5, 11) and the common difference 6. These can be computed with
         * the Lucas chain 1, 2, 3, 5, 6, 7, 11 at the cost of 4 dadd (=1M)
         * and 2 dbl (1S). If w=30, we could use 1,2,3,4,6,7,11,13 which
         * has 4 dadd and 3 dbl.
         */
        Residue ap1_0(m), ap1_1(m), ap5_0(m), ap5_1(m), X2(m), X6(m);

        /* ap1_0 = V_1(X), ap1_1 = V_7(X), ap5_0 = V_5(X), ap5_1 = V_11(X)
         * and X6 = V_6(X)
         */
        m.set(ap1_0, X);                 /* ap1_0 = V_1(X) */
        m.V_dbl(X2, X, two);             /* X2 = V_2(X) */
        m.V_dadd(X6, X2, X, X);          /* V_3(X) = V_2(X)*V_1(X) - V_1(X) */
        m.V_dadd(ap5_0, X6, X2, X);      /* V_5(X) = V_3(X)*V_2(X) - V_1(X) */
        m.V_dbl(X6, X6, two);            /* V_6(X) = V_3(X)*V_3(X) - 2 */
        m.V_dadd(ap1_1, X6, X, ap5_0);   /* V_7(X) = V_6(X)*V_1(X) - V_5(X) */
        m.V_dadd(ap5_1, X6, ap5_0, X);   /* V_11(X) = V_6(X)*V_5(X) - V_1(X) */

#ifdef PP1_STAGE2_DEBUG
        fmt::print("# {}: V_1 (X)={}\n", __func__, m.get(ap1_0));
        fmt::print("# {}: V_5 (X)={}\n", __func__, m.get(ap5_0));
        fmt::print("# {}: V_7 (X)={}\n", __func__, m.get(ap1_1));
        fmt::print("# {}: V_11 (X)={}\n", __func__, m.get(ap5_1));
#endif

        /* Now we generate all the V_u (X) for u in U. We treat the first
         * two manually because those might correspond to ap1_0 = V_1 (X)
         * and ap5_0 = V_5 (X) */
        unsigned int k = 0;
        if (plan.U_len > k && plan.U[k] == 1)
            m.set(Xu[k++], ap1_0);
        if (plan.U_len > k && plan.U[k] == 5)
            m.set(Xu[k++], ap5_0);

        unsigned int u_1mod6 = 7, u_5mod6 = 11;
        while (k < plan.U_len) {
            if (plan.U[k] == u_1mod6) {
                m.set(Xu[k++], ap1_1);
                continue;
            }
            if (plan.U[k] == u_5mod6) {
                m.set(Xu[k++], ap5_1);
                continue;
            }

            m.V_dadd(t, ap1_1, X6, ap1_0);
            m.set(ap1_0, ap1_1);
            m.set(ap1_1, t);
            u_1mod6 += 6;

            m.V_dadd(t, ap5_1, X6, ap5_0);
            m.set(ap5_0, ap5_1);
            m.set(ap5_1, t);
            u_5mod6 += 6;
        }

        /* Compute Xw = V_w (X) */
        if (plan.w == 6) {
            m.set(Xw, X6);
        } else if (plan.w == 12) {
            m.V_dbl(Xw, X6, two);
        } else if (plan.w % 12 == 0) {
            m.V_dadd(t, ap1_1, X6, ap1_0);
            m.V_dadd(Xw, t, ap5_1, X2);
        } else { /* plan.w % 12 == 6 */
            m.V_dbl(X2, X2, two); /* We need V_4 (X) for difference */
            m.V_dadd(Xw, ap1_1, ap5_1, X2);
        }

#ifdef PP1_STAGE2_DEBUG
        for (unsigned int k = 0; k < plan.U_len; k++)
            fmt::print("# {}: V_{} (X)={}\n", __func__, plan.U[k], m.get(Xu[k]));
        fmt::print("# {}: V_{} (X)={}\n", __func__, plan.w, m.get(Xw));
#endif
    }

    { /***************************** giant step *************************/
        /* Compute V_{v*w} (X) for vmin <= v <= vmax */
        typename layer::Integer const vmin(plan.vmin);
        if (plan.vmin == plan.vmax) {
            /* only V_{vmin*w} (X) is needed */
            m.V(Xvw[0], nullptr, Xw, vmin);
        } else {
            m.V(Xvw[0], &Xvw[1], Xw, vmin);
            unsigned int k = 2;
            for (unsigned int v = plan.vmin + 2; v <= plan.vmax; k++, v++) {
                if (v % 2 == 0 && v / 2 >= plan.vmin)
                    m.V_dbl(Xvw[k], Xvw[v / 2 - plan.vmin], two);
                else
                    m.V_dadd(Xvw[k], Xvw[k - 1], Xw, Xvw[k - 2]);
            }
        }

#ifdef PP1_STAGE2_DEBUG
        for (unsigned int k = 0; k < plan.vmax - plan.vmin + 1; k++)
            fmt::print("# {}: V_{}w (X)={}\n", __func__, plan.vmin + k, m.get(Xvw[k]));
#endif
    }

    { /***************************** product ****************************/
        /* r = prod_{u,v in plan.pairs}{V_{v*w} (X) - V_u (X)} */
        Residue r_bak(m); /* Backup value of r, in case we get r == 0 */

        m.set1(r);
        m.set(r_bak, r);

        unsigned int const * u_ptr = plan.pairs;
        for (unsigned int v_idx = 0;; v_idx++, u_ptr++) {
            ASSERT(v_idx <= plan.vmax - plan.vmin);
#if PP1_BACKTRACKING
            m.set(r_bak, r);
#endif
            for (; *u_ptr != PAIR_END && *u_ptr != PAIR_INCR_V; u_ptr++) {
                m.sub(t, Xvw[v_idx], Xu[*u_ptr]);
                m.mul(r, r, t);
            }

#if PP1_BACKTRACKING
            /* If we got r == 0, restore the previous value of r and end
             * stage 2. Let's hope not all factors were found since the
             * last increase of v.
             */
            if (m.is0(r)) {
                m.set(r, r_bak);
                bt = 1;
                break;
            }
#endif

            if (*u_ptr == PAIR_END)
                break;
        }

#ifdef PP1_STAGE2_DEBUG
        fmt::print("# {}: r={}\n", __func__, m.get(r));
#endif
    }

    return bt;
}

#endif /* CADO_SIEVE_ECM_PP1_STAGE2_HPP */
