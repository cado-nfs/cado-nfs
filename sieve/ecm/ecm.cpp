#include "cado.h" // IWYU pragma: keep

/* define to print number of operations of ECM (/!\ not thread-safe) */
//#define ECM_COUNT_OPS

//#define ECM_DEBUG /* define to print debug information */
//#define ECM_STAGE2_DEBUG /* define to print debug information for stage 2 */

#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <vector>

#include "arithxx/modredc64.hpp"
#include "bytecode.h"
#include "ec_arith_common.hpp"
#include "ec_arith_Edwards.hpp"
#include "ec_arith_Montgomery.hpp"
#include "ec_arith_Weierstrass_new.hpp"
#include "ec_parameterization.hpp"
#include "ecm.hpp"
#include "facul_arithxx.hpp"
#include "facul_ecm.h"
#include "macros.h"
#include "stage2.h"

/* Do we want backtracking when processing factors of 2 in E? */
#ifndef ECM_BACKTRACKING
/* Default is "yes." Set to 0 for "no." */
#define ECM_BACKTRACKING 1
#endif

#ifdef ECM_COUNT_OPS
static unsigned int _count_stage2_common_z, _count_stage2_product;
#define ECM_COUNT_OPS_STAGE1_TOTAL_M EDWARDS_COUNT_OPS_M+MONTGOMERY_COUNT_OPS_M
#define ECM_COUNT_OPS_STAGE2_TOTAL_M MONTGOMERY_COUNT_OPS_M \
                                   + _count_stage2_common_z \
                                   + _count_stage2_product
#define ECM_COUNT_OPS_RESET() do {                             \
      EDWARDS_COUNT_OPS_RESET(); MONTGOMERY_COUNT_OPS_RESET(); \
      _count_stage2_common_z = _count_stage2_product = 0;      \
    } while (0)
#endif

/******************************************************************************/
/*********************** bytecode interpreter functions ***********************/
/******************************************************************************/

template<typename layer>
using ec_points = std::vector<ec_point<layer>>;

template<typename layer>
static ec_points<layer> ec_points_alloc(size_t n, typename layer::Modulus const & m)
{
    ec_points<layer> R;
    R.reserve(n);
    for (size_t i = 0; i < n; i++)
        R.emplace_back(m);
    return R;
}

/**** Interpret the double-base chain bytecode for Twisted Edwards curves ****/
/* Return a pointer to the last parsed byte.
 * Assume that the bytecode is correct (for example, assume that the size
 * of the array is correct).
 */
template<typename layer>
static bytecode_const
bytecode_dbchain_interpret_edwards_internal(bytecode_const bc,
        ec_point<layer> * R, typename layer::Modulus const & m)
{
    for (;;) {
        uint8_t op, f, s, n, pow2 = 0, pow3 = 0;
        bytecode_elt_split_2_1_1_4(&op, &f, &s, &n, *bc);

        if (op == DBCHAIN_OP_TPLADD || op == DBCHAIN_OP_TPLDBLADD) {
            bc++;
            pow3 = bytecode_elt_to_uint8(*bc);
        }
        if (op == DBCHAIN_OP_DBLADD || op == DBCHAIN_OP_TPLDBLADD) {
            bc++;
            pow2 = bytecode_elt_to_uint8(*bc);
        }

        for (uint8_t i = 0; i < pow3; i++)
            edwards_tpl(R[0], R[0], m, (i + 1 == pow3 + pow2) ?
                    TWISTED_EDWARDS_ext : TWISTED_EDWARDS_proj);
        for (uint8_t i = 0; i < pow2; i++)
            edwards_dbl(R[0], R[0], m, (i + 1 == pow2) ?
                    TWISTED_EDWARDS_ext : TWISTED_EDWARDS_proj);

        ec_point_coord_type_t output;
        if (f) {
            uint8_t t;
            bytecode_elt_split_4_4(&t, nullptr, bc[1]);
            if (bc[1] == MISHMASH_FINAL || t == MISHMASH_PRAC_BLOCK)
                output = MONTGOMERY_xz;
            else
                output = TWISTED_EDWARDS_ext;
        } else {
            output = TWISTED_EDWARDS_proj;
        }

        if (s == 0)
            edwards_add(R[f], R[0], R[n], m, output);
        else
            edwards_sub(R[f], R[0], R[n], m, output);

        if (f) /* is it finished ? */
            break;
        bc++; /* go to next byte */
    }

    return bc;
}

/********* Interpret the precomp bytecode for Twisted Edwards curves **********/
/* Return a pointer to the last parsed byte.
 * Assume that the bytecode is correct (for example, assume that the size
 * of the array is correct).
 */
template<typename layer>
static bytecode_const
bytecode_precomp_interpret_edwards_internal(bytecode_const bc,
        ec_point<layer> * R, typename layer::Modulus const & m)
{
    while (*bc != PRECOMP_FINAL) {
        uint8_t op, a, s, i, j, k, pow2, pow3;
        bytecode_elt_split_2_1_1_4(&op, &a, &s, &k, *bc);
        ec_point_coord_type_t const out = a ? TWISTED_EDWARDS_ext : TWISTED_EDWARDS_proj;

        switch (op) {
            case PRECOMP_OP_ADD:
                bc++;
                bytecode_elt_split_4_4(&i, &j, *bc);
                if (s == 0)
                    edwards_add(R[k], R[i], R[j], m, out);
                else
                    edwards_sub(R[k], R[i], R[j], m, out);
                break;
            case PRECOMP_OP_DBL:
                bc++;
                pow2 = bytecode_elt_to_uint8(*bc);
                for (uint8_t i = 0; i < pow2; i++)
                    edwards_dbl(R[0], R[0], m, (i + 1 == pow2 && a) ?
                            TWISTED_EDWARDS_ext : TWISTED_EDWARDS_proj);
                if (k > 0)
                    ec_point_set(R[k], R[0], m, out);
                break;
            case PRECOMP_OP_TPL:
                bc++;
                pow3 = bytecode_elt_to_uint8(*bc);
                for (uint8_t i = 0; i < pow3; i++)
                    edwards_tpl(R[0], R[0], m, (i + 1 == pow3 && a) ?
                            TWISTED_EDWARDS_ext : TWISTED_EDWARDS_proj);
                if (k > 0)
                    ec_point_set(R[k], R[0], m, out);
                break;
            default:
                printf("Fatal error in %s at %s:%d -- unknown bytecode 0x%02x\n",
                        __func__, __FILE__, __LINE__, *bc);
                abort();
        }

        bc++; /* go to next byte */
    }

    return bc;
}

/************* Interpret the PRAC bytecode for Montgomery curves *************/
/* Return a pointer to the last parsed byte.
 * Assume that the bytecode is correct (for example, assume that the size
 * of the array R is >= 5).
 */
template<typename layer>
static bytecode_const
bytecode_prac_interpret_montgomery_internal(bytecode_const bc,
        ec_point<layer> * R, typename layer::Modulus const & m,
        typename layer::Residue const & b)
{
    for (bool finished = false; !finished; ) {
        switch (*bc) {
            case PRAC_SWAP: /* [ = 's' ] Swap R[0], R[1] */
                ec_point_swap(R[0], R[1], MONTGOMERY_xz);
                break;
            case PRAC_SUBBLOCK_INIT: /* [ = 'i' ] Start of a sub-block */
                ec_point_set(R[1], R[0], m, MONTGOMERY_xz);
                ec_point_set(R[2], R[0], m, MONTGOMERY_xz);
                montgomery_dbl(R[0], R[0], m, b);
                break;
            case PRAC_SUBBLOCK_FINAL: /* [ = 'f' ] End of a sub-block */
                montgomery_dadd(R[0], R[0], R[1], R[2], b, m);
                break;
            case PRAC_BLOCK_FINAL: /* [ = 'F' ] End of the block */
                montgomery_dadd(R[1], R[0], R[1], R[2], b, m);
                finished = true;
                break;
            case 1:
                montgomery_dadd(R[3], R[0], R[1], R[2], b, m);
                montgomery_dadd(R[4], R[3], R[0], R[1], b, m);
                montgomery_dadd(R[1], R[1], R[3], R[0], b, m);
                ec_point_set(R[0], R[4], m, MONTGOMERY_xz);
                break;
            case 2:
                montgomery_dadd(R[1], R[0], R[1], R[2], b, m);
                montgomery_dbl(R[0], R[0], m, b);
                break;
            case 3:
                montgomery_dadd(R[2], R[1], R[0], R[2], b, m);
                ec_point_swap(R[1], R[2], MONTGOMERY_xz);
                break;
            case 4:
                montgomery_dadd(R[1], R[1], R[0], R[2], b, m);
                montgomery_dbl(R[0], R[0], m, b);
                break;
            case 5:
                montgomery_dadd(R[2], R[2], R[0], R[1], b, m);
                montgomery_dbl(R[0], R[0], m, b);
                break;
            case 6:
                montgomery_dbl(R[3], R[0], m, b);
                montgomery_dadd(R[4], R[0], R[1], R[2], b, m);
                montgomery_dadd(R[0], R[3], R[0], R[0], b, m);
                montgomery_dadd(R[2], R[3], R[4], R[2], b, m);
                ec_point_swap(R[1], R[2], MONTGOMERY_xz);
                break;
            case 7:
                montgomery_dadd(R[3], R[0], R[1], R[2], b, m);
                montgomery_dadd(R[1], R[3], R[0], R[1], b, m);
                montgomery_dbl(R[3], R[0], m, b);
                montgomery_dadd(R[0], R[0], R[3], R[0], b, m);
                break;
            case 8:
                montgomery_dadd(R[3], R[0], R[1], R[2], b, m);
                montgomery_dadd(R[2], R[2], R[0], R[1], b, m);
                ec_point_swap(R[1], R[3], MONTGOMERY_xz);
                montgomery_dbl(R[3], R[0], m, b);
                montgomery_dadd(R[0], R[0], R[3], R[0], b, m);
                break;
            case 9:
                montgomery_dadd(R[2], R[2], R[1], R[0], b, m);
                montgomery_dbl(R[1], R[1], m, b);
                break;
            case 10:
                /* Combined final add of old subchain and init of new
                 * subchain [=fi] */
                montgomery_dadd(R[1], R[0], R[1], R[2], b, m);
                ec_point_set(R[2], R[1], m, MONTGOMERY_xz);
                montgomery_dbl(R[0], R[1], m, b);
                break;
            case 11:
                /* Combined rule 3 and rule 0 [=\x3s] */
                montgomery_dadd(R[2], R[1], R[0], R[2], b, m);
                /* (R[1],R[2],R[0]) := (R[0],R[1],R[2])  */
                ec_point_swap(R[1], R[2], MONTGOMERY_xz);
                ec_point_swap(R[0], R[1], MONTGOMERY_xz);
                break;
            case 12:
                /* Combined rule 3, then subchain end/start [=\x3fi] */
                montgomery_dadd(R[3], R[1], R[0], R[2], b, m);
                montgomery_dadd(R[2], R[0], R[3], R[1], b, m);
                ec_point_set(R[1], R[2], m, MONTGOMERY_xz);
                montgomery_dbl(R[0], R[2], m, b);
                break;
            case 13:
                /* Combined rule 3, swap, rule 3 and swap, merged a bit
                 * [=\x3s\x3s] */
                ec_point_set(R[3], R[1], m, MONTGOMERY_xz);
                montgomery_dadd(R[1], R[1], R[0], R[2], b, m);
                ec_point_set(R[2], R[0], m, MONTGOMERY_xz);
                montgomery_dadd(R[0], R[0], R[1], R[3], b, m);
                break;
            default:
                printf("Fatal error in %s at %s:%d -- unknown bytecode 0x%02x\n",
                        __func__, __FILE__, __LINE__, *bc);
                abort();
        }

        if (!finished)
            bc++; /* go to next byte */
    }

    return bc;
}

template<typename layer>
static void bytecode_prac_interpret_montgomery(ec_point<layer> & P,
        bytecode_const bc, typename layer::Modulus const & m,
        typename layer::Residue const & b)
{
    if (bc == nullptr)
        return;

    /* we need 5 points: 3 for PRAC + 2 temporary points */
    auto R = ec_points_alloc<layer>(5, m);

    /* current point (here starting point) go into R[0] at init */
    ec_point_set(R[0], P, m, MONTGOMERY_xz);

    bytecode_prac_interpret_montgomery_internal<layer>(bc, R.data(), m, b);

    /* output is in R[1] */
    ec_point_set(P, R[1], m, MONTGOMERY_xz);
}

/** Interpret the MISHMASH bytecode on Twisted Edwards and Montgomery curves **/
template<typename layer>
static void bytecode_mishmash_interpret_mixed_repr(ec_point<layer> & P,
        bytecode_const bc, typename layer::Modulus const & m,
        typename layer::Residue const & b)
{
    uint8_t n;

    /* parse first byte to alloc R and go to the next byte */
    bytecode_elt_split_4_4(nullptr, &n, *bc);
    bc++;

    auto R = ec_points_alloc<layer>(n + 2, m);

    /* starting point go into R[1] at init */
    ec_point_set(R[1], P, m, TWISTED_EDWARDS_ext);

    while (*bc != MISHMASH_FINAL) {
        uint8_t t, n;
        bytecode_elt_split_4_4(&t, &n, *bc);

        if (n != 0)
            ec_point_set(R[0], R[n], m, TWISTED_EDWARDS_ext);

        if (t == MISHMASH_DBCHAIN_BLOCK) {
            bc = bytecode_dbchain_interpret_edwards_internal<layer>(++bc, R.data(), m);
        } else if (t == MISHMASH_PRECOMP_BLOCK) {
            bc = bytecode_precomp_interpret_edwards_internal<layer>(++bc, R.data(), m);
        } else if (t == MISHMASH_PRAC_BLOCK) {
            bc = bytecode_prac_interpret_montgomery_internal<layer>(++bc, R.data(), m, b);
        } else { /* unexpected bytecode */
            printf("Fatal error in %s at %s:%d -- unexpected bytecode 0x%02x\n",
                    __func__, __FILE__, __LINE__, *bc);
            abort();
        }

        bc++; /* go to next byte */
    }

    /* output (always in Montgomery form) is in R[1] */
    ec_point_set(P, R[1], m, MONTGOMERY_xz);
}

/******************************************************************************/
/********************************* stage 2 ************************************/
/******************************************************************************/

/* Multiplies x[1] by z[2]*z[3]*z[4]...*z[n],
   x[2] by z[1]*z[3]*z[4]...*z[n] etc., generally
   x[i] by \prod_{1\leq j \leq n, j\neq i} z[j]
   Requires n > 1. Uses 4n-6 multiplications.
 * If n <= 1 there is nothing to do.
 */
template<typename layer>
static void ATTRIBUTE((__noinline__))
common_z(ec_point<layer> * L1, unsigned int const n1,
        ec_point<layer> * L2, unsigned int const n2,
        typename layer::Modulus const & m)
{
    unsigned int const n = n1 + n2;

    if (n <= 1) /* nothing to do in this case */
        return;

#ifdef ECM_COUNT_OPS
    _count_stage2_common_z += 4*n - 6;
#endif

    typename layer::Residue p(m);
    std::vector<typename layer::Residue> t;
    t.reserve(n - 1);
    for (unsigned int i = 0; i < n - 1; i++)
        t.emplace_back(m);

    /* Set t[i] = z_0 * z_1 * ... * z_i, where the z_i are the
     * Z-coordinates taken from the two lists L1 and L2.
     */
    if (n1)
        m.set(t[0], L1[0].z);
    else
        m.set(t[0], L2[0].z); /* L1 is empty */

    for (unsigned int i = 1; i < n - 1; i++) {
        if (i < n1)
            m.mul(t[i], t[i - 1], L1[i].z);
        else
            m.mul(t[i], t[i - 1], L2[i - n1].z);
    }
    /* cost: n-2 mul */

    if (n2) {
        m.mul(L2[n2 - 1].x, L2[n2 - 1].x, t[n - 2]);
        m.set(p, L2[n2 - 1].z);
    } else {
        m.mul(L1[n1 - 1].x, L1[n1 - 1].x, t[n - 2]);
        m.set(p, L1[n1 - 1].z);
    }
    /* cost: 1 mul */

    for (unsigned int i = n - 2; i > 0; i--) {
        if (i < n1) {
            m.mul(L1[i].x, L1[i].x, p);
            m.mul(L1[i].x, L1[i].x, t[i - 1]);
            m.mul(p, p, L1[i].z);
        } else {
            m.mul(L2[i - n1].x, L2[i - n1].x, p);
            m.mul(L2[i - n1].x, L2[i - n1].x, t[i - 1]);
            m.mul(p, p, L2[i - n1].z);
        }
    }
    /* cost: 3*n-6 mul */

    if (n1)
        m.mul(L1[0].x, L1[0].x, p);
    else
        m.mul(L2[0].x, L2[0].x, p);
    /* cost: 1 mul */
}

template<typename layer>
static int ATTRIBUTE((__noinline__))
ecm_stage2(typename layer::Residue & r, ec_point<layer> const & P,
        stage2_plan_t const & plan, typename layer::Residue const & b,
        typename layer::Modulus const & m)
{
    ASSERT(plan.w % 6 == 0); /* see stage2_make_plan */
    ASSERT(plan.vmin <= plan.vmax);

    ec_point<layer> wP(m), Pt(m); /* w*P and a temp */
    /* saved v*w*P, vmin <= v <= vmax, and uP, u in U. */
    auto vwP = ec_points_alloc<layer>(plan.vmax - plan.vmin + 1, m);
    auto uP = ec_points_alloc<layer>(plan.U_len, m);
    int bt = 0;

#ifdef ECM_STAGE2_DEBUG
    fmt::print("# {}: w={}\n", __func__, plan.w);
    fmt::print("# {}: input P=", __func__);
    ec_point_fprintf(stdout, P, MONTGOMERY_xz, m);
    fputc('\n', stdout);
#endif

    { /***************************** baby step **************************/
        /* Compute u*P for u in U. Compute all the u, 1 <= u < d/2,
         * gcd(u,w)=1 with two arithmetic progressions 1+6k and 5+6k
         * (this assumes 6|w). We need two values of each progression (1, 7
         * and 5, 11) and the common difference 6. These can be computed
         * with the Lucas chain 1, 2, 3, 5, 6, 7, 11 at the cost of 4
         * additions and 2 doublings. If w=30, we could use
         * 1,2,3,4,6,7,11,13 which has 4 additions and 3 doublings.
         */
        ec_point<layer> ap1_0(m), ap1_1(m), ap5_0(m), ap5_1(m), P2(m), P6(m);

        /* Init ap1_0 = 1P, ap1_1 = 7P, ap5_0 = 5P, ap5_1 = 11P and P6 = 6P */
        ec_point_set(ap1_0, P, m, MONTGOMERY_xz);   /* ap1_0 = 1*P */
        montgomery_dbl(P2, P, m, b);                /* P2 = 2*P */
        montgomery_dadd(P6, P2, P, P, b, m);        /* P6 = 3*P (for now) */
        montgomery_dadd(ap5_0, P6, P2, P, b, m);    /* 5*P = 3*P + 2*P */
        montgomery_dbl(P6, P6, m, b);               /* P6 = 6*P = 2*(3*P) */
        montgomery_dadd(ap1_1, P6, P, ap5_0, b, m); /* 7*P = 6*P + P */
        montgomery_dadd(ap5_1, P6, ap5_0, P, b, m); /* 11*P = 6*P + 5*P */

#ifdef ECM_STAGE2_DEBUG
        fmt::print("# {}: 1*P=", __func__);
        ec_point_fprintf(stdout, ap1_0, MONTGOMERY_xz, m);
        fmt::print("\n# {}: 5*P=", __func__);
        ec_point_fprintf(stdout, ap5_0, MONTGOMERY_xz, m);
        fmt::print("\n# {}: 7*P=", __func__);
        ec_point_fprintf(stdout, ap1_1, MONTGOMERY_xz, m);
        fmt::print("\n# {}: 11*P=", __func__);
        ec_point_fprintf(stdout, ap5_1, MONTGOMERY_xz, m);
        fputc('\n', stdout);
#endif

        /* Now we generate all the u*P for u in U. We treat the first two
         * manually because those might correspond to ap1_0 = 1*P and
         * ap5_0 = 5*P */
        unsigned int k = 0;
        if (k < plan.U_len && plan.U[k] == 1)
            ec_point_set(uP[k++], ap1_0, m, MONTGOMERY_xz);
        if (k < plan.U_len && plan.U[k] == 5)
            ec_point_set(uP[k++], ap5_0, m, MONTGOMERY_xz);

        unsigned int u_1mod6 = 7, u_5mod6 = 11;
        while (k < plan.U_len) {
            if (plan.U[k] == u_1mod6) {
                ec_point_set(uP[k++], ap1_1, m, MONTGOMERY_xz);
                continue;
            }
            if (plan.U[k] == u_5mod6) {
                ec_point_set(uP[k++], ap5_1, m, MONTGOMERY_xz);
                continue;
            }

            montgomery_dadd(Pt, ap1_1, P6, ap1_0, b, m);
            ec_point_set(ap1_0, ap1_1, m, MONTGOMERY_xz);
            ec_point_set(ap1_1, Pt, m, MONTGOMERY_xz);
            u_1mod6 += 6;

            montgomery_dadd(Pt, ap5_1, P6, ap5_0, b, m);
            ec_point_set(ap5_0, ap5_1, m, MONTGOMERY_xz);
            ec_point_set(ap5_1, Pt, m, MONTGOMERY_xz);
            u_5mod6 += 6;
        }

        /* Compute wP = w*P */
        if (plan.w == 6) {
            ec_point_set(wP, P6, m, MONTGOMERY_xz);
        } else if (plan.w == 12) {
            montgomery_dbl(wP, P6, m, b);
        } else if (plan.w % 12 == 0) {
            montgomery_dadd(Pt, ap1_1, P6, ap1_0, b, m);
            montgomery_dadd(wP, Pt, ap5_1, P2, b, m);
        } else { /* plan.w % 12 == 6 */
            montgomery_dbl(P2, P2, m, b); /* We need 4P for difference */
            montgomery_dadd(wP, ap1_1, ap5_1, P2, b, m);
        }

#ifdef ECM_STAGE2_DEBUG
        for (unsigned int k = 0; k < plan.U_len; k++) {
            fmt::print("{}# {}: {}*P=", k ? "\n" : "", __func__, plan.U[k]);
            ec_point_fprintf(stdout, uP[k], MONTGOMERY_xz, m);
        }
        fmt::print("\n# {}: {}*P=", __func__, plan.w);
        ec_point_fprintf(stdout, wP, MONTGOMERY_xz, m);
        fputc('\n', stdout);
#endif
    }

    { /**************************** giant step **************************/
        /* Compute vwP for vmin <= v < vmax */
        if (plan.vmin == plan.vmax) {
            /* If vmin == vmax, only vmin*w*P is needed */
            montgomery_smul_ul(vwP[0], (ec_point<layer> *) nullptr, wP, plan.vmin, m, b);
        } else {
            montgomery_smul_ul(vwP[0], &vwP[1], wP, plan.vmin, m, b);
            unsigned int k = 2;
            for (unsigned int v = plan.vmin + 2; v <= plan.vmax; k++, v++) {
                if (v % 2 == 0 && v / 2 >= plan.vmin)
                    montgomery_dbl(vwP[k], vwP[v / 2 - plan.vmin], m, b);
                else
                    montgomery_dadd(vwP[k], vwP[k - 1], wP, vwP[k - 2], b, m);
            }
        }

#ifdef ECM_STAGE2_DEBUG
        for (unsigned int k = 0; k < plan.vmax - plan.vmin + 1; k++) {
            fmt::print("{}# {}: {}*w*P=", k ? "\n" : "", __func__, plan.vmin + k);
            ec_point_fprintf(stdout, vwP[k], MONTGOMERY_xz, m);
        }
        fputc('\n', stdout);
#endif
    }

    { /***************************** common_z ***************************/
        /* Now we've computed all the points we need, so multiply each by
         * the Z-coordinates of all the others, using Zimmermann's two
         * product-lists trick. If vmin == 0, then vwP[0] is the point at
         * infinity (0::0), so we skip that one
         */
        unsigned int const skip = (plan.vmin == 0) ? 1 : 0;
        common_z<layer>(uP.data(), plan.U_len, vwP.data() + skip,
                plan.vmax - plan.vmin + 1 - skip, m);

#ifdef ECM_STAGE2_DEBUG
        for (unsigned int k = 0; k < plan.U_len; k++)
            fmt::print("# {}: after common_z, {}*P=({:#x} :: Z)\n", __func__,
                    plan.U[k], ec_residue_get_mpz<layer>(uP[k].x, m));
        for (unsigned int k = 0; k < plan.vmax - plan.vmin + 1; k++)
            fmt::print("# {}: after common_z, {}*w*P=({:#x} :: Z)\n", __func__,
                    plan.vmin + k, ec_residue_get_mpz<layer>(vwP[k].x, m));
#endif
    }

    { /***************************** product ****************************/
        /* Now compute
         *        r = prod_{u,v in plan.pairs}{X_vwP - X_uP}
         *
         * Initialize r with uP[0], which contains the product of the
         * Z-coordinates of all the precomputed points, except the
         * Z-coordinate of uP[0]. Note that uP[0] is equal to P, and we know
         * that the Z-coordinate of P is coprime to the modulus m. It cost
         * nothing to add this in the product and we could catch a factor
         * p if, by chance, one of the Z-coordinates is 0 modulo p.
         */
        typename layer::Residue t(m), r_bak(m); /* Backup value of r, in case we get r == 0 */

        m.set(r, uP[0].x);
        m.set(r_bak, r);

        unsigned int const * u_ptr = plan.pairs;
        for (unsigned int v_idx = 0;; v_idx++, u_ptr++) {
            ASSERT(v_idx <= plan.vmax - plan.vmin);
#if ECM_BACKTRACKING
            m.set(r_bak, r);
#endif

            for (; *u_ptr != PAIR_END && *u_ptr != PAIR_INCR_V; u_ptr++) {
                m.sub(t, vwP[v_idx].x, uP[*u_ptr].x);
                m.mul(r, r, t);
#ifdef ECM_COUNT_OPS
                _count_stage2_product++;
#endif
            }

#if ECM_BACKTRACKING
            /* See if we got r == 0. If yes, restore previous r value and
             * end stage 2. Let's hope not all factors were found since the
             * last v increase.
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

#ifdef ECM_STAGE2_DEBUG
        fmt::print("# {}: r={}\n", __func__, ec_residue_get_mpz<layer>(r, m));
#endif
    }

    return bt;
}

/******************************************************************************/
/*********************************** ECM **************************************/
/******************************************************************************/

/* Stores any factor found in f (1 if no factor found).
   If back-tracking was used, returns 1, otherwise returns 0. */
template<typename layer>
static int ecm_run(typename layer::Integer & f,
        typename layer::Modulus const & m,
        ecm_plan_t const & plan)
{
#ifdef ECM_DEBUG
    fmt::print("# {}: start with B1={} and B2={} for m = {}\n",
            __func__, plan.B1, plan.stage2.B2, m.getmod());
    fmt::print("# {}: using parameterization {:#x} with parameter {}\n",
            __func__, int(plan.parameterization), plan.parameter);
#endif

    typename layer::Residue u(m), b(m);
    ec_point<layer> P(m), Pt(m);
    ec_point_coord_type_t param_output_type;
    unsigned int i;
    int r, bt = 0;

    f = typename layer::Integer(uint64_t(1));

    /**************************** parameterization **************************/
    if (plan.parameterization & FULLMONTY) {
        param_output_type = MONTGOMERY_xz;
        if (plan.parameterization == BRENT12)
            r = ec_parameterization_Brent_Suyama<layer>(b, P, plan.parameter, m);
        else if (plan.parameterization == MONTY12)
            r = ec_parameterization_Montgomery12<layer>(b, P, plan.parameter, m);
        else /* if (plan.parameterization == MONTY16) */
            r = ec_parameterization_Montgomery16<layer>(b, P, plan.parameter, m);
    } else if (plan.parameterization == MONTYTWED12) {
        param_output_type = plan.bc[1] & 0x80 ? MONTGOMERY_xz : TWISTED_EDWARDS_ext;
        r = ec_parameterization_Z6<layer>(&b, P, plan.parameter, param_output_type, m);
    } else {
        fprintf(stderr, "%s: unknown parameterization\n", __func__);
        abort();
    }

    if (r == 0) {
        m.gcd(f, P.x);
#ifdef ECM_DEBUG
        fmt::print("# {}: during parameterization, found factor {}\n",
                __func__, f);
#endif
        return 0;
    }

#ifdef ECM_DEBUG
    {
        typename layer::Residue A(m);
        ec_point<layer> PM(m);
        montgomery_A_from_b<layer>(A, b, m);
        fmt::print("# {}: starting values:\n", __func__);

        if (param_output_type == MONTGOMERY_xz) {
            montgomery_curve_fprintf(stdout, "# ecm:   ", A, &P, m);
            fmt::print("# {}:                     = ", __func__);
            montgomery_point_fprintf_affine(stdout, P, m);
            fputc('\n', stdout);
        } else {
            typename layer::Residue d(m);
            edwards_d_from_montgomery_A<layer>(d, A, m);
            edwards_ext_curve_fprintf(stdout, "# ecm:   ", d, &P, m);

            fmt::print("# {}:   Equivalent to Montgomery curve with:\n", __func__);
            montgomery_point_from_edwards_point(PM, P, 1, m);
            montgomery_curve_fprintf(stdout, "# ecm:     ", A, &PM, m);
        }
    }
#endif

    /******************************** stage 1 *******************************/
#ifdef ECM_COUNT_OPS
    ECM_COUNT_OPS_RESET();
#endif

    /* output is always in Montgomery form */
    if (plan.parameterization & FULLMONTY)
        bytecode_prac_interpret_montgomery<layer>(P, plan.bc, m, b);
    else if (plan.parameterization & FULLMONTYTWED)
        bytecode_mishmash_interpret_mixed_repr<layer>(P, plan.bc, m, b);

#ifdef ECM_DEBUG
    fmt::print("# {}: after stage 1 (without the powers of 2):\n", __func__);
    fmt::print("# {}:   output (X::Z) = ", __func__);
    ec_point_fprintf(stdout, P, MONTGOMERY_xz, m);
    fmt::print("\n# {}:                 = ", __func__);
    montgomery_point_fprintf_affine(stdout, P, m);
    fputc('\n', stdout);
#endif

    /* Add prime 2 in the desired power. If a zero residue for the
     * Z-coordinate is encountered, we backtrack to previous point and stop.
     * NOTE: This is not as effective as I hoped. It prevents trivial
     * factorizations only if after processing the odd part of the stage 1
     * multiplier, the resulting point has power-of-2 order on E_p for all
     * p|N. If that were to happen, the point probably had that (presumably
     * small on most E_p) power-of-2 order during the last couple of primes
     * processed in the precomputed Lucas chain, and then quite likely the
     * Lucas chain incorrectly used an addition of identical points, causing
     * the Z-coordinate to become zero, leading to 0 (mod N) before we even
     * get here. For example, using 10^6 composites from an RSA155 sieving
     * experiment, without backtracking we get N as the factor 456 times,
     * with backtracking still 360 times.
     * TODO: this could probably be fixed by treating 3 separately, too,
     * instead of putting it in the precomputed Lucas chain. Then the
     * probability that a point of very small order on all E_p is
     * encountered during the Lucas chain is reduced, and so the probability
     * of using curve addition erroneously.
     *
     * The following code assume P is in Montgomery form.
     */
    ec_point_set(Pt, P, m, MONTGOMERY_xz);
    for (i = 0; i < plan.exp2; i++) {
        montgomery_dbl(P, P, m, b);
#if ECM_BACKTRACKING
        if (m.is0(P.z)) {
            ec_point_set(P, Pt, m, MONTGOMERY_xz);
            bt = 1;
            break;
        }
        ec_point_set(Pt, P, m, MONTGOMERY_xz);
#endif
    }
    m.gcd(f, P.z);

#ifdef ECM_DEBUG
    fmt::print("# {}: after stage 1 and {} power(s) of 2 (out of {}):\n",
            __func__, i, plan.exp2);
    fmt::print("# {}:   output (X::Z) = ", __func__);
    ec_point_fprintf(stdout, P, MONTGOMERY_xz, m);
    fmt::print("\n# {}:                 = ", __func__);
    montgomery_point_fprintf_affine(stdout, P, m);
#if ECM_BACKTRACKING
    fmt::print("\n# {}:   backtracking was{} used", __func__, bt ? "" : " not");
#endif
    fmt::print("\n# {}:   gcd = {}\n", __func__, f);
#endif

#ifdef ECM_COUNT_OPS
    unsigned int tot_stage1_M = ECM_COUNT_OPS_STAGE1_TOTAL_M;
    fmt::print("# {}: count ops: for stage 1 with B1 = {}\n"
            "# {}: count ops:         Edwards: {} DBL {} DBLext {} TPL {} TPLext "
                                              "{} ADD {} ADDext {} ADDmont\n"
            "# {}: count ops:      Montgomery: {} DBL {} dADD\n"
            "# {}: count ops:   stage 1 total: {} M\n", __func__, plan.B1,
            __func__, _count_edwards_dbl, _count_edwards_dblext,
            _count_edwards_tpl, _count_edwards_tplext, _count_edwards_add,
            _count_edwards_addext, _count_edwards_addmont, __func__,
            _count_montgomery_dbl, _count_montgomery_dadd, __func__,
            tot_stage1_M);
#endif

    /******************************** stage 2 *******************************/
#ifdef ECM_COUNT_OPS
    ECM_COUNT_OPS_RESET();
#endif

    if (bt == 0 && f == uint64_t(1) && plan.B1 < plan.stage2.B2) {
        bt = ecm_stage2<layer>(u, P, plan.stage2, b, m);
        m.gcd(f, u);
#ifdef ECM_DEBUG
        fmt::print("# {}: after stage 2, gcd={}\n", __func__, f);
#endif
    }
#ifdef ECM_DEBUG
    else
        fmt::print("# {}: stage 2 not done\n", __func__);
#endif

#ifdef ECM_COUNT_OPS
    unsigned int tot_stage2_M = ECM_COUNT_OPS_STAGE2_TOTAL_M;
    fmt::print("# {}: count ops: for stage 2 with B2 = {} [ with w = {} ]\n"
            "# {}: count ops:      Montgomery: {} DBL {} dADD\n"
            "# {}: count ops:        common_z: {} M\n"
            "# {}: count ops:         product: {} M\n"
            "# {}: count ops:   stage 2 total: {} M\n"
            "# {}: count ops: ECM total: {:5} M\n", __func__,
            plan.stage2.B2, plan.stage2.w, __func__, _count_montgomery_dbl,
            _count_montgomery_dadd, __func__, _count_stage2_common_z, __func__,
            _count_stage2_product, __func__, tot_stage2_M, __func__,
            tot_stage1_M + tot_stage2_M);
#endif

    return bt;
}

#define ECM_ENTRY_POINT(modint, modulus)                                \
int ecm(modint f, const modulus m, const ecm_plan_t * plan)             \
{                                                                       \
    return facul_arithxx_run(f, m,                                      \
            [&]<typename layer>(auto & g, auto const & mm) {            \
                return ecm_run<layer>(g, mm, *plan); });                \
}

ECM_ENTRY_POINT(modintredcul_t, modulusredcul_t)
ECM_ENTRY_POINT(modintredc15ul_t, modulusredc15ul_t)
ECM_ENTRY_POINT(modintredc2ul2_t, modulusredc2ul2_t)
ECM_ENTRY_POINT(modintmpz_t, modulusmpz_t)

/******************************************************************************/
/******************* parameters, point and curve orders ***********************/
/******************************************************************************/

int ec_parameter_is_valid(ec_parameterization_t const parameterization,
        unsigned long const parameter)
{
    switch (parameterization) {
        case BRENT12:
            return ec_parameterization_Brent_Suyama_is_valid(parameter);
        case MONTY12:
            return ec_parameterization_Montgomery12_is_valid(parameter);
        case MONTY16:
            return ec_parameterization_Montgomery16_is_valid(parameter);
        case MONTYTWED12:
            return ec_parameterization_Z6_is_valid(parameter);
        // TODO MONTYTWED16 ???
        default:
            printf("Fatal error in %s at %s:%d -- unknown parameterization %u\n",
                    __func__, __FILE__, __LINE__, parameterization);
            abort();
    }
}

/* return a valid parameter for the given parameterization, based on any
 * nonnegative integer. We make it so that:
 *  - the exceptional values, if any, come first.
 *  - all singular values are skipped.
 * this will eventually yield a finite number of exceptions, so a
 * finite-length code.
 */
unsigned long
ec_valid_parameter_from_sequence(ec_parameterization_t const parameterization,
        unsigned long const sequence_value)
{
    switch (parameterization) {
        case BRENT12:
            return ec_parameterization_Brent_Suyama_valid_parameter_from_sequence(sequence_value);
        case MONTY12:
            return ec_parameterization_Montgomery12_valid_parameter_from_sequence(sequence_value);
        case MONTY16:
            return ec_parameterization_Montgomery16_valid_parameter_from_sequence(sequence_value);
        case MONTYTWED12:
            return ec_parameterization_Z6_valid_parameter_from_sequence(sequence_value);
        // TODO MONTYTWED16 ???
        default:
            printf("Fatal error in %s at %s:%d -- unknown parameterization %u\n",
                    __func__, __FILE__, __LINE__, parameterization);
            abort();
    }
}

/* The curve given by parameterization and parameter, in Montgomery form:
 * its coefficient A and its starting point P. Returns 0 if an inversion
 * failed. */
template<typename layer>
static int ec_parameterization_Montgomery_curve(typename layer::Residue & A,
        ec_point<layer> & P, ec_parameterization_t const parameterization,
        unsigned long const parameter, typename layer::Modulus const & m)
{
    typename layer::Residue b(m);
    int r;

    if (parameterization == BRENT12)
        r = ec_parameterization_Brent_Suyama<layer>(b, P, parameter, m);
    else if (parameterization == MONTY12)
        r = ec_parameterization_Montgomery12<layer>(b, P, parameter, m);
    else if (parameterization == MONTY16)
        r = ec_parameterization_Montgomery16<layer>(b, P, parameter, m);
    else if (parameterization == MONTYTWED12)
        r = ec_parameterization_Z6<layer>(&b, P, parameter, MONTGOMERY_xz, m);
    else {
        fprintf(stderr, "%s: Unknown parameterization\n", __func__);
        abort();
    }

    if (r)
        montgomery_A_from_b<layer>(A, b, m);
    return r;
}

/* The variables "parameterization" and "parameter" are used to compute a
 * curve E and a point P over the finite field of size p. This function
 * returns the order of the point P.
 */
uint64_t ec_parameterization_point_order(ec_parameterization_t const parameterization,
        unsigned long const parameter, uint64_t const known_m,
        uint64_t const known_r, uint64_t const p, int const verbose)
{
    using layer = arithxx_modredc64;
    layer::Modulus const m(p);
    layer::Residue A(m), a(m), xw(m), yw(m);
    ec_point<layer> P(m);

    if (!ec_parameterization_Montgomery_curve<layer>(A, P, parameterization, parameter, m))
        return 0;

    if (verbose >= 2)
        printf("%s: Montgomery curve: B*Y^2 = X^3 + A*X^2*Z + X*Z^2, "
                "A = %" PRIu64 ", with point (X::Z) = (%" PRIu64 " :: %" PRIu64 ")\n",
                __func__,
                uint64_t(m.get(A)), uint64_t(m.get(P.x)), uint64_t(m.get(P.z)));

    if (!weierstrass_aff_from_montgomery<layer>(a, xw, yw, A, P, m))
        return 0;

    ECWeierstrass<layer> const E(m, a);
    typename ECWeierstrass<layer>::AffinePoint const Pw(E, xw, yw);
    return Pw.point_order(known_m, known_r, verbose);
}

/* Return the order of the curve defined by a parameterization and a
 * parameter, modulo p. Return 0 if a inversion failed during the
 * computation (which indicates that the curve is not defined modulo p).
 */
uint64_t ec_parameterization_curve_order(ec_parameterization_t const parameterization,
        unsigned long const parameter, uint64_t const p)
{
    using layer = arithxx_modredc64;
    layer::Modulus const m(p);
    layer::Residue A(m);
    ec_point<layer> P(m);

    if (!ec_parameterization_Montgomery_curve<layer>(A, P, parameterization, parameter, m))
        return 0;

    return montgomery_curve_order<layer>(A, P, m);
}
