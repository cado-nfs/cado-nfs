/* This file is not compilable as-is, it must be included by another file. */
#ifndef mod_init
  #error "One of the mod*_default.h headers must be included before this file"
#endif

#include "cado.h" // IWYU pragma: keep
#include <math.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include "facul_ecm.h"

/* define to print number of operations of ECM (/!\ not thread-safe) */
//#define ECM_COUNT_OPS

//#define ECM_DEBUG /* define to print debug information */
//#define ECM_STAGE2_DEBUG /* define to print debug information for stage 2 */

/* Define to 1 to make ellM_add() test if the two points are identical,
   and call ellM_double() if they are */
#ifndef ELLM_SAFE_ADD
#define ELLM_SAFE_ADD 0
#endif

/* Exported functions. */
#define ecm MOD_APPEND_TYPE(ecm)

/* Each static function are rename just before its declaration (there is no need
 * to rename them, but it's nice to have these functions distinguishable e.g. in
 * profiler output).
 */

#include "ec_arith_common.h"
#include "ec_arith_Edwards.h"
#include "ec_arith_Montgomery.h"
#include "ec_arith_Weierstrass.h"
#include "ec_parameterization.h"
#include "macros.h"     // CADO_STRINGIZE

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

/**** Interpret the double-base chain bytecode for Twisted Edwards curves ****/
/* Internal function.
 * Return a pointer to the last parsed byte.
 * Assume that the bytecode is correct (for example, assume that the size of the
 * array is correct).
 */
#define bytecode_dbchain_interpret_edwards_internal MOD_APPEND_TYPE(bytecode_dbchain_interpret_edwards_internal)
static bytecode_const
bytecode_dbchain_interpret_edwards_internal (bytecode_const bc, ec_point_t *R,
                                             const modulus_t m)
{
  while (1)
  {
    uint8_t op, f, s, n, pow2 = 0, pow3 = 0;
    bytecode_elt_split_2_1_1_4 (&op, &f, &s, &n, *bc);

    if (op == DBCHAIN_OP_TPLADD || op == DBCHAIN_OP_TPLDBLADD)
    {
      bc++;
      pow3 = bytecode_elt_to_uint8 (*bc);
    }
    if (op == DBCHAIN_OP_DBLADD || op == DBCHAIN_OP_TPLDBLADD)
    {
      bc++;
      pow2 = bytecode_elt_to_uint8 (*bc);
    }

    for (uint8_t i = 0; i < pow3; i++)
      edwards_tpl (R[0], R[0], m, (i+1==pow3+pow2) ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
    for (uint8_t i = 0; i < pow2; i++)
      edwards_dbl (R[0], R[0], m, (i+1==pow2) ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);

    ec_point_coord_type_t output;
    if (f)
    {
      uint8_t t;
      bytecode_elt_split_4_4 (&t, NULL, bc[1]);
      if (bc[1] == MISHMASH_FINAL || t == MISHMASH_PRAC_BLOCK)
        output = MONTGOMERY_xz;
      else
        output = TWISTED_EDWARDS_ext;
    }
    else
      output = TWISTED_EDWARDS_proj;

    if (s == 0)
      edwards_add (R[f], R[0], R[n], m, output);
    else
      edwards_sub (R[f], R[0], R[n], m, output);

    if (f) /* is it finished ? */
      break;
    else
      bc++; /* go to next byte */
  }

  return bc;
}

/********* Interpret the precomp bytecode for Twisted Edwards curves **********/
/* Internal function.
 * Return a pointer to the last parsed byte.
 * Assume that the bytecode is correct (for example, assume that the size of the
 * array is correct).
 */
#define bytecode_precomp_interpret_edwards_internal MOD_APPEND_TYPE(bytecode_precomp_interpret_edwards_internal)
static bytecode_const
bytecode_precomp_interpret_edwards_internal (bytecode_const bc, ec_point_t *R,
                                             const modulus_t m)
{
  while (*bc != PRECOMP_FINAL)
  {
    uint8_t op, a, s, i, j, k, pow2, pow3;
    bytecode_elt_split_2_1_1_4 (&op, &a, &s, &k, *bc);

    switch (op)
    {
      case PRECOMP_OP_ADD:
        bc++;
        bytecode_elt_split_4_4 (&i, &j, *bc);
        if (s == 0)
          edwards_add (R[k], R[i], R[j], m, a ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
        else
          edwards_sub (R[k], R[i], R[j], m, a ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
        break;
      case PRECOMP_OP_DBL:
        bc++;
        pow2 = bytecode_elt_to_uint8 (*bc);
        for (uint8_t i = 0; i < pow2; i++)
          edwards_dbl (R[0], R[0], m, (i+1==pow2 && a) ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
        if (k > 0)
          ec_point_set (R[k], R[0], m, a ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
        break;
      case PRECOMP_OP_TPL:
        bc++;
        pow3 = bytecode_elt_to_uint8 (*bc);
        for (uint8_t i = 0; i < pow3; i++)
          edwards_tpl (R[0], R[0], m, (i+1==pow3 && a) ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
        if (k > 0)
          ec_point_set (R[k], R[0], m, a ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
        break;
      default:
        printf ("Fatal error in %s at %s:%d -- unknown bytecode 0x%02x\n",
                __func__, __FILE__, __LINE__, *bc);
        abort ();
    }

    bc++; /* go to next byte */
  }

  return bc;
}

/************* Interpret the PRAC bytecode for Montgomery curves *************/
/* Internal function.
 * Return a pointer to the last parsed byte.
 * Assume that the bytecode is correct (for example, assume that the size of the
 * array R is >= 5).
 */
#define bytecode_prac_interpret_montgomery_internal MOD_APPEND_TYPE(bytecode_prac_interpret_montgomery_internal)
static bytecode_const
bytecode_prac_interpret_montgomery_internal (bytecode_const bc, ec_point_t *R,
                                             const modulus_t m,
                                             const residue_t b)
{
  while (1)
  {
    int finished = 0;
    switch (*bc)
    {
      case PRAC_SWAP: /* [ = 's' ] Swap R[0], R[1] */
        ec_point_swap (R[0], R[1], m, MONTGOMERY_xz);
        break;
      case PRAC_SUBBLOCK_INIT: /* [ = 'i' ] Start of a sub-block */
        ec_point_set (R[1], R[0], m, MONTGOMERY_xz);
        ec_point_set (R[2], R[0], m, MONTGOMERY_xz);
        montgomery_dbl (R[0], R[0], m, b);
        break;
      case PRAC_SUBBLOCK_FINAL: /* [ = 'f' ] End of a sub-block */
        montgomery_dadd (R[0], R[0], R[1], R[2], b, m);
        break;
      case PRAC_BLOCK_FINAL: /* [ = 'F' ] End of the block */
        montgomery_dadd (R[1], R[0], R[1], R[2], b, m);
        finished = 1;
        break;
      case 1:
        montgomery_dadd (R[3], R[0], R[1], R[2], b, m);
        montgomery_dadd (R[4], R[3], R[0], R[1], b, m);
        montgomery_dadd (R[1], R[1], R[3], R[0], b, m);
        ec_point_set (R[0], R[4], m, MONTGOMERY_xz);
        break;
      case 2:
        montgomery_dadd (R[1], R[0], R[1], R[2], b, m);
        montgomery_dbl (R[0], R[0], m, b);
        break;
      case 3:
        montgomery_dadd (R[2], R[1], R[0], R[2], b, m);
        ec_point_swap (R[1], R[2], m, MONTGOMERY_xz);
        break;
      case 4:
        montgomery_dadd (R[1], R[1], R[0], R[2], b, m);
        montgomery_dbl (R[0], R[0], m, b);
        break;
      case 5:
        montgomery_dadd (R[2], R[2], R[0], R[1], b, m);
        montgomery_dbl (R[0], R[0], m, b);
        break;
      case 6:
        montgomery_dbl (R[3], R[0], m, b);
        montgomery_dadd (R[4], R[0], R[1], R[2], b, m);
        montgomery_dadd (R[0], R[3], R[0], R[0], b, m);
        montgomery_dadd (R[2], R[3], R[4], R[2], b, m);
        ec_point_swap (R[1], R[2], m, MONTGOMERY_xz);
        break;
      case 7:
        montgomery_dadd (R[3], R[0], R[1], R[2], b, m);
        montgomery_dadd (R[1], R[3], R[0], R[1], b, m);
        montgomery_dbl (R[3], R[0], m, b);
        montgomery_dadd (R[0], R[0], R[3], R[0], b, m);
        break;
      case 8:
        montgomery_dadd (R[3], R[0], R[1], R[2], b, m);
        montgomery_dadd (R[2], R[2], R[0], R[1], b, m);
        ec_point_swap (R[1], R[3], m, MONTGOMERY_xz);
        montgomery_dbl (R[3], R[0], m, b);
        montgomery_dadd (R[0], R[0], R[3], R[0], b, m);
        break;
      case 9:
        montgomery_dadd (R[2], R[2], R[1], R[0], b, m);
        montgomery_dbl (R[1], R[1], m, b);
        break;
      case 10:
        /* Combined final add of old subchain and init of new subchain [=fi] */
        montgomery_dadd (R[1], R[0], R[1], R[2], b, m);
        ec_point_set (R[2], R[1], m, MONTGOMERY_xz);
        montgomery_dbl (R[0], R[1], m, b);
        break;
      case 11:
        /* Combined rule 3 and rule 0 [=\x3s] */
        montgomery_dadd (R[2], R[1], R[0], R[2], b, m);
        /* (R[1],R[2],R[0]) := (R[0],R[1],R[2])  */
        ec_point_swap (R[1], R[2], m, MONTGOMERY_xz);
        ec_point_swap (R[0], R[1], m, MONTGOMERY_xz);
        break;
      case 12:
        /* Combined rule 3, then subchain end/start [=\x3fi] */
        montgomery_dadd (R[3], R[1], R[0], R[2], b, m);
        montgomery_dadd (R[2], R[0], R[3], R[1], b, m);
        ec_point_set (R[1], R[2], m, MONTGOMERY_xz);
        montgomery_dbl (R[0], R[2], m, b);
        break;
      case 13:
        /* Combined rule 3, swap, rule 3 and swap, merged a bit [=\x3s\x3s] */
        ec_point_set (R[3], R[1], m, MONTGOMERY_xz);
        montgomery_dadd (R[1], R[1], R[0], R[2], b, m);
        ec_point_set (R[2], R[0], m, MONTGOMERY_xz);
        montgomery_dadd (R[0], R[0], R[1], R[3], b, m);
        break;
      default:
        printf ("Fatal error in %s at %s:%d -- unknown bytecode 0x%02x\n",
                __func__, __FILE__, __LINE__, *bc);
        abort ();
    }

    if (finished) /* is it finished ? */
      break;
    else
      bc++; /* go to next byte */
  }

  return bc;
}

#define bytecode_prac_interpret_montgomery MOD_APPEND_TYPE(bytecode_prac_interpret_montgomery)
static void
bytecode_prac_interpret_montgomery (ec_point_t P, bytecode_const bc,
                                    const modulus_t m, const residue_t b)
{
  ec_point_t *R = NULL;
  unsigned int R_nalloc;

  if (bc == NULL)
    return;

  R_nalloc = 5; /* we need 5 points: 3 for PRAC + 2 temporary points */
  R = (ec_point_t *) malloc (R_nalloc * sizeof (ec_point_t));
  FATAL_ERROR_CHECK (R == NULL, "could not malloc R");
  for (unsigned int i = 0; i < R_nalloc; i++)
    ec_point_init (R[i], m);

  /* current point (here starting point) go into R[0] at init */
  ec_point_set (R[0], P, m, MONTGOMERY_xz);

  bytecode_prac_interpret_montgomery_internal (bc, R, m, b);

  /* output is in R[1] */
  ec_point_set (P, R[1], m, MONTGOMERY_xz);

  for (unsigned int i = 0; i < R_nalloc; i++)
    ec_point_clear (R[i], m);
  free (R);
}

/** Interpret the MISHMASH bytecode on Twisted Edwards and Montgomery curves **/
#define bytecode_mishmash_interpret_mixed_repr MOD_APPEND_TYPE(bytecode_mishmash_interpret_mixed_repr)
static void
bytecode_mishmash_interpret_mixed_repr (ec_point_t P, bytecode_const bc,
                                        const modulus_t m, const residue_t b)
{
  ec_point_t *R = NULL;
  unsigned int R_nalloc;
  uint8_t n;

  /* parse first byte to alloc R and go to the next byte */
  bytecode_elt_split_4_4 (NULL, &n, *bc);
  bc++;

  R_nalloc = n + 2;
  R = (ec_point_t *) malloc (R_nalloc * sizeof (ec_point_t));
  FATAL_ERROR_CHECK (R == NULL, "could not malloc R");
  for (unsigned int i = 0; i < R_nalloc; i++)
    ec_point_init (R[i], m);

  /* starting point go into R[1] at init */
  ec_point_set (R[1], P, m, TWISTED_EDWARDS_ext);

  while (*bc != MISHMASH_FINAL)
  {
    uint8_t t, n;
    bytecode_elt_split_4_4 (&t, &n, *bc);

    if (n != 0)
      ec_point_set (R[0], R[n], m, TWISTED_EDWARDS_ext);

    if (t == MISHMASH_DBCHAIN_BLOCK)
      bc = bytecode_dbchain_interpret_edwards_internal (++bc, R, m);
    else if (t == MISHMASH_PRECOMP_BLOCK)
      bc = bytecode_precomp_interpret_edwards_internal (++bc, R, m);
    else if (t == MISHMASH_PRAC_BLOCK)
      bc = bytecode_prac_interpret_montgomery_internal (++bc, R, m, b);
    else /* unexpected bytecode */
    {
      printf ("Fatal error in %s at %s:%d -- unexpected bytecode 0x%02x\n",
              __func__, __FILE__, __LINE__, *bc);
      abort ();
    }

    bc++; /* go to next byte */
  }

  /* output (always in Montgomery form) is in R[1] */
  ec_point_set (P, R[1], m, MONTGOMERY_xz);

  for (unsigned int i = 0; i < R_nalloc; i++)
    ec_point_clear (R[i], m);
  free (R);
}


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/* Multiplies x[1] by z[2]*z[3]*z[4]...*z[n],
   x[2] by z[1]*z[3]*z[4]...*z[n] etc., generally
   x[i] by \prod_{1\leq j \leq n, j\neq i} z[j]
   Requires n > 1. Uses 4n-6 multiplications.
 * If n <= 1 there is nothing to do.
 */

#define common_z MOD_APPEND_TYPE(common_z)
static void ATTRIBUTE((__noinline__))
common_z (ec_point_t *L1, const unsigned int n1, ec_point_t *L2,
          const unsigned int n2, const modulus_t m)
{
  const unsigned int n = n1 + n2;
  residue_t *t = NULL, p;

  if (n <= 1) /* nothing to do in this case */
    return;

#ifdef ECM_COUNT_OPS
  _count_stage2_common_z += 4*n - 6;
#endif

  mod_init (p, m);
  t = (residue_t *) malloc ((n-1) * sizeof (residue_t));
  ASSERT_ALWAYS (t != NULL);
  for (unsigned int i = 0; i < n-1; i++)
    mod_init (t[i], m);

  /* Set t[i] = z_0 * z_1 * ... * z_i, where the z_i are the Z-coordinates taken
   * from the two lists L1 and L2.
   */
  if (n1)
    mod_set (t[0], L1[0]->z, m);
  else
    mod_set (t[0], L2[0]->z, m); /* L1 is empty */
  for (unsigned int i = 1; i < n-1; i++)
  {
    if (i < n1)
      mod_mul (t[i], t[i-1], L1[i]->z, m);
    else
      mod_mul (t[i], t[i-1], L2[i-n1]->z, m);
  }
  /* cost: n-2 mul */

  if (n2)
  {
    mod_mul (L2[n2-1]->x, L2[n2-1]->x, t[n-2], m);
    mod_set (p, L2[n2-1]->z, m);
  }
  else
  {
    mod_mul (L1[n1-1]->x, L1[n1-1]->x, t[n-2], m);
    mod_set (p, L1[n1-1]->z, m);
  }
  /* cost: 1 mul */

  for (unsigned int i = n-2; i > 0; i--)
  {
    if (i < n1)
    {
      mod_mul (L1[i]->x, L1[i]->x, p, m);
      mod_mul (L1[i]->x, L1[i]->x, t[i-1], m);
      mod_mul (p, p, L1[i]->z, m);
    }
    else
    {
      mod_mul (L2[i-n1]->x, L2[i-n1]->x, p, m);
      mod_mul (L2[i-n1]->x, L2[i-n1]->x, t[i-1], m);
      mod_mul (p, p, L2[i-n1]->z, m);
    }
  }
  /* cost: 3*n-6 mul */

  if (n1)
    mod_mul (L1[0]->x, L1[0]->x, p, m);
  else
    mod_mul (L2[0]->x, L2[0]->x, p, m);
  /* cost: 1 mul */

  /* free memory */
  mod_clear (p, m);
  for (unsigned int i = 0; i < n-1; i++)
    mod_clear (t[i], m);
  free (t);
}


#define ecm_stage2 MOD_APPEND_TYPE(ecm_stage_2)
static int ATTRIBUTE((__noinline__))
ecm_stage2 (residue_t r, const ec_point_t P, const stage2_plan_t *plan,
            const residue_t b, const modulus_t m)
{
  ASSERT (plan->w % 6 == 0); /* see stage2_make_plan */

  ec_point_t wP, Pt; /* w*P and a temp */
  ec_point_t *vwP, *uP; /* saved v*w*P, vmin <= v <= vmax, and uP, u in U. */
  int bt = 0;

  ec_point_init (wP, m);
  ec_point_init (Pt, m);

  ASSERT (plan->vmin <= plan->vmax);
  vwP = (ec_point_t *) malloc ((plan->vmax-plan->vmin+1) * sizeof (ec_point_t));
  ASSERT (vwP != NULL);
  uP = (ec_point_t *) malloc (plan->U_len * sizeof (ec_point_t));
  ASSERT (uP != NULL);

  for (unsigned int i = 0; i < plan->vmax - plan->vmin + 1; i++)
    ec_point_init_noset0 (vwP[i], m);
  for (unsigned int i = 0; i < plan->U_len; i++)
    ec_point_init_noset0 (uP[i], m);

#ifdef ECM_STAGE2_DEBUG
  modint_t _int;
  mod_intinit (_int);
  printf ("# %s: w=%u\n", __func__, plan->w);
  printf ("# %s: input P=", __func__);
  ec_point_fprintf (stdout, P, MONTGOMERY_xz, m);
  fputc ('\n', stdout);
#endif

  { /***************************** baby step **********************************/
    /* Compute u*P for u in U. Compute all the u, 1 <= u < d/2, gcd(u,w)=1
     * with two arithmetic progressions 1+6k and 5+6k (this assumes 6|w).
     * We need two values of each progression (1, 7 and 5, 11) and the
     * common difference 6. These can be computed with the Lucas chain
     * 1, 2, 3, 5, 6, 7, 11 at the cost of 4 additions and 2 doublings.
     * If w=30, we could use 1,2,3,4,6,7,11,13 which has 4 additions and 3
     * doublings.
     */
    ec_point_t ap1_0, ap1_1, ap5_0, ap5_1, P2, P6;
    ec_point_init_noset0 (ap1_0, m);
    ec_point_init_noset0 (ap1_1, m);
    ec_point_init_noset0 (ap5_0, m);
    ec_point_init_noset0 (ap5_1, m);
    ec_point_init_noset0 (P6, m);
    ec_point_init_noset0 (P2, m);

    /* Init ap1_0 = 1P, ap1_1 = 7P, ap5_0 = 5P, ap5_1 = 11P and P6 = 6P */
    ec_point_set (ap1_0, P, m, MONTGOMERY_xz);    /* ap1_0 = 1*P */
    montgomery_dbl (P2, P, m, b);                 /* P2 = 2*P */
    montgomery_dadd (P6, P2, P, P, b, m);         /* P6 = 3*P (for now) */
    montgomery_dadd (ap5_0, P6, P2, P, b, m);     /* 5*P = 3*P + 2*P */
    montgomery_dbl (P6, P6, m, b);                /* P6 = 6*P = 2*(3*P) */
    montgomery_dadd (ap1_1, P6, P, ap5_0, b, m);  /* 7*P = 6*P + P */
    montgomery_dadd (ap5_1, P6, ap5_0, P, b, m);  /* 11*P = 6*P + 5*P */

#ifdef ECM_STAGE2_DEBUG
    printf ("# %s: 1*P=", __func__);
    ec_point_fprintf (stdout, ap1_0, MONTGOMERY_xz, m);
    printf ("\n# %s: 5*P=", __func__);
    ec_point_fprintf (stdout, ap5_0, MONTGOMERY_xz, m);
    printf ("\n# %s: 7*P=", __func__);
    ec_point_fprintf (stdout, ap1_1, MONTGOMERY_xz, m);
    printf ("\n# %s: 11*P=", __func__);
    ec_point_fprintf (stdout, ap5_1, MONTGOMERY_xz, m);
    fputc ('\n', stdout);
#endif

    /* Now we generate all the u*P for u in U */
    /* We treat the first two manually because those might correspond to
     * ap1_0 = 1*P and ap5_0 = 5*P */
    unsigned int k = 0;
    if (k < plan->U_len && plan->U[k] == 1)
    {
      ec_point_set (uP[k], ap1_0, m, MONTGOMERY_xz);
      k++;
    }
    if (k < plan->U_len && plan->U[k] == 5)
    {
      ec_point_set (uP[k], ap5_0, m, MONTGOMERY_xz);
      k++;
    }

    unsigned int u_1mod6 = 7, u_5mod6 = 11;
    while (k < plan->U_len)
    {
      if (plan->U[k] == u_1mod6)
      {
        ec_point_set (uP[k], ap1_1, m, MONTGOMERY_xz);
        k++;
        continue;
      }
      if (plan->U[k] == u_5mod6)
      {
        ec_point_set (uP[k], ap5_1, m, MONTGOMERY_xz);
        k++;
        continue;
      }

      montgomery_dadd (Pt, ap1_1, P6, ap1_0, b, m);
      ec_point_set (ap1_0, ap1_1, m, MONTGOMERY_xz);
      ec_point_set (ap1_1, Pt, m, MONTGOMERY_xz);
      u_1mod6 += 6;

      montgomery_dadd (Pt, ap5_1, P6, ap5_0, b, m);
      ec_point_set (ap5_0, ap5_1, m, MONTGOMERY_xz);
      ec_point_set (ap5_1, Pt, m, MONTGOMERY_xz);
      u_5mod6 += 6;
    }

    /* Compute compute wP = w*P */
    if (plan->w == 6)
      ec_point_set (wP, P6, m, MONTGOMERY_xz);
    else if (plan->w == 12)
      montgomery_dbl (wP, P6, m, b);
    else if (plan->w % 12 == 0)
    {
      montgomery_dadd (Pt, ap1_1, P6, ap1_0, b, m);
      montgomery_dadd (wP, Pt, ap5_1, P2, b, m);
    }
    else /* plan->w % 12 == 6 */
    {
      montgomery_dbl (P2, P2, m, b); /* We need 4P for difference */
      montgomery_dadd (wP, ap1_1, ap5_1, P2, b, m);
    }

#ifdef ECM_STAGE2_DEBUG
    for (unsigned int k = 0; k < plan->U_len; k++)
    {
      printf ("%s# %s: %u*P=", k ? "\n" : "", __func__, plan->U[k]);
      ec_point_fprintf (stdout, uP[k], MONTGOMERY_xz, m);
    }
    printf ("\n# %s: %u*P=", __func__, plan->w);
    ec_point_fprintf (stdout, wP, MONTGOMERY_xz, m);
    fputc ('\n', stdout);
#endif

    ec_point_clear (ap1_0, m);
    ec_point_clear (ap1_1, m);
    ec_point_clear (ap5_0, m);
    ec_point_clear (ap5_1, m);
    ec_point_clear (P6, m);
    ec_point_clear (P2, m);
  } /* end baby step */

  { /**************************** giant step **********************************/
    /* Compute vwP for vmin <= v < vmax */
    if (plan->vmin == plan->vmax) /* If vmin == vmax, only vmin*w*P is needed */
      montgomery_smul_ul (vwP[0], NULL, wP, plan->vmin, m, b);
    else
    {
      montgomery_smul_ul (vwP[0], vwP[1], wP, plan->vmin, m, b);

      unsigned int k = 2, v = plan->vmin+2;

      for ( ; v <= plan->vmax; k++, v++)
      {
        if (v % 2 == 0 && v/2 >= plan->vmin)
          montgomery_dbl (vwP[k], vwP[v/2-plan->vmin], m, b);
        else
          montgomery_dadd (vwP[k], vwP[k-1], wP, vwP[k-2], b, m);
      }
    }

#ifdef ECM_STAGE2_DEBUG
    for (unsigned int k = 0; k < plan->vmax - plan->vmin + 1; k++)
    {
      printf ("%s# %s: %u*w*P=", k ? "\n" : "", __func__, plan->vmin+k);
      ec_point_fprintf (stdout, vwP[k], MONTGOMERY_xz, m);
    }
    fputc ('\n', stdout);
#endif
  } /* end giant step */


  { /***************************** common_z ***********************************/
    /* Now we've computed all the points we need, so multiply each by the
     * Z-coordinates of all the others, using Zimmermann's two product-lists
     * trick. If vmin == 0, then vwP[0] is the point at infinity (0::0), so we
     * skip that one
     */
    int skip = (plan->vmin == 0) ? 1 : 0;
    common_z (uP, plan->U_len, vwP+skip, plan->vmax-plan->vmin+1-skip, m);

#ifdef ECM_STAGE2_DEBUG
    for (unsigned int k = 0; k < plan->U_len; k++)
    {
      mod_get_int (_int, uP[k]->x, m);
      mod_printf ("# %s: after common_z, %u*P=(0x%" PRIMODx " :: Z)\n",
                  __func__, plan->U[k], MOD_PRINT_INT (_int));
    }
    for (unsigned int k = 0; k < plan->vmax - plan->vmin + 1; k++)
    {
      mod_get_int (_int, vwP[k]->x, m);
      mod_printf ("# %s: after common_z, %u*w*P=(0x%" PRIMODx " :: Z)\n",
                  __func__, plan->vmin+k, MOD_PRINT_INT (_int));
    }
#endif
  } /* end of common_z */

  { /***************************** product ************************************/
    /* Now compute
     *        r = prod_{u,v in plan->pairs}{X_vwP - X_uP}
     *
     * Initialize r with uP[0], which contains the product of the Z-coordinates
     * of all the precomputed points, except the Z-coordinate of uP[0]. Note
     * that uP[0] is equal to P, and we know that the Z-coordinate of P is
     * coprime to the modulus m. It cost nothing to add this in the product and
     * we could catch a factor p if, by chance, one of the Z-coordinates is 0
     * modulo p.
     */
    residue_t t, r_bak; /* Backup value of a, in case we get r == 0 */
    mod_init_noset0 (t, m);
    mod_init_noset0 (r_bak, m);

    mod_set (r, uP[0]->x, m);
    mod_set (r_bak, r, m);

    for (unsigned int v_idx = 0, *u_ptr = plan->pairs; ; v_idx++, u_ptr++)
    {
      ASSERT (v_idx <= plan->vmax - plan->vmin);
#if ECM_BACKTRACKING
      mod_set (r_bak, r, m);
#endif

      for (; *u_ptr != PAIR_END && *u_ptr != PAIR_INCR_V; u_ptr++)
      {
        mod_sub (t, vwP[v_idx]->x, uP[*u_ptr]->x, m);
        mod_mul (r, r, t, m);
#ifdef ECM_COUNT_OPS
        _count_stage2_product++;
#endif
      }

#if ECM_BACKTRACKING
      /* See if we got r == 0. If yes, restore previous r value and end stage 2.
       * Let's hope not all factors were found since the last v increase.
       */
      if (mod_is0 (r, m))
      {
        mod_set (r, r_bak, m);
        bt = 1;
        break;
      }
#endif

      if (*u_ptr == PAIR_END)
        break;
    }

#ifdef ECM_STAGE2_DEBUG
    mod_get_int (_int, r, m);
    mod_printf ("# %s: r=%" PRIMODu "\n", __func__, MOD_PRINT_INT (_int));
    mod_intclear (_int);
#endif

    mod_clear (t, m);
    mod_clear (r_bak, m);
  }

  /* Clear everything */
  for (unsigned int i = 0; i < plan->U_len; i++)
    ec_point_clear (uP[i], m);
  for (unsigned int i = 0; i < plan->vmax - plan->vmin + 1; i++)
    ec_point_clear (vwP[i], m);

  free (uP);
  free (vwP);

  ec_point_clear (wP, m);
  ec_point_clear (Pt, m);
  return bt;
}

/* Stores any factor found in f_out (1 if no factor found).
   If back-tracking was used, returns 1, otherwise returns 0. */

int
ecm (modint_t f, const modulus_t m, const ecm_plan_t *plan)
{
#ifdef ECM_DEBUG
    mod_printf ("# %s: start with B1=%u and B2=%u for m = %" PRIMODu "\n",
                __func__, plan->B1, plan->stage2.B2, MOD_PRINT_MODULUS (m));
    printf ("# %s: using parameterization 0x%x with parameter %lu\n",
            __func__, plan->parameterization, plan->parameter);
#endif

  residue_t u, b;
  ec_point_t P, Pt;
  ec_point_coord_type_t param_output_type;

  unsigned int i;
  int r, bt = 0;

  ec_point_init (P, m);
  mod_init_noset0 (b, m);

  mod_intset_ul (f, 1UL);

  /**************************** parameterization ******************************/
  if (plan->parameterization & FULLMONTY)
  {
    param_output_type = MONTGOMERY_xz;
    if (plan->parameterization == BRENT12)
      r = ec_parameterization_Brent_Suyama (b, P, plan->parameter, m);
    else if (plan->parameterization == MONTY12)
      r = ec_parameterization_Montgomery12 (b, P, plan->parameter, m);
    else /* if (plan->parameterization == MONTY16) */
      r = ec_parameterization_Montgomery16 (b, P, plan->parameter, m);
  }
  else if (plan->parameterization == MONTYTWED12)
  {
    param_output_type = plan->bc[1] & 0x80 ? MONTGOMERY_xz: TWISTED_EDWARDS_ext;
    r = ec_parameterization_Z6 (b, P, plan->parameter, param_output_type, m);
  }
  else
  {
    fprintf (stderr, "%s: unknown parameterization\n", __func__);
    abort();
  }

  if (r == 0)
  {
    mod_gcd (f, P->x, m);
#ifdef ECM_DEBUG
    mod_printf ("# %s: during parameterization, found factor %" PRIMODu "\n",
                MOD_PRINT_INT (f));
#endif
    mod_clear (b, m);
    ec_point_clear (P, m);
    return 0;
  }

#ifdef ECM_DEBUG
  residue_t A;
  ec_point_t PM;
  ec_point_init (PM, m);
  mod_init_noset0 (A, m);
  montgomery_A_from_b (A, b, m);
  printf ("# %s: starting values:\n", __func__);

  if (param_output_type == MONTGOMERY_xz)
  {
    montgomery_curve_fprintf (stdout, "# " CADO_STRINGIZE(ecm) ":   ", A, P, m);
    printf ("# %s:                     = ", __func__);
    montgomery_point_fprintf_affine (stdout, P, m);
    fputc ('\n', stdout);
  }
  else
  {
    residue_t d;
    mod_init_noset0 (d, m);
    edwards_d_from_montgomery_A (d, A, m);
    edwards_ext_curve_fprintf (stdout, "# " CADO_STRINGIZE(ecm) ":   ", d, P, m);
    mod_clear (d, m);

    printf ("# %s:   Equivalent to Montgomery curve with:\n", __func__);
    montgomery_point_from_edwards_point (PM, P, 1, m);
    montgomery_curve_fprintf (stdout, "# " CADO_STRINGIZE(ecm) ":     ", A, PM, m);
  }
  mod_clear (A, m);
  ec_point_clear (PM, m);
#endif

  /******************************** stage 1 ***********************************/
#ifdef ECM_COUNT_OPS
  ECM_COUNT_OPS_RESET();
#endif

  /* output is always in Montgomery form */
  if (plan->parameterization & FULLMONTY)
    bytecode_prac_interpret_montgomery (P, plan->bc, m, b);
  else if (plan->parameterization & FULLMONTYTWED)
    bytecode_mishmash_interpret_mixed_repr (P, plan->bc, m, b);

#ifdef ECM_DEBUG
  printf ("# %s: after stage 1 (without the powers of 2):\n", __func__);
  printf ("# %s:   output (X::Z) = ", __func__);
  ec_point_fprintf (stdout, P, MONTGOMERY_xz, m);
  printf ("\n# %s:                 = ", __func__);
  montgomery_point_fprintf_affine (stdout, P, m);
  fputc ('\n', stdout);
#endif

  /* Add prime 2 in the desired power. If a zero residue for the Z-coordinate is
   * encountered, we backtrack to previous point and stop.
   * NOTE: This is not as effective as I hoped. It prevents trivial
   * factorizations only if after processing the odd part of the stage 1
   * multiplier, the resulting point has power-of-2 order on E_p for all p|N.
   * If that were to happen, the point probably had that (presumably small on
   * most E_p) power-of-2 order during the last couple of primes processed in
   * the precomputed Lucas chain, and then quite likely the Lucas chain
   * incorrectly used an addition of identical points, causing the Z-coordinate
   * to become zero, leading to 0 (mod N) before we even get here.  For example,
   * using 10^6 composites from an RSA155 sieving experiment, without
   * backtracking we get N as the factor 456 times, with backtracking still 360
   * times.
   * TODO: this could probably be fixed by treating 3 separately, too, instead
   * of putting it in the precomputed Lucas chain. Then the probability that a
   * point of very small order on all E_p is encountered during the Lucas chain
   * is reduced, and so the probability of using curve addition erroneously.
   *
   * The following code assume P is in Montgomery form.
   */
  ec_point_init (Pt, m);
  ec_point_set (Pt, P, m, MONTGOMERY_xz);
  for (i = 0; i < plan->exp2; i++)
  {
    montgomery_dbl (P, P, m, b);
#if ECM_BACKTRACKING
    if (mod_is0 (P->z, m))
    {
      ec_point_set (P, Pt, m, MONTGOMERY_xz);
      bt = 1;
      break;
    }
    ec_point_set (Pt, P, m, MONTGOMERY_xz);
#endif
  }
  mod_gcd (f, P->z, m);

#ifdef ECM_DEBUG
  printf ("# %s: after stage 1 and %u power(s) of 2 (out of %u):\n", __func__,
          i, plan->exp2);
  printf ("# %s:   output (X::Z) = ", __func__);
  ec_point_fprintf (stdout, P, MONTGOMERY_xz, m);
  printf ("\n# %s:                 = ", __func__);
  montgomery_point_fprintf_affine (stdout, P, m);
#if ECM_BACKTRACKING
  mod_printf ("\n# %s:   backtracking was%s used", __func__, bt ? "" :" not");
#endif
  mod_printf ("\n# %s:   gcd = %" PRIMODu "\n", __func__, MOD_PRINT_INT (f));
#endif

#ifdef ECM_COUNT_OPS
  unsigned int tot_stage1_M = ECM_COUNT_OPS_STAGE1_TOTAL_M;
  printf ("# %s: count ops: for stage 1 with B1 = %u\n"
          "# %s: count ops:         Edwards: %u DBL %u DBLext %u TPL %u TPLext "
                                            "%u ADD %u ADDext %u ADDmont\n"
          "# %s: count ops:      Montgomery: %u DBL %u dADD\n"
          "# %s: count ops:   stage 1 total: %u M\n", __func__, plan->B1,
          __func__, _count_edwards_dbl, _count_edwards_dblext,
          _count_edwards_tpl, _count_edwards_tplext, _count_edwards_add,
          _count_edwards_addext, _count_edwards_addmont, __func__,
          _count_montgomery_dbl, _count_montgomery_dadd, __func__,
          tot_stage1_M);
#endif

  /******************************** stage 2 ***********************************/
#ifdef ECM_COUNT_OPS
  ECM_COUNT_OPS_RESET();
#endif

  mod_init (u, m);
  if (bt == 0 && mod_intcmp_ul (f, 1UL) == 0 && plan->B1 < plan->stage2.B2)
  {
    bt = ecm_stage2 (u, P, &(plan->stage2), b, m);
    mod_gcd (f, u, m);
#ifdef ECM_DEBUG
    mod_printf ("# %s: after stage 2, gcd=%" PRIMODu "\n", __func__,
                                                           MOD_PRINT_INT(f));
#endif
  }
#ifdef ECM_DEBUG
  else
    printf ("# %s: stage 2 not done\n", __func__);
#endif

#ifdef ECM_COUNT_OPS
  unsigned int tot_stage2_M = ECM_COUNT_OPS_STAGE2_TOTAL_M;
  printf ("# %s: count ops: for stage 2 with B2 = %u [ with w = %u ]\n"
          "# %s: count ops:      Montgomery: %u DBL %u dADD\n"
          "# %s: count ops:        common_z: %u M\n"
          "# %s: count ops:         product: %u M\n"
          "# %s: count ops:   stage 2 total: %u M\n"
          "# %s: count ops: ECM total: %5u M\n", __func__,
          plan->stage2.B2, plan->stage2.w, __func__, _count_montgomery_dbl,
          _count_montgomery_dadd, __func__, _count_stage2_common_z, __func__,
          _count_stage2_product, __func__, tot_stage2_M, __func__,
          tot_stage1_M + tot_stage2_M);
#endif

  mod_clear (u, m);
  mod_clear (b, m);
  ec_point_clear (P, m);
  ec_point_clear (Pt, m);

  return bt;
}
