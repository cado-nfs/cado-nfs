#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include "pp1.h"
#include "lucas_V_mod.h"
#include "macros.h"

//#define PP1_STAGE2_DEBUG /* define to print debug information for stage 2 */

/* Do we want backtracking when processing factors of 2 in E? */
#ifndef PP1_BACKTRACKING
/* Default is "yes." Set to 0 for "no." */
#define PP1_BACKTRACKING 1
#endif

/* Interpret the PRAC bytecode for P+1 */
static void
pp1_stage1 (residue_t X, bytecode_const bc, const residue_t two,
            const modulus_t m)
{
  residue_t *R = NULL;
  unsigned int R_nalloc;

  R_nalloc = 5; /* we need 5 points: 3 for PRAC + 2 temporary points */
  R = (residue_t *) malloc (R_nalloc * sizeof (residue_t));
  FATAL_ERROR_CHECK (R == NULL, "could not malloc R");
  for (unsigned int i = 0; i < R_nalloc; i++)
    mod_init (R[i], m);

  /* current point (here starting point) go into R[0] at init */
  mod_set (R[0], X, m);

  while (1)
  {
    int finished = 0;
    switch (*bc)
    {
      case PRAC_SWAP: /* [ = 's' ] Swap R[0], R[1] */
        mod_swap (R[0], R[1], m);
        break;
      case PRAC_SUBBLOCK_INIT: /* [ = 'i' ] Start of a sub-block */
        mod_set (R[1], R[0], m);
        mod_set (R[2], R[0], m);
        mod_V_dbl (R[0], R[0], two, m);
        break;
      case PRAC_SUBBLOCK_FINAL: /* [ = 'f' ] End of a sub-block */
        mod_V_dadd (R[0], R[0], R[1], R[2], m);
        break;
      case PRAC_BLOCK_FINAL: /* [ = 'F' ] End of the block */
        mod_V_dadd (R[1], R[0], R[1], R[2], m);
        finished = 1;
        break;
      case 1:
        mod_V_dadd (R[3], R[0], R[1], R[2], m);
        mod_V_dadd (R[4], R[3], R[0], R[1], m);
        mod_V_dadd (R[1], R[1], R[3], R[0], m);
        mod_set (R[0], R[4], m);
        break;
      case 2:
        mod_V_dadd (R[1], R[0], R[1], R[2], m);
        mod_V_dbl (R[0], R[0], two, m);
        break;
      case 3:
        mod_V_dadd (R[3], R[1], R[0], R[2], m);
        mod_set (R[2], R[1], m);
        mod_set (R[1], R[3], m);
        break;
      case 4:
        mod_V_dadd (R[1], R[1], R[0], R[2], m);
        mod_V_dbl (R[0], R[0], two, m);
        break;
      case 5:
        mod_V_dadd (R[2], R[2], R[0], R[1], m);
        mod_V_dbl (R[0], R[0], two, m);
        break;
      case 6:
        mod_V_dbl (R[3], R[0], two, m);
        mod_V_dadd (R[4], R[0], R[1], R[2], m);
        mod_V_dadd (R[4], R[3], R[4], R[2], m);
        mod_set (R[2], R[4], m);
        mod_V_dadd (R[4], R[3], R[0], R[0], m);
        mod_set (R[0], R[4], m);
        mod_swap (R[1], R[2], m);
        break;
      case 7:
        mod_V_dadd (R[3], R[0], R[1], R[2], m);
        mod_V_dadd (R[4], R[3], R[0], R[1], m);
        mod_set (R[1], R[4], m);
        mod_V_dbl (R[3], R[0], two, m);
        mod_set (R[4], R[0], m);
        mod_V_dadd (R[0], R[0], R[3], R[4], m);
        break;
      case 8:
        mod_V_dadd (R[3], R[0], R[1], R[2], m);
        mod_V_dadd (R[2], R[2], R[0], R[1], m);
        mod_swap (R[1], R[3], m);
        mod_V_dbl (R[3], R[0], two, m);
        mod_V_dadd (R[4], R[0], R[3], R[0], m);
        mod_set (R[0], R[4], m);
        break;
      case 9:
        mod_V_dadd (R[2], R[2], R[1], R[0], m);
        mod_V_dbl (R[1], R[1], two, m);
        break;
      case 10:
        /* Combined final add of old subchain and init of new subchain [=fi] */
        mod_V_dadd (R[1], R[0], R[1], R[2], m);
        mod_set (R[2], R[1], m);
        mod_V_dbl (R[0], R[1], two, m);
        break;
      case 11:
        /* Combined rule 3 and rule 0 [=\x3s] */
        mod_set (R[3], R[0], m);
        mod_V_dadd (R[0], R[1], R[0], R[2], m);
        mod_set (R[2], R[1], m);
        mod_set (R[1], R[3], m);
        break;
      case 12:
        /* Combined rule 3, then subchain end/start [=\x3fi] */
        mod_V_dadd (R[3], R[1], R[0], R[2], m);
        mod_V_dadd (R[2], R[0], R[3], R[1], m);
        mod_set (R[1], R[2], m);
        mod_V_dbl (R[0], R[2], two, m);
        break;
      case 13:
        /* Combined rule 3, swap, rule 3 and swap, merged a bit [=\x3s\x3s] */
        mod_set (R[3], R[1], m);
        mod_V_dadd (R[1], R[1], R[0], R[2], m);
        mod_set (R[2], R[0], m);
        mod_V_dadd (R[0], R[0], R[1], R[3], m);
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

  mod_set (X, R[1], m);

  for (unsigned int i = 0; i < R_nalloc; i++)
    mod_clear (R[i], m);
  free (R);
}

int
pp1_stage2 (residue_t r, const residue_t X, const stage2_plan_t *plan,
            const residue_t two, const modulus_t m)
{
  ASSERT (plan->w % 6 == 0); /* see stage2_make_plan */

  residue_t Xw, t; /* V_w (X) and a temp */
  residue_t *Xvw, *Xu; /* saved V_{v*w} (X), vmin <= v <= vmax,
                        * and V_u (X), u in U */
  int bt = 0;

  mod_init_noset0 (Xw, m);
  mod_init_noset0 (t, m);

  ASSERT (plan->vmin <= plan->vmax);
  Xvw = (residue_t *) malloc ((plan->vmax-plan->vmin+1) * sizeof(residue_t));
  ASSERT (Xvw != NULL);
  Xu = malloc (plan->U_len * sizeof(residue_t));
  ASSERT (Xu != NULL);

  for (unsigned int i = 0; i < plan->vmax - plan->vmin + 1; i++)
    mod_init_noset0 (Xvw[i], m);
  for (unsigned int i = 0; i < plan->U_len; i++)
    mod_init_noset0 (Xu[i], m);

#ifdef PP1_STAGE2_DEBUG
  modint_t _int;
  mod_intinit (_int);
  mod_get_int (_int, X, m);
  printf ("# %s: w=%u\n", __func__, plan->w);
  mod_printf ("# %s: input X=%" PRIMODu "\n", __func__, MOD_PRINT_INT(_int));
#endif

  { /***************************** baby step **********************************/
    /* Compute V_u (X) for u in U. Compute all the u, 1 <= u < d/2, gcd(u,w)=1
     * with two arithmetic progressions 1+6k and 5+6k (this assumes 6|w).
     * We need two values of each progression (1, 7 and 5, 11) and the
     * common difference 6. These can be computed with the Lucas chain
     * 1, 2, 3, 5, 6, 7, 11 at the cost of 4 dadd (=1M) and 2 dbl (1S).
     * If w=30, we could use 1,2,3,4,6,7,11,13 which has 4 dadd and 3 dbl.
     */
    residue_t ap1_0, ap1_1, ap5_0, ap5_1, X2, X6;
    mod_init_noset0 (ap1_0, m);
    mod_init_noset0 (ap1_1, m);
    mod_init_noset0 (ap5_0, m);
    mod_init_noset0 (ap5_1, m);
    mod_init_noset0 (X6, m);
    mod_init_noset0 (X2, m);

    /* Init ap1_0 = V_1(X), ap1_1 = V_7(X), ap5_0 = V_5(X), ap5_1 = V_11(X)
     * and X6 = V_6(X)
     */
    mod_set (ap1_0, X, m);                /* ap1_0 = V_1(X) */
    mod_V_dbl (X2, X, two, m);            /* X2 = V_2(X) */
    mod_V_dadd (X6, X2, X, X, m);         /* V_3(X) = V_2(X)*V_1(X) - V_1(X) */
    mod_V_dadd (ap5_0, X6, X2, X, m);     /* V_5(X) = V_3(X)*V_2(X) - V_1(X) */
    mod_V_dbl (X6, X6, two, m);           /* V_6(X) = V_3(X)*V_3(X) - 2 */
    mod_V_dadd (ap1_1, X6, X, ap5_0, m);  /* V_7(X) = V_6(X)*V_1(X) - V_5(X) */
    mod_V_dadd (ap5_1, X6, ap5_0, X, m);  /* V_11(X) = V_6(X)*V_5(X) - V_1(X) */

#ifdef PP1_STAGE2_DEBUG
    mod_get_int (_int, ap1_0, m);
    mod_printf ("# %s: V_1 (X)=%" PRIMODu "\n", __func__, MOD_PRINT_INT(_int));
    mod_get_int (_int, ap5_0, m);
    mod_printf ("# %s: V_5 (X)=%" PRIMODu "\n", __func__, MOD_PRINT_INT(_int));
    mod_get_int (_int, ap1_1, m);
    mod_printf ("# %s: V_7 (X)=%" PRIMODu "\n", __func__, MOD_PRINT_INT(_int));
    mod_get_int (_int, ap5_1, m);
    mod_printf ("# %s: V_11 (X)=%" PRIMODu "\n", __func__, MOD_PRINT_INT(_int));
#endif

    /* Now we generate all the V_u (X) for u in U */
    /* We treat the first two manually because those might correspond to
     * ap1_0 = V_1 (X) and ap5_0 = V_5 (X) */
    unsigned int k = 0;
    if (plan->U_len > k && plan->U[k] == 1)
    {
      mod_set (Xu[k], ap1_0, m);
      k++;
    }
    if (plan->U_len > k && plan->U[k] == 5)
    {
      mod_set (Xu[k], ap5_0, m);
      k++;
    }

    unsigned int u_1mod6 = 7, u_5mod6 = 11;
    while (k < plan->U_len)
    {
      if (plan->U[k] == u_1mod6)
      {
        mod_set (Xu[k], ap1_1, m);
        k++;
        continue;
      }
      if (plan->U[k] == u_5mod6)
      {
        mod_set (Xu[k], ap5_1, m);
        k++;
        continue;
      }

      mod_V_dadd (t, ap1_1, X6, ap1_0, m);
      mod_set (ap1_0, ap1_1, m);
      mod_set (ap1_1, t, m);
      u_1mod6 += 6;

      mod_V_dadd (t, ap5_1, X6, ap5_0, m);
      mod_set (ap5_0, ap5_1, m);
      mod_set (ap5_1, t, m);
      u_5mod6 += 6;
    }

    /* Compute compute Xw = V_w (X) */
    if (plan->w == 6)
      mod_set (Xw, X6, m);
    else if (plan->w == 12)
      mod_V_dbl (Xw, X6, two, m);
    else if (plan->w % 12 == 0)
    {
      mod_V_dadd (t, ap1_1, X6, ap1_0, m);
      mod_V_dadd (Xw, t, ap5_1, X2, m);
    }
    else /* plan->w % 12 == 6 */
    {
      mod_V_dbl (X2, X2, two, m); /* We need V_4 (X) for difference */
      mod_V_dadd (Xw, ap1_1, ap5_1, X2, m);
    }

#ifdef PP1_STAGE2_DEBUG
    for (unsigned int k = 0; k < plan->U_len; k++)
    {
      mod_get_int (_int, Xu[k], m);
      mod_printf ("# %s: V_%u (X)=%" PRIMODu "\n", __func__, plan->U[k],
                                                           MOD_PRINT_INT(_int));
    }
    mod_get_int (_int, Xw, m);
    mod_printf ("# %s: V_%u (X)=%" PRIMODu "\n", __func__, plan->w,
                                                         MOD_PRINT_INT(_int));
#endif

    mod_clear (ap1_0, m);
    mod_clear (ap1_1, m);
    mod_clear (ap5_0, m);
    mod_clear (ap5_1, m);
    mod_clear (X6, m);
    mod_clear (X2, m);
  } /* end baby step */


  { /**************************** giant step **********************************/
    /* Compute V_{v*w} (X) for vmin <= v < vmax */
    if (plan->vmin == plan->vmax) /* If vmin == vmax, only V_{vmin*w} (X) is
                                   * needed */
      mod_V_eval_ul (Xvw[0], NULL, Xw, plan->vmin, m);
    else
    {
      mod_V_eval_ul (Xvw[0], Xvw[1], Xw, plan->vmin, m);

      unsigned int k = 2, v = plan->vmin+2;

      for ( ; v <= plan->vmax; k++, v++)
      {
        if (v % 2 == 0 && v/2 >= plan->vmin)
          mod_V_dbl (Xvw[k], Xvw[v/2-plan->vmin], two, m);
        else
          mod_V_dadd (Xvw[k], Xvw[k-1], Xw, Xvw[k-2], m);
      }
    }

#ifdef PP1_STAGE2_DEBUG
    for (unsigned int k = 0; k < plan->vmax - plan->vmin + 1; k++)
    {
      mod_get_int (_int, Xvw[k], m);
      mod_printf ("# %s: V_%uw (X)=%" PRIMODu "\n", __func__, plan->vmin+k,
                                                    MOD_PRINT_INT(_int));
    }
#endif
  } /* end giant step */

  { /***************************** product ************************************/
    /* Now compute
     *        r = prod_{u,v in plan->pairs}{V_{v*w} (X) - V_u (X)}
     */
    residue_t r_bak; /* Backup value of r, in case we get r == 0 */
    mod_init_noset0 (r_bak, m);

    mod_set1 (r, m);
    mod_set (r_bak, r, m);

    for (unsigned int v_idx = 0, *u_ptr = plan->pairs; ; v_idx++, u_ptr++)
    {
      ASSERT (v_idx <= plan->vmax - plan->vmin);
#if PP1_BACKTRACKING
      mod_set (r_bak, r, m);
#endif

      for (; *u_ptr != PAIR_END && *u_ptr != PAIR_INCR_V; u_ptr++)
      {
        mod_sub (t, Xvw[v_idx], Xu[*u_ptr], m);
        mod_mul (r, r, t, m);
      }

#if PP1_BACKTRACKING
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

#ifdef PP1_STAGE2_DEBUG
    mod_get_int (_int, r, m);
    mod_printf ("# %s: r=%" PRIMODu "\n", __func__, MOD_PRINT_INT(_int));
    mod_intclear (_int);
#endif

    mod_clear (r_bak, m);
  }

  /* Clear everything */
  for (unsigned int i = 0; i < plan->U_len; i++)
    mod_clear (Xu[i], m);
  for (unsigned int i = 0; i < plan->vmax - plan->vmin + 1; i++)
    mod_clear (Xvw[i], m);

  free (Xu);
  free (Xvw);

  mod_clear (Xw, m);
  mod_clear (t, m);
  return bt;
}

static inline int
pp1 (modint_t f, residue_t b, const residue_t two, const modulus_t m,
     const pp1_plan_t *plan)
{
  residue_t t;
  int bt = 0;

  mod_init_noset0 (t, m);

  /* stage1 */
  pp1_stage1 (b, plan->bc, two, m);

  /* Backtracking for the 2's in the exponent */
  mod_set (t, b, m);
  for (unsigned int i = 0; i < plan->exp2; i++)
  {
    mod_V_dbl (b, b, two, m);
#if PP1_BACKTRACKING
    if (mod_equal (b, two, m))
    {
      mod_set (b, t, m);
      bt = 1;
      break;
    }
    mod_set (t, b, m);
#endif
  }
  mod_sub (t, b, two, m);
  mod_gcd (f, t, m);

  /* stage 2 */
  if (bt == 0 && mod_intcmp_ul (f, 1UL) == 0 && plan->B1 < plan->stage2.B2)
  {
    bt = pp1_stage2 (t, b, &(plan->stage2), two, m);
    mod_gcd (f, t, m);
  }

  mod_clear (t, m);
  return bt;
}

int
pp1_27 (modint_t f, const modulus_t m, const pp1_plan_t *plan)
{
  residue_t b, two;

  mod_init_noset0 (b, m);
  mod_init_noset0 (two, m);
  mod_set1 (two, m);
  mod_add (two, two, two, m);

  /* Compute 2/7 (mod N) */
  mod_set (b, two, m);
  mod_div7 (b, b, m);

  int ret = pp1 (f, b, two, m, plan);

  mod_clear (b, m);
  mod_clear (two, m);
  return ret;
}

int
pp1_65 (modint_t f, const modulus_t m, const pp1_plan_t *plan)
{
  residue_t b, two;

  mod_init_noset0 (b, m);
  mod_init_noset0 (two, m);
  mod_set1 (two, m);
  mod_add (two, two, two, m);

  /* Compute 6/5 (mod N) */
  mod_set (b, two, m);
  mod_add (b, b, two, m);
  mod_add (b, b, two, m);
  mod_div5 (b, b, m);

  int ret = pp1 (f, b, two, m, plan);

  mod_clear (b, m);
  mod_clear (two, m);
  return ret;
}
