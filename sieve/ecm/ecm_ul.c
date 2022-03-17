#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h> // abort
#include "modredc_ul.h"
#include "modredc_ul_default.h"
#include "ecm.c"
#include "facul_ecm.h"             // for BRENT12, MONTY12, MONTY16, MONTYTW...
#include "ec_parameterization.h"   // ec_parameterization_Brent_Suyama_ul
#include "ec_arith_common.h"       // ec_point_clear
#include "ec_arith_Weierstrass.h"  // weierstrass_aff_from_montgomery_ul


/* Exported functions (specific to _ul arithmetic). */
#define ec_parameterization_point_order MOD_APPEND_TYPE(ec_parameterization_point_order)
#define ec_parameterization_curve_order MOD_APPEND_TYPE(ec_parameterization_curve_order)

/* This function do not depend on modredc arithmetic, we define it only here
 * (not in ecm_15ul.c, ecm_2ul2.c nor ecm_mpz.c) and we do not append _ul to its
 * name */
int
ec_parameter_is_valid (ec_parameterization_t parameterization,
                       const unsigned long parameter)
{
  switch (parameterization)
  {
    case BRENT12:
      return ec_parameterization_Brent_Suyama_is_valid (parameter);
    case MONTY12:
      return ec_parameterization_Montgomery12_is_valid (parameter);
    case MONTY16:
      return ec_parameterization_Montgomery16_is_valid (parameter);
    case MONTYTWED12:
      return ec_parameterization_Z6_is_valid (parameter);
    // TODO MONTYTWED16 ???
    default:
      printf ("Fatal error in %s at %s:%d -- unknown parameterization %u\n",
              __func__, __FILE__, __LINE__, parameterization);
      abort ();
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
ec_valid_parameter_from_sequence (ec_parameterization_t parameterization,
                       const unsigned long sequence_value)
{
  switch (parameterization)
  {
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
      printf ("Fatal error in %s at %s:%d -- unknown parameterization %u\n",
              __func__, __FILE__, __LINE__, parameterization);
      abort ();
  }
}

/* The variables "parameterization" and "parameter" are used to computed a curve
 * E and a point P over the finite field of size m.
 * This function returns the order of the point P.
 * If the curve order is known to be == r (mod m), this can be supplied in
 * the variables "known_r" and" known_m".
 * Return 0 if a inversion failed during the computation (which indicates that
 * the curve is not defined modulo m).
 * Only the _ul version of this function is defined, because it depends on
 * weierstrass_aff_point_order.
 */
unsigned long
ec_parameterization_point_order (const ec_parameterization_t parameterization,
                                 const unsigned long parameter,
                                 const unsigned long known_m,
                                 const unsigned long known_r, const modulus_t m,
                                 const int verbose)
{
  residue_t a;
  unsigned long order = 0UL;
  int r = 1;
  ec_point_t P;

  mod_init (a, m);
  ec_point_init (P, m);

  if (parameterization & FULLMONTY || parameterization == MONTYTWED12)
  {
    residue_t A, b;

    mod_init (b, m);
    mod_init (A, m);

    if (parameterization == BRENT12)
      r = ec_parameterization_Brent_Suyama (b, P, parameter, m);
    else if (parameterization == MONTY12)
      r = ec_parameterization_Montgomery12 (b, P, parameter, m);
    else if (parameterization == MONTY16)
      r = ec_parameterization_Montgomery16 (b, P, parameter, m);
    else /* if (parameterization == MONTYTWED12) */
      r = ec_parameterization_Z6 (b, P, parameter, MONTGOMERY_xz, m);

    if (r)
    {
      montgomery_A_from_b (A, b, m);

      if (verbose >= 2)
      {
        printf ("%s: ", __func__);
        montgomery_curve_fprintf (stdout, NULL, A, P, m);
      }

      r = weierstrass_aff_from_montgomery (a, P, A, P, m);
    }

    mod_clear (b, m);
    mod_clear (A, m);
  }
  else
  {
    fprintf (stderr, "%s: Unknown parameterization\n", __func__);
    abort();
  }

  if (r)
    order = weierstrass_aff_point_order (a, P, known_m, known_r, m, verbose);

  mod_clear (a, m);
  ec_point_clear (P, m);
  return order;
}

/* Return the order of the curve define by a parameterization and a parameter.
 * Return 0 if a inversion failed during the computation (which indicates that
 * the curve is not defined modulo m).
 * Only the _ul version of this function is defined, because it depends on
 * montgomery_curve_order.
 */
unsigned long
ec_parameterization_curve_order (const ec_parameterization_t parameterization,
                                 const unsigned long parameter,
                                 const modulus_t m)
{
  residue_t A;
  unsigned long order = 0UL;
  int r = 1;
  ec_point_t P;

  mod_init_noset0 (A, m);
  ec_point_init (P, m);

  if (parameterization & FULLMONTY || parameterization == MONTYTWED12)
  {
    residue_t b;

    mod_init_noset0 (b, m);

    if (parameterization == BRENT12)
      r = ec_parameterization_Brent_Suyama (b, P, parameter, m);
    else if (parameterization == MONTY12)
      r = ec_parameterization_Montgomery12 (b, P, parameter, m);
    else if (parameterization == MONTY16)
      r = ec_parameterization_Montgomery16 (b, P, parameter, m);
    else /* if (parameterization == MONTYTWED12) */
      r = ec_parameterization_Z6 (b, P, parameter, MONTGOMERY_xz, m);

    if (r)
      montgomery_A_from_b (A, b, m);
    mod_clear (b, m);
  }
  else
  {
    fprintf (stderr, "%s: Unknown parameterization\n", __func__);
    abort();
  }

  if (r)
    order = montgomery_curve_order (A, P, m);

  mod_clear (A, m);
  ec_point_clear (P, m);
  return order;
}
