/**
   These files (mpz_vector.*) implement vectors whose coefficients are
   multiprecision integers (using mpz_t from GNU MP).
   We use them in twoquadractics and twocubics.
*/

#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mpz_poly.h"
#include "mpz_vector.h"
#include "gmp_aux.h"       // mpz_set_uint64
#include "macros.h" // ASSERT_ALWAYS


void
mpz_vector_init (mpz_vector_ptr v, unsigned int d)
{
  ASSERT_ALWAYS (d > 0);
  v->dim = d;
  v->c = (mpz_t *) malloc (d * sizeof(mpz_t));
  ASSERT_ALWAYS (v->c != NULL);
  for (unsigned int i = 0; i < d; i++)
    mpz_init (v->c[i]);
}

void
mpz_vector_clear(mpz_vector_ptr v)
{
  for (unsigned int i = 0; i < v->dim; i++)
    mpz_clear (v->c[i]);
  free (v->c);
}

void
mpz_vector_swap (mpz_vector_ptr v1, mpz_vector_ptr v2)
{
  ASSERT_ALWAYS (v1->dim == v2->dim);
  mpz_t *tmp = v1->c;
  v1->c = v2->c;
  v2->c = tmp;
}

void
mpz_vector_set (mpz_vector_ptr v, mpz_vector_srcptr s)
{
  ASSERT_ALWAYS (v->dim == s->dim);
  for (unsigned int i = 0; i < v->dim; i++)
    mpz_set (v->c[i], s->c[i]);
}

void
mpz_vector_setcoordinate (mpz_vector_ptr v, unsigned int i, mpz_srcptr z)
{
  ASSERT_ALWAYS (i < v->dim);
  mpz_set(v->c[i], z);
}

void
mpz_vector_setcoordinate_ui (mpz_vector_ptr v, unsigned int i, unsigned int z)
{
  ASSERT_ALWAYS (i < v->dim);
  mpz_set_ui (v->c[i], z);
}

void
mpz_vector_setcoordinate_si (mpz_vector_ptr v, unsigned int i, int z)
{
  ASSERT_ALWAYS (i < v->dim);
  mpz_set_si (v->c[i], z);
}

void
mpz_vector_setcoordinate_uint64 (mpz_vector_ptr v, unsigned int i, uint64_t z)
{
  ASSERT_ALWAYS (i < v->dim);
  mpz_set_uint64 (v->c[i], z);
}

void
mpz_vector_setcoordinate_int64 (mpz_vector_ptr v, unsigned int i, int64_t z)
{
  ASSERT_ALWAYS (i < v->dim);
  mpz_set_int64 (v->c[i], z);
}

int
mpz_vector_is_coordinate_zero (mpz_vector_srcptr v, unsigned int i)
{
  ASSERT_ALWAYS (i < v->dim);
  return (mpz_sgn (v->c[i]) == 0);
}

void
mpz_vector_get_mpz_poly (mpz_poly_ptr p, mpz_vector_srcptr v)
{
  for (unsigned int i = 0; i < v->dim; i++)
    mpz_poly_setcoeff (p, i, v->c[i]);
}

int mpz_vector_cmp (mpz_vector_srcptr a, mpz_vector_srcptr b)
{
  int r = (a->dim > b->dim) - (b->dim > a->dim);
  if (r) return r;
  for(int d = a->dim - 1; d >= 0 ; d--) {
    r = mpz_cmp(a->c[d], b->c[d]);
    if (r) return r;
  }
  return 0;
}

void mpz_vector_fprintf(FILE * file, mpz_vector_srcptr v)
{
  fprintf(file, "[");
  for (unsigned int i = 0; i < v->dim - 1; i++)
    gmp_fprintf(file, "%Zd, ", v->c[i]);
  gmp_fprintf(file, "%Zd]\n", v->c[v->dim - 1]);
}

void
mpz_vector_dot_product (mpz_ptr res, mpz_vector_srcptr a, mpz_vector_srcptr b)
{
  ASSERT_ALWAYS (a->dim == b->dim);
  mpz_mul(res, a->c[0], b->c[0]);
  for (unsigned int i = 1; i < a->dim; i++)
    mpz_addmul(res, a->c[i], b->c[i]);
}

void
mpz_vector_norm (mpz_ptr res, mpz_vector_srcptr a)
{
  mpz_vector_dot_product (res, a, a);
}

void
mpz_vector_skew_dot_product (mpz_ptr res, mpz_vector_srcptr a, mpz_vector_srcptr b,
                             mpz_srcptr skewness)
{
  ASSERT_ALWAYS (a->dim == b->dim);
  mpz_t tmp, s, s2;
  mpz_init (tmp);
  mpz_init (s2);
  mpz_init_set_ui (s, 1);
  mpz_mul (s2, skewness, skewness);

  mpz_mul(res, a->c[0], b->c[0]);
  for (unsigned int i = 1; i < a->dim; i++)
  {
    mpz_mul(tmp, a->c[i], b->c[i]);
    mpz_mul(s, s, s2);
    mpz_addmul(res, tmp, s);
  }

  mpz_clear(tmp);
  mpz_clear(s2);
  mpz_clear(s);
}

void
mpz_vector_skew_norm (mpz_ptr res, mpz_vector_srcptr a, mpz_srcptr skewness)
{
  mpz_vector_skew_dot_product (res, a, a, skewness);
}


/* compute r <- r - q*v */
void
mpz_vector_submul (mpz_vector_ptr r, mpz_srcptr q, mpz_vector_srcptr v)
{
  ASSERT_ALWAYS (r->dim == v->dim);
  for (unsigned int i = 0; i < r->dim; i++)
    mpz_submul(r->c[i], q, v->c[i]);
}

/* Lagrange Algorithm.
   Reduce lattice of rank 2.
   Reduce the lattice spanned by u and v and set a and b to the reduced basis.
   The obtained basis is always the two shortest vector.
   See, Nguyen And Vallee, "The LLL algorithm: Survey and applications" p.41
   If skewness parameter is lesser or equal to 1, then the non-skew version of
   mpz_vector_norm and mpz_vector_dot_product is used.
   */
void
mpz_vector_Lagrange (mpz_vector_ptr a, mpz_vector_ptr b,
                     mpz_vector_srcptr u, mpz_vector_srcptr v, mpz_srcptr skewness)
{
  mpz_t norm_a, norm_b, q, r, tmp;
  mpz_init (norm_a);
  mpz_init (norm_b);
  mpz_init (q);
  mpz_init (r);
  mpz_init (tmp);

  mpz_vector_set (a, u);
  mpz_vector_set (b, v);

  int skew = (mpz_cmp_ui (skewness, 1) > 0);

  if (skew)
  {
    mpz_vector_skew_norm (norm_a, a, skewness);
    mpz_vector_skew_norm (norm_b, b, skewness);
  }
  else
  {
    mpz_vector_norm (norm_a, a);
    mpz_vector_norm (norm_b, b);
  }

  if (mpz_cmp (norm_a, norm_b) < 0)
  {
    mpz_vector_swap (a, b);
    mpz_swap (norm_a, norm_b);
  }
  while (mpz_cmp (norm_a, norm_b) > 0)
  {
    if (skew)
      mpz_vector_skew_dot_product (tmp, a, b, skewness);
    else
      mpz_vector_dot_product (tmp, a, b);
    mpz_ndiv_q (q, tmp, norm_b);
    mpz_vector_submul (a, q, b);
    mpz_vector_swap (a, b);
    mpz_set (norm_a, norm_b);
    if (skew)
      mpz_vector_skew_norm (norm_b, b, skewness);
    else
      mpz_vector_norm (norm_b, b);
  }

  mpz_clear (norm_a);
  mpz_clear (norm_b);
  mpz_clear (q);
  mpz_clear (r);
  mpz_clear (tmp);
}


/* Find the maximun value of skewness for which mpz_vector_Lagrange returns 2
  vectors with non-zero d-th coordinates (i.e. 2 corresponding polynomials of
  degree d).
  This function assume that it exists S such that for skew <= S, the two
  vectors have non-zero d-th coordinates and for skew > S, this no more
  the case.

  Output: reduced_a, reduced_b and skew_used.
  Input: a, b, max_skewness and d.
*/

void
mpz_vector_reduce_with_max_skew (mpz_vector_ptr reduced_a, mpz_vector_ptr reduced_b,
                                 mpz_ptr skew_used, mpz_vector_srcptr a,
                                 mpz_vector_srcptr b, mpz_srcptr max_skewness, int d)
{
  mpz_t max_skew, min_skew, tmp;

  mpz_init_set (max_skew, max_skewness);
  mpz_init_set_ui (min_skew, 1);
  mpz_init (tmp);

  mpz_set (skew_used, max_skew);

  do
  {
    // Perform Lagrange algorithm a and b
    mpz_vector_Lagrange (reduced_a, reduced_b, a, b, skew_used);

    if (mpz_vector_is_coordinate_zero (reduced_a, d)
        || mpz_vector_is_coordinate_zero (reduced_b, d))
      mpz_set (max_skew, skew_used);
    else
    {
      mpz_set (min_skew, skew_used);
      mpz_sub (tmp, max_skew, min_skew);
      if (mpz_cmp_ui (tmp, 1) == 0) /* if max_skew - min_skew == 1, then stop */
        mpz_set (max_skew, skew_used);
    }
    mpz_add (skew_used, min_skew, max_skew);
    mpz_tdiv_q_2exp (skew_used, skew_used, 1);
  } while (mpz_cmp(min_skew, max_skew) < 0); // min_skew < max_skew

  mpz_clear (max_skew);
  mpz_clear (min_skew);
  mpz_clear (tmp);
}
