#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <gmp.h>
#include "utils_lattice.h"
#include "ideal.h"
#include "utils_mpz.h"
#include "sieving_bound.h"
#include "mat_double.h"
#include "double_vector.h"
#include "list_int64_vector.h"
#include "list_double_vector.h"
#include "list_int64_vector_index.h"
#include "gcd.h"
#include "utils_int64.h"
#include "macros.h"

#define swap(x, y) { unsigned int _tmp = (x); (x) = (y); (y) = _tmp; }
#define int64_swap_n(x, y) { int64_t *_tmp = (x); (x) = (y); (y) = _tmp; }

/*
 * Compute angle between two vectors.
 *
 * v0: first vector.
 * v1: second vector.
 */
static double angle_2_coordinate(int64_vector_srcptr v0, int64_vector_srcptr v1)
{
  ASSERT(v0->dim == v1->dim);
  ASSERT(v0->dim >= 2);

#ifndef NDEBUG
  for (unsigned int i = 2; i < v0->dim; i++) {
    ASSERT(v0->c[i] == v1->c[i]);
  } 
#endif // NDEBUG

  double dot = 0;
  for (unsigned int i = 0; i < 2; i++) {
    dot += (double)(v0->c[i] * v1->c[i]);
  }
  dot = dot / (sqrt((double)(v0->c[0] * v0->c[0] + v0->c[1] * v0->c[1])));
  dot = dot / (sqrt((double)(v1->c[0] * v1->c[0] + v1->c[1] * v1->c[1])));
  return acos(dot);
}

/*
 * Compute the determinant of a 2 * 2 matrix.
 */
MAYBE_UNUSED static int determinant2(mat_int64_srcptr matrix)
{
  return (int) matrix->coeff[1][1] * matrix->coeff[2][2] - matrix->coeff[2][1] *
    matrix->coeff[1][2];
}

int gauss_reduction(int64_vector_ptr v0, int64_vector_ptr v1,
    mat_int64_ptr U, int64_vector_srcptr v0_root, int64_vector_srcptr v1_root)
{
  int64_vector_set(v0, v0_root);
  int64_vector_set(v1, v1_root);

  //Verify that the vector are in the same plane.
  ASSERT(v0_root->dim == v1_root->dim);
  ASSERT(v0_root->dim >= 2);
#ifndef NDEBUG
  for (unsigned int i = 3; i < v1->dim; i++) {
    ASSERT(v0_root->c[i] == v1_root->c[i]);
  }

  if (U != NULL) {
    ASSERT(U->NumRows == U->NumCols);
    ASSERT(U->NumRows == 2);
  }

  int64_vector_t v0_tmp;
  int64_vector_init(v0_tmp, v0->dim);
  int64_vector_set(v0_tmp, v0);

  int64_vector_t v1_tmp;
  int64_vector_init(v1_tmp, v1->dim);
  int64_vector_set(v1_tmp, v1);
#endif // NDEBUG
  ASSERT(v0->dim == v1->dim);
  ASSERT(v0_root->dim >= v1->dim);

  //determinant of U, 0 if U is NULL.
  int det = 0;
  
  double B0 = int64_vector_norml2sqr(v0);
  double m = (double)int64_vector_dot_product(v0, v1) / B0;
  for (unsigned int i = 0; i < 2; i++) {
    v1->c[i] = v1->c[i] - (int64_t)round(m) * v0->c[i];
  }

  if (U != NULL) {
    U->coeff[1][1] = 1;
    U->coeff[1][2] = -(int64_t)round(m);
    U->coeff[2][1] = 0;
    U->coeff[2][2] = 1;
    det = 1;
  }

  double B1 = int64_vector_norml2sqr(v1);

  while (B1 < B0) {
    int64_vector_swap(v0, v1);
    B0 = B1;
    m = (double)int64_vector_dot_product(v0, v1) / B0;
    for (unsigned int i = 0; i < 2; i++) {
      v1->c[i] = v1->c[i] - (int64_t)round(m) * v0->c[i];
    }

    //Update U by U = U * transformation_matrix.
    if (U != NULL) {
      int64_t tmp = U->coeff[1][2];
      U->coeff[1][2] = U->coeff[1][1] - (int64_t)round(m) * U->coeff[1][2];
      U->coeff[1][1] = tmp;
      tmp = U->coeff[2][2];
      U->coeff[2][2] = U->coeff[2][1] - (int64_t)round(m) * U->coeff[2][2];
      U->coeff[2][1] = tmp;
      det = -det;
    }

    B1 = int64_vector_norml2sqr(v1);

    //To have an acute basis.
  }

  if (M_PI_2 < angle_2_coordinate(v0, v1))
  {
    for (unsigned int i = 0; i < 2; i++) {
      v1->c[i] = -v1->c[i];
    }

    if (U != NULL) {
      U->coeff[1][2] = -U->coeff[1][2];
      U->coeff[2][2] = -U->coeff[2][2];
      det = -det;
    }
  }

#ifndef NDEBUG
  if (U != NULL) {
    ASSERT(v0_tmp->c[0] * U->coeff[1][1] + v1_tmp->c[0] * U->coeff[2][1]
        == v0->c[0]);
    ASSERT(v0_tmp->c[1] * U->coeff[1][1] + v1_tmp->c[1] * U->coeff[2][1]
        == v0->c[1]);
    ASSERT(v0_tmp->c[0] * U->coeff[1][2] + v1_tmp->c[0] * U->coeff[2][2]
        == v1->c[0]);
    ASSERT(v0_tmp->c[1] * U->coeff[1][2] + v1_tmp->c[1] * U->coeff[2][2]
        == v1->c[1]);
    ASSERT(determinant2(U) == det);
    ASSERT(ABS(det) == 1);
  }

  int64_vector_clear(v0_tmp);
  int64_vector_clear(v1_tmp);
#endif // NDEBUG

  return det;
}

int gauss_reduction_zero(int64_vector_ptr v0, int64_vector_ptr v1,
    mat_int64_ptr U, int64_vector_srcptr v0_root, int64_vector_srcptr v1_root)
{
  int det = gauss_reduction(v0, v1, U, v0_root, v1_root);

  for (unsigned int i = 3; i < v1->dim; i++) {
    v0->c[i] = 0;
    v1->c[i] = 0;
  }

  return det; 
}

#ifndef SLLL_SAFE
void skew_LLL(mat_int64_ptr MSLLL, mat_int64_srcptr M_root,
    int64_vector_srcptr skewness, FILE * errstd)
{
  ASSERT(skewness->dim == M_root->NumCols);
  ASSERT(M_root->NumRows == M_root->NumCols);
  ASSERT(MSLLL->NumRows == M_root->NumCols);
  ASSERT(MSLLL->NumRows == M_root->NumCols);


#if 0
  mat_int64_t I_s;
  mat_int64_init(I_s, M_root->NumRows, M_root->NumCols);
  mat_int64_set_diag(I_s, skewness);

  mat_int64_t M;
  mat_int64_init(M, M_root->NumRows, M_root->NumCols);
  mat_int64_mul_mat_int64(M, I_s, M_root);


  mat_int64_t U;
  mat_int64_init(U, M_root->NumRows, M_root->NumCols);
  lll_Mqr_unimodular(U, M, errstd);

  mat_int64_mul_mat_int64(MSLLL, M_root, U);

  mat_int64_clear(U);
  mat_int64_clear(I_s);
  mat_int64_clear(M);
#else // 1
  int64_vector_t skew;
  int64_vector_init(skew, skewness->dim);
  int overflow = 0;
  for (unsigned int i = 0; i < skewness->dim; i++) {
    if (__builtin_smull_overflow(skewness->c[i], skewness->c[i], &(skew->c[i])))
    {
      overflow = 1;
      break;
    }
  }
  if (!overflow) {
    slll_Mqr(MSLLL, M_root, skew, errstd);
  } else {
    fprintf(errstd, "# Can not use slll_Mqr, use safe skew LLL instead.\n");

    mat_int64_t I_s;
    mat_int64_init(I_s, M_root->NumRows, M_root->NumCols);
    mat_int64_set_diag(I_s, skewness);

    mat_int64_t M;
    mat_int64_init(M, M_root->NumRows, M_root->NumCols);
    mat_int64_mul_mat_int64(M, I_s, M_root);


    mat_int64_t U;
    mat_int64_init(U, M_root->NumRows, M_root->NumCols);
    mat_int64_LLL_unimodular_transpose(U, M);

    mat_int64_mul_mat_int64(MSLLL, M_root, U);

    mat_int64_clear(U);
    mat_int64_clear(I_s);
    mat_int64_clear(M);
  }
  int64_vector_clear(skew);
#endif // 1

}
#else // SLLL_SAFE
void skew_LLL_safe(mat_int64_ptr MSLLL, mat_int64_srcptr M_root,
    int64_vector_srcptr skewness)
{
  ASSERT(skewness->dim == M_root->NumCols);
  ASSERT(M_root->NumRows == M_root->NumCols);
  ASSERT(MSLLL->NumRows == M_root->NumCols);
  ASSERT(MSLLL->NumRows == M_root->NumCols);

  mat_int64_t I_s;
  mat_int64_init(I_s, M_root->NumRows, M_root->NumCols);
  mat_int64_set_diag(I_s, skewness);

  mat_int64_t M;
  mat_int64_init(M, M_root->NumRows, M_root->NumCols);
  mat_int64_mul_mat_int64(M, I_s, M_root);


  mat_int64_t U;
  mat_int64_init(U, M_root->NumRows, M_root->NumCols);
  mat_int64_LLL_unimodular_transpose(U, M);

  mat_int64_mul_mat_int64(MSLLL, M_root, U);

  mat_int64_clear(U);
  mat_int64_clear(I_s);
  mat_int64_clear(M);
}
#endif // SLLL_SAFE

int reduce_qlattice(int64_vector_ptr v0, int64_vector_ptr v1,
    int64_vector_srcptr v0_root, int64_vector_srcptr v1_root, int64_t I)
{
  ASSERT(v0_root->c[1] == 0);
  ASSERT(v1_root->c[1] == 1);
  ASSERT(v0_root->dim == v1_root->dim);
  ASSERT(v0->dim == v1->dim);
  ASSERT(v0_root->dim == v1->dim);
  ASSERT(I > 0);
#ifndef NDEBUG
  for (unsigned int i = 3; i < v1->dim; i++) {
    ASSERT(v0_root->c[i] == v1_root->c[i]);
  }
#endif // NDEBUG

  int64_vector_set(v0, v0_root);
  int64_vector_set(v1, v1_root);

  int64_t a0 = -v0->c[0], b0 = v1->c[0], a1 = 0, b1 = 1, k;

  const int64_t hI = I;
  const int64_t mhI = -hI;

  while (b0 >= hI) {
    k = a0 / b0; a0 %= b0; a1 -= k * b1;

    if (a0 > mhI) {
      break;
    }
    k = b0 / a0; b0 %= a0; b1 -= k * a1;

    if (b0 < hI) {
      break;
    }
    k = a0 / b0; a0 %= b0; a1 -= k * b1;

    if (a0 > mhI) {
      break;
    }
    k = b0 / a0; b0 %= a0; b1 -= k * a1;

  }
  k = b0 - hI - a0;
  if (b0 > -a0) {
    if (!a0) {
      return 0;
    }
    k /= a0; b0 -= k * a0; b1 -= k * a1;

  } else {
    if (!b0) {
      return 0;
    }
    k /= b0; a0 += k * b0; a1 += k * b1;

  }

  ASSERT(a0 > mhI);
  ASSERT(0 >= a0);
  ASSERT(0 <= b0);
  ASSERT(b0 < hI);
  ASSERT(a1 > 0);
  ASSERT(b1 > 0);

  v0->c[0] = a0;
  v0->c[1] = a1;
  v1->c[0] = b0;
  v1->c[1] = b1;

  return 1;
}

int reduce_qlattice_zero(int64_vector_ptr v0, int64_vector_ptr v1,
    int64_vector_srcptr v0_root, int64_vector_srcptr v1_root, int64_t I)
{
  int res = reduce_qlattice(v0, v1, v0_root, v1_root, I);

  for (unsigned int i = 3; i < v1->dim; i++) {
    v0->c[i] = 0;
    v1->c[i] = 0;
  }

  return res;
}

/*
 * For instance, it is SV3.
 */
void SV4(list_int64_vector_ptr SV, int64_vector_srcptr v0_root,
    int64_vector_srcptr v1_root, int64_vector_srcptr v2)
{
  ASSERT(v0_root->dim == v1_root->dim);
#ifndef NDEBUG
  for (unsigned int j = 3; j < v0_root->dim; j++) {
    ASSERT(v0_root->c[j] == v1_root->c[j]);
    ASSERT(v0_root->c[j] == 0);
  }
#endif // NDEBUG
  ASSERT(v2->c[1] == 0);
  unsigned int v2_assert = 0;
  for (unsigned int i = 1; i < v2->dim; i++) {
    if (v2->c[i] == 1) {
      v2_assert++;
    } else {
      ASSERT(v2_assert < 2);
      ASSERT(v2->c[i] == 0);
    }
  }
  ASSERT(v2_assert == 1);

  int64_vector_t u;
  int64_vector_init(u, v0_root->dim);

  int64_vector_t v0;
  int64_vector_t v1;
  int64_vector_init(v0, v0_root->dim);
  int64_vector_init(v1, v1_root->dim);
  int64_vector_set(v0, v0_root);
  int64_vector_set(v1, v1_root);

  mat_int64_t U;
  mat_int64_init(U, 2, 2);

  int det = gauss_reduction(v0, v1, U, v0, v1);
  ASSERT(det != 0);

  double target0 = (double)v2->c[0] / (double)v0_root->c[0];
  //U^-1 * [target0, 0] = (v2_new_base_x, v2_new_base_y).
  double v2_new_base_x = (double)det * (double)U->coeff[2][2] * target0;
  double v2_new_base_y = -(double)det * (double)U->coeff[2][1] * target0;
  int64_t a = (int64_t)round(v2_new_base_x);
  int64_t b = (int64_t)round(v2_new_base_y);

  for (unsigned int i = 0; i < 2; i++) {
    u->c[i] = a * v0->c[i] + b * v1->c[i];
  }
  for (unsigned int i = 2; i < u->dim; i++) {
    u->c[i] = 0;
  }
  list_int64_vector_add_int64_vector(SV, u);

  //Build a triangle around the projection of v2 in the plane z = 0.
  if (v2_new_base_x - (double)a < 0) {
    for (unsigned int i = 0; i < 2; i++) {
      u->c[i] = SV->v[0]->c[i] - v0->c[i];
    }
    for (unsigned int i = 2; i < u->dim; i++) {
      u->c[i] = 0;
    }
    list_int64_vector_add_int64_vector(SV, u);
  } else {
    for (unsigned int i = 0; i < 2; i++) {
      u->c[i] = SV->v[0]->c[i] + v0->c[i];
    }
    for (unsigned int i = 2; i < u->dim; i++) {
      u->c[i] = 0;
    }
    list_int64_vector_add_int64_vector(SV, u);
  }

  if (v2_new_base_y - (double)b < 0) {
    for (unsigned int i = 0; i < 2; i++) {
      u->c[i] = SV->v[0]->c[i] - v1->c[i];
    }
    for (unsigned int i = 2; i < u->dim; i++) {
      u->c[i] = 0;
    }
    list_int64_vector_add_int64_vector(SV, u);
  } else {
    for (unsigned int i = 0; i < 2; i++) {
      u->c[i] = SV->v[0]->c[i] + v1->c[i];
    }
    for (unsigned int i = 2; i < u->dim; i++) {
      u->c[i] = 0;
    }
    list_int64_vector_add_int64_vector(SV, u);
  }

#ifndef NDEBUG
  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, v2->dim);
  int64_vector_set(v_tmp, v2);
  for (unsigned int i = 2; i < v_tmp->dim; i++) {
    v_tmp->c[i] = 0;
  }
  ASSERT(int64_vector_in_polytop_list_int64_vector(v_tmp, SV) == 1);
  int64_vector_clear(v_tmp);
#endif // NDEBUG

  for (unsigned int i = 0; i < 3; i++) {
    int64_vector_sub(SV->v[i], v2, SV->v[i]);
  }

  int64_vector_clear(v0);
  int64_vector_clear(v1);
  int64_vector_clear(u);
  mat_int64_clear(U);
}

/*
 * Do SV4 on Mqr. The vectors of Mqr must be written in columns.
 */
void SV4_Mqr(list_int64_vector_ptr SV, mat_int64_srcptr Mqr)
{
  ASSERT(Mqr->NumRows == Mqr->NumCols);
  ASSERT(Mqr->NumRows == 3);

  int64_vector_t * v = malloc(sizeof(int64_vector_t) * Mqr->NumRows);
  for (unsigned int i = 0; i < Mqr->NumCols; i++) {
    int64_vector_init(v[i], Mqr->NumRows);
    mat_int64_extract_vector(v[i], Mqr, i);
  }
  SV4(SV, v[0], v[1], v[2]);

  for (unsigned int i = 0; i < Mqr->NumCols; i++) {
    int64_vector_clear(v[i]);
  }
  free(v);
}

/*
 * Return the gap between the x coordinate and the closest border of the sieving
 * region defined by the sieving bound H.
 *
 * v: current vector.
 * H: sieving bound.
 */
/*static uint64_t difference_bound_x(int64_vector_srcptr v,*/
    /*sieving_bound_srcptr H)*/
/*{*/
  /*return (uint64_t) MIN(ABS(-(int64_t)H->h[0] - v->c[0]),*/
      /*ABS((int64_t)H->h[0] - 1 - v->c[0]));*/
/*}*/

/*
 * Sub task, compute how many times we need to add the Franke-Kleinjung vector.
 *
 * x: x coordinate of the current starting point.
 * xfk: x coordinate of the Franke-Kleinjung vector.
 * H0: bound of the sieving region on x.
 */
static int64_t compute_k(int64_t x, int64_t xfk, int64_t H0)
{
  return (int64_t)ceil((double)(-H0 - x) / (double)xfk);
}

void add_FK_vector(int64_vector_ptr v, list_int64_vector_srcptr list,
    int64_vector_srcptr e0, int64_vector_srcptr e1, sieving_bound_srcptr H)
{
  ASSERT(v->dim == H->t);

  int64_vector_t v_use;
  int64_vector_init(v_use, e0->dim);

  //Select the Franke-Kleinjung vector which have the larger x coordinate.
  if (-e0->c[0] > e1->c[0]) {
    ASSERT(v->dim == e0->dim);
    //Set v_use = -e0;
    for (unsigned int i = 0; i < e0->dim; i++) {
      v_use->c[i] = -e0->c[i];
    }
  } else {
    int64_vector_set(v_use, e1);
  }

  ASSERT(v_use->c[0] > 0);

  //Compute the distance.
  unsigned int pos = 0;
  int64_t k = compute_k(list->v[0]->c[0], v_use->c[0], (int64_t)H->h[0]);

  ASSERT(list->v[0]->c[0] + k * v_use->c[0] < (int64_t)H->h[0]);
  ASSERT(list->v[0]->c[0] + k * v_use->c[0] >= -(int64_t)H->h[0]);

  uint64_t min_y = (uint64_t)ABS(list->v[0]->c[1] + k * v_use->c[1]);

  for (unsigned int i = 1; i < list->length; i++) {
    int64_t k_tmp = compute_k(list->v[i]->c[0], v_use->c[0], (int64_t)H->h[0]);
   
    ASSERT(list->v[i]->c[0] + k_tmp * v_use->c[0] < (int64_t)H->h[0]);
    ASSERT(list->v[i]->c[0] + k_tmp * v_use->c[0] >= -(int64_t)H->h[0]);

    uint64_t tmp = (uint64_t)ABS(list->v[i]->c[1] + k_tmp * v_use->c[1]);
    if (tmp < min_y) {
      min_y = tmp;
      pos = i;
      k = k_tmp;
    }
  }

  int64_vector_set(v, list->v[pos]);

  for (unsigned int i = 0; i < 2; i++) {
    v->c[i] = v->c[i] + k * v_use->c[i];
  }

  ASSERT(v->c[0] < (int64_t)H->h[0]);
  ASSERT(v->c[0] >= -(int64_t)H->h[0]);

  int64_vector_clear(v_use);
}

unsigned int enum_pos_with_FK(int64_vector_ptr v, int64_vector_srcptr v_old,
    int64_vector_srcptr v0, int64_vector_srcptr v1, int64_t A, int64_t I)
{
  //v0 = (alpha, beta) and v1 = (gamma, delta)
  ASSERT(I > 0);
  ASSERT(v0->c[0] > -I);
  ASSERT(v0->c[0] <= 0);
  ASSERT(v0->c[1] > 0);
  ASSERT(v1->c[0] < I);
  ASSERT(v1->c[0] >= 0);
  ASSERT(v1->c[1] > 0);
  ASSERT(v1->c[0] - v0->c[0] >= I);
#ifndef NDEBUG
  for (unsigned int j = 3; j < v0->dim; j++) {
    ASSERT(v0->c[j] == 0);
    ASSERT(v1->c[j] == 0);
  }
#endif // NDEBUG
  ASSERT(v_old->c[0] < A + I);
  ASSERT(v_old->c[0] >= A);

  if (v_old->c[0] >= A - v0->c[0]) {
    int64_vector_add(v, v_old, v0);
    return 0;
  } else if (v_old->c[0] < A + I - v1->c[0]) {
    int64_vector_add(v, v_old, v1);
    return 1;
  } else {
    int64_vector_add(v, v_old, v0);
    int64_vector_add(v, v, v1);
    return 2;
  }
}

unsigned int enum_neg_with_FK(int64_vector_ptr v, int64_vector_srcptr v_old,
    int64_vector_srcptr v0, int64_vector_srcptr v1, int64_t A, int64_t I)
{
  //v0 = (alpha, beta) and v1 = (gamma, delta)
  ASSERT(I > 0);
  ASSERT(v0->c[0] > -I);
  ASSERT(v0->c[0] <= 0);
  ASSERT(v0->c[1] > 0);
  ASSERT(v1->c[0] < I);
  ASSERT(v1->c[0] >= 0);
  ASSERT(v1->c[1] > 0);
  ASSERT(v1->c[0] - v0->c[0] >= I);
#ifndef NDEBUG
  for (unsigned int j = 3; j < v0->dim; j++) {
    ASSERT(v0->c[j] == 0);
    ASSERT(v1->c[j] == 0);
  }
#endif // NDEBUG
  ASSERT(v_old->c[0] > -A - I);
  ASSERT(v_old->c[0] <= -A);

  if (v_old->c[0] <= -A + v0->c[0]) {
    int64_vector_sub(v, v_old, v0);
    return 0;
  } else if (v_old->c[0] > -A - I + v1->c[0]) {
    int64_vector_sub(v, v_old, v1);
    return 1;
  } else {
    int64_vector_sub(v, v_old, v0);
    int64_vector_sub(v, v, v1);
    return 2;
  }
}

void coordinate_FK_vector(uint64_t * coord_v0, uint64_t * coord_v1,
    int64_vector_srcptr v0, int64_vector_srcptr v1, sieving_bound_srcptr H,
    uint64_t number_element)
{
  ASSERT(2*H->h[0] > 0);
  ASSERT(v0->c[0] > -2*(int64_t)H->h[0]);
  ASSERT(v0->c[0] <= 0);
  ASSERT(v0->c[1] > 0);
  ASSERT(v1->c[0] < 2*(int64_t)H->h[0]);
  ASSERT(v1->c[0] >= 0);
  ASSERT(v1->c[1] > 0);
  ASSERT(v1->c[0] - v0->c[0] >= 2*(int64_t)H->h[0]);
#ifndef NDEBUG
  for (unsigned int j = 3; j < v0->dim; j++) {
    ASSERT(v0->c[j] == 0);
    ASSERT(v1->c[j] == 0);
  }
#endif // NDEBUG

  //TODO: there is maybe a better way to achieve this.
  * coord_v0 = index_vector(v0, H, number_element);
  * coord_v1 = index_vector(v1, H, number_element);
}

//TODO: is it a good idea to retun 0 if fail?
uint64_t index_vector(int64_vector_srcptr v, sieving_bound_srcptr H,
    uint64_t number_element)
{
  uint64_t index = 0;

  int64_vector_t v_start;
  int64_vector_init(v_start, v->dim);
  int64_vector_set_zero(v_start);
  int64_vector_t v_end;
  int64_vector_init(v_end, v->dim);
  int64_vector_set(v_end, v);
  for (unsigned int i = 0; i < v->dim - 1; i++) {
    if (v->c[i] >= 0) {
      v_start->c[i] = -(int64_t)H->h[i];
    } else {
      v_start->c[i] = (int64_t)(H->h[i] - 1);
    }
    v_end->c[i] = v_start->c[i] + v->c[i];
  }
  if (int64_vector_in_sieving_region(v_end, H)) {
    uint64_t index_end = array_int64_vector_index(v_end, H, number_element);
    uint64_t index_start =
      array_int64_vector_index(v_start, H, number_element);
    if (index_end > index_start) {
      index = index_end - index_start;
    } else {
      index = index_start - index_end;
    }
  } else {
    index = 0;
  }

  ASSERT(index < number_element);
  int64_vector_clear(v_start);
  int64_vector_clear(v_end);

  return index;
}

/*
 * Return the position of the vector with a x coordinate minimal, according to
 *  the classification defined by the stamp array and the value of the
 *  classification val_stamp.
 *
 * SV: list of vector.
 * H: sieving bound.
 * stamp: array with length equal to SV, gives the classification of the
 *  vectors.
 * val_stamp: value we want for the stamp of the vector.
 */
static unsigned int find_min_y(list_int64_vector_srcptr SV,
    unsigned char * stamp, unsigned char val_stamp)
{
  int64_t y = -1; 
  unsigned int pos = 0;

  for (unsigned int i = 0; i < SV->length; i++) {
    if (stamp[i] == val_stamp) {
      //Find the closest y to 0.
      //TODO: warning, with this, we forget the problem in H1 - 1.
      int64_t tmp = ABS(SV->v[i]->c[1]);
      if (y == -1) {
        y = tmp;
      }
      if (y >= tmp) {
        y = tmp;
        pos = i;
      }
    }
  }
  return pos;
}

void plane_sieve_next_plane(int64_vector_ptr vs, list_int64_vector_srcptr SV,
    int64_vector_srcptr e0, int64_vector_srcptr e1, sieving_bound_srcptr H,
    int boolean)
{
  // Contain all the possible vectors to go from z=d to z=d+1.
  list_int64_vector_t list;
  list_int64_vector_init(list, SV->v[0]->dim);
  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, SV->v[0]->dim);
  /*
   * 0: in the sieving region.
   * 1: in [-H_0, H_0[
   * 2: outside
   */
  unsigned char * assert = malloc(sizeof(unsigned char) * SV->length);
  unsigned char found = 2;

  for (unsigned int i = 0; i < SV->length; i++) {
    int64_vector_add(v_tmp, vs, SV->v[i]);
    list_int64_vector_add_int64_vector(list, v_tmp);

    if (-(int64_t) H->h[0] <= v_tmp->c[0] && (int64_t)H->h[0] > v_tmp->c[0]) {
      if (-(int64_t) H->h[1] <= v_tmp->c[1] && (int64_t)H->h[1] > v_tmp->c[1])
      {
        assert[i] = 0;
        found = 0;
      } else {
        assert[i] = 1;
        if (found > 1) {
          found = 1;
        }
      }
    } else {
      assert[i] = 2;
    }
  }
  int64_vector_clear(v_tmp);

  if (found < 2) {
    unsigned pos = find_min_y(list, assert, found);
    int64_vector_set(vs, list->v[pos]);
  } else {
    ASSERT(found == 2);
    
    if (boolean) {
      add_FK_vector(vs, list, e0, e1, H);
      ASSERT(vs->c[0] < (int64_t)H->h[0]);
      ASSERT(vs->c[0] >= -(int64_t)H->h[0]);
    } else {
      ASSERT(boolean == 0);

      unsigned pos = find_min_y(list, assert, 2);
      int64_vector_set(vs, list->v[pos]);
      //TODO: not the closest point to [-H0, H0[.
      while (vs->c[0] < -(int64_t)H->h[0]) {
        int64_vector_add(vs, vs, e0);
      }
      while (vs->c[0] > (int64_t)H->h[0] - 1) {
        int64_vector_sub(vs, vs, e0);
      }
    }
  }
  
  free(assert);
  
  list_int64_vector_clear(list);
}

void double_vector_gram_schmidt(list_double_vector_ptr list_new,
    mat_double_ptr m, list_double_vector_srcptr list_old)
{
  for (unsigned int i = 0; i < list_old->length; i++) {
    list_double_vector_add_double_vector(list_new, list_old->v[i]);
  }

  double_vector_t v_tmp;
  double_vector_init(v_tmp, list_new->v[0]->dim);
  for (unsigned int i = 0; i < list_new->length; i++) {
    m->coeff[i + 1][i + 1] = 1;
  }

  for (unsigned int i = 0; i < list_new->length; i++) {
    for (unsigned int j = 0; j < i; j++) {
      m->coeff[i + 1][j + 1] = double_vector_orthogonal_projection(v_tmp,
          list_new->v[j], list_new->v[i]); 
      double_vector_sub(list_new->v[i], list_new->v[i], v_tmp);
    }
  }

  double_vector_clear(v_tmp);
}

#ifdef SPACE_SIEVE_REDUCE_QLATTICE
/*
 * Perform reduce q lattice if a Franke-Kleinjung vector is known. list contains
 *  two times the same vector. r is the r of the ideal we sieve, and I is equal
 *  to 2*H0.
 */
static int reduce_qlattice_output(list_int64_vector_ptr list, int64_t r,
    int64_t I)
{
  //list->v[0]->c[0] <= 0, list->v[1]->c[0] >= 0
  //index: where is the vector of Franke-Kleinjung.
  unsigned int index = 0;
  unsigned int index_new = 1;
  if (list->v[0]->c[0] > 0) {
    index = 1;
    index_new = 0;
  }

  ASSERT(list->v[index]->c[1] >= 0);
  ASSERT(I > 0);
#ifndef NDEBUG
  for (unsigned int i = 3; i < list->v[index]->dim; i++) {
    ASSERT(0 == list->v[index]->c[i]);
  }
#endif // NDEBUG

  //TODO: not sure that only this case can fail.
  if (list->v[index]->c[0] == 0 || list->v[index]->c[1] == 0) {
    return 0;
  }

  int64_t x = 0;
  int64_t y = 0;

  if (index == 0) {
    ASSERT(list->v[0]->c[0] < 0);

    int64_t g = 0;
    int64_xgcd(&g, &x, &y, list->v[0]->c[1], -list->v[0]->c[0]);

    if (g == -1) {
      g = 1;
      x = -x;
      y = -y;
    }

    ASSERT(g == 1);
  } else {
    ASSERT(list->v[0]->c[0] > 0);

    int64_t g = 0;
    int64_xgcd(&g, &x, &y, -list->v[0]->c[1], list->v[0]->c[0]);

    if (g != 1) {
      g = -g;
      x = -x;
      y = -y;
    }

    ASSERT(g == 1);
  }

  //Beware of overflow
  x = x * r;
  y = y * r;

  int64_t k = 0;
  if (index == 0) {
    //k = ceil((I - x) / v[0])
    k = siceildiv(I - x, list->v[0]->c[0]);
    /*k = (int64_t) ceil(((double)I - (double)x) / (double)list->v[0]->c[0]);*/
    if (I - x == k * list->v[0]->c[0]) {
      k++;
    }
  } else {
    //k = ceil((-I - x) / v[0])
    k = siceildiv(-I - x, list->v[1]->c[0]);
    /*k = (int64_t) ceil((-(double)I - (double)x) / (double)list->v[1]->c[0]);*/
    if (-I - x == k * list->v[0]->c[0]) {
      k++;
    }
  }
  list->v[index_new]->c[0] = x + k * list->v[index]->c[0];
  list->v[index_new]->c[1] = y + k * list->v[index]->c[1];

#ifndef NDEBUG
  ASSERT(list->v[0]->c[0] * list->v[1]->c[1] - list->v[0]->c[1]
    * list->v[1]->c[0] == -r);
#endif // NDEBUG

  return 1;
}
#endif // SPACE_SIEVE_REDUCE_QLATTICE

/*
 * Verify if a vector can be interesting for the space sieve, ie the t-1 first
 *  coordinate must be in [-2*Hi, 2*Hi].
 */
static int space_sieve_good_vector(int64_vector_srcptr v,
    sieving_bound_srcptr H)
{
  ASSERT(v->dim == H->t);

  for (unsigned int i = 0; i < v->dim - 1; i++) {
    if ((-2 * (int64_t)H->h[i]) >= v->c[i] ||
        (2 * (int64_t)H->h[i]) <= v->c[i]) {
      return 0;
    }
  }
  if (0 > v->c[v->dim - 1] || (int64_t)H->h[v->dim - 1] <= v->c[v->dim - 1]) {
    return 0;
  }
  return 1; 
}

/*
 * Verify if a vector is in a list in which vector has the last coordinate equal
 *  to 0.
 */
static int int64_vector_in_list_zero(int64_vector_srcptr v_tmp,
    list_int64_vector_index_srcptr list_zero)
{
  for (unsigned int i = 0; i < list_zero->length; i++) {
    for (unsigned int j = 0; j < v_tmp->dim - 1; j++) {
      //TODO: ABS is no longer required.
      if (ABS(v_tmp->c[j]) == ABS(list_zero->v[i]->vec->c[j])) {
        return 1;
      }
    }
  }
  return 0;
}

/*
 * Store v in the good list (list_zero has vector with 0 for last coordinate).
 * Return 1 if a vector in list has a z coordinate equal to 1.
 */
static unsigned int good_vector_in_list(list_int64_vector_index_ptr list,
    list_int64_vector_index_ptr list_zero, int64_vector_srcptr v,
    sieving_bound_srcptr H, int64_t number_element)
{
  unsigned int vector_1 = 0;
  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, v->dim);
  int64_vector_set(v_tmp, v);
  int64_vector_reduce(v_tmp, v_tmp);
  if (space_sieve_good_vector(v_tmp, H)) {
    if (v_tmp->c[2] == 0) {
      //v_tmp != 0.
      if (v_tmp->c[1] != 0 || v_tmp->c[0] != 0) {
        if (v_tmp->c[1] < 0) {
          v_tmp->c[0] = -v_tmp->c[0];
          v_tmp->c[1] = -v_tmp->c[1];
        }
        if (!int64_vector_in_list_zero(v_tmp, list_zero)) {
          list_int64_vector_index_add_int64_vector_index(list_zero,
              v_tmp, index_vector(v_tmp, H, number_element));
        }
      }
    } else {
      list_int64_vector_index_add_int64_vector_index(list, v_tmp, 0);
      if (v_tmp->c[2] == 1) {
        vector_1 = 1;
      }
    }
  }
  int64_vector_clear(v_tmp);
  return vector_1;
}

//TODO: change that.
/*
 * Do some linear combination of the vector of MSLLL to initialize space sieve.
 * Return 1 if a vector with a z coordinate equal to 1 is in list.
 */
static unsigned int space_sieve_linear_combination(
    list_int64_vector_index_ptr list, list_int64_vector_index_ptr list_zero,
    mat_int64_srcptr MSLLL, sieving_bound_srcptr H, uint64_t number_element)
{
  unsigned int vector_1 = 0;

  list_int64_vector_t list_tmp;
  list_int64_vector_init(list_tmp, 3);
  list_int64_vector_extract_mat_int64(list_tmp, MSLLL);

  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, MSLLL->NumRows);
  int64_vector_t v_tmp_n;
  int64_vector_init(v_tmp_n, MSLLL->NumRows);
  int64_vector_set_zero(v_tmp);
  int64_t l = -1;
  int64_t m = -1;
  int64_t n = -2;
  //v_tmp = l * v[0] + n * v[1] + 1 * v[2], l in [-1, 2[, m in [-1, 2[, n in
  //[-1, 2[.
  int64_vector_sub(v_tmp, v_tmp, list_tmp->v[0]);
  int64_vector_sub(v_tmp, v_tmp, list_tmp->v[1]);
  int64_vector_addmul(v_tmp, v_tmp, list_tmp->v[2], -2);

  //14 = 3^3 / 2 + 1;
  for (unsigned int i = 0; i < 14; i++) {
    if (n < 1) {
      n++;
      int64_vector_add(v_tmp, v_tmp, list_tmp->v[2]);
    } else {
      n = -1;
      int64_vector_addmul(v_tmp, v_tmp, list_tmp->v[2], -2);
      if (m < 1) {
        m++;
        int64_vector_add(v_tmp, v_tmp, list_tmp->v[1]);
      } else {
        m = -1;
        int64_vector_addmul(v_tmp, v_tmp, list_tmp->v[1], -2);
        l = l + 1;
        int64_vector_add(v_tmp, v_tmp, list_tmp->v[0]);
      }
    }

#ifndef NDEBUG
  int64_vector_t v_tmp1;
  int64_vector_init(v_tmp1, MSLLL->NumRows);
  int64_vector_set_zero(v_tmp1);
  int64_vector_addmul(v_tmp1, v_tmp1, list_tmp->v[0], l); 
  int64_vector_addmul(v_tmp1, v_tmp1, list_tmp->v[1], m); 
  int64_vector_addmul(v_tmp1, v_tmp1, list_tmp->v[2], n);
  ASSERT(int64_vector_equal(v_tmp, v_tmp1));
  int64_vector_clear(v_tmp1);
#endif // NDEBUG

    if (v_tmp->c[2] < 0) {
      for (unsigned int i = 0; i < v_tmp->dim; i++) {
        v_tmp_n->c[i] = -v_tmp->c[i];
      }
      vector_1 = good_vector_in_list(list, list_zero, v_tmp_n, H,
          number_element);
    } else {
      vector_1 = good_vector_in_list(list, list_zero, v_tmp, H,
          number_element);
    }
  }

  for (unsigned int i = 0; i < list_zero->length; i++) {
    ASSERT(int64_vector_gcd(list_zero->v[i]->vec) == 1);
  }

  int64_vector_clear(v_tmp);
  int64_vector_clear(v_tmp_n);
  list_int64_vector_clear(list_tmp);

  return vector_1;
}

#ifdef SPACE_SIEVE_ENTROPY
/*
 * When a new vector is added to list, do some linear combination of the last
 *  vector and the previous vectors to find news vectors.
 */
static void space_sieve_generate_new_vectors(
    list_int64_vector_index_ptr list, list_int64_vector_index_ptr list_zero,
    sieving_bound_srcptr H, uint64_t number_element)
{
  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, list->vector_dim);
  unsigned int ind = list->length - 1;
  for (unsigned int i = 0; i < ind; i++) {
    //v_tmp = l * v[0] + n * v[1] + m * v[2]
    int64_vector_set(v_tmp, list->v[ind]->vec);
    int64_vector_add(v_tmp, v_tmp, list->v[i]->vec);
    good_vector_in_list(list, list_zero, v_tmp, H, number_element);
    int64_vector_set(v_tmp, list->v[ind]->vec);
    int64_vector_sub(v_tmp, v_tmp, list->v[i]->vec);
    good_vector_in_list(list, list_zero, v_tmp, H, number_element);
  }

#ifndef NDEBUG
  for (unsigned int i = 0; i < list_zero->length; i++) {
    ASSERT(int64_vector_gcd(list_zero->v[i]->vec) == 1);
  }
#endif // NDEBUG

  int64_vector_clear(v_tmp);
}
#endif // SPACE_SIEVE_ENTROPY

/*
 * Return 1 if a vector is in the sieving region except the last coordinate, 0
 *  otherwise.
 */
int int64_vector_in_sieving_region_dim(int64_vector_srcptr v,
    sieving_bound_srcptr H) {
  ASSERT(v->dim == H->t);

  for (unsigned int i = 0; i < v->dim - 1; i++) {
    if (-(int64_t)H->h[i] > v->c[i] || ((int64_t)H->h[i]) <= v->c[i]) {
      return 0;
    }
  }
  return 1; 
}

void plane_sieve_1_enum_plane_incomplete(int64_vector_ptr v_in,
    int64_vector_srcptr e0, int64_vector_srcptr e1, sieving_bound_srcptr H)
{
  int64_vector_t v;
  int64_vector_init(v, v_in->dim);
  int64_vector_set(v, v_in);
  //Perform Franke-Kleinjung enumeration.
  //x increases.
  while (v->c[1] < (int64_t)H->h[1]) {
    if (v->c[1] >= -(int64_t)H->h[1]) {
      int64_vector_set(v_in, v);
      int64_vector_clear(v);
      return;
    }
    enum_pos_with_FK(v, v, e0, e1, -(int64_t)H->h[0], 2 * (int64_t)
        H->h[0]);
  }

  //x decreases.
  int64_vector_set(v, v_in);
  enum_neg_with_FK(v, v, e0, e1, -(int64_t)H->h[0] + 1, 2 * (int64_t)H->h[0]);
  while (v->c[1] >= -(int64_t)H->h[1]) {
    if (v->c[1] < (int64_t)H->h[1]) {
      int64_vector_set(v_in, v);
      int64_vector_clear(v);
      return;
    }
    enum_neg_with_FK(v, v, e0, e1, -(int64_t)H->h[0] + 1, 2 *
        (int64_t) H->h[0]);
  }
  int64_vector_clear(v);
}

void plane_sieve_1_incomplete(int64_vector_ptr s_out, int64_vector_srcptr s,
    MAYBE_UNUSED mat_int64_srcptr Mqr, sieving_bound_srcptr H,
    list_int64_vector_srcptr list_FK, list_int64_vector_srcptr list_SV)
{
  ASSERT(Mqr->NumRows == Mqr->NumCols);
  ASSERT(Mqr->NumRows == 3);

  int64_vector_set(s_out, s);

  plane_sieve_next_plane(s_out, list_SV, list_FK->v[0], list_FK->v[1], H, 1);
  //Enumerate the element of the sieving region.
  for (unsigned int d = (unsigned int) s->c[2] + 1; d < H->h[2]; d++) {
    plane_sieve_1_enum_plane_incomplete
      (s_out, list_FK->v[0], list_FK->v[1], H);
    if (int64_vector_in_sieving_region(s_out, H)) {
      return;
    }

    //Jump in the next plane.
    plane_sieve_next_plane(s_out, list_SV, list_FK->v[0], list_FK->v[1], H, 1);
  }
}

unsigned int space_sieve_1_init(list_int64_vector_index_ptr list_vec,
    list_int64_vector_index_ptr list_vec_zero, ideal_1_srcptr r,
    mat_int64_srcptr Mqr, sieving_bound_srcptr H, uint64_t number_element,
    MAYBE_UNUSED unsigned int * skew_lll_fail,
    MAYBE_UNUSED FILE * file_space_sieve_stat, MAYBE_UNUSED FILE * errstd)
{
  int64_vector_t skewness;
  int64_vector_init(skewness, 3);
#ifdef SKEWNESS_TRUE
  skewness->c[0] = (int64_t)r->ideal->r;
  skewness->c[1] = (int64_t) (H->h[0] / H->h[1]) * (int64_t)r->ideal->r;
#ifdef SKEWNESS
  skewness->c[2] = (int64_t)(SKEWNESS * H->h[0] * H->h[0] * H->h[1]);
  //Theoretical bound.
  /*skewness->c[2] = (int64_t)(8 * H->h[0] * H->h[0] * H->h[1]);*/
#else //SKEWNESS
  skewness->c[2] = (int64_t)(2 * H->h[0] * H->h[0] * H->h[1]);
#endif //SKEWNESS
#else // SKEWNESS_TRUE
  skewness->c[0] = 1;
  skewness->c[1] = (int64_t) (H->h[0] / H->h[1]);
#ifdef SKEWNESS
  skewness->c[2] = (int64_t)(SKEWNESS * H->h[0] * H->h[0] * H->h[1]) / (int64_t)r->ideal->r;
  //Theoretical bound.
  /*skewness->c[2] = (int64_t)(8 * H->h[0] * H->h[0] * H->h[1]) / (int64_t)r->ideal->r;*/
#else //SKEWNESS
  skewness->c[2] = (int64_t)(2 * H->h[0] * H->h[0] * H->h[1]) / (int64_t)r->ideal->r;
#endif // SKEWNESS
#endif // SKEWNESS_TRUE

  for (unsigned int i = 0; i < skewness->dim; i++) {
    ASSERT(skewness->c[i] >= 0);
    if (skewness->c[i] == 0) {
      skewness->c[i] = 1;
    }
  }

  mat_int64_t MSLLL;
  mat_int64_init(MSLLL, Mqr->NumRows, Mqr->NumCols);
#ifdef SLLL_SAFE
  skew_LLL_safe(MSLLL, Mqr, skewness);
#else // SLLL_SAFE
  skew_LLL(MSLLL, Mqr, skewness, errstd);
#endif // SLLL_SAFE
  int64_vector_clear(skewness);

  for (unsigned int col = 1; col <= MSLLL->NumCols; col++) {
    if (MSLLL->coeff[3][col] < 0) {
      for (unsigned int row = 1; row <= MSLLL->NumRows; row++) {
        MSLLL->coeff[row][col] = - MSLLL->coeff[row][col];
      }
    }
  }

#ifdef SPACE_SIEVE_STAT
  int64_t * target = (int64_t *) malloc(sizeof(int64_t) * 2);
  target[0] = (int64_t) (2 * H->h[0]);
  target[1] = (int64_t) (2 * H->h[1]);
  * skew_lll_fail = 0;
  unsigned int out = 0;
  unsigned int col = 1;
  while (col <= MSLLL->NumCols) {
    //Just verify the 2 first coefficients of a vector.
    for (unsigned int row = 1; row < MSLLL->NumRows; row++) {
      if (ABS(MSLLL->coeff[row][col]) > target[row - 1]) {
        * skew_lll_fail =  1;
        out++;
        break;
      }
    }
    col++;
  }

  fprintf(file_space_sieve_stat, "Mqr =\n");
  mat_int64_fprintf(file_space_sieve_stat, Mqr);
  fprintf(file_space_sieve_stat, "target = [%" PRId64 ", %" PRId64 ", ?]\n",
      target[0], target[1]);
  fprintf(file_space_sieve_stat, "MSLLL =\n");
  mat_int64_fprintf(file_space_sieve_stat, MSLLL);
  if (out > 1) {
    fprintf(file_space_sieve_stat, "%u vectors out\n", out);
  }

  free(target);
#endif // SPACE_SIEVE_STAT
  
  //TODO: we must generate all the zero vectors!
  unsigned int vector_1 = space_sieve_linear_combination(
      list_vec, list_vec_zero, MSLLL, H, number_element);
  mat_int64_clear(MSLLL);

  list_int64_vector_index_sort_last(list_vec);

#ifdef SPACE_SIEVE_STAT
  fprintf(file_space_sieve_stat, "Generate %u vectors with z == 0 and %u \
vectors with z != 0.\n", list_vec_zero->length, list_vec->length);
#endif // SPACE_SIEVE_STAT

  return vector_1;
}

int space_sieve_1_plane_sieve_init(list_int64_vector_ptr list_SV,
    list_int64_vector_ptr list_FK, list_int64_vector_index_ptr list_vec,
    list_int64_vector_index_ptr list_vec_zero, MAYBE_UNUSED ideal_1_srcptr r,
    sieving_bound_srcptr H, mat_int64_srcptr Mqr,
    unsigned int vector_1, uint64_t number_element, unsigned int * new_vec)
{
  int64_vector_t * vec = malloc(sizeof(int64_vector_t) * Mqr->NumRows);
  for (unsigned int i = 0; i < Mqr->NumCols; i++) {
    int64_vector_init(vec[i], Mqr->NumRows);
    mat_int64_extract_vector(vec[i], Mqr, i);
  }
  if (vector_1 && list_vec->length != 0) {
    unsigned int cpt = 0;
    while (cpt < list_vec->length && list_vec->v[cpt]->vec->c[2] < 2) {
      if (list_vec->v[cpt]->vec->c[2] == 1) {
        list_int64_vector_add_int64_vector(list_SV,
            list_vec->v[cpt]->vec);
      }
      cpt++;
    }
  } else {
    ASSERT(vector_1 == 0);
    SV4(list_SV, vec[0], vec[1], vec[2]);
  }
  int boolean = 0;
#ifdef SPACE_SIEVE_REDUCE_QLATTICE
  if (list_vec_zero->length == 0) {
    //Stupid add, just to declare this two vector.
    list_int64_vector_add_int64_vector(list_FK, vec[0]);
    list_int64_vector_add_int64_vector(list_FK, vec[1]);
    boolean = reduce_qlattice(list_FK->v[0], list_FK->v[1], vec[0],
        vec[1], (int64_t)(2 * H->h[0]));
  } else if (list_vec_zero->length == 1) {
    list_int64_vector_add_int64_vector(list_FK,
        list_vec_zero->v[0]->vec);
    list_int64_vector_add_int64_vector(list_FK,
        list_vec_zero->v[0]->vec);
    boolean = reduce_qlattice_output(list_FK, (int64_t)r->ideal->r,
        (int64_t)(2 * H->h[0]));
  } else {
    ASSERT(list_vec_zero->length == 2);
    if (0 >= list_vec_zero->v[0]->vec->c[0]) {
      list_int64_vector_add_int64_vector(list_FK,
          list_vec_zero->v[0]->vec);
      list_int64_vector_add_int64_vector(list_FK,
          list_vec_zero->v[1]->vec);
    } else {
      list_int64_vector_add_int64_vector(list_FK,
          list_vec_zero->v[1]->vec);
      list_int64_vector_add_int64_vector(list_FK,
          list_vec_zero->v[0]->vec);
    }
    boolean = 1;
  }
#else
  list_int64_vector_add_int64_vector(list_FK, vec[0]);
  list_int64_vector_add_int64_vector(list_FK, vec[1]);
  boolean = reduce_qlattice(list_FK->v[0], list_FK->v[1], vec[0],
      vec[1], (int64_t)(2 * H->h[0]));
#endif

  for (unsigned int i = 0; i < Mqr->NumCols; i++) {
    int64_vector_clear(vec[i]);
  }
  free(vec);


  if (!boolean) {
    ASSERT(boolean == 0);
    ASSERT(list_FK->v[0]->c[1] == 0);
    ASSERT(list_FK->v[1]->c[0] == 0);

    return 0;
  }

  ASSERT(list_FK->v[0]->c[0] > -(int64_t)(2 * H->h[0]));
  ASSERT(0 >= list_FK->v[0]->c[0]);
  ASSERT(0 <= list_FK->v[1]->c[0]);
  ASSERT(list_FK->v[1]->c[0] < (int64_t)(2 * H->h[0]));
  ASSERT(list_FK->v[0]->c[1] > 0);
  ASSERT(list_FK->v[1]->c[1] > 0);
  //TODO: Go up, and just when we need.

  if (space_sieve_good_vector(list_FK->v[0], H)) {
    if (!int64_vector_in_list_zero(list_FK->v[0], list_vec_zero)) {
      list_int64_vector_index_add_int64_vector_index(list_vec_zero,
          list_FK->v[0],
          index_vector(list_FK->v[0], H, number_element));
      * new_vec = 0;
    }
  }
  if (space_sieve_good_vector(list_FK->v[1], H)) {
    if (!int64_vector_in_list_zero(list_FK->v[1], list_vec_zero)) {
      list_int64_vector_index_add_int64_vector_index(list_vec_zero,
          list_FK->v[1],
          index_vector(list_FK->v[1], H, number_element));
      * new_vec = 0;
    }
  }

  return 1;
}

unsigned int space_sieve_1_next_plane_seek(int64_vector_ptr s_tmp,
    unsigned int * index_vec, unsigned int * s_change,
    list_int64_vector_srcptr list_s, list_int64_vector_index_srcptr list_vec,
    sieving_bound_srcptr H, MAYBE_UNUSED int64_vector_srcptr s)
{
  ASSERT(* s_change == 0);
  ASSERT(int64_vector_equal(list_s->v[0], s));
  ASSERT(s->dim == s_tmp->dim);

  unsigned int hit = 0;
  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, s_tmp->dim);

  for (unsigned int i = 0; i < s_tmp->dim; i++) {
    s_tmp->c[i] = (int64_t)H->h[i];
  }
  for (unsigned int i = 0; i < list_s->length; i++) {
    ASSERT(list_s->v[i]->c[2] == s->c[2]);

    for (unsigned int j = 0; j < list_vec->length; j++) {
      int64_vector_add(v_tmp, list_s->v[i], list_vec->v[j]->vec);
      if (int64_vector_in_sieving_region_dim(v_tmp, H)) {
        if (v_tmp->c[2] < s_tmp->c[2]) {
          int64_vector_set(s_tmp, v_tmp);
          * index_vec = j;
          if (i != 0) {
            * s_change = 1;
          }
        }
        hit = 1;
#ifndef SPACE_SIEVE_SEEK_ENUM_ALL
        int64_vector_clear(v_tmp);
        return hit;
#endif // SPACE_SIEVE_SEEK_ENUM_ALL
      }
    }
  }

  int64_vector_clear(v_tmp);
  return hit;
}

#ifndef SLLL_SAFE
/* This uses Paul Zimmermann implementation.

   LLL using exact multiprecision arithmetic.
   Translated from NTL 4.1a <http://www.shoup.net/>
   into GMP <http://www.swox.se/gmp/> 
   by Paul Zimmermann, July 2000.

   Revised April 4, 2002 (bug found by Jens Franke <franke (at) math
   (dot) uni-bonn (dot) de>).

   This program is open-source software distributed under the terms 
   of the GNU General Public License <http://www.fsf.org/copyleft/gpl.html>.

   This follows the implementation in utils/lll.c, but for mat_int64_t.
 */

// Return 1 if not overflow, 0 if overflow.
// Return 1 if conversion is possible, 0 otherwise.
static int int128_convert_64(int64_t * x_64, __int128_t x)
{
  __int128_t min = INT64_MIN;
  __int128_t max = INT64_MAX;

  * x_64 = (int64_t) x;

  int ret = (min <= x && x <= max);
  ASSERT(ret == 0 || ret == 1);

  return ret;
}

static int skewinnerproduct_square(int64_t * x, int64_t * a, int64_t * b,
    unsigned int n, FILE * errstd, int64_vector_srcptr skewness)
{
  int64_t tmp_skew = 0;
  * x = 0;
  if (__builtin_smull_overflow(b[1], skewness->c[0], &tmp_skew) ||
      __builtin_smull_overflow(a[1], tmp_skew, x)) {
    fprintf(errstd, "# ----- LLL int64 -----\n");
    fprintf(errstd, "# Overflow during skewinnerproduct_square, ");
    return 0;
  }
  int64_t tmp = 0;
  for (unsigned int i = 2; i <= n; i++) {
    /*x = x + a[i] * b[i];*/
    if (__builtin_smull_overflow(b[i], skewness->c[i - 1], &tmp_skew) ||
        __builtin_smull_overflow(a[i], tmp_skew, &tmp) ||
        __builtin_saddl_overflow(* x, tmp, x)) {
      fprintf(errstd, "# ----- LLL int64 -----\n");
      fprintf(errstd, "# Overflow during skewinnerproduct_square, ");
      return 0;
    }
  }

#ifdef ASSERT_LLL
  mpz_t a_Z, b_Z, skew, res, ret;
  mpz_init(a_Z);
  mpz_init(b_Z);
  mpz_init(skew);
  mpz_init(res);
  mpz_init(ret);
  mpz_set_int64(ret, * x);
  mpz_set_int64(a_Z, a[1]);
  mpz_set_int64(b_Z, b[1]);
  mpz_set_int64(skew, skewness->c[0]);

  mpz_mul (skew, skew, skew);
  mpz_mul (b_Z, b_Z, skew);
  mpz_mul (res, a_Z, b_Z);
  for (unsigned int i = 2; i <= n; i++) {
    mpz_set_int64(a_Z, a[i]);
    mpz_set_int64(b_Z, b[i]);
    mpz_set_int64(skew, skewness->c[i - 1]);
    mpz_mul (skew, skew, skew);
    mpz_mul (b_Z, b_Z, skew);
    mpz_addmul(res, a_Z, b_Z);
  }
  ASSERT_ALWAYS(mpz_cmp(res, ret) == 0);
  mpz_clear(a_Z);
  mpz_clear(b_Z);
  mpz_clear(skew);
  mpz_clear(res);
  mpz_clear(ret);
#endif // ASSERT_LLL

  return 1;
}

// Return 0 if overflow, 1 otherwise.
static int skewinnerproduct(int64_t * x, int64_t * a, int64_t * b,
    unsigned int n, FILE * errstd, int64_vector_srcptr skewness)
{
  int64_t tmp_skew = 0;
  * x = 0;
  if (__builtin_smull_overflow(skewness->c[0], skewness->c[0], &tmp_skew) ||
      __builtin_smull_overflow(b[1], tmp_skew, &tmp_skew) ||
      __builtin_smull_overflow(a[1], tmp_skew, x)) {
    fprintf(errstd, "# ----- LLL int64 -----\n");
    fprintf(errstd, "# Overflow during skewinnerproduct, ");
    return 0;
  }
  int64_t tmp = 0;
  for (unsigned int i = 2; i <= n; i++) {
    /*x = x + a[i] * b[i];*/
    if (__builtin_smull_overflow(skewness->c[i - 1], skewness->c[i - 1], &tmp_skew) ||
        __builtin_smull_overflow(b[i], tmp_skew, &tmp_skew) ||
        __builtin_smull_overflow(a[i], tmp_skew, &tmp) ||
        __builtin_saddl_overflow(* x, tmp, x)) {
      fprintf(errstd, "# ----- LLL int64 -----\n");
      fprintf(errstd, "# Overflow during skewinnerproduct, ");
      return 0;
    }
  }

#ifdef ASSERT_LLL
  mpz_t a_Z, b_Z, skew, res, ret;
  mpz_init(a_Z);
  mpz_init(b_Z);
  mpz_init(skew);
  mpz_init(res);
  mpz_init(ret);
  mpz_set_int64(ret, * x);
  mpz_set_int64(a_Z, a[1]);
  mpz_set_int64(b_Z, b[1]);
  mpz_set_int64(skew, skewness->c[0]);

  mpz_mul (skew, skew, skew);
  mpz_mul (b_Z, b_Z, skew);
  mpz_mul (res, a_Z, b_Z);
  for (unsigned int i = 2; i <= n; i++) {
    mpz_set_int64(a_Z, a[i]);
    mpz_set_int64(b_Z, b[i]);
    mpz_set_int64(skew, skewness->c[i - 1]);
    mpz_mul (skew, skew, skew);
    mpz_mul (b_Z, b_Z, skew);
    mpz_addmul(res, a_Z, b_Z);
  }
  ASSERT_ALWAYS(mpz_cmp(res, ret) == 0);
  mpz_clear(a_Z);
  mpz_clear(b_Z);
  mpz_clear(skew);
  mpz_clear(res);
  mpz_clear(ret);
#endif // ASSERT_LLL

  return 1;
}

// Return 0 if overflow, 1 otherwise.
static int innerproduct(int64_t * x, int64_t * a, int64_t * b,
    unsigned int n, FILE * errstd)
{
  * x = 0;
  if (__builtin_smull_overflow(a[1], b[1], x)) {
    fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
    fprintf(errstd, "# Overflow during innerproduct, ");
    return 0;
  }
  int64_t tmp = 0;
  for (unsigned int i = 2; i <= n; i++) {
    /*x = x + a[i] * b[i];*/
    if (__builtin_smull_overflow(a[i], b[i], &tmp) ||
        __builtin_saddl_overflow(* x, tmp, x)) {
      fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
      fprintf(errstd, "# Overflow during innerproduct, ");
      return 0;
    }
  }

#ifdef ASSERT_LLL
  mpz_t a_Z, b_Z, res, ret;
  mpz_init(a_Z);
  mpz_init(b_Z);
  mpz_init(res);
  mpz_init(ret);
  mpz_set_int64(ret, * x);
  mpz_set_int64(a_Z, a[1]);
  mpz_set_int64(b_Z, b[1]);

  mpz_mul (res, a_Z, b_Z);
  for (unsigned int i = 2; i <= n; i++) {
    mpz_set_int64(a_Z, a[i]);
    mpz_set_int64(b_Z, b[i]);
    mpz_addmul(res, a_Z, b_Z);
  }
  ASSERT_ALWAYS(mpz_cmp(res, ret) == 0);
  mpz_clear(a_Z);
  mpz_clear(b_Z);
  mpz_clear(res);
  mpz_clear(ret);
#endif // ASSERT_LLL

  return 1;
}

/* c = (x*c1 + y*c2)/z */
static int muladddiv(int64_t * tmp, int64_t c1, int64_t c2, int64_t x,
    int64_t y, int64_t z, FILE * errstd)
{
#ifdef OLD_LLL
  ASSERT(z != 0);

  double d_tmp0 = (double)x * (double)c1;
  double d_tmp1 = (double)y * (double)c2;
  

  if (ceil(log2(fabs(d_tmp0))) < 52.0 && ceil(log2(fabs(d_tmp1))) < 52.0) {

    tmp = ((int64_t) d_tmp0 + (int64_t) d_tmp1) / z;
  } else if (ceil(log2(fabs(d_tmp0))) < 62.0
      && ceil(log2(fabs(d_tmp1))) < 62.0) {

    return (x * c1 + y * c2) / z;
  } else {
    tmp = 0;
    mpz_t c1_Z, c2_Z, x_Z, y_Z, z_Z, c, t1;
    mpz_init(c1_Z);
    mpz_init(c2_Z);
    mpz_init(x_Z);
    mpz_init(y_Z);
    mpz_init(z_Z);
    mpz_init(c);
    mpz_init(t1);

    mpz_set_int64(c1_Z, c1);
    mpz_set_int64(c2_Z, c2);
    mpz_set_int64(x_Z, x);
    mpz_set_int64(y_Z, y);
    mpz_set_int64(z_Z, z);

    mpz_mul(t1, x_Z, c1_Z);
    mpz_addmul(t1, y_Z, c2_Z);
    mpz_divexact(c, t1, z_Z);

    mpz_clear(c1_Z);
    mpz_clear(c2_Z);
    mpz_clear(x_Z);
    mpz_clear(y_Z);
    mpz_clear(z_Z);
    mpz_clear(t1);

    ASSERT(mpz_sizeinbase(c, 2) < 63);
    tmp = mpz_get_int64(c);
    mpz_clear(c);
  }
#else // OLD_LLL
  int64_t tmpb = 0;
  if (__builtin_smull_overflow(c1, x, tmp) ||
      __builtin_smull_overflow(c2, y, &tmpb) ||
      __builtin_saddl_overflow(* tmp, tmpb, tmp)) {
    __int128_t c1_128, c2_128, x_128, y_128, z_128;
    c1_128 = (__int128_t) c1;
    c2_128 = (__int128_t) c2;
    x_128 = (__int128_t) x;
    y_128 = (__int128_t) y;
    z_128 = (__int128_t) z;

    if(!int128_convert_64(tmp, (x_128 * c1_128 + y_128 * c2_128) /
        z_128)) {
      fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
      fprintf(errstd, "# Overflow during muladddiv, ");
      return 0;
    }
  } else {
    * tmp = * tmp / z;
  }
#endif // OLD_LLL

#ifdef ASSERT_LLL
  mpz_t c1_Z, c2_Z, x_Z, y_Z, z_Z, c, t1;
  mpz_init(c1_Z);
  mpz_init(c2_Z);
  mpz_init(x_Z);
  mpz_init(y_Z);
  mpz_init(z_Z);
  mpz_init(c);
  mpz_init(t1);

  mpz_set_int64(c1_Z, c1);
  mpz_set_int64(c2_Z, c2);
  mpz_set_int64(x_Z, x);
  mpz_set_int64(y_Z, y);
  mpz_set_int64(z_Z, z);

  mpz_mul(t1, x_Z, c1_Z);
  mpz_addmul(t1, y_Z, c2_Z);
  mpz_divexact(c, t1, z_Z);

  mpz_clear(c1_Z);
  mpz_clear(c2_Z);
  mpz_clear(x_Z);
  mpz_clear(y_Z);
  mpz_clear(z_Z);

  mpz_set_int64(t1, * tmp);
  ASSERT_ALWAYS(mpz_cmp(t1, c) == 0);

  mpz_clear(t1);
  mpz_clear(c);
#endif // ASSERT_LLL

  return 1;
}

/* c = (x*c1 - y*c2)/z */
static int mulsubdiv(int64_t * tmp, int64_t c1, int64_t c2, int64_t x,
    int64_t y, int64_t z, FILE * errstd)
{
  ASSERT(z != 0);

#ifdef OLD_LLL
  double d_tmp0 = (double)x * (double)c1;
  double d_tmp1 = (double)y * (double)c2;

  if (ceil(log2(fabs(d_tmp0))) < 52.0 && ceil(log2(fabs(d_tmp1))) < 52.0) {
    tmp = ((int64_t) d_tmp0 - (int64_t) d_tmp1) / z;
  } else if (ceil(log2(fabs(d_tmp0))) < 62.0 &&
      ceil(log2(fabs(d_tmp1))) < 62.0) {
    tmp = (x * c1 - y * c2) / z;
  } else {
    mpz_t c1_Z, c2_Z, x_Z, y_Z, z_Z, c, t1;
    mpz_init(c1_Z);
    mpz_init(c2_Z);
    mpz_init(x_Z);
    mpz_init(y_Z);
    mpz_init(z_Z);
    mpz_init(c);
    mpz_init(t1);

    mpz_set_int64(c1_Z, c1);
    mpz_set_int64(c2_Z, c2);
    mpz_set_int64(x_Z, x);
    mpz_set_int64(y_Z, y);
    mpz_set_int64(z_Z, z);

    mpz_mul(t1, x_Z, c1_Z);
    mpz_submul(t1, y_Z, c2_Z);
    mpz_divexact(c, t1, z_Z);

    mpz_clear(c1_Z);
    mpz_clear(c2_Z);
    mpz_clear(x_Z);
    mpz_clear(y_Z);
    mpz_clear(z_Z);
    mpz_clear(t1);

    ASSERT(mpz_sizeinbase(c, 2) < 63);
    tmp = mpz_get_int64(c);
    mpz_clear(c);
  }
#else // OLD_LLL
  int64_t tmpb = 0;
  if (__builtin_smull_overflow(c1, x, tmp) ||
      __builtin_smull_overflow(c2, y, &tmpb) ||
      __builtin_ssubl_overflow(* tmp, tmpb, tmp)) {
    __int128_t c1_128, c2_128, x_128, y_128, z_128;
    c1_128 = (__int128_t) c1;
    c2_128 = (__int128_t) c2;
    x_128 = (__int128_t) x;
    y_128 = (__int128_t) y;
    z_128 = (__int128_t) z;

    if (!int128_convert_64(tmp, (x_128 * c1_128 - y_128 * c2_128) /
        z_128)) {
      fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
      fprintf(errstd, "# Overflow during mulsubdiv, ");
      return 0;
    }
  } else {
    * tmp = * tmp / z;
  }
#endif // OLD_LLL

#ifdef ASSERT_LLL
  mpz_t c1_Z, c2_Z, x_Z, y_Z, z_Z, c, t1;
  mpz_init(c1_Z);
  mpz_init(c2_Z);
  mpz_init(x_Z);
  mpz_init(y_Z);
  mpz_init(z_Z);
  mpz_init(c);
  mpz_init(t1);

  mpz_set_int64(c1_Z, c1);
  mpz_set_int64(c2_Z, c2);
  mpz_set_int64(x_Z, x);
  mpz_set_int64(y_Z, y);
  mpz_set_int64(z_Z, z);

  mpz_mul(t1, x_Z, c1_Z);
  mpz_submul(t1, y_Z, c2_Z);
  mpz_divexact(c, t1, z_Z);

  mpz_clear(c1_Z);
  mpz_clear(c2_Z);
  mpz_clear(x_Z);
  mpz_clear(y_Z);
  mpz_clear(z_Z);

  mpz_set_int64(t1, * tmp);
  ASSERT_ALWAYS(mpz_cmp(t1, c) == 0);

  mpz_clear(t1);
  mpz_clear(c);
#endif // ASSERT_LLL

  return 1;
}

static int skewincrementalgs_square(mat_int64_srcptr B, unsigned int * P,
    int64_t * D, int64_t ** lam, unsigned int * s, unsigned int k,
    int64_vector_srcptr skewness, FILE * errstd)
{
  unsigned int n = B->NumCols;
  int64_t u = 0;

  for (unsigned int j = 1; j <= k - 1; j++) {
    unsigned int posj = P[j];
    if (posj == 0) {
      continue;
    }

    if (!skewinnerproduct_square(&u, B->coeff[k], B->coeff[j], n, errstd,
          skewness)) {
      fprintf(errstd, "skewincrementalgs_square, ");
      return 0;
    }
    for (unsigned int i = 1; i <= posj - 1; i++) {
      ASSERT(D[i - 1] != 0);

      /*u = (D[i] * u - lam[k][i] * lam[j][i]) / D[i - 1];*/
      if (!mulsubdiv(&u, u, lam[j][i], D[i], lam[k][i], D[i - 1], errstd)) {
        fprintf(errstd, "skewincrementalgs_square, ");
        return 0;
      }
    }

    lam[k][posj] = u;
  }

  if (!skewinnerproduct_square(&u, B->coeff[k], B->coeff[k], n, errstd, skewness)) {
    fprintf(errstd, "skewincrementalgs_square, ");
    return 0;
  }

  for (unsigned int i = 1; i <= * s; i++) {
    ASSERT(D[i - 1] != 0);

    /*u = (D[i] * u - lam[k][i] * lam[k][i]) / D[i - 1];*/
    if (!mulsubdiv(&u, u, lam[k][i], D[i], lam[k][i], D[i - 1], errstd)) {
      fprintf(errstd, "skewincrementalgs_square, ");
      return 0;
    }
  }

  if (u == 0) {
    P[k] = 0;
  } else {
    * s = * s + 1;
    P[k] = * s;
    D[* s] = u;
  }

  return 1;
}

static int skewincrementalgs(mat_int64_srcptr B, unsigned int * P, int64_t * D,
    int64_t ** lam, unsigned int * s, unsigned int k,
    int64_vector_srcptr skewness, FILE * errstd)
{
  unsigned int n = B->NumCols;
  int64_t u = 0;

  for (unsigned int j = 1; j <= k - 1; j++) {
    unsigned int posj = P[j];
    if (posj == 0) {
      continue;
    }

    if (!skewinnerproduct(&u, B->coeff[k], B->coeff[j], n, errstd, skewness)) {
      fprintf(errstd, "skewincrementalgs, ");
      return 0;
    }
    for (unsigned int i = 1; i <= posj - 1; i++) {
      ASSERT(D[i - 1] != 0);

      /*u = (D[i] * u - lam[k][i] * lam[j][i]) / D[i - 1];*/
      if (!mulsubdiv(&u, u, lam[j][i], D[i], lam[k][i], D[i - 1], errstd)) {
        fprintf(errstd, "skewincrementalgs, ");
        return 0;
      }
    }

    lam[k][posj] = u;
  }

  if (!skewinnerproduct(&u, B->coeff[k], B->coeff[k], n, errstd, skewness)) {
    fprintf(errstd, "skewincrementalgs, ");
    return 0;
  }

  for (unsigned int i = 1; i <= * s; i++) {
    ASSERT(D[i - 1] != 0);

    /*u = (D[i] * u - lam[k][i] * lam[k][i]) / D[i - 1];*/
    if (!mulsubdiv(&u, u, lam[k][i], D[i], lam[k][i], D[i - 1], errstd)) {
      fprintf(errstd, "skewincrementalgs, ");
      return 0;
    }
  }

  if (u == 0) {
    P[k] = 0;
  } else {
    * s = * s + 1;
    P[k] = * s;
    D[* s] = u;
  }

  return 1;
}

static int incrementalgs(mat_int64_srcptr B, unsigned int * P, int64_t * D,
    int64_t ** lam, unsigned int * s, unsigned int k, FILE * errstd)
{
  unsigned int n = B->NumCols;
  int64_t u = 0;

  for (unsigned int j = 1; j <= k - 1; j++) {
    unsigned int posj = P[j];
    if (posj == 0) {
      continue;
    }

    if (!innerproduct(&u, B->coeff[k], B->coeff[j], n, errstd)) {
      fprintf(errstd, "incrementalgs, ");
      return 0;
    }
    for (unsigned int i = 1; i <= posj - 1; i++) {
      ASSERT(D[i - 1] != 0);

      /*u = (D[i] * u - lam[k][i] * lam[j][i]) / D[i - 1];*/
      if (!mulsubdiv(&u, u, lam[j][i], D[i], lam[k][i], D[i - 1], errstd)) {
        fprintf(errstd, "incrementalgs, ");
        return 0;
      }
    }

    lam[k][posj] = u;
  }

  if (!innerproduct(&u, B->coeff[k], B->coeff[k], n, errstd)) {
    fprintf(errstd, "incrementalgs, ");
    return 0;
  }

  for (unsigned int i = 1; i <= * s; i++) {
    ASSERT(D[i - 1] != 0);

    /*u = (D[i] * u - lam[k][i] * lam[k][i]) / D[i - 1];*/
    if (!mulsubdiv(&u, u, lam[k][i], D[i], lam[k][i], D[i - 1], errstd)) {
      fprintf(errstd, "incrementalgs, ");
      return 0;
    }
  }

  if (u == 0) {
    P[k] = 0;
  } else {
    * s = * s + 1;
    P[k] = * s;
    D[* s] = u;
  }

  return 1;
}

/*  rounds a/d to nearest integer, breaking ties
    by rounding towards zero.  Assumes d > 0. */
static int64_t baldiv(int64_t a, int64_t d)
{
  int64_t q = 0;
  int64_t r = 0;
  int64_fdiv_qr(&q, &r, a, d);
  r = r * 2;

  if (r > d || (r == d && q < 0)) {
    q = q + 1;
  }

  return q;
}

/* c0 = c0 - x*c1 */
static int mulsubn (int64_t * c0, int64_t * c1, int64_t x, unsigned int n,
    FILE * errstd)
{
#ifdef ASSERT_LLL
  mpz_t * c = (mpz_t *) malloc(sizeof(mpz_t) * n);
  mpz_t * c2 = (mpz_t *) malloc(sizeof(mpz_t) * n);
  for (unsigned int i = 1; i <= n; i++) {
    mpz_init(c[i - 1]);
    mpz_set_int64(c[i - 1], c0[i]);
    mpz_init(c2[i - 1]);
    mpz_set_int64(c2[i - 1], c1[i]);
  }
  mpz_t x_Z;
  mpz_init(x_Z);
  mpz_set_int64(x_Z, x);
  mpz_t tmp_Z;
  mpz_init(tmp_Z);
#endif // ASSERT_LLL

#ifdef OLD_LLL
  for (unsigned int i = 1; i <= n; i++) {
    ASSERT(log2(fabs((double)c0[i] - (double)c1[i] * (double)x)) < 63.0);
    c0[i] = c0[i] - c1[i] * x;
  }
#else // OLD_LLL
  int64_t tmp = 0;
  for (unsigned int i = 1; i <= n; i++) {
    if (__builtin_smull_overflow(x, c1[i], &tmp) ||
        __builtin_ssubl_overflow(c0[i], tmp, c0 + i)) {
      fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
      fprintf(errstd, "# Overflow during mulsubn, ");
      return 0;
    }
  }
#endif // OLD_LLL

#ifdef ASSERT_LLL
  unsigned int i;
  signed long int x0;

  x0 = mpz_get_si (x_Z);
  if (mpz_cmp_si (x_Z, x0) == 0 && 
      x0 != ((signed long int) 1 << (mp_bits_per_limb - 1))) {
    if (x0 > 0) {
      for (i = 0; i < n; i++) {
        mpz_mul_ui (tmp_Z, c2[i], x0);
        mpz_sub (c[i], c[i], tmp_Z);
      }
    } else if (x0 < 0) {
      x0 = -x0;
      for (i = 0; i < n; i++)
        mpz_addmul_ui (c[i], c2[i], x0);
    }
  } else {
    for (i = 0; i < n; i++) {
      mpz_mul (tmp_Z, c2[i], x_Z);
      mpz_sub (c[i], c[i], tmp_Z);
    }
  }

  for (i = 1; i <= n; i++) {
    mpz_set_int64(tmp_Z, c0[i]);
    ASSERT_ALWAYS(mpz_cmp(tmp_Z, c[i - 1]) == 0);
  }

  for (unsigned int i = 0; i < n; i++) {
    mpz_clear(c[i]);
    mpz_clear(c2[i]);
  }
  mpz_clear(x_Z);
  mpz_clear(tmp_Z);
  free(c);
  free(c2);
#endif // ASSERT_LLL

  return 1;
}

static int reduce(unsigned int k, unsigned int l, mat_int64_ptr B,
    unsigned int * P, int64_t * D, int64_t ** lam, mat_int64_ptr U,
    FILE * errstd)
{
  if (P[l] == 0) {
    return 1;
  }

  int64_t t1 = lam[k][P[l]] * 2;
  int64_t r = 0;
  t1 = ABS(t1);
  if (t1 <= D[P[l]]) {
    return 1;
  }

  r = baldiv(lam[k][P[l]], D[P[l]]);
  if (!mulsubn(B->coeff[k], B->coeff[l], r, B->NumCols, errstd)) {
    fprintf(errstd, "reduce, ");
    return 0;
  }

  if (U != NULL) {
    if (!mulsubn(U->coeff[k], U->coeff[l], r, B->NumRows, errstd)) {
      fprintf(errstd, "reduce, ");
      return 0;
    }
  }

  for (unsigned int j = 1; j <= l - 1; j++) {
    if (P[j] != 0) {
#ifdef ASSERT_LLL
      mpz_t lamkPj;
      mpz_init(lamkPj);
      mpz_set_int64(lamkPj, lam[k][P[j]]);
      mpz_t lamlPj;
      mpz_init(lamlPj);
      mpz_set_int64(lamlPj, lam[l][P[j]]);
      mpz_t r_Z;
      mpz_init(r_Z);
      mpz_set_int64(r_Z, r);
#endif // ASSERT_LLL
#ifdef OLD_LLL
      lam[k][P[j]] = lam[k][P[j]] - lam[l][P[j]] * r;
#else // OLD_LLL
      int64_t tmp = 0;
      if (__builtin_smull_overflow(r, lam[l][P[j]], &tmp) ||
          __builtin_ssubl_overflow(lam[k][P[j]] , tmp, lam[k] + P[j])) {
        fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
        fprintf(errstd, "# Overflow during reduce, ");
        return 0;
      }
#endif // OLD_LLL
#ifdef ASSERT_LLL
      mpz_submul(lamkPj, lamlPj, r_Z);
      mpz_set_int64(r_Z, lam[k][P[j]]);
      ASSERT_ALWAYS(mpz_cmp(r_Z, lamkPj) == 0);
      mpz_clear(lamkPj);
      mpz_clear(lamlPj);
      mpz_clear(r_Z);
#endif // ASSERT_LLL
    }
  }

#ifdef ASSERT_LLL
  mpz_t lamkPl;
  mpz_init(lamkPl);
  mpz_set_int64(lamkPl, lam[k][P[l]]);
  mpz_t DPl;
  mpz_init(DPl);
  mpz_set_int64(DPl, D[P[l]]);
  mpz_t r_Z;
  mpz_init(r_Z);
  mpz_set_int64(r_Z, r);
#endif // ASSERT_LLL
#ifdef OLD_LLL
  lam[k][P[l]] = lam[k][P[l]] - D[P[l]] * r;
#else // OLD_LLL
  int64_t tmp = 0;
  if (__builtin_smull_overflow(r, D[P[l]], &tmp) ||
      __builtin_ssubl_overflow(lam[k][P[l]] , tmp, lam[k] + P[l])) {
    fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
    fprintf(errstd, "# Overflow during reduce, ");
    return 0;
  }
#endif // OLD_LLL
#ifdef ASSERT_LLL
  mpz_submul(lamkPl, DPl, r_Z);
  mpz_set_int64(r_Z, lam[k][P[l]]);
  ASSERT_ALWAYS(mpz_cmp(r_Z, lamkPl) == 0);
  mpz_clear(lamkPl);
  mpz_clear(DPl);
  mpz_clear(r_Z);
#endif // ASSERT_LLL

  return 1;
}

/* test if a*d1^2 > b*(d0*d2 + lam^2)
   t1 and t2 are temporary variables */
static int swaptest(int64_t d0, int64_t d1, int64_t d2, int64_t lam,
    int64_t a, int64_t b)
{
  double t1 = ((double)d0 * (double)d2 + (double)lam * (double)lam) *
    (double)b;
  double t2 = (double)d1 * (double)d1 * (double)a;

#ifdef ASSERT_LLL
  mpz_t d0_Z;
  mpz_init(d0_Z);
  mpz_set_int64(d0_Z, d0);
  mpz_t d1_Z;
  mpz_init(d1_Z);
  mpz_set_int64(d1_Z, d1);
  mpz_t d2_Z;
  mpz_init(d2_Z);
  mpz_set_int64(d2_Z, d2);
  mpz_t lam_Z;
  mpz_init(lam_Z);
  mpz_set_int64(lam_Z, lam);
  mpz_t a_Z;
  mpz_init(a_Z);
  mpz_set_int64(a_Z, a);
  mpz_t b_Z;
  mpz_init(b_Z);
  mpz_set_int64(b_Z, b);
  mpz_t t1_Z;
  mpz_init(t1_Z);
  mpz_t t2_Z;
  mpz_init(t2_Z);

  mpz_mul (t1_Z, d0_Z, d2_Z);
  mpz_mul (t2_Z, lam_Z, lam_Z);
  mpz_add (t1_Z, t1_Z, t2_Z);
  mpz_mul (t1_Z, t1_Z, b_Z);

  mpz_mul (t2_Z, d1_Z, d1_Z);
  mpz_mul (t2_Z, t2_Z, a_Z);

  ASSERT_ALWAYS((mpz_cmp (t2_Z, t1_Z) > 0) == (t2 > t1));

  mpz_clear(d0_Z);
  mpz_clear(d1_Z);
  mpz_clear(d2_Z);
  mpz_clear(lam_Z);
  mpz_clear(a_Z);
  mpz_clear(b_Z);
  mpz_clear(t1_Z);
  mpz_clear(t2_Z);
#endif // ASSERT_LLL

  return (t2 > t1);
}

/* (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2) */
static int rowtransform(int64_t * c1, int64_t * c2, int64_t x, int64_t y,
    int64_t u, int64_t v, FILE * errstd)
{
#ifdef ASSERT_LLL
  mpz_t t1_Z, t2_Z;
  mpz_t c1_Z, c2_Z;
  mpz_t x_Z, y_Z;
  mpz_t u_Z, v_Z;

  mpz_init(c1_Z);
  mpz_init(c2_Z);
  mpz_init(x_Z);
  mpz_init(y_Z);
  mpz_init(u_Z);
  mpz_init(v_Z);
  mpz_set_int64(c1_Z, * c1);
  mpz_set_int64(c2_Z, * c2);
  mpz_set_int64(x_Z, x);
  mpz_set_int64(y_Z, y);
  mpz_set_int64(u_Z, u);
  mpz_set_int64(v_Z, v);

  mpz_init(t1_Z);
  mpz_init(t2_Z);

  mpz_mul(t1_Z, x_Z, c1_Z);
  mpz_mul(t2_Z, y_Z, c2_Z);
  mpz_add(t1_Z, t1_Z, t2_Z);

  mpz_mul(t2_Z, u_Z, c1_Z);
  mpz_set(c1_Z, t1_Z);
  mpz_mul(t1_Z, v_Z, c2_Z);
  mpz_add(c2_Z, t1_Z, t2_Z);

  mpz_clear(t2_Z);
  mpz_clear(x_Z);
  mpz_clear(y_Z);
  mpz_clear(u_Z);
  mpz_clear(v_Z);
#endif // ASSERT_LLL

#ifdef OLD_LLL
  double d_tmp0 = (double)x * (double) (* c1);
  double d_tmp1 = (double)y * (double) (* c2);

  if (ceil(log2(fabs(d_tmp0))) > 62.0 &&
      ceil(log2(fabs(d_tmp1))) > 62.0) {
    printf("Overflow.\n");
    ASSERT_ALWAYS(0);
  }

  d_tmp0 = (double)u * (double) (* c1);
  d_tmp1 = (double)v * (double) (* c2);

  if (ceil(log2(fabs(d_tmp0))) > 62.0 &&
      ceil(log2(fabs(d_tmp1))) > 62.0) {
    printf("Overflow.\n");
    ASSERT_ALWAYS(0);
  }

  int64_t t1 = x * * c1 + y * * c2;
  int64_t t2 = u * * c1;
  * c1 = t1;
  t1 = v * * c2;
  * c2 = t1 + t2;
#else // OLD_LLL
  int64_t xc1 = 0;
  int64_t yc2 = 0;
  if (__builtin_smull_overflow(x, * c1, &xc1) ||
      __builtin_smull_overflow(y, * c2, &yc2)) {
    fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
    fprintf(errstd, "# Overflow during rowtransform, ");
    return 0;
  }

  int64_t uc1 = 0;
  int64_t vc2 = 0;
  if (__builtin_smull_overflow(u, * c1, &uc1) ||
      __builtin_smull_overflow(v, * c2, &vc2)) {
    fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
    fprintf(errstd, "# Overflow during rowtransform, ");
    return 0;
  }

  int64_t t1 = 0;
  if (__builtin_saddl_overflow(xc1, yc2, &t1)) {
    fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
    fprintf(errstd, "# Overflow during rowtransform, ");
    return 0;
  }
  int64_t t2 = uc1;
  * c1 = t1;
  t1 = vc2;
  if (__builtin_saddl_overflow(t1, t2, c2)) {
    fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
    fprintf(errstd, "# Overflow during rowtransform, ");
    return 0;
  }
#endif // OLD_LLL

#ifdef ASSERT_LLL
  mpz_set_int64(t1_Z, *c1);
  ASSERT(mpz_cmp(t1_Z, c1_Z) == 0);
  mpz_set_int64(t1_Z, *c2);
  ASSERT(mpz_cmp(t1_Z, c2_Z) == 0);
  mpz_clear(c1_Z);
  mpz_clear(c2_Z);
  mpz_clear(t1_Z);
#endif // ASSERT_LLL

  return 1;
}

/* (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2) */
static int rowtransformn(int64_t * c1, int64_t * c2, int64_t x, int64_t y,
    int64_t u, int64_t v, unsigned int n, FILE * errstd)
{
  for (unsigned int i = 1; i <= n; i++) {
#ifdef ASSERT_LLL
    mpz_t t1_Z, t2_Z;
    mpz_t c1_Z, c2_Z, x_Z, y_Z, u_Z, v_Z;
    mpz_init(t1_Z);
    mpz_init(t2_Z);
    mpz_init(c1_Z);
    mpz_init(c2_Z);
    mpz_init(x_Z);
    mpz_init(y_Z);
    mpz_init(u_Z);
    mpz_init(v_Z);
    mpz_set_int64(c1_Z, c1[i]);
    mpz_set_int64(c2_Z, c2[i]);
    mpz_set_int64(x_Z, x);
    mpz_set_int64(y_Z, y);
    mpz_set_int64(u_Z, u);
    mpz_set_int64(v_Z, v);
    mpz_mul(t1_Z, x_Z, c1_Z);
    mpz_mul(t2_Z, y_Z, c2_Z);
    mpz_add(t1_Z, t1_Z, t2_Z);

    mpz_mul(t2_Z, u_Z, c1_Z);
    mpz_set(c1_Z, t1_Z);
    mpz_mul(t1_Z, v_Z, c2_Z);
    mpz_add(c2_Z, t1_Z, t2_Z);
#endif // ASSERT_LLL

#ifdef OLD_LLL
    double d_tmp0 = (double)x * (double) (c1[i]);
    double d_tmp1 = (double)y * (double) (c2[i]);

    if (ceil(log2(fabs(d_tmp0))) > 62.0 &&
        ceil(log2(fabs(d_tmp1))) > 62.0) {
      printf("Overflow.\n");
      ASSERT_ALWAYS(0);
    }

    d_tmp0 = (double)u * (double) (c1[i]);
    d_tmp1 = (double)v * (double) (c2[i]);

    if (ceil(log2(fabs(d_tmp0))) > 62.0 &&
        ceil(log2(fabs(d_tmp1))) > 62.0) {
      printf("Overflow.\n");
      ASSERT_ALWAYS(0);
    }
    int64_t t1 = x * c1[i] + y * c2[i];
    int64_t t2 = u * c1[i];
    c1[i] = t1;
    t1 = v * c2[i];
    c2[i] = t1 + t2;
#else // OLD_LLL
    int64_t xc1 = 0;
    int64_t yc2 = 0;
    if (__builtin_smull_overflow(x, c1[i], &xc1) ||
        __builtin_smull_overflow(y, c2[i], &yc2)) {
      fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
      fprintf(errstd, "# Overflow during rowtransformn, ");
      return 0;
    }

    int64_t uc1 = 0;
    int64_t vc2 = 0;
    if (__builtin_smull_overflow(u, c1[i], &uc1) ||
        __builtin_smull_overflow(v, c2[i], &vc2)) {
      fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
      fprintf(errstd, "# Overflow during rowtransformn, ");
      return 0;
    }

    int64_t t1 = 0;
    if (__builtin_saddl_overflow(xc1, yc2, &t1)) {
      fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
      fprintf(errstd, "# Overflow during rowtransformn, ");
      return 0;
    }
    int64_t t2 = uc1;
    c1[i] = t1;
    t1 = vc2;
    if (__builtin_saddl_overflow(t1, t2, c2 + i)) {
      fprintf(errstd, "# ----- LLL int64, fall back to mpz LLL -----\n");
      fprintf(errstd, "# Overflow during rowtransformn, ");
      return 0;
    }
#endif // OLD_LLL

#ifdef ASSERT_LLL
    mpz_clear(t2_Z);
    mpz_clear(x_Z);
    mpz_clear(y_Z);
    mpz_clear(u_Z);
    mpz_clear(v_Z);
    mpz_set_int64(t1_Z, c1[i]);
    ASSERT(mpz_cmp(t1_Z, c1_Z) == 0);
    mpz_set_int64(t1_Z, c2[i]);
    ASSERT(mpz_cmp(t1_Z, c2_Z) == 0);
    mpz_clear(c1_Z);
    mpz_clear(c2_Z);
    mpz_clear(t1_Z);
#endif // ASSERT_LLL
  }

  return 1;
}

/* swaps vectors k-1 and k;  assumes P(k-1) != 0 */
static int swaplll (unsigned int k, mat_int64_ptr B, unsigned int * P,
    int64_t * D, int64_t ** lam, mat_int64_ptr U, unsigned int m, FILE * errstd)
{
  unsigned int i, j;
  int64_t t1, t2, t3, e, x, y;

  if (P[k] != 0) {
    int64_swap_n(B->coeff[k - 1], B->coeff[k]);
    if (U != NULL) {
      int64_swap_n(U->coeff[k - 1], U->coeff[k]);
    }

    for (unsigned int j = 1; j <= k - 2; j++) {
      if (P[j] != 0) {
        swap_int64(&(lam[k - 1][P[j]]), &(lam[k][P[j]]));
      }
    }

    for (unsigned int i = k + 1; i <= m; i++) {
      if (!muladddiv(&t1, lam[i][P[k] - 1], lam[i][P[k]],
          lam[k][P[k] - 1], D[P[k] - 2], D[P[k] - 1], errstd)) {
        fprintf(errstd, "swaplll, ");
        return 0;
      }
      if (!mulsubdiv(lam[i] + P[k], lam[i][P[k] - 1], lam[i][P[k]], 
          D[P[k]], lam[k][P[k] - 1], D[P[k] - 1], errstd)) {
        fprintf(errstd, "swaplll, ");
        return 0;
      }
      lam[i][P[k] - 1] = t1;
    }

    if (!muladddiv(D + (P[k] - 1), D[P[k]], lam[k][P[k] - 1],
        D[P[k] - 2], lam[k][P[k] - 1], D[P[k] - 1], errstd)) {
      fprintf(errstd, "swaplll, ");
      return 0;
    }
  } else if (lam[k][P[k - 1]] != 0) {
    int64_gcdext(&e, &x, &y, lam[k][P[k - 1]], D[P[k - 1]]);

    t1 = lam[k][P[k - 1]] / e;
    t2 = D[P[k - 1]] / e;

    t3 = t2;
    t2 = -t2;
    if (!rowtransformn(B->coeff[k - 1], B->coeff[k], t1, t2, y, x, B->NumCols,
          errstd))
    {
      fprintf(errstd, "swaplll, ");
      return 0;
    }
    if (U != NULL) {
      if (!rowtransformn(U->coeff[k - 1], U->coeff[k], t1, t2, y, x,
            B->NumCols, errstd)) {
        fprintf(errstd, "swaplll, ");
        return 0;
      }
    }
    for (unsigned j = 1; j <= k - 2; j++) {
      if (P[j] != 0) {
        if (!rowtransform(&(lam[k - 1][P[j]]), &(lam[k][P[j]]), t1, t2, y, x,
              errstd)) {
          fprintf(errstd, "swaplll, ");
          return 0;
        }
      }
    }

    if (__builtin_smull_overflow(t2, t2, &t2)) {
      return 0;
    }

    D[P[k - 1]] = D[P[k - 1]] / t2;

    for (unsigned i = k + 1; i <= m; i++) {
      if (P[i] != 0) {
        D[P[i]] = D[P[i]] / t2;
        for (j = i + 1; j <= m; j++) {
          lam[j][P[i]] = lam[j][P[i]] / t2;
        }
      }
    }
    for (i = k + 1; i <= m; i++) {
      lam[i][P[k - 1]] = lam[i][P[k - 1]] / t3;
    }

    swap(P[k - 1], P[k]);
  } else {
    int64_swap_n(B->coeff[k - 1], B->coeff[k]);
    if (U != NULL) {
      int64_swap_n(U->coeff[k - 1], U->coeff[k]);
    }

    for (j = 1; j <= k - 2; j++) {
      if (P[j] != 0) {
        swap_int64(&(lam[k - 1][P[j]]), &(lam[k][P[j]]));
      }
    }

    swap(P[k - 1], P[k]);
  }

  return 1;
}

int slll_square(unsigned int * s, int64_t * det, mat_int64_ptr B,
    mat_int64_ptr U, int64_t a, int64_t b, int64_vector_srcptr skewness,
    FILE * errstd)
{
  unsigned int m = B->NumRows;
  unsigned int n = B->NumCols;
  ASSERT_ALWAYS(n >= m);

  unsigned int * P = (unsigned int *) malloc((m + 1) * sizeof(unsigned int));

  int64_t * D = (int64_t *) malloc((m + 1) * sizeof(int64_t));
  for (unsigned int j = 0; j <= m; j++) {
    D[j] = (j == 0);
  }

  int64_t ** lam = (int64_t **) malloc((m + 1) * sizeof(int64_t *));
  for (unsigned int j = 0; j <= m; j++) {
    lam[j] = (int64_t *) malloc((m + 1) * sizeof(int64_t));
    for (unsigned int k = 0; k <= m; k++) {
      lam[j][k] = 0;
    }
  }

  if (U != NULL) {
    ASSERT(U->NumRows == m);
    ASSERT(U->NumCols == m);

    mat_int64_set_identity(U);
  }

  * s = 0;

  unsigned int k = 1;
  unsigned int max_k = 0;

  while (k <= m) {
    if (k > max_k) {
      if (!skewincrementalgs_square(B, P, D, lam, s, k, skewness, errstd)) {
        fprintf(errstd, "slll_square.\n");
        free(D);
        for (unsigned int j = 0; j <= m; j++) {
          free (lam[j]);
        }
        free (lam);

        free(P);
        return 0;
      }
      max_k = k;
    }

    if (k == 1) {
      k++;
      continue;
    }

    if (!reduce(k, k - 1, B, P, D, lam, U, errstd)) {
      fprintf(errstd, "slll_square.\n");
      free(D);
      for (unsigned int j = 0; j <= m; j++) {
        free (lam[j]);
      }
      free (lam);

      free(P);
      return 0;
    }

    if (P[k - 1] != 0 && (P[k] == 0 || 
          swaptest(D[P[k]], D[P[k] - 1], D[P[k] - 2],
            lam[k][P[k] - 1], a, b))) {
      if (!swaplll(k, B, P, D, lam, U, max_k, errstd)) {
        fprintf(errstd, "slll_square.\n");
        free(D);
        for (unsigned int j = 0; j <= m; j++) {
          free (lam[j]);
        }
        free (lam);

        free(P);
        return 0;
      }
      k--;
    } else {
      for (unsigned int j = k - 2; j >= 1; j--) {
        if (!reduce(k, j, B, P, D, lam, U, errstd)) {
          fprintf(errstd, "slll_square.\n");
          free(D);
          for (unsigned int j = 0; j <= m; j++) {
            free (lam[j]);
          }
          free (lam);

          free(P);
          return 0;
        }
      }
      k++;
    }
  }

  * det = D[* s];
  free(D);
  for (unsigned int j = 0; j <= m; j++) {
    free (lam[j]);
  }
  free (lam);

  free(P);

  return 1;
}

int slll(unsigned int * s, int64_t * det, mat_int64_ptr B,
    mat_int64_ptr U, int64_t a, int64_t b, int64_vector_srcptr skewness,
    FILE * errstd)
{
  unsigned int m = B->NumRows;
  unsigned int n = B->NumCols;
  ASSERT_ALWAYS(n >= m);

  unsigned int * P = (unsigned int *) malloc((m + 1) * sizeof(unsigned int));

  int64_t * D = (int64_t *) malloc((m + 1) * sizeof(int64_t));
  for (unsigned int j = 0; j <= m; j++) {
    D[j] = (j == 0);
  }

  int64_t ** lam = (int64_t **) malloc((m + 1) * sizeof(int64_t *));
  for (unsigned int j = 0; j <= m; j++) {
    lam[j] = (int64_t *) malloc((m + 1) * sizeof(int64_t));
    for (unsigned int k = 0; k <= m; k++) {
      lam[j][k] = 0;
    }
  }

  if (U != NULL) {
    ASSERT(U->NumRows == m);
    ASSERT(U->NumCols == m);

    mat_int64_set_identity(U);
  }

  * s = 0;

  unsigned int k = 1;
  unsigned int max_k = 0;

  while (k <= m) {
    if (k > max_k) {
      if (!skewincrementalgs(B, P, D, lam, s, k, skewness, errstd)) {
        fprintf(errstd, "slll.\n");
        free(D);
        for (unsigned int j = 0; j <= m; j++) {
          free (lam[j]);
        }
        free (lam);

        free(P);
        return 0;
      }
      max_k = k;
    }

    if (k == 1) {
      k++;
      continue;
    }

    if (!reduce(k, k - 1, B, P, D, lam, U, errstd)) {
      fprintf(errstd, "slll\n");
      free(D);
      for (unsigned int j = 0; j <= m; j++) {
        free (lam[j]);
      }
      free (lam);

      free(P);
      return 0;
    }

    if (P[k - 1] != 0 && (P[k] == 0 || 
          swaptest(D[P[k]], D[P[k] - 1], D[P[k] - 2],
            lam[k][P[k] - 1], a, b))) {
      if (!swaplll(k, B, P, D, lam, U, max_k, errstd)) {
        fprintf(errstd, "slll.\n");
        free(D);
        for (unsigned int j = 0; j <= m; j++) {
          free (lam[j]);
        }
        free (lam);

        free(P);
        return 0;
      }
      k--;
    } else {
      for (unsigned int j = k - 2; j >= 1; j--) {
        if (!reduce(k, j, B, P, D, lam, U, errstd)) {
          fprintf(errstd, "slll.\n");
          free(D);
          for (unsigned int j = 0; j <= m; j++) {
            free (lam[j]);
          }
          free (lam);

          free(P);
          return 0;
        }
      }
      k++;
    }
  }

  * det = D[* s];
  free(D);
  for (unsigned int j = 0; j <= m; j++) {
    free (lam[j]);
  }
  free (lam);

  free(P);

  return 1;
}

/* LLL-reduce the matrix B (whose rows represent vectors, with indices
   starting at 1):
 * det (output) is the determinant
 * U (output) is the transformation matrix (NULL if not needed)
 * a, b are parameters (delta = a/b = 3/4 classically, we must have
   1/4 < delta < 1, the closer delta is from 1, the better is the reduction)
 m is the number of vectors (i.e., number of rows)
 n is the number of columns (i.e., length of each vector)
 */
int lll(unsigned int * s, int64_t * det, mat_int64_ptr B,
    mat_int64_ptr U, int64_t a, int64_t b, FILE * errstd)
{
  unsigned int m = B->NumRows;
  unsigned int n = B->NumCols;
  ASSERT_ALWAYS(n >= m);

  unsigned int * P = (unsigned int *) malloc((m + 1) * sizeof(unsigned int));

  int64_t * D = (int64_t *) malloc((m + 1) * sizeof(int64_t));
  for (unsigned int j = 0; j <= m; j++) {
    D[j] = (j == 0);
  }

  int64_t ** lam = (int64_t **) malloc((m + 1) * sizeof(int64_t *));
  for (unsigned int j = 0; j <= m; j++) {
    lam[j] = (int64_t *) malloc((m + 1) * sizeof(int64_t));
    for (unsigned int k = 0; k <= m; k++) {
      lam[j][k] = 0;
    }
  }

  if (U != NULL) {
    ASSERT(U->NumRows == m);
    ASSERT(U->NumCols == m);

    mat_int64_set_identity(U);
  }

  * s = 0;

  unsigned int k = 1;
  unsigned int max_k = 0;

  while (k <= m) {
    if (k > max_k) {
      if (!incrementalgs(B, P, D, lam, s, k, errstd)) {
        fprintf(errstd, "lll.\n");
        free(D);
        for (unsigned int j = 0; j <= m; j++) {
          free (lam[j]);
        }
        free (lam);

        free(P);
        return 0;
      }
      max_k = k;
    }

    if (k == 1) {
      k++;
      continue;
    }

    if (!reduce(k, k - 1, B, P, D, lam, U, errstd)) {
      fprintf(errstd, "lll.\n");
      free(D);
      for (unsigned int j = 0; j <= m; j++) {
        free (lam[j]);
      }
      free (lam);

      free(P);
      return 0;
    }

    if (P[k - 1] != 0 && (P[k] == 0 || 
          swaptest(D[P[k]], D[P[k] - 1], D[P[k] - 2],
            lam[k][P[k] - 1], a, b))) {
      if (!swaplll(k, B, P, D, lam, U, max_k, errstd)) {
        fprintf(errstd, "lll.\n");
        free(D);
        for (unsigned int j = 0; j <= m; j++) {
          free (lam[j]);
        }
        free (lam);

        free(P);
        return 0;
      }
      k--;
    } else {
      for (unsigned int j = k - 2; j >= 1; j--) {
        if (!reduce(k, j, B, P, D, lam, U, errstd)) {
          fprintf(errstd, "lll.\n");
          free(D);
          for (unsigned int j = 0; j <= m; j++) {
            free (lam[j]);
          }
          free (lam);

          free(P);
          return 0;
        }
      }
      k++;
    }
  }

  * det = D[* s];
  free(D);
  for (unsigned int j = 0; j <= m; j++) {
    free (lam[j]);
  }
  free (lam);

  free(P);

  return 1;
}

void slll_Mqr(mat_int64_ptr C, mat_int64_srcptr A,
    int64_vector_srcptr skewness, FILE * errstd)
{
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(A->NumCols == C->NumCols);

  mat_int64_transpose(C, A);
  int64_t a = 3;
  int64_t b = 4;
  int64_t det = 0;
  unsigned int s = 0;

  if (!slll_square(&s, &det, C, NULL, a, b, skewness, errstd)) {
    fprintf(errstd, "# Mqr=\n");
    mat_int64_fprintf_comment(errstd, A);
    fprintf(errstd, "# -----\n");

    int64_vector_t skew;
    int64_vector_init(skew, skewness->dim);
    for (unsigned int i = 0; i < skewness->dim; i++) {
      skew->c[i] = (int64_t) sqrt((double)skewness->c[i]);
    }

    mat_int64_t I_s;
    mat_int64_init(I_s, A->NumRows, A->NumCols);
    mat_int64_set_diag(I_s, skew);

    mat_int64_t M;
    mat_int64_init(M, A->NumRows, A->NumCols);
    mat_int64_mul_mat_int64(M, I_s, A);


    mat_int64_t U;
    mat_int64_init(U, A->NumRows, A->NumCols);
    mat_int64_LLL_unimodular_transpose(U, M);

    mat_int64_mul_mat_int64(C, A, U);

    mat_int64_clear(U);
    mat_int64_clear(I_s);
    mat_int64_clear(M);
  } else {
    mat_int64_transpose(C, C);
  }


/*#ifndef NDEBUG*/
  /*mat_int64_t C_tmp;*/
  /*mat_int64_init(C_tmp, C->NumRows, C->NumCols);*/
  /*mat_int64_set_zero(C_tmp);*/
  /*mat_int64_LLL_transpose(C_tmp, A);*/
  /*ASSERT(mat_int64_equal(C_tmp, C));*/
  /*mat_int64_clear(C_tmp);*/
/*#endif // NDEBUG*/
}

void lll_Mqr(mat_int64_ptr C, mat_int64_srcptr A, FILE * errstd)
{
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(A->NumCols == C->NumCols);

  mat_int64_transpose(C, A);
  int64_t a = 3;
  int64_t b = 4;
  int64_t det = 0;
  unsigned int s = 0;

  if (!lll(&s, &det, C, NULL, a, b, errstd)) {
    fprintf(errstd, "# Mqr=\n");
    mat_int64_fprintf_comment(errstd, A);
    fprintf(errstd, "# -----\n");
    mat_int64_LLL_transpose(C, A);
  } else {
    mat_int64_transpose(C, C);
  }


#ifndef NDEBUG
  mat_int64_t C_tmp;
  mat_int64_init(C_tmp, C->NumRows, C->NumCols);
  mat_int64_set_zero(C_tmp);
  mat_int64_LLL_transpose(C_tmp, A);
  ASSERT(mat_int64_equal(C_tmp, C));
  mat_int64_clear(C_tmp);
#endif // NDEBUG
}

void lll_Mqr_unimodular(mat_int64_ptr U, mat_int64_srcptr A, FILE * errstd)
{
  ASSERT(A->NumRows == U->NumRows);
  ASSERT(A->NumCols == U->NumCols);

  mat_int64_t C;
  mat_int64_init(C, A->NumRows, A->NumCols);
  mat_int64_transpose(C, A);
  int64_t a = 3;
  int64_t b = 4;
  int64_t det = 0;
  unsigned int s = 0;

  if (!lll(&s, &det, C, U, a, b, errstd)) {
    fprintf(errstd, "# Mqr =\n");
    mat_int64_fprintf_comment(errstd, A);
    fprintf(errstd, "# -----\n");
    mat_int64_LLL_unimodular_transpose(U, A);
  } else {
    mat_int64_transpose(U, U);
  }

  mat_int64_clear(C);

#ifndef NDEBUG
  mat_int64_t U_tmp;
  mat_int64_init(U_tmp, U->NumRows, U->NumCols);
  mat_int64_set_zero(U_tmp);
  mat_int64_LLL_unimodular_transpose(U_tmp, A);
  ASSERT(mat_int64_equal(U_tmp, U));
  mat_int64_clear(U_tmp);
#endif // NDEBUG
}
#endif // SLLL_SAFE

#ifdef MAIN_LLL_INT64
int main(int argc, char ** argv)
{
  if (argc < 3) {
    fprintf(stderr, "The matrix is defined as matrix([[argc], [0, 1, …, 0], …,"
      "[0, …, 0, 1]]).\n");
    fprintf(stderr, "Need at least two argument for lll_int64.\n");
    return EXIT_FAILURE;
  }
  mat_int64_t M;
  mat_int64_init(M, (unsigned int) argc - 1, (unsigned int) argc - 1);
  mat_int64_set_zero(M);
  for (unsigned int i = 1; i <= (unsigned int) argc - 1; i++) {
    sscanf(argv[i], "%" PRId64 "", &M->coeff[1][i]);
    if (i != 1) {
      M->coeff[i][i] = 1;
    }
  }

  printf("M = ");
  mat_int64_fprintf(stdout, M);

  mat_int64_t MLLL;
  mat_int64_init(MLLL, (unsigned int) argc - 1, (unsigned int) argc - 1);
  mat_int64_t MLLL_safe;
  mat_int64_init(MLLL_safe, (unsigned int) argc - 1, (unsigned int) argc - 1);

  lll_Mqr(MLLL, M, stderr);

  mat_int64_LLL_transpose(MLLL_safe, M);
  printf("MLLL = ");
  mat_int64_fprintf(stdout, MLLL);

  ASSERT(mat_int64_equal(MLLL, MLLL_safe));

  lll_Mqr_unimodular(MLLL, M, stderr);

  mat_int64_LLL_unimodular_transpose(MLLL_safe, M);

  ASSERT(mat_int64_equal(MLLL, MLLL_safe));

  printf("U = ");
  mat_int64_fprintf(stdout, MLLL);

  mat_int64_clear(M);
  mat_int64_clear(MLLL);
  mat_int64_clear(MLLL_safe);

  return EXIT_SUCCESS;
}
#endif // MAIN_LLL_INT64
