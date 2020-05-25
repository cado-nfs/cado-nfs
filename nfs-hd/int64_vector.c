#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "macros.h"
#include "sieving_bound.h"
#include "int64_vector.h"
#include <inttypes.h>
#include "gcd.h"

void int64_vector_init(int64_vector_ptr v, unsigned int d)  
{
  ASSERT(d > 0);

  v->dim = d;
  v->c = (int64_t *) malloc(sizeof(int64_t) * d);

  ASSERT_ALWAYS (v->c != NULL);
}

void int64_vector_clear(int64_vector_ptr v)
{
  free(v->c);
}

void int64_vector_swap(int64_vector_ptr v0, int64_vector_ptr v1)
{
  ASSERT(v0->dim == v1->dim);

  int64_t * tmp = v0->c;
  v0->c = v1->c;
  v1->c = tmp;
}

void int64_vector_set(int64_vector_ptr v, int64_vector_srcptr s)
{
  ASSERT(v->dim == s->dim);

  for (unsigned int i = 0; i < v->dim; i++) {
    v->c[i] = s->c[i];
  }
}

void int64_vector_set_zero(int64_vector_ptr v)
{
  memset(v->c, 0, v->dim * sizeof(int64_t));
}

int int64_vector_equal(int64_vector_srcptr a, int64_vector_srcptr b)
{
  int r = (a->dim > b->dim) - (b->dim > a->dim);
  if (r) return 0;
  for(int d = (int)(a->dim - 1); d >= 0 ; d--) {
    r = a->c[d] - b->c[d];
    if (r) return 0;
  }
  return 1;
}

void int64_vector_add(int64_vector_ptr a, int64_vector_srcptr b,
    int64_vector_srcptr c)
{
  ASSERT(a->dim == b->dim);
  ASSERT(c->dim == b->dim);

  for (unsigned int i = 0; i < a->dim; i++) {
    a->c[i] = b->c[i] + c->c[i];
  }
}

void int64_vector_sub(int64_vector_ptr a, int64_vector_srcptr b,
    int64_vector_srcptr c)
{
  ASSERT(a->dim == b->dim);
  ASSERT(c->dim == b->dim);

  for (unsigned int i = 0; i < a->dim; i++) {
    a->c[i] = b->c[i] - c->c[i];
  }
}

void int64_vector_mul(int64_vector_ptr a, int64_vector_srcptr b,
    int64_t c)
{
  ASSERT(a->dim == b->dim);

  for (unsigned int i = 0; i < a->dim; i++) {
    a->c[i] = c * b->c[i];
  }
}

void int64_vector_div(int64_vector_ptr a, int64_vector_srcptr b,
    int64_t c)
{
  ASSERT(a->dim == b->dim);

  for (unsigned int i = 0; i < a->dim; i++) {
    a->c[i] = b->c[i] / c;
  }
}

void int64_vector_addmul(int64_vector_ptr a, int64_vector_srcptr b,
    int64_vector_srcptr c, int64_t d)
{
  ASSERT(a->dim == b->dim);
  ASSERT(c->dim == b->dim);

  for (unsigned int i = 0; i < a->dim; i++) {
    a->c[i] = b->c[i] + d * c->c[i];
  }
}

void int64_vector_fprintf(FILE * file, int64_vector_srcptr v)
{
  fprintf(file, "[");
  for (unsigned int i = 0; i < v->dim - 1; i++)
    fprintf(file, "%" PRId64 ", ", v->c[i]);
  if (v->dim != 0) {
    fprintf(file, "%" PRId64 "]\n", v->c[v->dim - 1]);
  } else {
    fprintf(file, "]\n");
  }
}

void int64_vector_add_one(int64_vector_ptr v, sieving_bound_srcptr H)
{
  int64_vector_add_one_i(v, 0, H);
}

unsigned int int64_vector_add_one_i(int64_vector_ptr v, unsigned int i,
    sieving_bound_srcptr H)
{
  ASSERT(v->dim == H->t);
  ASSERT(i < v->dim);

#ifndef NDEBUG
  int64_t tmp = 0;
  int64_t tmp2 = 0;
  for (unsigned int j = 0; j < v->dim; j++) {
    tmp = tmp + v->c[j];
    tmp2 = tmp2 + H->h[j];
  }
  ASSERT(tmp <= tmp2);
#endif

  unsigned int k = i;
  while(k < v->dim) {
    if (v->c[k] == H->h[k] - 1) {
      int64_vector_setcoordinate(v, k, -(int64_t)H->h[k]);
      k++;
    } else {
      break;
    }
  }
  if (k < v->dim) {
    v->c[k]++;
  }

#ifndef NDEBUG
  for (unsigned int j = 0; j < v->dim - 1; j++) {
    tmp = v->c[j];
    ASSERT(tmp >= -(int64_t)H->h[j]);
    ASSERT(tmp < (int64_t)H->h[j]);
  }
  tmp = v->c[v->dim -1];
  ASSERT(tmp >= 0);
  ASSERT(tmp < H->h[v->dim - 1]);
#endif

  return k;
}

void int64_vector_to_mpz_vector(mpz_vector_ptr a, int64_vector_srcptr b)
{
  ASSERT(a->dim == b->dim);  

  for (unsigned int i = 0; i < b->dim; i++) {
    mpz_vector_setcoordinate_int64 (a, i, b->c[i]);
  }
}

uint64_t int64_vector_norml2sqr(int64_vector_srcptr a)
{
  ASSERT(a->dim >= 2);

  uint64_t norm = 0;
  for (unsigned int i = 0; i < a->dim; i++) {
    ASSERT(a->c[i] * a->c[i] >= 0);

    norm += (uint64_t)(a->c[i] * a->c[i]);
  }
  return norm;
}

double int64_vector_norml2(int64_vector_srcptr a)
{
  ASSERT(a->dim >= 2);

  return sqrt((double)int64_vector_norml2sqr(a));
}

int64_t int64_vector_dot_product(int64_vector_srcptr v0, int64_vector_srcptr v1)
{
  ASSERT(v0->dim == v1->dim);

  int64_t dot = 0;
  for (unsigned int i = 0; i < v0->dim; i++) {
    dot += v0->c[i] * v1->c[i];
  }
  return dot;
}

int int64_vector_in_sieving_region(int64_vector_srcptr v,
    sieving_bound_srcptr H)
{
  ASSERT(v->dim == H->t);

  for (unsigned int i = 0; i < v->dim - 1; i++) {
    if (-(int64_t)H->h[i] > v->c[i] || (int64_t)H->h[i] <= v->c[i]) {
      return 0;
    }
  }
  if (0 > v->c[v->dim - 1] || (int64_t)H->h[v->dim - 1] <= v->c[v->dim - 1]) {
    return 0;
  }
  return 1;
}

void double_vector_in_int64_vector(int64_vector_ptr v_i,
    double_vector_srcptr v_d)
{
  ASSERT(v_i->dim == v_d->dim);

  for (unsigned int i = 0; i < v_i->dim; i++) {
    v_d->c[i] = (int64_t) v_i->c[i];
  }
}

int64_t int64_vector_gcd(int64_vector_srcptr v)
{
  ASSERT(v->dim >= 2);

  int64_t gcd = gcd_int64(v->c[0], v->c[1]);
  for (unsigned int i = 2; i < v->dim; i++) {
    gcd = gcd_int64(gcd, v->c[i]);
  }
  return gcd;
}

void int64_vector_reduce(int64_vector_ptr u, int64_vector_srcptr v)
{
  ASSERT(v->dim >= 2);
  ASSERT(u->dim == v->dim);
  
  int64_t gcd = int64_vector_gcd(v);
  if (gcd != 0) {
    int64_vector_div(u, v, gcd);
  }
}

/*int64_t int64_vector_orthogonal_projection(int64_vector_ptr res,*/
    /*int64_vector_srcptr u, int64_vector_srcptr v)*/
/*{*/
  /*ASSERT(res->dim == u->dim);*/
  /*ASSERT(u->dim == v->dim);*/
  
  /*int64 coeff =*/
    /*int64_vector_dot_product(u, v);*/
  /*for (unsigned int i = 0; i < res->dim; i++) {*/
    /*res->c[i] = coeff * u->c[i];*/
  /*}*/
  /*return coeff;*/
/*}*/
