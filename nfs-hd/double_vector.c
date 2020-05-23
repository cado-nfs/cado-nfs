#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include "macros.h"
#include "sieving_bound.h"
#include "double_vector.h"

void double_vector_init(double_vector_ptr v, unsigned int d)  
{
  ASSERT(d > 0);

  v->dim = d;
  v->c = (double *) malloc(sizeof(double) * d);

  ASSERT_ALWAYS (v->c != NULL);
}

void double_vector_clear(double_vector_ptr v)
{
  free(v->c);
}

void double_vector_swap(double_vector_ptr v0, double_vector_ptr v1)
{
  ASSERT(v0->dim == v1->dim);

  double * tmp = v0->c;
  v0->c = v1->c;
  v1->c = tmp;
}

void double_vector_set(double_vector_ptr v, double_vector_srcptr s)
{
  ASSERT(v->dim == s->dim);

  for (unsigned int i = 0; i < v->dim; i++) {
    v->c[i] = s->c[i];
  }
}

void double_vector_set_zero(double_vector_ptr v)
{
  for (unsigned int i = 0; i < v->dim; i++) {
    v->c[i] = 0;
  }
}

int double_vector_equal(double_vector_srcptr a, double_vector_srcptr b)
{
  int r = (a->dim > b->dim) - (b->dim > a->dim);
  if (r) return 0;
  for(int d = a->dim; d >= 0 ; d--) {
    r = a->c[d] - b->c[d];
    if (r) return 0;
  }
  return 1;
}

void double_vector_add(double_vector_ptr a, double_vector_srcptr b,
    double_vector_srcptr c)
{
  ASSERT(a->dim == b->dim);
  ASSERT(c->dim == b->dim);

  for (unsigned int i = 0; i < a->dim; i++) {
    a->c[i] = b->c[i] + c->c[i];
  }
}

void double_vector_sub(double_vector_ptr a, double_vector_srcptr b,
    double_vector_srcptr c)
{
  ASSERT(a->dim == b->dim);
  ASSERT(c->dim == b->dim);

  for (unsigned int i = 0; i < a->dim; i++) {
    a->c[i] = b->c[i] - c->c[i];
  }
}


void double_vector_fprintf(FILE * file, double_vector_srcptr v)
{
  fprintf(file, "[");
  for (unsigned int i = 0; i < v->dim - 1; i++)
    fprintf(file, "%f, ", v->c[i]);
  if (v->dim != 0) {
    fprintf(file, "%f]\n", v->c[v->dim - 1]);
  } else {
    fprintf(file, "]\n");
  }
}

double double_vector_norml2sqr(double_vector_srcptr a)
{
  ASSERT(a->dim >= 2);

  return double_vector_dot_product(a, a);
}

double double_vector_norml2(double_vector_srcptr a)
{
  ASSERT(a->dim >= 2);

  return sqrt((double)double_vector_norml2sqr(a));
}

double double_vector_dot_product(double_vector_srcptr v0,
    double_vector_srcptr v1)
{
  ASSERT(v0->dim == v1->dim);

  double dot = 0;
  for (unsigned int i = 0; i < v0->dim; i++) {
    dot += v0->c[i] * v1->c[i];
  }
  return dot;
}

int double_vector_in_sieving_region(double_vector_srcptr v,
    sieving_bound_srcptr H)
{
  ASSERT(v->dim == H->t);

  for (unsigned int i = 0; i < v->dim - 1; i++) {
    if (-(double)H->h[i] > v->c[i] || (double)H->h[i] <= v->c[i]) {
      return 0;
    }
  }
  if (0 > v->c[v->dim - 1] || (double)H->h[v->dim - 1] <= v->c[v->dim - 1]) {
    return 0;
  }
  return 1;
}

void int64_vector_to_double_vector(double_vector_ptr v_d,
    int64_vector_srcptr v_i)
{
  ASSERT(v_d->dim == v_i->dim);

  for (unsigned int i = 0; i < v_d->dim; i++) {
    v_d->c[i] = (double)v_i->c[i];
  } 
}

double double_vector_orthogonal_projection(double_vector_ptr res,
    double_vector_srcptr u, double_vector_srcptr v)
{
  ASSERT(res->dim == u->dim);
  ASSERT(u->dim == v->dim);

  double coeff =
    double_vector_dot_product(u, v) / double_vector_norml2sqr(u);
  for (unsigned int i = 0; i < res->dim; i++) {
    res->c[i] = coeff * u->c[i];
  }
  return coeff;
}
