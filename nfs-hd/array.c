#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include "sieving_bound.h"
#include "array.h"
#include "int64_vector.h"

void array_init(array_ptr array, uint64_t number_element)
{
  ASSERT(number_element != 0);

  array->number_element = number_element;
  array->array = malloc(sizeof(unsigned char) * number_element);

  ASSERT_ALWAYS(array->array != NULL);
}

void array_clear(array_ptr array)
{
  free(array->array);
  array->number_element = 0;
}

void array_fprintf(FILE * filew, array_srcptr array)
{
  ASSERT(array->number_element != 0);

  fprintf(filew, "[");
  for (uint64_t i = 0; i < array->number_element - 1; i++) {
    fprintf(filew, "%u, ", array->array[i]);
  }
  fprintf(filew, "%u]\n", array->array[array->number_element -1]);
}

void array_index_mpz_vector(mpz_vector_ptr v, uint64_t index,
    sieving_bound_srcptr H, uint64_t number_element)
{
  ASSERT(v->dim == H->t);

  unsigned int k = v->dim - 1;
  uint64_t prod = number_element / ((uint64_t)H->h[k]);
  uint64_t res = index / prod;
  mpz_vector_setcoordinate_si(v, k, res);
  index = index - res * prod;
  k = k - 1;
  prod = prod / (2 * H->h[k]);
  for ( ; k > 0; k--) {
    res = index / prod;
    mpz_vector_setcoordinate_si(v, k, res - H->h[k]);
    index = index - res * prod;
    prod = prod / (2 * H->h[k - 1]);
  }
  res = index / prod;
  mpz_vector_setcoordinate_si(v, k, res - H->h[k]);

#ifndef NDEBUG
  mpz_t tmp;
  mpz_init(tmp);
  for (unsigned int i = 0; i < v->dim - 1; i++) {
    mpz_set(tmp, v->c[i]);
    ASSERT(mpz_cmp_si(tmp, -(int64_t)H->h[i]) >= 0);
    ASSERT(mpz_cmp_ui(tmp, H->h[i]) < 0);
  }
  mpz_set(tmp, v->c[v->dim - 1]);
  ASSERT(mpz_cmp_ui(tmp, 0) >= 0);
  ASSERT(mpz_cmp_ui(tmp, H->h[H->t - 1]) < 0);
  mpz_clear(tmp);
#endif
}

uint64_t array_mpz_vector_index(mpz_vector_srcptr v,
    sieving_bound_srcptr H, MAYBE_UNUSED uint64_t number_element)
{
#ifndef NDEBUG
  ASSERT(v->dim == H->t);
  mpz_t tmp;
  mpz_init(tmp);
  for (unsigned int i = 0; i < v->dim - 1; i++) {
    mpz_set(tmp, v->c[i]);
    ASSERT(mpz_cmp_si(tmp, -(int64_t)H->h[i]) >= 0);
    ASSERT(mpz_cmp_ui(tmp, H->h[i]) < 0);
  }
  mpz_set(tmp, v->c[v->dim - 1]);
  ASSERT(mpz_cmp_ui(tmp, 0) >= 0);
  ASSERT(mpz_cmp_ui(tmp, H->h[H->t - 1]) < 0);
  mpz_clear(tmp);
#else
  mpz_t tmp;
#endif

  mpz_init(tmp);
  uint64_t prod = 1;
  unsigned int k = 0;
  mpz_set(tmp, v->c[k]);
  uint64_t index = (uint64_t)(mpz_get_si(tmp) + H->h[k]);

  for (k = 1; k < v->dim - 1; k++) {
    prod = prod * (2 * H->h[k - 1]);
    mpz_set(tmp, v->c[k]);
    index = index + (uint64_t)(mpz_get_si(tmp) + H->h[k]) * prod;
  }
  prod = prod * (2 * H->h[k - 1]);
  mpz_set(tmp, v->c[k]);
  index = index + (uint64_t)(mpz_get_si(tmp)) * prod;
  mpz_clear(tmp);

  ASSERT(index < number_element);

  return index;
}

void array_index_int64_vector(int64_vector_ptr v, uint64_t index,
    sieving_bound_srcptr H, uint64_t number_element)
{
  ASSERT(v->dim == H->t);

  unsigned int k = v->dim - 1;
  uint64_t prod = number_element / ((uint64_t)H->h[k]);
  uint64_t res = index / prod;
  v->c[k] = (int64_t)res;
  index = index - res * prod;
  k = k - 1;
  prod = prod / (2 * H->h[k]);
  for ( ; k > 0; k--) {
    res = index / prod;
    v->c[k] = (int64_t)res - (int64_t)H->h[k];
    index = index - res * prod;
    prod = prod / (2 * H->h[k - 1]);
  }
  res = index / prod;
  v->c[k] = (int64_t)res - (int64_t)H->h[k];

#ifndef NDEBUG
  int64_t tmp = 0;
  for (unsigned int i = 0; i < v->dim - 1; i++) {
    tmp = v->c[i];
    ASSERT(tmp >= -(int64_t)H->h[i]);
    ASSERT(tmp < (int64_t)H->h[i]);
  }
  tmp = v->c[v->dim - 1];
  ASSERT(tmp >= 0);
  ASSERT(tmp < (int64_t)H->h[H->t - 1]);
#endif
}

uint64_t array_int64_vector_index(int64_vector_srcptr v,
    sieving_bound_srcptr H, MAYBE_UNUSED uint64_t number_element)
{
#ifndef NDEBUG
  ASSERT(v->dim == H->t);
  int64_t tmp = 0;
  for (unsigned int i = 0; i < v->dim - 1; i++) {
    tmp = v->c[i];
    ASSERT(tmp >= -(int64_t)H->h[i]);
    ASSERT(tmp < (int64_t)H->h[i]);
  }
  tmp = v->c[v->dim - 1];
  ASSERT(tmp >= 0);
  ASSERT(tmp < (int64_t)H->h[H->t - 1]);
#endif

  uint64_t prod = 1;
  unsigned int k = 0;
  uint64_t index = (uint64_t)v->c[k] + H->h[k];

  for (k = 1; k < v->dim - 1; k++) {
    prod = prod * (2 * H->h[k - 1]);
    index = index + ((uint64_t)v->c[k] + H->h[k]) * prod;
  }
  prod = prod * (2 * H->h[k - 1]);
  index = index + (uint64_t)(v->c[k]) * prod;

  ASSERT(index < number_element);

  return index;
}

void array_set_at(array_ptr array, int * vec, unsigned char val,
    sieving_bound_srcptr H)
{
  uint64_t prod = 1;
  unsigned int k = 0;
  uint64_t index = (uint64_t) vec[k] + (uint64_t) H->h[k];

  for (k = 1; k < H->t - 1; k++) {
    prod = prod * (uint64_t)(2 * H->h[k - 1]);
    index = index + ((uint64_t)vec[k] + (uint64_t)H->h[k]) * prod;
  }
  prod = prod * (2 * H->h[k - 1]);
  index = index + (uint64_t)(vec[k]) * prod;

  array_set(array, index, val);
}

unsigned char array_get_at(array_ptr array, int * vec,
    sieving_bound_srcptr H)
{
  uint64_t prod = 1;
  unsigned int k = 0;
  uint64_t index = (uint64_t) vec[k] + (uint64_t) H->h[k];

  for (k = 1; k < H->t - 1; k++) {
    prod = prod * (uint64_t)(2 * H->h[k - 1]);
    index = index + ((uint64_t)vec[k] + (uint64_t)H->h[k]) * prod;
  }
  prod = prod * (2 * H->h[k - 1]);
  index = index + (uint64_t)(vec[k]) * prod;

  return array_get(array, index);
}
