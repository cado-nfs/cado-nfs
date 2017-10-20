#include "cado.h"
#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>
#include <limits.h>
#include "sieving_bound.h"
#include "array.h"
#include "utils.h"
#include "macros.h"
#include "int64_vector.h"

int main()
{
  sieving_bound_t H;
  unsigned int t = 4;

  sieving_bound_init(H, t);
  for (unsigned int i = 0; i < t; i++) {
    sieving_bound_set_hi(H, i, 2);
  }

  uint64_t nb = sieving_bound_number_element(H);

  array_t array;
  array_init(array, nb);
  array_set_all_elements(array, 0);
  for (uint64_t i = 0; i < nb; i++) {
    ASSERT_ALWAYS(array_get(array, i) == 0);
  }

  array_set_all_elements(array, UCHAR_MAX);
  for (uint64_t i = 0; i < nb; i++) {
    ASSERT_ALWAYS(array_get(array, i) == UCHAR_MAX);
  }

  int64_vector_t vector;
  int64_vector_init(vector, t);

  mpz_vector_t v;
  mpz_vector_init(v, t);

  for (unsigned int i = 0; i < t - 1; i++) {
    int64_vector_setcoordinate(vector, i, -(int64_t)H->h[i]);
  }
  int64_vector_setcoordinate(vector, t - 1, 0);

  array_set(array, 0, 42);
  ASSERT_ALWAYS(array_get(array, 0) == 42);

  int * vec = (int *)malloc(sizeof(int) * t);
  for (unsigned int i = 0; i < t - 1; i++) {
    vec[i] = -(int)H->h[i];
  }
  vec[t - 1] = 0;

  array_set_at(array, vec, 17, H);
  ASSERT_ALWAYS(array_get_at(array, vec, H) == 17);

  uint64_t index;
  uint64_t index1;

  index = array_int64_vector_index(vector, H, nb);
  ASSERT_ALWAYS(index == 0);
  for (unsigned int i = 1; i < nb; i++) {
    int64_vector_add_one(vector, H);
    index = array_int64_vector_index(vector, H, nb);
    ASSERT_ALWAYS(index == i);
    array_index_mpz_vector(v, index, H, nb);
    for (unsigned int j = 0; j < t; j++) {
      ASSERT_ALWAYS(mpz_cmp_si(v->c[j], vector->c[j]) == 0);
    }
    index1 = array_mpz_vector_index(v, H, nb);
    ASSERT_ALWAYS(index == index1);
  }



  mpz_vector_clear(v);
  int64_vector_clear(vector);
  sieving_bound_clear(H);
  array_clear(array);

  return 0;
}
