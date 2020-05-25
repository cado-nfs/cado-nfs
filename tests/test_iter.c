#include "cado.h" // IWYU pragma: keep
#include <stdint.h>
#ifdef TESTDRIVE
#include <stdio.h>
#endif
#include "test_iter.h"
#include "macros.h"

void
test_iter_init (test_iter_t iter, const unsigned long max,
                test_iter_next_func_ptr next_func)
{
  iter->next = 0;
  iter->max = max;
  iter->next_func = next_func;
}

int
test_iter_isdone (test_iter_t iter)
{
  ASSERT(iter->next <= iter->max);
  return iter->next == iter->max;
}

void
test_iter_next (void *output, test_iter_t iter)
{
  ASSERT(!test_iter_isdone(iter));
  iter->next_func(output, iter);
  iter->next++;
}

void
test_iter_int64_next(void * const output, test_iter_t iter)
{
  int64_t * const v_out = (int64_t *)output;
  unsigned long limit1 = (iter->max + 3) / 4;
  unsigned long limit2 = 2 * (iter->max - limit1) / 3;

  unsigned long i = iter->next;
  if (i < limit1) {
    *v_out = INT64_MIN + (int64_t) iter->next;
    return;
  }
  i -= limit1;
  if (i < limit2) {
    *v_out = -(int64_t) ((iter->max - limit1) / 3) + (int64_t) i;
    return;
  }
  i -= limit2;
  *v_out = INT64_MAX - (int64_t) i;
}

void
test_iter_uint64_next(void *output, test_iter_t iter)
{
  uint64_t * const v_out = (uint64_t *)output;
  unsigned long limit1 = (iter->max + 1) / 2;

  unsigned long i = iter->next;
  if (i < limit1) {
    *v_out = (uint64_t) i;
    return;
  }
  i -= limit1;
  *v_out = UINT64_MAX - (uint64_t) i;
}


#ifdef TESTDRIVE
int main(int argc, char **argv)
{
  test_iter_t iter;
  test_iter_init(iter, atoi(argv[1]), test_iter_int64_next);
  while (!test_iter_isdone(iter)) {
    int64_t v;
    test_iter_next(&v, iter);
    printf ("int64_t %lld\n", (long long) v);
  }

  test_iter_init(iter, atoi(argv[1]), test_iter_uint64_next);
  while (!test_iter_isdone(iter)) {
    uint64_t v;
    test_iter_next(&v, iter);
    printf ("uint64_t %llu\n", (unsigned long long) v);
  }

  return 0;
}
#endif
