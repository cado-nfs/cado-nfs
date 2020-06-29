#ifndef TEST_ITER_H_
#define TEST_ITER_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Iterators that produce the smallest and largest possible values uint64_t,
   and smallest, largest, and closest to 0 for int64_t. */

typedef struct test_iter_s * test_iter_ptr;
typedef const struct test_iter_s * test_iter_srcptr;

typedef void (*test_iter_next_func_ptr)(void *, test_iter_ptr);

struct test_iter_s {
  unsigned long next, max;
  test_iter_next_func_ptr next_func;
};

typedef struct test_iter_s test_iter_t[1];

void test_iter_init (test_iter_t, unsigned long, test_iter_next_func_ptr);
void test_iter_next (void *, test_iter_t);
int test_iter_isdone (test_iter_t);
void test_iter_int64_next(void *, test_iter_t);
void test_iter_uint64_next(void *, test_iter_t);

#ifdef __cplusplus
}
#endif

#endif	/* TEST_ITER_H_ */
