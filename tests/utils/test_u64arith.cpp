#include "cado.h"
#include <cstdint>
#include <cinttypes>
#include "macros.h"
#include "tests_common.h"
#include "test_iter.h"
#include "u64arith.h"

void
test_one_u64arith_gt_2_2(const uint64_t a1, const uint64_t a2,
  const uint64_t b1, const uint64_t b2, const int v) {
  const int r = u64arith_gt_2_2(a1, a2, b1, b2);
  if (v != r) {
    printf("%s: Error, got result %d, but expected %d\n", __func__, v, r);
    abort();
  }
}

void test_u64arith_gt_2_2() {
  test_one_u64arith_gt_2_2(0,0,0,0,0);
  test_one_u64arith_gt_2_2(1,0,0,0,1);
  test_one_u64arith_gt_2_2(0,1,0,0,1);
  test_one_u64arith_gt_2_2(1,1,0,0,1);
  test_one_u64arith_gt_2_2(0,0,1,0,0);
  test_one_u64arith_gt_2_2(1,0,1,0,0);
  test_one_u64arith_gt_2_2(0,1,1,0,1);
  test_one_u64arith_gt_2_2(1,1,1,0,1);
  test_one_u64arith_gt_2_2(0,0,0,1,0);
  test_one_u64arith_gt_2_2(1,0,0,1,0);
  test_one_u64arith_gt_2_2(0,1,0,1,0);
  test_one_u64arith_gt_2_2(1,1,0,1,1);
  test_one_u64arith_gt_2_2(0,0,1,1,0);
  test_one_u64arith_gt_2_2(1,0,1,1,0);
  test_one_u64arith_gt_2_2(0,1,1,1,0);
  test_one_u64arith_gt_2_2(1,1,1,1,0);
}

void
test_one_u64arith_add_1_2(uint64_t r1, uint64_t r2,
    const uint64_t a, const uint64_t v1, const uint64_t v2) {
  u64arith_add_1_2(&r1, &r2, a);
  if (r1 != v1 || r2 != v2) {
    printf("%s: Error, got result %" PRIu64 ":%" PRIu64 ", but expected %" PRIu64 ":%" PRIu64 "\n",
        __func__, r2, r1, v2, v1);
    exit(EXIT_FAILURE);
  }
}

void
test_u64arith_add_1_2() {
  test_one_u64arith_add_1_2(1, 4, 1, 2, 4);
  test_one_u64arith_add_1_2(UINT64_MAX, 1, 1, 0, 2);
  test_one_u64arith_add_1_2(UINT64_MAX, UINT64_MAX, 1, 0, 0);
}


void
test_one_u64arith_add_2_2(uint64_t r1, uint64_t r2,
    const uint64_t a1, const uint64_t a2, const uint64_t v1, const uint64_t v2) {
  u64arith_add_2_2(&r1, &r2, a1, a2);
  if (r1 != v1 || r2 != v2) {
    printf("%s: Error, got result %" PRIu64 ":%" PRIu64 ", but expected %" PRIu64 ":%" PRIu64 "\n",
        __func__, r2, r1, v2, v1);
    exit(EXIT_FAILURE);
  }
}


void
test_u64arith_add_2_2() {
  test_one_u64arith_add_2_2(1, 4, 1, 0, 2, 4);
  test_one_u64arith_add_2_2(UINT64_MAX, 1, 1, 0, 0, 2);
  test_one_u64arith_add_2_2(UINT64_MAX, UINT64_MAX, 1, 0, 0, 0);
  test_one_u64arith_add_2_2(1, 4, 0, 1, 1, 5);
}


void
test_one_u64arith_add_2_2_cy(uint64_t r1, uint64_t r2,
    const uint64_t a1, const uint64_t a2, const uint64_t v1, const uint64_t v2,
    const char vcy) {
  unsigned char cy = u64arith_add_2_2_cy(&r1, &r2, a1, a2);
  if (r1 != v1 || r2 != v2 || vcy != cy) {
    printf("%s: Error, got result %hhu:%" PRIu64 ":%" PRIu64 ", but expected %hhu:%" PRIu64 ":%" PRIu64 "\n",
        __func__, cy, r2, r1, vcy, v2, v1);
    exit(EXIT_FAILURE);
  }
}


void
test_u64arith_add_2_2_cy() {
  test_one_u64arith_add_2_2_cy(1, 4, 1, 0, 2, 4, 0);
  test_one_u64arith_add_2_2_cy(UINT64_MAX, 1, 1, 0, 0, 2, 0);
  test_one_u64arith_add_2_2_cy(UINT64_MAX, UINT64_MAX, 1, 0, 0, 0, 1);
  test_one_u64arith_add_2_2_cy(1, 4, 0, 1, 1, 5, 0);
}

static void check_result(const uint64_t r, const uint64_t v, const char *func) {
  if (r != v) {
    printf("%s: Error, got result %" PRIu64 ", but expected %" PRIu64 "\n",
        func, r, v);
    exit(EXIT_FAILURE);
  }
}

static void
test_one_u64arith_addmod_1_1(const uint64_t a, const uint64_t b,
    const uint64_t m, const uint64_t v) {
  uint64_t r, ca = a, cb = b;
  u64arith_addmod_1_1 (&r, a, b, m);
  check_result(r, v, __func__);
  /* Test with aliased output and input operands */
  u64arith_addmod_1_1 (&ca, ca, b, m);
  check_result(ca, v, __func__);
  u64arith_addmod_1_1 (&cb, a, cb, m);
  check_result(cb, v, __func__);
}


static void
test_u64arith_addmod_1_1() {
  test_one_u64arith_addmod_1_1(1, 1, 5, 2);
  test_one_u64arith_addmod_1_1(2, 3, 5, 0);
  test_one_u64arith_addmod_1_1(2, 4, 5, 1);
  test_one_u64arith_addmod_1_1(1, UINT64_MAX - 1, UINT64_MAX, 0);
  test_one_u64arith_addmod_1_1(2, UINT64_MAX - 1, UINT64_MAX, 1);
  /* TODO: Write better tests */
  /* u64arith_addmod_1_1() is specified to be correct for second argument
     equal to modulus. */
  test_one_u64arith_addmod_1_1(2, 5, 5, 2);
}

/* TODO: add tests for u64arith_sub_1_2() */
/* TODO: add tests for u64arith_sub_2_2() */
/* TODO: add tests for u64arith_sub_2_2_cy() */
/* TODO: add tests for u64arith_sub_1_1_ge() */
/* TODO: add tests for u64arith_sub_2_2_ge() */
/* TODO: add tests for u64arith_submod_1_1() */

void
test_one_u64arith_mul_1_1_2(const uint64_t a, const uint64_t b,
    const uint64_t v1, const uint64_t v2) {
  uint64_t r1, r2;
  u64arith_mul_1_1_2 (&r1, &r2, a, b);
  if (r1 != v1 || r2 != v2) {
    printf("%s: Error, got result %" PRIu64 ":%" PRIu64 ", but expected %" PRIu64 ":%" PRIu64 "\n",
        __func__, r2, r1, v2, v1);
    exit(EXIT_FAILURE);
  }
}

void
test_u64arith_mul_1_1_2() {
  test_one_u64arith_mul_1_1_2(1, 2, 2, 0);
  uint64_t a = 1, b = 1;
  a <<= 31; b <<= 33;
  test_one_u64arith_mul_1_1_2(a, b, 0, 1);
  test_one_u64arith_mul_1_1_2(UINT64_MAX, UINT64_MAX - 1, 2, UINT64_MAX - 2);
  test_one_u64arith_mul_1_1_2(UINT64_MAX, UINT64_MAX, 1, UINT64_MAX - 1);
}


void
test_one_u64arith_sqr_1_2(const uint64_t a, const uint64_t v1, const uint64_t v2) {
  uint64_t r1, r2;
  u64arith_sqr_1_2 (&r1, &r2, a);
  if (r1 != v1 || r2 != v2) {
    printf("%s: Error, got result %" PRIu64 ":%" PRIu64 ", but expected %" PRIu64 ":%" PRIu64 "\n",
        __func__, r2, r1, v2, v1);
    exit(EXIT_FAILURE);
  }
}


void
test_u64arith_sqr_1_2() {
  test_one_u64arith_sqr_1_2(2, 4, 0);
  uint64_t a = 1;
  a <<= 32;
  test_one_u64arith_sqr_1_2(a, 0, 1);
  test_one_u64arith_sqr_1_2(UINT64_MAX, 1, UINT64_MAX - 1);
}


void
test_one_u64arith_reciprocal_for_div(const uint64_t d) {
  uint64_t q, r;
  const uint64_t v = u64arith_reciprocal_for_div(d);
  u64arith_divqr_2_1_1(&q, &r, UINT64_MAX, UINT64_MAX - d, d);
  if (q != v) {
    printf("%s(%" PRIu64 "): Error, got result %" PRIu64 ", but expected %" PRIu64 "\n",
        __func__, d, v, q);
    exit(EXIT_FAILURE);
  }
  /* printf("u64arith_reciprocal_for_div(%" PRIu64 ") = %" PRIu64 "\n", d, v); */
}

void
test_u64arith_reciprocal_for_div() {
  const uint64_t one = 1;
  unsigned long i, iter = 100;
  tests_common_get_iter(&iter);
  for (i = 0; i < iter; i++) {
    test_one_u64arith_reciprocal_for_div((one << 63) + i);
    test_one_u64arith_reciprocal_for_div(UINT64_MAX - 1);
    test_one_u64arith_reciprocal_for_div((one << 63) + UINT64_MAX / 2 / iter);
    test_one_u64arith_reciprocal_for_div((one << 63) | random_uint64());
  }
}


void
test_one_u64arith_divqr_2_1_1(const uint64_t b, const uint64_t cq, const uint64_t cr)
{
  uint64_t a1, a2, q, r;
  u64arith_mul_1_1_2(&a1, &a2, cq, b);
  u64arith_add_2_2(&a1, &a2, cr, 0);
  u64arith_divqr_2_1_1(&q, &r, a1, a2, b);
  if (q != cq || r != cr) {
    printf("%s(%" PRIu64 ", %" PRIu64 ", %" PRIu64 "): Error, got result q=%"
        PRIu64 ", r=%" PRIu64 ", but expected q=%" PRIu64 ", r=%" PRIu64 "\n",
        __func__, a1, a2, b, q, r, cq, cr);
    exit(EXIT_FAILURE);
  }
}


void
test_u64arith_divqr_2_1_1()
{
  unsigned long i, iter = 100;

  tests_common_get_iter(&iter);
  test_one_u64arith_divqr_2_1_1(123, 1, 0);
  test_one_u64arith_divqr_2_1_1(123, 1, 1);
  test_one_u64arith_divqr_2_1_1(123, UINT64_MAX, 122);
  test_one_u64arith_divqr_2_1_1(UINT64_MAX, UINT64_MAX, UINT64_MAX - 1);

  for (i = 0; i < iter; i++) {
    uint64_t b = random_uint64(),
      q = random_uint64(),
      r = random_uint64() % b;
      test_one_u64arith_divqr_2_1_1(b, q, r);
  }
}

/* TODO: add tests for u64arith_divr_2_1_1() */
/* TODO: add tests for u64arith_shrd() */
/* TODO: add tests for u64arith_shld() */
/* TODO: add tests for u64arith_ctz() */
/* TODO: add tests for u64arith_clz() */
/* TODO: add tests for u64arith_invmod() */
/* TODO: add tests for u64arith_div2mod() */
/* TODO: add tests for u64arith_sqrt() */
/* TODO: add tests for u64arith_post_process_inverse() */
/* TODO: add tests for u64arith_redc() */

int
main (int argc, const char *argv[])
{
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  test_u64arith_gt_2_2();
  test_u64arith_add_1_2();
  test_u64arith_add_2_2();
  test_u64arith_add_2_2_cy();
  test_u64arith_addmod_1_1();
  test_u64arith_mul_1_1_2();
  test_u64arith_sqr_1_2();
  test_u64arith_divqr_2_1_1();
  test_u64arith_reciprocal_for_div();
  tests_common_clear();
  exit (EXIT_SUCCESS);
}
