#include "cado.h" // IWYU pragma: keep
#include <cstdint>
#include <cinttypes>
#include <stdio.h>
#include <stdlib.h>
#include "macros.h"
#include "tests_common.h"
#include "utils/u64arith.h"
#include "misc.h"

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
  const uint64_t msb = UINT64_C(1) << 63;
  unsigned long i, iter = 100;
  tests_common_get_iter(&iter);
  for (i = 0; i < iter; i++) {
    test_one_u64arith_reciprocal_for_div(msb + i);
    test_one_u64arith_reciprocal_for_div(UINT64_MAX - 1);
    test_one_u64arith_reciprocal_for_div(msb + UINT64_MAX / 2 / iter);
    test_one_u64arith_reciprocal_for_div(msb | u64_random(state));
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
  test_one_u64arith_divqr_2_1_1(123, 1, 0);
  test_one_u64arith_divqr_2_1_1(123, 1, 1);
  test_one_u64arith_divqr_2_1_1(123, UINT64_MAX, 122);
  test_one_u64arith_divqr_2_1_1(UINT64_MAX, UINT64_MAX, UINT64_MAX - 1);
}

void
test_one_u64arith_reciprocal_for_div_3by2(const uint64_t d0, const uint64_t d1) {
    uint64_t t, p0, p1, p2;
    const uint64_t v = u64arith_reciprocal_for_div_3by2(d0, d1);
    /* Let D = d1*2^64+d0
     * v + beta = floor((beta^3 - 1) / D)
     * (beta^3 - 1) / D - 1 < v + beta <= (beta^3 - 1) / D
     * beta^3 - 1 - D < v*D + beta*D <= beta^3 - 1
     */
    u64arith_mul_1_1_2(&p0, &p1, v, d0);
    u64arith_mul_1_1_2(&t, &p2, v, d1);
    bool cy = u64arith_add_2_2_cy(&p1, &p2, t, 0);
    if (cy) {
        printf("%s(%" PRIu64 ":%"  PRIu64 "): Error, overflow occurred, v = %" PRIu64 "\n",
               __func__, d1, d0, v);
        exit(EXIT_FAILURE);
    }
    /* Now p2:p1:p0 = v*D */
    cy = u64arith_add_2_2_cy(&p1, &p2, d0, d1);
    /* Now cy:p2:p1:p0 = v*D + beta*D <= beta^3 - 1 */
    if (cy) {
        printf("%s(%" PRIu64 ":%"  PRIu64"): Error, overflow occurred, v = %" PRIu64 "\n",
               __func__, d1, d0, v);
        exit(EXIT_FAILURE);
    }
    /* Now p2:p1:p0 = v*D + beta*D <= beta^3 - 1
     * We still need to check beta^3 - 1 - D < p2:p1:p0  <=>
     * beta^3 - 1 - p2:p1:p0 < D
     * Compute beta^3 - 1 - (p2:p1:p0) */
    p0 = ~p0;
    p1 = ~p1;
    p2 = ~p2;

    if (p2 != 0 || u64arith_ge_2_2(p0, p1, d0, d1)) {
        printf("%s(%" PRIu64 ":%"  PRIu64 "): Error, v = %" PRIu64 " incorrect\n",
               __func__, d0, d1, v);
        exit(EXIT_FAILURE);
    }
    /* printf("u64arith_reciprocal_for_div_3by2(" PRIu64 ":%" PRIu64 ") = %" PRIu64 "\n", __func__, d0, d1, v); */
}

void test_u64arith_reciprocal_for_div_3by2(const unsigned long iter) {
    const uint64_t msb = UINT64_C(1) << 63;
    test_one_u64arith_reciprocal_for_div_3by2(0, msb);
    test_one_u64arith_reciprocal_for_div_3by2(0, UINT64_MAX);
    test_one_u64arith_reciprocal_for_div_3by2(UINT64_MAX, UINT64_MAX);
    for (unsigned long i = 0; i < iter; i++) {
        test_one_u64arith_reciprocal_for_div_3by2(u64_random(state), msb);
        test_one_u64arith_reciprocal_for_div_3by2(u64_random(state), UINT64_MAX);
        test_one_u64arith_reciprocal_for_div_3by2(0, u64_random(state) | msb);
        test_one_u64arith_reciprocal_for_div_3by2(UINT64_MAX, u64_random(state) | msb);
        test_one_u64arith_reciprocal_for_div_3by2(u64_random(state), u64_random(state) | msb);
    }
}

void test_one_u64arith_divqr_3_2_1_recip_precomp(
    const uint64_t q, const uint64_t d0, const uint64_t d1, const uint64_t r0, const uint64_t r1) 
{
    /* D = d0 + 2^64*d1, R = r0 + 2^64*r1 */
    ASSERT_ALWAYS(d1 != 0);
    ASSERT_ALWAYS(u64arith_lt_2_2(r0, r1, d0, d1));
    uint64_t u0, u1, u2, t;
    char cy;
    u64arith_mul_1_1_2(&u0, &u1, d0, q);
    u64arith_mul_1_1_2(&t, &u2, d1, q);
    u64arith_add_1_2(&u1, &u2, t); /* U = u0 + 2^64*u1 + 2^128*u2 = D*q */
    cy = u64arith_add_2_2_cy(&u0, &u1, r0, r1);
    ASSERT_ALWAYS(u2 <= UINT64_MAX - cy);
    u2 += cy; /* U = D*q + R */
    /* D <= 2^128-1, q <= 2^64-1, R <= D-1 <= 2^128-2
       D*q + R <= 2^192-2^64-1 */
    
    uint64_t q_new, r0_new, r1_new;
    u64arith_divqr_3_2_1 (&q_new, &r0_new, &r1_new, u0, u1, u2, d0, d1);

    if (q_new != q || r0_new != r0 || r1_new != r1) {
        printf("test_one_u64arith_divqr_3_2_1_recip_precomp("
               "q=%" PRIu64 ", d0=%" PRIu64 ", d1=%" PRIu64 ", r0=%" PRIu64 ", r1=%" PRIu64 "): "
               "u0 = %" PRIu64 ", u1 = %" PRIu64 ", u2 = %" PRIu64
               " Error, got result q=%" PRIu64 ", r0=%" PRIu64 ", r1=%" PRIu64 "\n",
               q, d0, d1, r0, r1, u0, u1, u2, q_new, r0_new, r1_new);
        exit(EXIT_FAILURE);
    }
}

void test_u64arith_divqr_3_2_1_recip_precomp(const unsigned long iter) {
    test_one_u64arith_divqr_3_2_1_recip_precomp(1, 1, 1, 0, 1);
    for (unsigned long i = 0; i < iter; i++) {
        const uint64_t q = u64_random(state),
                       d0 = u64_random(state);
        uint64_t d1 = u64_random(state),
                 r0 = u64_random(state),
                 r1 = u64_random(state);
        if (d1 == 0)
            d1 = 1;
        r1 %= d1;
        test_one_u64arith_divqr_3_2_1_recip_precomp(q, d0, d1, r0, r1);
    }
}


static void test_one_u64arith_div2mod(const uint64_t a, const uint64_t m) {
    uint64_t r, t;
    r = u64arith_div2mod(a, m);
    u64arith_addmod_1_1(&t, r, r, m);
    if (t != a) {
        printf("u64arith_div2mod(%" PRIu64 ", %" PRIu64 ") = %" PRIu64 " wrong\n", a, m, r);
        exit(EXIT_FAILURE);
    }
}

void test_u64arith_div2mod(unsigned long iter) {
    for (unsigned long i = 0; i < iter; i++) {
        const uint64_t a = u64_random(state);
        const uint64_t m = u64_random(state) | 1;
        test_one_u64arith_div2mod(a % m, m);
    }
}

/* TODO: add tests for u64arith_divr_2_1_1() */
/* TODO: add tests for u64arith_shrd() */
/* TODO: add tests for u64arith_shld() */
/* TODO: add tests for u64arith_ctz() */
/* TODO: add tests for u64arith_clz() */
/* TODO: add tests for u64arith_invmod() */
/* TODO: add tests for u64arith_sqrt() */
/* TODO: add tests for u64arith_post_process_inverse() */
/* TODO: add tests for u64arith_redc() */

// coverity[root_function]
int
main (int argc, const char *argv[])
{
  unsigned long iter = 100;
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);
  test_u64arith_gt_2_2();
  test_u64arith_add_1_2();
  test_u64arith_add_2_2();
  test_u64arith_add_2_2_cy();
  test_u64arith_addmod_1_1();
  test_u64arith_mul_1_1_2();
  test_u64arith_sqr_1_2();
  test_u64arith_divqr_2_1_1();
  test_u64arith_reciprocal_for_div();
  test_u64arith_reciprocal_for_div_3by2(iter);
  test_u64arith_divqr_3_2_1_recip_precomp(iter);
  test_u64arith_div2mod(iter);
  tests_common_clear();
  exit (EXIT_SUCCESS);
}
