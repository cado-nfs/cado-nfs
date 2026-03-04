#include "cado.h" // IWYU pragma: keep

#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <memory>

#include "fmt/base.h"
#include <gmp.h>

#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "misc.h"
#include "portability.h"
#include "relation-tools.h"
#include "relation.hpp"
#include "tests_common.h"
#include "typedefs.h"

static unsigned long mpz_compute_r(cxx_mpz const & a, cxx_mpz & b,
                                   cxx_mpz const & p)
{
    if (mpz_invert(b, b, p) == 0)
        return mpz_get_ui(p);
    else {
        mpz_mul(b, a, b);
        mpz_mod(b, b, p);

        return mpz_get_ui(b);
    }
}

static int test_compute_r(unsigned int nb)
{
    int err = 0;
    cxx_mpz tp, ta, tb;

    for (unsigned int i = 0; i < nb; i++) {
        int64_t a;
        uint64_t b;
        unsigned long p;

        /* 5% of tests are for the case where p = 2^k */
        if (i % 20 == 0) {
            unsigned int exp = gmp_urandomb_ui(state, 5);
            exp = (exp == 0) ? 1 : exp;
            p = 1UL << exp;
            mpz_set_ui(tp, p);
        } else {
            mpz_set_ui(tp, gmp_urandomb_ui(state, 31));
            mpz_nextprime(tp, tp);
            p = mpz_get_ui(tp);
        }

        a = i64_random(state);
        /* 5% of tests are for the case where b = 0 mod p (with b > 0)
         * We do not need to test for free relations as they never go throught
         * relation_compute_r
         */
        if (i < (nb / 20))
            b = gmp_urandomb_ui(state, 31) * p;
        else {
            b = u64_random(state);
            b += !b;
        }
        mpz_set_int64(ta, a);
        mpz_set_uint64(tb, b);

        unsigned long const r = relation_compute_r(a, b, p);

        unsigned long const r2 = mpz_compute_r(ta, tb, tp);
        if (r != r2) {
            gmp_fprintf(stderr,
                        "ERROR: a=%" PRId64 " b=%" PRIu64 " p=%" PRpr "\n"
                        "Got r=%" PRpr " instead of %" PRpr "\n",
                        a, b, p, r, r2);
            err++;
        }
    }
    return err;
}

static int test_compute_r_large_ab(unsigned int nb)
{
    int err = 0;
    cxx_mpz a, b, tp;

    for (unsigned int i = 0; i < nb; i++) {
        unsigned long p;

        /* 5% of tests are for the case where p = 2^k */
        if (i % 20 == 0) {
            unsigned int exp = gmp_urandomb_ui(state, 5);
            exp = (exp == 0) ? 1 : exp;
            p = 1UL << exp;
            mpz_set_ui(tp, p);
        } else {
            mpz_set_ui(tp, gmp_urandomb_ui(state, 31));
            mpz_nextprime(tp, tp);
            p = mpz_get_ui(tp);
        }

        mpz_urandomb(a, state, 256);
        if (gmp_urandomb_ui(state, 1)) {
            mpz_neg(a, a);
        }
        /* 5% of tests are for the case where b = 0 mod p (with b > 0)
         * We do not need to test for free relations as they never go throught
         * relation_compute_r
         */
        mpz_urandomb(b, state, 192);
        mpz_add_ui(b, b, 1u);
        if (i < (nb / 20)) {
            mpz_mul_ui(b, b, p);
        }

        unsigned long const r = relation_compute_r(a, b, p);

        unsigned long const r2 = mpz_compute_r(a, b, tp);
        if (r != r2) {
            fmt::print(stderr,
                       "ERROR: a={} b={} p={}\nGot r={} instead of {}\n",
                       a, b, p, r, r2);
            err++;
        }
    }
    return err;
}

static int test_compute_all_r(unsigned int nb)
{
    int err = 0;
    cxx_mpz tp, ta, tb;

    for (unsigned int i = 0; i < nb; i++) {
        int64_t const a = i64_random(state);
        uint64_t b = u64_random(state);
        b += !b;
        relation t1(a, b);

        for (unsigned long k = 0; k <= gmp_urandomm_ui(state, 5); k++) {
            mpz_set_ui(tp, gmp_urandomb_ui(state, 31));
            mpz_nextprime(tp, tp);
            t1.add(0, tp, cxx_mpz(0));
        }
        for (unsigned long k = 0; k <= gmp_urandomm_ui(state, 5); k++) {
            mpz_set_ui(tp, gmp_urandomb_ui(state, 31));
            mpz_nextprime(tp, tp);
            t1.add(1, tp, cxx_mpz(0));
        }

        relation t2 = t1;
        t1.fixup_r();

        for (size_t k = 0; k < t2.sides[1].size(); k++) {
            mpz_set_int64(ta, t2.a);
            mpz_set_uint64(tb, t2.b);
            mpz_set(tp, t2.sides[1][k].p);
            unsigned long const r = mpz_compute_r(ta, tb, tp);
            if (r != mpz_get_ui(t1.sides[1][k].r)) {
                gmp_fprintf(stderr,
                            "ERROR: a=%" PRId64 " b=%" PRIu64 " p=%" PRpr "\n"
                            "Got r=%" PRpr " instead of %" PRpr "\n",
                            t2.a, t2.b, (mpz_srcptr)t2.sides[1][k].p,
                            (mpz_srcptr)t1.sides[1][k].r, r);
                err++;
            }
        }
    }

    return err;
}

static int check_str_err(char const * s1, char const * s2, cxx_mpz const & t)
{
    if (strcmp(s1, s2) != 0) {
        fmt::print(stderr,
                   "ERROR with integer {} got \"{}\" instead of "
                   "\"{}\"\n",
                   t, s2, s1);
        return 1;
    } else
        return 0;
}

static int test_conversion(unsigned int nb)
{
    int err = 0;

    const std::unique_ptr<char[]> s1(new char[25]);
    const std::unique_ptr<char[]> s2(new char[25]);

    const std::unique_ptr<char[]> l1(new char[100]);
    const std::unique_ptr<char[]> l2(new char[100]);

    char * tmp;
    cxx_mpz t;

    for (unsigned int i = 0; i < nb; i++) {
        int64_t const a = i64_random(state);
        uint64_t const b = u64_random(state);

        mpz_set_int64(t, a);

        mpz_get_str(s1.get(), 10, t);
        tmp = d64toa10(s2.get(), a);
        *tmp = '\0';
        err += check_str_err(s1.get(), s2.get(), t);

        mpz_get_str(s1.get(), 16, t);
        tmp = d64toa16(s2.get(), a);
        *tmp = '\0';
        err += check_str_err(s1.get(), s2.get(), t);

        mpz_set_uint64(t, b);

        mpz_get_str(s1.get(), 10, t);
        tmp = u64toa10(s2.get(), b);
        *tmp = '\0';
        err += check_str_err(s1.get(), s2.get(), t);

        mpz_get_str(s1.get(), 16, t);
        tmp = u64toa16(s2.get(), b);
        *tmp = '\0';
        err += check_str_err(s1.get(), s2.get(), t);

        /* for large a and b */
        {
            cxx_mpz a, b;
            mpz_urandomb(a, state, 256);
            if (gmp_urandomb_ui(state, 1)) {
                mpz_neg(a, a);
            }
            mpz_urandomb(b, state, 192);

            mpz_get_str(l1.get(), 16, a);
            tmp = d64toa16(l2.get(), a);
            *tmp = '\0';
            err += check_str_err(l1.get(), l2.get(), t);

            mpz_get_str(l1.get(), 16, b);
            tmp = u64toa16(l2.get(), b);
            *tmp = '\0';
            err += check_str_err(l1.get(), l2.get(), t);
        }
    }

    return err;
}

int main(int argc, char const * argv[])
{
    int err = 0;
    unsigned long iter = 10000;

    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
    tests_common_get_iter(&iter);

    err += test_compute_r(iter);
    err += test_compute_r_large_ab(iter);
    err += test_compute_all_r(iter / 10);
    err += test_conversion(iter);

    if (err)
        fprintf(stderr, "# %d erro%s found\n", err, (err == 1) ? "r" : "rs");
    tests_common_clear();
    return (err) ? EXIT_FAILURE : EXIT_SUCCESS;
}
