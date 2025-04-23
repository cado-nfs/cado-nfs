#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdint>

#include <iostream>         // for std::cout
#include <vector>

#include "tests_common.h"   // for tests_common_cmdline, tests_common_clear, ...
#include "utils/galois_action.hpp"

static bool
test_galois_apply_one(galois_action const & G,
                      bool (*is_fixed_point) (unsigned long r, unsigned long p)){
    bool ret = true;

    std::vector<unsigned long> const P { 2, 3, 5, 11, 23, 101 };

    for (unsigned long const p: P) {
        if ((p == 2 && (G.get_order() == 4 || G.get_order() == 6))
                || (p == 3 && G.get_order() == 6))
            /* skip 2 for action of order 4 and 6, skip 3 for action
             * of order 6: very special case
             */
            continue;
        for (unsigned long r = 0; r <= p; r++) {
            unsigned long sigma_r = G.apply(r, p);
            bool b;
            if (is_fixed_point(r, p)) {
                b = sigma_r == r;
                if (!b) {
                    std::cout << "error in " << __func__ << ": with " << G
                              << ": G.apply(" << r << ", " << p << ") should "
                              << "be a fixed point but got " << sigma_r
                              << "\n";
                }
            } else {
                /* Upper bound the nb of iter by 1000 to avoid endless loop that
                 * could occurs with some errors.
                 */
                unsigned int n = 1;
                while (sigma_r != r && n < 1000) {
                    sigma_r = G.apply(sigma_r, p);
                    n++;
                }
                b = n == G.get_order();
                if (!b) {
                    std::cout << "error in " << __func__ << ": with " << G
                              << ": orbit of (" << r << ", " << p << ") should "
                              << "of length " << G.get_order() << " but has "
                              << "length " << n << "\n";
                }
            }
            ret &= b;
        }
    }

    return ret;
}

static bool
test_galois_apply()
{
    bool ret = true;

    /* galois action identity: every r is a fixed point */
    galois_action const Gnone("none");
    auto is_fixed_point_none = [](unsigned long, unsigned long) {
        return true;
    };
    ret &= test_galois_apply_one(Gnone, is_fixed_point_none);

    /* galois action x->-x: only fixed points are 0, p (=oo) and (r=1,p=2) */
    galois_action const Gneg("_y");
    auto is_fixed_point_neg = [](unsigned long r, unsigned long p) {
            return p == 2 || r == 0 || r == p;
    };
    ret &= test_galois_apply_one(Gneg, is_fixed_point_neg);

    /* galois action x->1/x: only fixed points are 1 and -1 */
    galois_action const Ginv("1/y");
    auto is_fixed_point_inv = [](unsigned long r, unsigned long p) {
            return r == 1 || r == p-1;
    };
    ret &= test_galois_apply_one(Ginv, is_fixed_point_inv);

    /* galois action x->1-1/x: fixed points are such that r^2-r+1 mod p == 0 */
    galois_action const G31("autom3.1");
    auto is_fixed_point_31 = [](unsigned long r, unsigned long p) {
            /* only use with small values to avoid overflow */
            return ((r*r-r+1) % p) == 0;
    };
    ret &= test_galois_apply_one(G31, is_fixed_point_31);

    /* galois action x->-1-1/x: fixed points are such that r^2+r+1 mod p == 0 */
    galois_action const G32("autom3.2");
    auto is_fixed_point_32 = [](unsigned long r, unsigned long p) {
            /* only use with small values to avoid overflow */
            return ((r*r+r+1) % p) == 0;
    };
    ret &= test_galois_apply_one(G32, is_fixed_point_32);

    /* galois action x->-(x+1)/(x-1): fixed points are s.t. r^2+1 mod p == 0 */
    galois_action const G41("autom4.1");
    auto is_fixed_point_41 = [](unsigned long r, unsigned long p) {
            /* only use with small values to avoid overflow */
            return ((r*r+1) % p) == 0;
    };
    ret &= test_galois_apply_one(G41, is_fixed_point_41);

    /* galois action x->-(2*x+1)/(x-1): fixed points are s.t. r^2+r+1 mod p == 0 */
    galois_action const G61("autom6.1");
    auto is_fixed_point_61 = [](unsigned long r, unsigned long p) {
            /* only use with small values to avoid overflow */
            return ((r*r+r+1) % p) == 0;
    };
    ret &= test_galois_apply_one(G61, is_fixed_point_61);
    return ret;
}

static bool
test_galois_hash()
{
    const uint64_t CA = UINT64_C(314159265358979323);
    const uint64_t CB = UINT64_C(271828182845904523);

    bool ret = true;

#define TEST_HASH_INNER(g, a0, b0, a1, b1, op, inv_op_str) do {               \
            const uint64_t h0 = g.hash_ab(INT64_C(a0), UINT64_C(b0), CA, CB);       \
            const uint64_t h1 = g.hash_ab(INT64_C(a1), UINT64_C(b1), CA, CB);       \
            if (!(h0 op h1)) {                                                \
                std::cout << "error in " << __func__ << ": with " << g << ":" \
                          << " hash(" #a0 ", " #b0 ") " inv_op_str            \
                          << " hash(" #a1 ", " #b1 ")" << "\n";          \
                ret = false;                                                  \
            }                                                                 \
        }while(0)
#define TEST_HASH_EQ(g, a0, b0, a1, b1) \
        TEST_HASH_INNER(g, a0, b0, a1, b1, ==, "!=")
#define TEST_HASH_NEQ(g, a0, b0, a1, b1) \
        TEST_HASH_INNER(g, a0, b0, a1, b1, !=, "==")

    /*
     * galois action identity
     */
    galois_action const Gnone("none");
    TEST_HASH_NEQ(Gnone, -42, 17, -17, 42); /* h(a, b) != h(-b, a) */
    TEST_HASH_NEQ(Gnone, 42, 17, 17, 42);   /* h(a, b) != h(b, a) */
    TEST_HASH_NEQ(Gnone, -42, 17, 42, 17);  /* h(a, b) != h(-a, b) */

    /*
     * galois action x->-x
     */
    galois_action const Gneg("_y");
    TEST_HASH_NEQ(Gneg, -42, 17, -17, 42);  /* h(a, b) != h(-b, a) */
    TEST_HASH_NEQ(Gneg, 42, 17, 17, 42);    /* h(a, b) != h(b, a) */
    TEST_HASH_EQ(Gneg, -42, 17, 42, 17);    /* h(a, b) == h(-a, b) */

    /*
     * galois action x->1/x
     */
    galois_action const Ginv("1/y");
    TEST_HASH_EQ(Ginv, -42, 17, -17, 42);   /* h(a, b) == h(-b, a) */
    TEST_HASH_EQ(Ginv, 42, 17, 17, 42);     /* h(a, b) == h(b, a) */
    TEST_HASH_EQ(Ginv, -42, 17, -17, 42);   /* h(a, b) == h(-b, -a) */
    TEST_HASH_NEQ(Ginv, -42, 17, 42, 17);   /* h(a, b) != h(-a, b) */

    /*
     * galois action x->1-1/x
     */
    galois_action const G31("autom3.1");
    /* test case: a < 0 < b */
    TEST_HASH_EQ(G31, -17, 42, 42, 59);     /* h(a, b) == h(b, b-a) */
    TEST_HASH_EQ(G31, -17, 42, 59, 17);     /* h(a, b) == h(b-a, -a) */
    TEST_HASH_EQ(G31, -42, 17, 17, 59);     /* h(a, b) == h(b, b-a) */
    TEST_HASH_EQ(G31, -42, 17, 59, 42);     /* h(a, b) == h(b-a, -a) */
    TEST_HASH_EQ(G31, -1, 1, 1, 2);         /* h(a, b) == h(b, b-a) */
    TEST_HASH_EQ(G31, -1, 1, 2, 1);         /* h(a, b) == h(b-a, -a) */
    /* test case: 0 < a < b */
    TEST_HASH_EQ(G31, 17, 42, 42, 25);      /* h(a, b) == h(b, b-a) */
    TEST_HASH_EQ(G31, 17, 42, -25, 17);     /* h(a, b) == h(a-b, a) */
    /* test case: 0 < b < a */
    TEST_HASH_EQ(G31, 42, 17, -17, 25);     /* h(a, b) == h(-b, a-b) */
    TEST_HASH_EQ(G31, 42, 17, 25, 42);      /* h(a, b) == h(a-b, a) */

    TEST_HASH_NEQ(G31, -42, 17, 42, 17);    /* h(a, b) != h(-a, b) */
    TEST_HASH_NEQ(G31, 42, 17, 17, 42);     /* h(a, b) != h(b, a) */
    TEST_HASH_NEQ(G31, -42, 17, -17, 42);   /* h(a, b) != h(-b, -a) */

    /*
     * galois action x->-1-1/x
     */
    galois_action const G32("autom3.2");
    /* test case: 0 < a */
    TEST_HASH_EQ(G32, 17, 42, -42, 59);     /* h(a, b) == h(-b, a+b) */
    TEST_HASH_EQ(G32, 17, 42, -59, 17);     /* h(a, b) == h(-a-b, a) */
    TEST_HASH_EQ(G32, 42, 17, -17, 59);     /* h(a, b) == h(-b, a+b) */
    TEST_HASH_EQ(G32, 42, 17, -59, 42);     /* h(a, b) == h(-a-b, a) */
    TEST_HASH_EQ(G32, 1, 1, -1, 2);         /* h(a, b) == h(-b, a+b) */
    TEST_HASH_EQ(G32, 1, 1, -2, 1);         /* h(a, b) == h(-a-b, a) */
    /* test case -b < a < 0 */
    TEST_HASH_EQ(G32, -17, 42, -42, 25);    /* h(a, b) == h(-b, a+b) */
    TEST_HASH_EQ(G32, -17, 42, 25, 17);     /* h(a, b) == h(a+b, -a) */
    /* test case: a < -b < 0 */
    TEST_HASH_EQ(G32, -42, 17, 17, 25);     /* h(a, b) == h(b, -a-b) */
    TEST_HASH_EQ(G32, -42, 17, -25, 42);    /* h(a, b) == h(a+b, -a) */

    TEST_HASH_NEQ(G32, -42, 17, 42, 17);    /* h(a, b) != h(-a, b) */
    TEST_HASH_NEQ(G32, 42, 17, 17, 42);     /* h(a, b) != h(b, a) */
    TEST_HASH_NEQ(G32, -42, 17, -17, 42);   /* h(a, b) != h(-b, -a) */

    /*
     * galois action x->-(x+1)/(x-1)
     */
    galois_action const G41("autom4.1");
    /* test case: a < -b < 0 */
    TEST_HASH_EQ(G41, -42, 17, 59, 25);     /* h(a, b) == h(b-a, -(a+b)) */
    TEST_HASH_EQ(G41, -42, 17, 17, 42);     /* h(a, b) == h(b, -a) */
    TEST_HASH_EQ(G41, -42, 17, -25, 59);    /* h(a, b) == h(a+b, b-a) */
    /* test case: -b < a < 0 */
    TEST_HASH_EQ(G41, -17, 42, -59, 25);    /* h(a, b) == h(a-b, a+b) */
    TEST_HASH_EQ(G41, -17, 42, 42, 17);     /* h(a, b) == h(b, -a) */
    TEST_HASH_EQ(G41, -17, 42, 25, 59);     /* h(a, b) == h(a+b, b-a) */
    /* test case: 0 < a < b */
    TEST_HASH_EQ(G41, 17, 42, -25, 59);     /* h(a, b) == h(a-b, a+b) */
    TEST_HASH_EQ(G41, 17, 42, -42, 17);     /* h(a, b) == h(-b, a) */
    TEST_HASH_EQ(G41, 17, 42, 59, 25);      /* h(a, b) == h(a+b, b-a) */
    /* test case: 0 < b < a */
    TEST_HASH_EQ(G41, 42, 17, 25, 59);      /* h(a, b) == h(a-b, a+b) */
    TEST_HASH_EQ(G41, 42, 17, -17, 42);     /* h(a, b) == h(-b, a) */
    TEST_HASH_EQ(G41, 42, 17, -59, 25);     /* h(a, b) == h(-(a+b), a-b) */

    TEST_HASH_NEQ(G41, -42, 17, 42, 17);    /* h(a, b) != h(-a, b) */
    TEST_HASH_NEQ(G41, 42, 17, 17, 42);     /* h(a, b) != h(b, a) */
    TEST_HASH_NEQ(G41, -42, 17, -17, 42);   /* h(a, b) != h(-b, -a) */

    /*
     * galois action x->-(2*x+1)/(x-1)
     */
    galois_action const G61("autom6.1");
    /* test case: 0 < b < a */
    TEST_HASH_EQ(G61, 42, 17, 25, 76);    /* h(a, b) == h(a-b, a+2*b) */
    TEST_HASH_EQ(G61, 42, 17, -17, 59);   /* h(a, b) == h(-b, a+b) */
    TEST_HASH_EQ(G61, 42, 17, -76, 101);  /* h(a, b) == h(-a-2*b, 2*a+b) */
    TEST_HASH_EQ(G61, 42, 17, -59, 42);   /* h(a, b) == h(-a-b, a) */
    TEST_HASH_EQ(G61, 42, 17, -101, 25);  /* h(a, b) == h(-2*a-b, a-b) */
    /* test case: 0 < a < b */
    TEST_HASH_EQ(G61, 17, 42, -25, 101);  /* h(a, b) == h(a-b, a+2*b) */
    TEST_HASH_EQ(G61, 17, 42, -42, 59);   /* h(a, b) == h(-b, a+b) */
    TEST_HASH_EQ(G61, 17, 42, -101, 76);  /* h(a, b) == h(-a-2*b, 2*a+b) */
    TEST_HASH_EQ(G61, 17, 42, -59, 17);   /* h(a, b) == h(-a-b, a) */
    TEST_HASH_EQ(G61, 17, 42, 76, 25);    /* h(a, b) == h(2*a+b, -a+b) */
    /* test case: -b/2 < a < 0 */
    TEST_HASH_EQ(G61, -17, 42, -59, 67);  /* h(a, b) == h(a-b, a+2*b) */
    TEST_HASH_EQ(G61, -17, 42, -42, 25);  /* h(a, b) == h(-b, a+b) */
    TEST_HASH_EQ(G61, -17, 42, -67, 8);   /* h(a, b) == h(-a-2*b, 2*a+b) */
    TEST_HASH_EQ(G61, -17, 42, 25, 17);   /* h(a, b) == h(a+b, -a) */
    TEST_HASH_EQ(G61, -17, 42, 8, 59);    /* h(a, b) == h(2*a+b, -a+b) */
    /* test case: -b < a < -b/2 < 0 */
    TEST_HASH_EQ(G61, -25, 42, -67, 59);  /* h(a, b) == h(a-b, a+2*b) */
    TEST_HASH_EQ(G61, -25, 42, -42, 17);  /* h(a, b) == h(-b, a+b) */
    TEST_HASH_EQ(G61, -25, 42, 59, 8);    /* h(a, b) == h(a+2*b, -2*a-b) */
    TEST_HASH_EQ(G61, -25, 42, 17, 25);   /* h(a, b) == h(a+b, -a) */
    TEST_HASH_EQ(G61, -25, 42, -8, 67);   /* h(a, b) == h(2*a+b, -a+b) */
    /* test case: -2b < a < -b < 0 */
    TEST_HASH_EQ(G61, -25, 17, -42, 9);   /* h(a, b) == h(a-b, a+2*b) */
    TEST_HASH_EQ(G61, -25, 17, 17, 8);    /* h(a, b) == h(b, -a-b) */
    TEST_HASH_EQ(G61, -25, 17, 9, 33);    /* h(a, b) == h(a+2*b, -2*a-b) */
    TEST_HASH_EQ(G61, -25, 17, -8, 25);   /* h(a, b) == h(a+b, -a) */
    TEST_HASH_EQ(G61, -25, 17, -33, 42);  /* h(a, b) == h(2*a+b, -a+b) */
    /* test case: a < -2b < 0 */
    TEST_HASH_EQ(G61, -42, 17, 59, 8);    /* h(a, b) == h(-a+b, -a-2*b) */
    TEST_HASH_EQ(G61, -42, 17, 17, 25);   /* h(a, b) == h(b, -a-b) */
    TEST_HASH_EQ(G61, -42, 17, -8, 67);   /* h(a, b) == h(a+2*b, -2*a-b) */
    TEST_HASH_EQ(G61, -42, 17, -25, 42);  /* h(a, b) == h(a+b, -a) */
    TEST_HASH_EQ(G61, -42, 17, -67, 59);  /* h(a, b) == h(2*a+b, -a+b) */

    TEST_HASH_NEQ(G61, -42, 17, 42, 17);    /* h(a, b) != h(-a, b) */
    TEST_HASH_NEQ(G61, 42, 17, 17, 42);     /* h(a, b) != h(b, a) */
    TEST_HASH_NEQ(G61, -42, 17, -17, 42);   /* h(a, b) != h(-b, -a) */

#undef TEST_HASH_INNER
#undef TEST_HASH_EQ
#undef TEST_HASH_NEQ
    return ret;
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    unsigned long iter = 100;
    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
    tests_common_get_iter(&iter);

    bool all_ok = true;

    all_ok &= test_galois_apply();
    all_ok &= test_galois_hash();

    tests_common_clear();
    return all_ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
