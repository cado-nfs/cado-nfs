#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdint>

#include <vector>

#include "fmt/base.h"

#include "gmp_aux.h"
#include "tests_common.h"
#include "misc.h"
#include "galois_action.hpp"
#include "special-q.hpp"
#include "arithxx/mod64.hpp"
#include "gcd.h"

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
            special_q q0 { p, r, 1 };
            special_q q1 = G.apply(q0);
            if (q1.r != sigma_r) {
                ret = false;
                fmt::print("error in {}:"
                        " G.apply({})={} while G.apply({},{})=({},{})\n",
                        __func__, q0, q1, r, p, sigma_r, p);
            }
            if (is_fixed_point(r, p)) {
                if (sigma_r != r) {
                    ret = false;
                    fmt::print("error in {}:"
                            " G.apply({},{}) should be a fixed point"
                            " but we got {}\n",
                            __func__, r, p, sigma_r);
                }
                if (q1.r != q0.r) {
                    ret = false;
                    fmt::print("error in {}:"
                            " G.apply({}) should be a fixed point"
                            " but we got {}\n",
                            __func__, q0, q1);
                }
            } else {
                /* Upper bound the nb of iter by 1000 to avoid endless loop that
                 * could occurs with some errors.
                 */
                unsigned int n = 1;
                while (sigma_r != r && q1.r != q0.r && n < 1000) {
                    sigma_r = G.apply(sigma_r, p);
                    q1 = G.apply(q1);
                    n++;
                }
                if (n != G.get_order()) {
                    ret = false;
                    fmt::print("error in {}:"
                            " orbit of {} should have length {}, not {}\n",
                            __func__, q0, G.get_order(), n);
                }
            }
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

static bool test_hash(galois_action const & g,
        int64_t a0, uint64_t b0,
        int64_t a1, uint64_t b1,
        bool eq = true)
{
    const uint64_t CA = UINT64_C(314159265358979323);
    const uint64_t CB = UINT64_C(271828182845904523);

    const uint64_t h0 = g.hash_ab(a0, b0, CA, CB);
    const uint64_t h1 = g.hash_ab(a1, b1, CA, CB);
    if ((h0 == h1) != eq) {
        fmt::print("error in test_hash with {}: hash({},{}) {} hash({},{})\n",
                g, a0, b0, eq ? "!=" : "==", a1, b1);
        return false;
    }
    return true;
}

/* form a random set of orbits,  make sure that hashes are consistent
 * along the way
 */
static bool test_orbits(galois_action const & g)
{
    const uint64_t CA = UINT64_C(314159265358979323);
    const uint64_t CB = UINT64_C(271828182845904523);

    for(int i = 0 ; i < 10 ; i++) {
        auto const n = g.get_order();

        int64_t a0;
        uint64_t b0;
        uint64_t p;

        for( ;; ) {
            a0 = i64_random(state);
            b0 = u64_random(state) >> 1;
            if (n == 3 || n == 4) {
                /* make sure we can't overflow. Note that each iteration
                 * may double (a,b), but we have several in a cyle! */
                a0 /= 4;
                b0 /= 4;
            } else if (n == 6) {
                /* make sure we can't overflow. Note that each iteration
                 * may trible (a,b), but we have several in a cyle! */
                a0 /= 27;
                b0 /= 27;
            }

            /* create a random evaluation point */
            p = u64_random(state);

#if ULONG_BITS == 32
            /* since galois_action::apply takes and returns unsigned
             * longs, we're bound to use unsigned longs. We could just as
             * well let p and r actually *be* of unsigned long type, but
             * then we would run into the problem that arithxx only has
             * fixed-width types.
             */
            p >>= 32;
#endif

#if 0
            /* for easier gdb debugging at first */
            p = p % 1000;
            a0 >>= 52;
            b0 >>= 52;
#endif

            if (-p < 60) continue;

            if (gcd_uint64(safe_abs64(a0), b0) != 1) continue;

            break;
        }
        p = uint64_nextprime(p);
        const uint64_t r = u64_random(state) % p;
        const uint64_t s = g.apply(r, p);
        /* TODO: also check what happens for a projective root ? */
        const arithxx_mod64::Modulus pp(p);
        arithxx_mod64::Residue rr(pp);
        arithxx_mod64::Residue ss(pp);
        pp.set(rr, r);
        pp.set(ss, s);

        auto h0 = g.hash_ab(a0, b0, CA, CB);
        auto a = a0;
        auto b = b0;
        /* the Galois identity, after (a,b) -> (a',b'), should be
         * (a-b*r') * c = m * (a' - b' * r)
         *
         * (which should also be ok if c == 0)
         */
        arithxx_mod64::Residue cc(pp);
        pp.set(cc, static_cast<uint64_t>(g.apply_ab_cofactor(r, p)));

        for(unsigned int i = 0 ; i < n ; i++) {
            arithxx_mod64::Residue prev(pp);
            arithxx_mod64::Residue aa(pp);
            arithxx_mod64::Residue bb(pp);
            pp.set(aa, a);
            pp.set(bb, b);
            pp.mul(prev, bb, ss);
            pp.sub(prev, aa, prev);
            pp.mul(prev, prev, cc);
            int64_t m = g.apply_ab(a, b);
            auto h = g.hash_ab(a, b, CA, CB);
            arithxx_mod64::Residue rhs(pp);
            arithxx_mod64::Residue mm(pp);
            pp.set(mm, m);
            pp.set(aa, a);
            pp.set(bb, b);
            pp.mul(rhs, bb, rr);
            pp.sub(rhs, aa, rhs);
            pp.mul(rhs, rhs, mm);
            ASSERT_ALWAYS(pp.equal(prev, rhs));
            ASSERT_ALWAYS(h == h0);
        }
        ASSERT_ALWAYS(a == a0);
        ASSERT_ALWAYS(b == b0);
    }
    return true;
}

static bool
test_galois_hash()
{
    bool ret = true;

    /*
     * galois action identity
     */
    galois_action const Gnone("none");
    ret &= test_hash(Gnone, -42, 17, -17, 42, false); /* h(a, b) != h(-b, a) */
    ret &= test_hash(Gnone, 42, 17, 17, 42, false);   /* h(a, b) != h(b, a) */
    ret &= test_hash(Gnone, -42, 17, 42, 17, false);  /* h(a, b) != h(-a, b) */

    /*
     * galois action x->-x
     */
    galois_action const Gneg("_y");
    ret &= test_hash(Gneg, -42, 17, -17, 42, false);  /* h(a, b) != h(-b, a) */
    ret &= test_hash(Gneg, 42, 17, 17, 42, false);    /* h(a, b) != h(b, a) */
    ret &= test_hash(Gneg, -42, 17, 42, 17);    /* h(a, b) == h(-a, b) */

    /*
     * galois action x->1/x
     */
    galois_action const Ginv("1/y");
    ret &= test_hash(Ginv, -42, 17, -17, 42);   /* h(a, b) == h(-b, a) */
    ret &= test_hash(Ginv, 42, 17, 17, 42);     /* h(a, b) == h(b, a) */
    ret &= test_hash(Ginv, -42, 17, -17, 42);   /* h(a, b) == h(-b, -a) */
    ret &= test_hash(Ginv, -42, 17, 42, 17, false);   /* h(a, b) != h(-a, b) */

    /*
     * galois action x->1-1/x
     */
    galois_action const G31("autom3.1");
    /* test case: a < 0 < b */
    ret &= test_hash(G31, -17, 42, 42, 59);     /* h(a, b) == h(b, b-a) */
    ret &= test_hash(G31, -17, 42, 59, 17);     /* h(a, b) == h(b-a, -a) */
    ret &= test_hash(G31, -42, 17, 17, 59);     /* h(a, b) == h(b, b-a) */
    ret &= test_hash(G31, -42, 17, 59, 42);     /* h(a, b) == h(b-a, -a) */
    ret &= test_hash(G31, -1, 1, 1, 2);         /* h(a, b) == h(b, b-a) */
    ret &= test_hash(G31, -1, 1, 2, 1);         /* h(a, b) == h(b-a, -a) */
    /* test case: 0 < a < b */
    ret &= test_hash(G31, 17, 42, 42, 25);      /* h(a, b) == h(b, b-a) */
    ret &= test_hash(G31, 17, 42, -25, 17);     /* h(a, b) == h(a-b, a) */
    /* test case: 0 < b < a */
    ret &= test_hash(G31, 42, 17, -17, 25);     /* h(a, b) == h(-b, a-b) */
    ret &= test_hash(G31, 42, 17, 25, 42);      /* h(a, b) == h(a-b, a) */

    ret &= test_hash(G31, -42, 17, 42, 17, false);    /* h(a, b) != h(-a, b) */
    ret &= test_hash(G31, 42, 17, 17, 42, false);     /* h(a, b) != h(b, a) */
    ret &= test_hash(G31, -42, 17, -17, 42, false);   /* h(a, b) != h(-b, -a) */

    /*
     * galois action x->-1-1/x
     */
    galois_action const G32("autom3.2");
    /* test case: 0 < a */
    ret &= test_hash(G32, 17, 42, -42, 59);     /* h(a, b) == h(-b, a+b) */
    ret &= test_hash(G32, 17, 42, -59, 17);     /* h(a, b) == h(-a-b, a) */
    ret &= test_hash(G32, 42, 17, -17, 59);     /* h(a, b) == h(-b, a+b) */
    ret &= test_hash(G32, 42, 17, -59, 42);     /* h(a, b) == h(-a-b, a) */
    ret &= test_hash(G32, 1, 1, -1, 2);         /* h(a, b) == h(-b, a+b) */
    ret &= test_hash(G32, 1, 1, -2, 1);         /* h(a, b) == h(-a-b, a) */
    /* test case -b < a < 0 */
    ret &= test_hash(G32, -17, 42, -42, 25);    /* h(a, b) == h(-b, a+b) */
    ret &= test_hash(G32, -17, 42, 25, 17);     /* h(a, b) == h(a+b, -a) */
    /* test case: a < -b < 0 */
    ret &= test_hash(G32, -42, 17, 17, 25);     /* h(a, b) == h(b, -a-b) */
    ret &= test_hash(G32, -42, 17, -25, 42);    /* h(a, b) == h(a+b, -a) */

    ret &= test_hash(G32, -42, 17, 42, 17, false);    /* h(a, b) != h(-a, b) */
    ret &= test_hash(G32, 42, 17, 17, 42, false);     /* h(a, b) != h(b, a) */
    ret &= test_hash(G32, -42, 17, -17, 42, false);   /* h(a, b) != h(-b, -a) */

    /*
     * galois action x->-(x+1)/(x-1)
     */
    galois_action const G41("autom4.1");
    /* test case: a < -b < 0 */
    ret &= test_hash(G41, -42, 17, 59, 25);     /* h(a, b) == h(b-a, -(a+b)) */
    ret &= test_hash(G41, -42, 17, 17, 42);     /* h(a, b) == h(b, -a) */
    ret &= test_hash(G41, -42, 17, -25, 59);    /* h(a, b) == h(a+b, b-a) */
    /* test case: -b < a < 0 */
    ret &= test_hash(G41, -17, 42, -59, 25);    /* h(a, b) == h(a-b, a+b) */
    ret &= test_hash(G41, -17, 42, 42, 17);     /* h(a, b) == h(b, -a) */
    ret &= test_hash(G41, -17, 42, 25, 59);     /* h(a, b) == h(a+b, b-a) */
    /* test case: 0 < a < b */
    ret &= test_hash(G41, 17, 42, -25, 59);     /* h(a, b) == h(a-b, a+b) */
    ret &= test_hash(G41, 17, 42, -42, 17);     /* h(a, b) == h(-b, a) */
    ret &= test_hash(G41, 17, 42, 59, 25);      /* h(a, b) == h(a+b, b-a) */
    /* test case: 0 < b < a */
    ret &= test_hash(G41, 42, 17, 25, 59);      /* h(a, b) == h(a-b, a+b) */
    ret &= test_hash(G41, 42, 17, -17, 42);     /* h(a, b) == h(-b, a) */
    ret &= test_hash(G41, 42, 17, -59, 25);     /* h(a, b) == h(-(a+b), a-b) */

    ret &= test_hash(G41, -42, 17, 42, 17, false);    /* h(a, b) != h(-a, b) */
    ret &= test_hash(G41, 42, 17, 17, 42, false);     /* h(a, b) != h(b, a) */
    ret &= test_hash(G41, -42, 17, -17, 42, false);   /* h(a, b) != h(-b, -a) */

    /*
     * galois action x->-(2*x+1)/(x-1)
     */
    galois_action const G61("autom6.1");
    /* test case: 0 < b < a */
    ret &= test_hash(G61, 42, 17, 25, 76);    /* h(a, b) == h(a-b, a+2*b) */
    ret &= test_hash(G61, 42, 17, -17, 59);   /* h(a, b) == h(-b, a+b) */
    ret &= test_hash(G61, 42, 17, -76, 101);  /* h(a, b) == h(-a-2*b, 2*a+b) */
    ret &= test_hash(G61, 42, 17, -59, 42);   /* h(a, b) == h(-a-b, a) */
    ret &= test_hash(G61, 42, 17, -101, 25);  /* h(a, b) == h(-2*a-b, a-b) */
    /* test case: 0 < a < b */
    ret &= test_hash(G61, 17, 42, -25, 101);  /* h(a, b) == h(a-b, a+2*b) */
    ret &= test_hash(G61, 17, 42, -42, 59);   /* h(a, b) == h(-b, a+b) */
    ret &= test_hash(G61, 17, 42, -101, 76);  /* h(a, b) == h(-a-2*b, 2*a+b) */
    ret &= test_hash(G61, 17, 42, -59, 17);   /* h(a, b) == h(-a-b, a) */
    ret &= test_hash(G61, 17, 42, 76, 25);    /* h(a, b) == h(2*a+b, -a+b) */
    /* test case: -b/2 < a < 0 */
    ret &= test_hash(G61, -17, 42, -59, 67);  /* h(a, b) == h(a-b, a+2*b) */
    ret &= test_hash(G61, -17, 42, -42, 25);  /* h(a, b) == h(-b, a+b) */
    ret &= test_hash(G61, -17, 42, -67, 8);   /* h(a, b) == h(-a-2*b, 2*a+b) */
    ret &= test_hash(G61, -17, 42, 25, 17);   /* h(a, b) == h(a+b, -a) */
    ret &= test_hash(G61, -17, 42, 8, 59);    /* h(a, b) == h(2*a+b, -a+b) */
    /* test case: -b < a < -b/2 < 0 */
    ret &= test_hash(G61, -25, 42, -67, 59);  /* h(a, b) == h(a-b, a+2*b) */
    ret &= test_hash(G61, -25, 42, -42, 17);  /* h(a, b) == h(-b, a+b) */
    ret &= test_hash(G61, -25, 42, 59, 8);    /* h(a, b) == h(a+2*b, -2*a-b) */
    ret &= test_hash(G61, -25, 42, 17, 25);   /* h(a, b) == h(a+b, -a) */
    ret &= test_hash(G61, -25, 42, -8, 67);   /* h(a, b) == h(2*a+b, -a+b) */
    /* test case: -2b < a < -b < 0 */
    ret &= test_hash(G61, -25, 17, -42, 9);   /* h(a, b) == h(a-b, a+2*b) */
    ret &= test_hash(G61, -25, 17, 17, 8);    /* h(a, b) == h(b, -a-b) */
    ret &= test_hash(G61, -25, 17, 9, 33);    /* h(a, b) == h(a+2*b, -2*a-b) */
    ret &= test_hash(G61, -25, 17, -8, 25);   /* h(a, b) == h(a+b, -a) */
    ret &= test_hash(G61, -25, 17, -33, 42);  /* h(a, b) == h(2*a+b, -a+b) */
    /* test case: a < -2b < 0 */
    ret &= test_hash(G61, -42, 17, 59, 8);    /* h(a, b) == h(-a+b, -a-2*b) */
    ret &= test_hash(G61, -42, 17, 17, 25);   /* h(a, b) == h(b, -a-b) */
    ret &= test_hash(G61, -42, 17, -8, 67);   /* h(a, b) == h(a+2*b, -2*a-b) */
    ret &= test_hash(G61, -42, 17, -25, 42);  /* h(a, b) == h(a+b, -a) */
    ret &= test_hash(G61, -42, 17, -67, 59);  /* h(a, b) == h(2*a+b, -a+b) */

    ret &= test_hash(G61, -42, 17, 42, 17, false);    /* h(a, b) != h(-a, b) */
    ret &= test_hash(G61, 42, 17, 17, 42, false);     /* h(a, b) != h(b, a) */
    ret &= test_hash(G61, -42, 17, -17, 42, false);   /* h(a, b) != h(-b, -a) */

    test_orbits(Gnone);
    test_orbits(Gneg);
    test_orbits(Ginv);
    test_orbits(G31);
    test_orbits(G32);
    test_orbits(G41);
    test_orbits(G61);

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
