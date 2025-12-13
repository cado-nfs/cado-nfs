#include "cado.h" // IWYU pragma: keep

#include <cstdint> /* AIX wants it first (it's a bug) */
#include <cstdlib>

#include <iostream>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>

#include "fmt/base.h"
#include <gmp.h>

#include "arithxx/mod64.hpp"
#include "arithxx/mod_mpz_new.hpp"
#include "arithxx/modint.hpp"
#include "arithxx/modredc126.hpp"
#include "arithxx/modredc96.hpp"
#include "arithxx/modredc64.hpp"
#include "cxx_mpz.hpp"
#include "macros.h"
#include "misc.h"
#include "tests_common.h"

#ifdef __GNUG__
#include <cxxabi.h>
#include <memory>

static std::string demangle(char const * name)
{

    int status = -4; // some arbitrary value to eliminate the compiler warning

    // enable c++11 by passing the flag -std=c++11 to g++
    std::unique_ptr<char, void (*)(void *)> const res {
        abi::__cxa_demangle(name, nullptr, nullptr, &status), std::free};

    return (status == 0) ? res.get() : name;
}

#else

// does nothing if not g++
static std::string demangle(char const * name)
{
    return name;
}

#endif

template <typename T> static std::string tname()
{
    return demangle(typeid(T).name());
}

template <typename T> static T randomInteger();

template <> Integer64 randomInteger<Integer64>()
{
    return Integer64(u64_random(state));
}

template <> Integer128 randomInteger<Integer128>()
{
    return { u64_random(state), u64_random(state) };
}

template <> cxx_mpz randomInteger<cxx_mpz>()
{
    uint64_t randomWords[10];
    size_t const len = u64_random(state) % 10 + 1;
    for (size_t i = 0; i < len; i++)
        randomWords[i] = u64_random(state);
    return {randomWords, len};
}


template <class layer> class Tests
{
  public:
    using Modulus = typename layer::Modulus;
    using Integer = typename layer::Integer;
    using Residue = typename layer::Residue;
    using ResidueOp = typename layer::Modulus::ResidueOp;

    Modulus randomModulus(bool odd = false) const;

    static cxx_mpz modToMpz(Modulus const & m)
    {
        Integer const n = m.getmod();
        return (cxx_mpz)n;
    }

    bool test_one_init(cxx_mpz const & n) const
    {
        if (!Modulus::valid(n))
            return true;
        Modulus const m {Integer {n}};
        Residue const r = m(1);

        Integer const one = m.get(r);
        if (one != 1) {
            std::cerr << tname<layer>()
                      << " Precomputed constant 1 wrong for modulus " << n
                      << "\n";
            std::cerr << one << "\n";
            return false;
        } else {
            return true;
        }
    }

    bool test_init() const
    {
        bool ok = true;
        ok &= test_one_init(727);

        cxx_mpz n = mpz_from_expression("2^64+1");
        for (uint64_t i = 0; i < 10; i++)
            ok &= test_one_init(n + 2 * i);

        if (ok)
            std::cout << "Tests<" << tname<layer>() << ">::tests_init() passed"
                      << "\n";
        return ok;
    }

    bool test_set(unsigned long const iter) const
    {
        for (unsigned long i = 0; i < iter; i++) {
            Modulus const m = randomModulus();
            Integer const a = randomInteger<Integer>();
            Residue const t = m(a);
            Integer const r = m.get(t);

            cxx_mpz N = modToMpz(m);
            auto const A = cxx_mpz(a);
            auto const R = cxx_mpz(r);
            cxx_mpz R1;
            mpz_mod(R1, A, N);
            if (R1 != R) {
                std::cerr << tname<layer>() << "::set(" << A
                          << ") wrong for modulus " << N << "\n";
                return false;
            }
        }
        return true;
    }

    bool test_neg() const
    {
        Modulus const m = randomModulus();
        Residue r(m);
        m.neg(r, r);
        Integer a = m.get(r);
        if (a != 0) {
            cxx_mpz const N = modToMpz(m);
            std::cerr << tname<layer>() << "::neg(0) wrong for modulus " << N
                      << "\n";
            return false;
        }
        Integer const M = m.getmod();

        r = m(1);
        m.neg(r, r);
        a = m.get(r);
        if (a != M - 1) {
            cxx_mpz const N = modToMpz(m);
            std::cerr << tname<layer>() << "::neg(1) wrong for modulus " << N
                      << "\n";
            return false;
        }
        return true;
    }

    bool test_one_sqr(Integer const & i, Integer const & n) const
    {
        Modulus const m(n);
        Residue a = m(i);
        m.sqr(a, a);
        Integer const r = m.get(a);
        auto const A = cxx_mpz(i);
        auto const N = cxx_mpz(n);
        auto const R = (A * A) % N;
        if (R != (cxx_mpz)r) {
            std::cerr << tname<layer>() << "::sqr(" << A << ") mod " << n
                      << " wrong result: " << r << "\n";
            return false;
        }
        return true;
    }

    bool test_one_mul(Integer const & i, Integer const & j,
                      Integer const & n) const
    {
        Modulus const m(n);
        Residue a = m(i), b = m(j);
        m.mul(a, a, b);
        Integer const r = m.get(a);
        auto const A = cxx_mpz(i);
        auto const B = cxx_mpz(j);
        auto const N = cxx_mpz(n);
        auto const R = (A * B) % N;
        if (R != (cxx_mpz)r) {
            std::cerr << tname<layer>() << "::mul(" << A << ", " << B
                      << ") mod " << N << " wrong result: " << r << "\n";
            return false;
        }
        return true;
    }

    bool test_one_mul(Integer const & n) const
    {
        if (!Modulus::valid(n))
            return true;

        bool ok = true;

        for (uint64_t i = 1; i < 10; i++) {
            ok &= test_one_sqr(Integer(i), n);
            ok &= test_one_sqr(n - Integer(i), n);
        }

        for (uint64_t i = 1; i < 10; i++) {
            for (uint64_t j = 1; j < 10; j++) {
                ok &= test_one_mul(Integer(i), Integer(j), n);
                ok &= test_one_mul(n - Integer(i), Integer(j), n);
                ok &= test_one_mul(Integer(i), n - Integer(j), n);
                ok &= test_one_mul(n - Integer(i), n - Integer(j), n);
            }
        }
        return ok;
    }

    bool test_one_add_sub(Modulus const & m) const
    {
        Integer const a = randomInteger<Integer>();
        Integer const b = randomInteger<Integer>();
        Residue const ar = m(a);
        Residue const br = m(b);
        Residue cr(m);
        {
            m.add(cr, ar, br);
            Integer const c = m.get(cr);
            if (cxx_mpz(c) != (cxx_mpz(a) + cxx_mpz(b)) % cxx_mpz(m.getmod()))
                return false;
        }
        {
            m.add1(cr, ar);
            Integer const c = m.get(cr);
            if (cxx_mpz(c) != (cxx_mpz(a) + 1) % cxx_mpz(m.getmod()))
                return false;
        }
        {
            m.sub(cr, ar, br);
            Integer const c = m.get(cr);
            cxx_mpz e = (cxx_mpz(a) - cxx_mpz(b)) % cxx_mpz(m.getmod());
            if (e < 0)
                e += cxx_mpz(m.getmod());
            if (cxx_mpz(c) != e)
                return false;
        }
        {
            m.sub1(cr, ar);
            Integer const c = m.get(cr);
            cxx_mpz e = (cxx_mpz(a) - 1) % cxx_mpz(m.getmod());
            if (e < 0)
                e += cxx_mpz(m.getmod());
            if (cxx_mpz(c) != e)
                return false;
        }
        return true;
    }
    bool test_add_sub(unsigned long iter = 10) const
    {
        bool ok = true;
        for (unsigned long i = 0; ok && i < iter; i++) {
            Modulus const m = randomModulus(true);
            ok = ok && test_one_add_sub(m);
            if (!ok) {
                std::cout << "Tests<" << tname<layer>() << ">::tests_add_sub() "
                          << "FAILED for m=" << m.getmod() << "\n";
            }
        }
        if (ok)
            std::cout << "Tests<" << tname<layer>()
                      << ">::tests_add_sub() passed\n";
        return ok;
    }

    bool test_mul() const
    {
        bool ok = true;
        ok &= test_one_mul(Integer(727));
        for (uint64_t i = 0; i < 10; i++) {
            ok &= test_one_mul(Integer(UINT64_MAX - i));
        }
        cxx_mpz const a = mpz_from_expression("2^64 + 1");
        if (Modulus::valid(a)) {
            auto const n = Integer(a);
            ok &= test_one_mul(n);
            for (uint64_t i = 0; i < 10; i++) {
                ok &= test_one_mul(n + i);
            }
        }
        if (ok) {
            std::cout << "Tests<" << tname<layer>() << ">::tests_mul() passed"
                      << "\n";
        }
        return ok;
    }

    bool cmp_one_divn(Residue const & r, Residue const & a, uint64_t const b,
                      Modulus const & m) const
    {
        /* Check that r * b == a */
        Residue const _b = m(b);
        Residue rb(m);
        m.mul(rb, r, _b);
        if (!m.equal(a, rb)) {
            Integer const & mi = m.getmod();
            Integer const ai = m.get(a);
            Integer const ri = m.get(r);
            std::cerr << tname<layer>() << "::div" << b << "(" << ai
                      << ") = " << ri << " wrong for modulus " << mi << "\n";
            return false;
        }
        return true;
    }

    bool test_one_divn(Residue const & a, Modulus const & m) const
    {
        bool ok = true;
        Residue r(m);

        if (m.div2(r, a))
            ok &= cmp_one_divn(r, a, 2, m);
        if (m.div3(r, a))
            ok &= cmp_one_divn(r, a, 3, m);
        if (m.div5(r, a))
            ok &= cmp_one_divn(r, a, 5, m);
        if (m.div7(r, a))
            ok &= cmp_one_divn(r, a, 7, m);
        if (m.div11(r, a))
            ok &= cmp_one_divn(r, a, 11, m);
        if (m.div13(r, a))
            ok &= cmp_one_divn(r, a, 13, m);

        return ok;
    }

    bool test_divn(unsigned long const iter) const
    {
        bool ok = true;
        for (unsigned long i = 0; i < iter; i++) {
            Modulus const m = randomModulus();
            Integer const ai = randomInteger<Integer>();
            Residue const a = m(ai);
            ok &= test_one_divn(a, m);
        }
        if (ok) {
            std::cout << "Tests<" << tname<layer>() << ">::test_divn() passed"
                      << "\n";
        }
        return ok;
    }

    cxx_mpz compute_mpz_pow(Integer const & base, Integer const & exponent,
                            Integer const & n) const
    {
        cxx_mpz B = (cxx_mpz)base, E = (cxx_mpz)exponent, N = (cxx_mpz)n, R;
        mpz_powm(R, B, E, N);
        return R;
    }

    cxx_mpz compute_mpz_pow(Integer const & base, uint64_t const * exponent,
                            size_t const len, Integer const & n) const
    {
        cxx_mpz R;
        mpz_powm(R, cxx_mpz(base), cxx_mpz(exponent, len), cxx_mpz(n));
        return R;
    }

    bool test_one_pow(Integer const & base, uint64_t const * exponent,
                      size_t const len, Integer const & n) const
    {
        if (!Modulus::valid(n))
            return true;
        Modulus const m(n);
        Residue p(m);
        Integer r, e;
        cxx_mpz const R = compute_mpz_pow(base, exponent, len, n);

        auto const ez = cxx_mpz(exponent, len);
        bool const eFitsInt = Modulus::valid(ez);
        if (eFitsInt)
            e = Integer(ez);

        Residue const b = m(base);
        m.pow(p, b, exponent, len);
        r = m.get(p);
        if (R != (cxx_mpz)r) {
            std::cerr << tname<layer>() << "::pow(" << base << ", " << exponent
                      << ", " << len << ") mod " << n << " wrong result: " << r
                      << "\n";
            return false;
        }
        if (eFitsInt) {
            m.pow(p, b, e);
            r = m.get(p);
            if (R != (cxx_mpz)r) {
                std::cerr << tname<layer>() << "::pow(" << base << ", " << e
                          << ") mod " << n << " wrong result: " << r << "\n";
                return false;
            }
        }
        if (len == 1) {
            m.pow(p, b, exponent[0]);
            r = m.get(p);
            if (R != (cxx_mpz)r) {
                std::cerr << tname<layer>() << "::pow(" << base << ", "
                          << exponent << ") mod " << n << " wrong result: " << r
                          << "\n";
                return false;
            }
        }
        if (base == Integer(2)) {
            m.pow2(p, exponent, len);
            r = m.get(p);
            if (R != (cxx_mpz)r) {
                std::cerr << tname<layer>() << "::pow2(" << exponent << ", "
                          << len << ") mod " << n << " wrong result: " << r
                          << "\n";
                return false;
            }
            if (eFitsInt) {
                m.pow2(p, e);
                r = m.get(p);
                if (R != (cxx_mpz)r) {
                    std::cerr << tname<layer>() << "::pow2(" << e << ") mod "
                              << n << " wrong result: " << r << "\n";
                    return false;
                }
            }
            if (len == 1) {
                m.pow2(p, exponent[0]);
                r = m.get(p);
                if (R != (cxx_mpz)r) {
                    std::cerr << tname<layer>() << "::pow2(" << exponent
                              << ") mod " << n << " wrong result: " << r
                              << "\n";
                    return false;
                }
            }
        }
        return true;
    }

    bool test_pow(unsigned long const iter) const
    {
        uint64_t e = 3;
        bool ok = true;
        ok &= test_one_pow(Integer(2), &e, 1, Integer(727));
        for (e = 0; e < 20; e++) {
            for (uint64_t b = 2; b < 20; b++) {
                ok &= test_one_pow(Integer(b), &e, 1, Integer(727));
            }
        }
        cxx_mpz const a = mpz_from_expression("2^64 + 1");
        if (Modulus::valid(a)) {
            auto const n = Integer(a);
            e = 3;
            ok &= test_one_pow(Integer(2), &e, 1, n);
        }

        for (unsigned long i = 0; i < iter; i++) {
            Modulus const m = randomModulus();
            Integer const b = randomInteger<Integer>();
            Integer const n = m.getmod();

            size_t const maxlen = 10;
            uint64_t e2[maxlen];
            size_t const len = u64_random(state) % (maxlen - 1) + 1;
            for (size_t i = 0; i < len; i++)
                e2[i] = u64_random(state);
            ok &= test_one_pow(Integer(2), e2, 1, n);
            ok &= test_one_pow(b, e2, 1, n);
            if (len >= 2) {
                ok &= test_one_pow(Integer(2), e2, 2, n);
                ok &= test_one_pow(b, e2, 2, n);
            }
            ok &= test_one_pow(Integer(2), e2, len, n);
            ok &= test_one_pow(b, e2, len, n);
        }

        if (ok) {
            std::cout << "Tests<" << tname<layer>() << ">::test_pow() passed"
                      << "\n";
        }
        return ok;
    }

    bool test_one_sprp(cxx_mpz const & n, bool const isPrime) const
    {
        char const * prime_str[2] = {"composite", "prime"};
        if (!Modulus::valid(n))
            return true;
        auto const N = Integer(n);
        Modulus const m(N);
        if (m.is_strong_pseudoprime_base2() != isPrime) {
            fmt::print(stderr,
                       "{} incorrectly declared {} by "
                       "{}::is_strong_pseudoprime_base2()\n",
                       n, prime_str[!isPrime], tname<layer>());
            return false;
        }
        return true;
    }

    bool test_one_lucas_sprp(cxx_mpz const & n, bool const isPrime) const
    {
        char const * prime_str[2] = {"composite", "prime"};
        if (!Modulus::valid(n))
            return true;
        Modulus const m {Integer {n}};
        if (m.is_strong_lucas_pseudoprime() != isPrime) {
            fmt::print(stderr,
                       "{} incorrectly declared {} by "
                       "{}::is_strong_lucas_pseudoprime()\n",
                       n, prime_str[!isPrime], tname<layer>());
            return false;
        }
        return true;
    }

    bool test_sprp(unsigned long const iter MAYBE_UNUSED) const
    {
        bool ok = true;
        ok &= test_one_sprp((15841), true);                       // psp2
        ok &= test_one_sprp((15859), true);                       // prime
        ok &= test_one_sprp((31631), false);                      // lpsp
        ok &= test_one_sprp((75350936251), true);                 // psp2
        ok &= test_one_sprp((75350936291), true);                 // prime
        ok &= test_one_sprp(mpz_from_expression("2^64+1"), true); // psp2
        if (ok) {
            std::cout << "Tests<" << tname<layer>() << ">::tests_sprp() passed"
                      << "\n";
        }
        return ok;
    }

    bool test_lucas_sprp(unsigned long const iter MAYBE_UNUSED) const
    {
        bool ok = true;
        ok &= test_one_lucas_sprp((15841), false);       // psp2
        ok &= test_one_lucas_sprp((15859), true);        // prime
        ok &= test_one_lucas_sprp((31631), true);        // lpsp
        ok &= test_one_lucas_sprp((75350936251), false); // psp2
        ok &= test_one_lucas_sprp((75350936291), true);  // prime
        ok &= test_one_lucas_sprp(mpz_from_expression("2^64+1"), false); // psp2
        if (ok) {
            std::cout << "Tests<" << tname<layer>()
                      << ">::tests_lucas_sprp() passed" << "\n";
        }
        return ok;
    }

    bool test_one_isprime(Modulus const & m, bool const isPrime) const
    {
        char const * prime_str[2] = {"composite", "prime"};
        if (m.is_prime() != isPrime) {
            std::cerr << modToMpz(m) << " incorrectly declared "
                      << prime_str[!isPrime] << " by " << tname<layer>()
                      << "::isprime()" << "\n";
            return false;
        }
        return true;
    }

    bool test_isprime(unsigned long const iter) const
    {
        bool ok = true;
        for (unsigned long i = 0; i < iter; i++) {
            Modulus const m = randomModulus();
            cxx_mpz M = modToMpz(m);
            bool const gmpIsPrime = mpz_probab_prime_p(M, 10);
            ok &= test_one_isprime(m, gmpIsPrime);
        }
        if (ok) {
            std::cout << "Tests<" << tname<layer>()
                      << ">::test_isprime() passed" << "\n";
        }
        return ok;
    }

    bool test_one_gcd(Integer const & i, Modulus const & m) const
    {
        bool ok = true;
        Integer g;
        Residue const r = m(i);
        m.gcd(g, r);

        cxx_mpz I = (cxx_mpz)i, M = modToMpz(m), G;
        mpz_gcd(G, I, M);

        if (G != (cxx_mpz)g) {
            std::cerr << tname<layer>() << "::gcd(" << I << ", " << M
                      << ") = " << g << " wrong" << "\n";
            ok = false;
        }

        return ok;
    }

    bool test_gcd(unsigned long const iter) const
    {
        bool ok = true;

        Modulus const m = randomModulus();
        ok &= test_one_gcd(Integer(0), m);
        ok &= test_one_gcd(Integer(1), m);
        Integer const a = m.getmod();
        ASSERT_ALWAYS(a > 0);
        ok &= test_one_gcd(a - 1, m);

        for (unsigned long i = 0; i < iter; i++) {
            Modulus const m = randomModulus();
            Integer const a = randomInteger<Integer>();
            ok &= test_one_gcd(a, m);
        }
        if (ok) {
            std::cout << "Tests<" << tname<layer>() << ">::test_gcd() passed"
                      << "\n";
        }
        return ok;
    }

    bool test_one_inv(Integer const & a, Modulus const & m) const
    {
        bool ok = true;

        cxx_mpz A = (cxx_mpz)a, M = modToMpz(m), I;
        bool const invExists0 = mpz_invert(I, A, M);

        Residue const ar = m(a);
        Residue ir(m);
        bool const invExists1 = m.inv(ir, ar);
        Integer const i = m.get(ir);

        if (invExists0 != invExists1) {
            std::cerr << tname<layer>() << "::inv(" << A << ", " << M
                      << ") wrongly thinks inverse "
                      << (invExists1 ? "exists" : "does not exist") << "\n";
            ok = false;
        } else if (invExists1 && I != (cxx_mpz)i) {
            std::cerr << tname<layer>() << "::inv(" << A << ", " << M
                      << ") = " << i << " wrong" << "\n";
            ok = false;
        }

        if (mpz_odd_p(A)) {
            bool const invExists2 = m.inv_odd(ir, ar);
            Integer const i = m.get(ir);
            if (invExists2 != invExists0) {
                std::cerr << tname<layer>() << "::inv_odd(" << A << ", " << M
                          << ") wrongly thinks inverse "
                          << (invExists2 ? "exists" : "does not exist") << "\n";
                ok = false;
            } else if (invExists0 && I != (cxx_mpz)i) {
                std::cerr << tname<layer>() << "::inv_odd(" << A << ", " << M
                          << ") = " << i << " wrong" << "\n";
                ok = false;
            }
        } else if (A != 0 && (A & (A - 1)) == 0) {
            bool const invExists3 = m.inv_powerof2(ir, ar);
            Integer const i = m.get(ir);
            if (invExists3 != invExists0) {
                std::cerr << tname<layer>() << "::inv_powerof2(" << A << ", "
                          << M << ") wrongly thinks inverse "
                          << (invExists3 ? "exists" : "does not exist") << "\n";
                ok = false;
            } else if (invExists0 && I != (cxx_mpz)i) {
                std::cerr << tname<layer>() << "::inv_powerof2(" << A << ", "
                          << M << ") = " << i << " wrong" << "\n";
                ok = false;
            }
        }

        {
            Integer i;
            bool const invExists5 = m.intinv(i, a);
            if (invExists5 != invExists0) {
                std::cerr << tname<layer>() << "::intinv(" << A << ", " << M
                          << ") wrongly thinks inverse "
                          << (invExists5 ? "exists" : "does not exist") << "\n";
                ok = false;
            } else if (invExists0 && I != (cxx_mpz)i) {
                std::cerr << tname<layer>() << "::intinv(" << A << ", " << M
                          << ") = " << i << " wrong" << "\n";
                ok = false;
            }
        }

        return ok;
    }

    bool test_inv(unsigned long const iter) const
    {
        bool ok = true;

        Modulus const m = randomModulus();
        ok &= test_one_inv(Integer(0), m);
        ok &= test_one_inv(Integer(1), m);
        Integer const & a = m.getmod();
        ASSERT_ALWAYS(a > 0);
        ok &= test_one_inv(a - 1, m);

        for (unsigned long i = 0; i < iter; i++) {
            Modulus const m = randomModulus();
            Integer const a = randomInteger<Integer>();
            ok &= test_one_inv(a, m);
        }
        if (ok) {
            std::cout << "Tests<" << tname<layer>() << ">::test_inv() passed"
                      << "\n";
        }
        return ok;
    }

    bool test_one_batchinv(Residue const * a, size_t const len,
                           Residue const * c, Modulus const & m) const
    {
        ASSERT_ALWAYS(len <= 10);
        auto r = m.template make_array<10>();
        Residue t(m);
        bool const batchinv_valid = m.batchinv(r.data(), a, len, c);
        bool only_trivial_gcds = true;

        Integer g;
        for (size_t i = 0; i < len; i++) {
            m.gcd(g, a[i]);
            if (g != 1) {
                only_trivial_gcds = false;
            }
        }

        if (batchinv_valid) {
            if (!only_trivial_gcds) {
                std::cerr << tname<layer>()
                          << "::batchinv() wrongly thinks inverse exists"
                          << "\n";
                return false;
            }
        } else if (only_trivial_gcds) {
            std::cerr << tname<layer>()
                      << "::batchinv() wrongly thinks inverse does not exist"
                      << "\n";
            return false;
        } else {
            /* batchinv reported no inverse exists and there is, in fact, a
             * non-trivial gcd. All is well. */
            return true;
        }

        for (size_t i = 0; i < len; i++) {
            m.inv(t, a[i]);
            if (c != nullptr) {
                m.mul(t, t, *c);
            }
            if (!m.equal(t, r[i])) {
                std::cerr << tname<layer>()
                          << "::batchinv() computed wrong inverse" << "\n";
                return false;
            }
        }
        return true;
    }

    bool test_batchinv(unsigned long const iter) const
    {
        bool ok = true;

        for (unsigned long i_test = 0; i_test < iter; i_test++) {
            Modulus const m = randomModulus();
            auto a = m.template make_array<10>();
            Residue c(m);
            for (auto & x: a)
                x = m(randomInteger<Integer>());
            c = m(randomInteger<Integer>());
            for (size_t i_size = 0; i_size <= 10; i_size++) {
                ok &= test_one_batchinv(a.data(), i_size, nullptr, m);
                ok &= test_one_batchinv(a.data(), i_size, &c, m);
            }
        }
        if (ok) {
            std::cout << "Tests<" << tname<layer>()
                      << ">::test_batchinv() passed" << "\n";
        }
        return ok;
    }

    bool test_one_batch_Q_to_Fp(Integer const & num, Integer const & den,
                                int k) const
    {
        if (!Modulus::valid(den))
            return true;

        static constexpr int const batch_size = 4;
        std::vector<uint64_t> M;
        M.reserve(batch_size);
        for (int i = 0; i < batch_size; i++) {
            cxx_mpz p;
            do {
                p = cxx_mpz(randomModulus(true).getmod());
            } while (!mpz_probab_prime_p(p, 2));
            M.emplace_back(p);
        }
        auto R = Modulus::batch_Q_to_Fp(Integer(num), Integer(den), k, M);
        if (R.empty())
            return true;
        bool ok = true;
        for (int i = 0; ok && i < batch_size; i++) {
            cxx_mpz t = cxx_mpz(R[i]) * cxx_mpz(den);
            mpz_mul_2exp(t, t, k);
            t -= cxx_mpz(num);
            t %= M[i];
            ok = t == 0;
            if (!ok) {
                fmt::print(stderr,
                           "Tests<{}>::test_one_batch_Q_to_Fp fails on "
                           "{}/{}/2^{} mod {}\n",
                           tname<layer>(), cxx_mpz(num), cxx_mpz(den), k, M[i]);
            }
        }
        return ok;
    }

    template <typename T, typename = void>
    struct has_batch_Q_to_Fp : std::false_type { };

    template <typename T>
    struct has_batch_Q_to_Fp<
        T, std::void_t<decltype(&T::Modulus::batch_Q_to_Fp)>>
        : std::true_type {
    };
    template <typename T>
    static constexpr bool has_batch_Q_to_Fp_v = has_batch_Q_to_Fp<T>::value;

#if 0
    template<typename L>
        bool
        test_batch_Q_to_Fp(const unsigned long) const
        requires(!std::is_same_v<L, arithxx_modredc64>
        { return true; }

    template<typename L>
        typename std::enable_if<std::is_same<L, arithxx_modredc64>::value, bool>::type
#else
    template <typename L>
    bool
    test_batch_Q_to_Fp(const unsigned long) const
    requires(!has_batch_Q_to_Fp_v<L>)
    {
        return true;
    }

    template <typename L>
        bool
#endif
    test_batch_Q_to_Fp(const unsigned long iter) const
    requires has_batch_Q_to_Fp_v<L>
    {
        bool ok = true;
        ok &= test_one_batch_Q_to_Fp(Integer(42), Integer(1009), 0);
        ok &= test_one_batch_Q_to_Fp(Integer(42), Integer(17), 0);
        ok &= test_one_batch_Q_to_Fp(Integer(42), Integer(17), 3);
        unsigned long spin = 0;
        for (unsigned long i_test = 0; ok && i_test < iter + spin; i_test++) {
            Integer const num = randomInteger<Integer>();
            Integer const den = randomModulus(true).getmod();
            if (!Modulus::valid(den)) {
                spin++;
                continue;
            }
            int const k = gmp_urandomm_ui(state, 10);
            ok &= test_one_batch_Q_to_Fp(num, den, k);
        }
        if (ok)
            std::cout << "Tests<" << tname<layer>()
                      << ">::test_batch_Q_to_Fp() passed" << "\n";
        return ok;
    }

    bool test_one_jacobi(Integer const & a, Modulus const & m) const
    {
        bool ok = true;
        Residue const ar = m(a);
        int const jacobi1 = m.jacobi(ar);

        cxx_mpz A = (cxx_mpz)a, M = modToMpz(m);
        int const jacobi2 = mpz_jacobi(A, M);

        if (jacobi1 != jacobi2) {
            std::cerr << tname<layer>() << "::jacobi(" << A << ", " << M
                      << ") = " << jacobi1 << " wrong" << "\n";
            ok = false;
        }

        return ok;
    }

    bool test_jacobi(unsigned long const iter) const
    {
        bool ok = true;

        if (Modulus::valid(Integer(1))) {
            Modulus const m(Integer(1));
            ok &= test_one_jacobi(Integer(0), m);
        }

        Modulus const m = randomModulus(true);
        ok &= test_one_jacobi(Integer(0), m);
        ok &= test_one_jacobi(Integer(1), m);
        Integer const a = m.getmod();
        ASSERT_ALWAYS(a > 0);
        ok &= test_one_jacobi(a - 1, m);

        for (unsigned long i = 0; i < iter; i++) {
            Modulus const m = randomModulus(true);
            Integer const a = randomInteger<Integer>();
            ok &= test_one_jacobi(a, m);
        }
        if (ok) {
            std::cout << "Tests<" << tname<layer>() << ">::test_jacobi() passed"
                      << "\n";
        }
        return ok;
    }

    bool test_modop() const
    {
        Integer n;
        {
            Modulus const m = randomModulus();
            ResidueOp r(m), s(m), t(m);
            n = randomInteger<Integer>();

            r = n;
            s = 2;
            t = r * s;
            if (t != r + r) {
                return false;
            }
            t += n;
            if (t != r + r + r) {
                return false;
            }
        }
        n = 727;
        if (Modulus::valid(n)) {
            Modulus const m(n);
            ResidueOp r(m);

            r = 2;
            r = r.pow(n - 1);
            if ((Integer)r != 1)
                return false;
            r = 5;
            r = 5 / r;
            if ((Integer)r != 1)
                return false;
        }
        std::cout << "Tests<" << tname<layer>() << ">::test_modop() passed"
                  << "\n";
        return true;
    }

    bool runTests(unsigned long const iter) const
    {
        bool ok = true;
        ok &= test_init();
        ok &= test_set(iter);
        ok &= test_neg();
        ok &= test_add_sub(iter);
        ok &= test_mul();
        ok &= test_divn(iter);
        ok &= test_pow(iter);
        ok &= test_sprp(iter);
        ok &= test_lucas_sprp(iter);
        ok &= test_isprime(iter);
        ok &= test_gcd(iter);
        ok &= test_inv(iter);
        ok &= test_jacobi(iter);
        ok &= test_batchinv(iter);
        ok &= test_batch_Q_to_Fp<layer>(iter);
        ok &= test_modop();
        return ok;
    }
};

template <>
arithxx_mod64::Modulus Tests<arithxx_mod64>::randomModulus(bool const odd) const
{
    return Modulus(u64_random(state) | (odd ? 1 : 0));
}

template <>
arithxx_modredc64::Modulus
Tests<arithxx_modredc64>::randomModulus(bool const odd MAYBE_UNUSED) const
{
    return Modulus(u64_random(state) | 1);
}

template <>
arithxx_modredc96::Modulus
Tests<arithxx_modredc96>::randomModulus(bool const odd MAYBE_UNUSED) const
{
    Integer128 const m(u64_random(state) | 1,
                       u64_random(state) & (UINT64_MAX >> 32));
    return Modulus(m);
}

template <>
arithxx_modredc126::Modulus
Tests<arithxx_modredc126>::randomModulus(bool const odd MAYBE_UNUSED) const
{
    Integer128 const m(u64_random(state) | 1,
                       u64_random(state) & (UINT64_MAX >> 2));
    return Modulus(m);
}

template <>
arithxx_mod_mpz_new::Modulus
Tests<arithxx_mod_mpz_new>::randomModulus(bool const odd) const
{
    Integer i = randomInteger<arithxx_mod_mpz_new::Integer>();
    if (odd && mpz_even_p(i))
        i += 1;
    return Modulus(i);
}

int main(int argc, char const * argv[])
{
    unsigned long iter = 100;
    bool ok = true;

    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
    tests_common_get_iter(&iter);

    Tests<arithxx_mod64> const test1;
    ok &= test1.runTests(iter);

    Tests<arithxx_modredc64> const test2;
    ok &= test2.runTests(iter);

    Tests<arithxx_modredc96> const test3a;
    ok &= test3a.runTests(iter);

    Tests<arithxx_modredc126> const test3;
    ok &= test3.runTests(iter);

    Tests<arithxx_mod_mpz_new> const test4;
    ok &= test4.runTests(iter);

    tests_common_clear();
    exit(ok ? EXIT_SUCCESS : EXIT_FAILURE);
}
