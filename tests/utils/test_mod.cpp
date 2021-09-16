#include "cado.h" // IWYU pragma: keep
#include <cstdint>     /* AIX wants it first (it's a bug) */
#include <cstdlib>
#include <gmp.h>
#include <iostream>
#include <typeinfo>
#include "cxx_mpz.hpp"
#include "tests_common.h"
#include "mod64.hpp"
#include "modredc64.hpp"
#include "modredc126.hpp"
#include "mod_mpz_new.hpp"
#include "modint.hpp"
#include "macros.h"
#include "misc.h"

template <typename T>
T randomInteger();

template<>
Integer64 randomInteger<Integer64>() {
    return Integer64(u64_random(state));
}

template<>
Integer128 randomInteger<Integer128>() {
    return Integer128(u64_random(state), u64_random(state));
}

template<>
cxx_mpz randomInteger<cxx_mpz>() {
    cxx_mpz r;
    uint64_t randomWords[10];
    size_t len = u64_random(state) % 10 + 1;
    for (size_t i = 0; i < len; i++)
        randomWords[i] = u64_random(state);
    const bool ok = r.set(randomWords, len);
    ASSERT_ALWAYS(ok);
    return r;
}

template <class T, typename = typename T::Residue>
class Tests {
public:
    typedef T Modulus;
    typedef typename T::Integer Integer;
    typedef typename T::Residue Residue;
    typedef typename T::ResidueOp ResidueOp;

    Modulus randomModulus(bool odd = false) const;
    
    static cxx_mpz modToMpz(const Modulus &m) {
        Integer n;
        m.getmod(n);
        return (cxx_mpz) n;
    }

    bool test_one_init(Integer n) const {
        if (!T::valid(n))
            return true;
        Modulus m(n);
        Residue r(m);
        
        m.set1(r);
        Integer one;
        m.get(one, r);
        if (one != 1){
            std::cerr << typeid(T).name() << " Precomputed constant 1 wrong for modulus " << n << std::endl;
            std::cerr << one << std::endl;
            return false;
        } else {
            return true;
        }
    }

    bool test_init() const {
        bool ok = true;
        ok &= test_one_init(Integer(727));

        const uint64_t a[] = {1, 1};
        Integer n;
        if (n.set(a, 2)) {
            for (uint64_t i = 0; i < 10; i++)
                ok &= test_one_init(n + 2*i);
        }
        if (ok) {
            std::cout << "tests<" << typeid(T).name() << ">::tests_init() passed" << std::endl;
        }
        return ok;
    }
    
    bool test_set(const unsigned long iter) const {
        for (unsigned long i = 0; i < iter; i++) {
            Modulus m = randomModulus();
            Integer a = randomInteger<Integer>();
            Residue t(m);
            m.set(t, a);
            Integer r;
            m.get(r, t);

            cxx_mpz N  = modToMpz(m);
            cxx_mpz A = (cxx_mpz) a;
            cxx_mpz R = (cxx_mpz) r;
            cxx_mpz R1;
            mpz_mod(R1, A, N);
            if (R1 != R) {
                std::cerr << typeid(T).name() << "::set(" << A << ") wrong for modulus " << N << std::endl;
                return false;
            }
        }
        return true;
    }
    
    bool test_neg() const {
        Modulus m = randomModulus();
        Integer a, M;
        Residue r(m);
        m.neg(r, r);
        m.get(a, r);
        if (a != 0) {
            cxx_mpz N  = modToMpz(m);
            std::cerr << typeid(T).name() << "::neg(0) wrong for modulus " << N << std::endl;
            return false;
        }
        m.getmod(M);
        a = 1;
        m.set(r, a);
        m.neg(r, r);
        m.get(a, r);
        if (a != M - 1) {
            cxx_mpz N  = modToMpz(m);
            std::cerr << typeid(T).name() << "::neg(1) wrong for modulus " << N << std::endl;
            return false;
        }
        return true;
    }

    bool test_one_sqr(Integer i, Integer &n) const {
        Modulus m(n);
        Residue a(m);
        Integer r;
        m.set(a, i);
        m.sqr(a, a);
        m.get(r, a);
        cxx_mpz A = (cxx_mpz) i, N = (cxx_mpz) n;
        cxx_mpz R;
        mpz_mul(R, A, A);
        mpz_mod(R, R, N);
        if (R != (cxx_mpz) r) {
            std::cerr << typeid(T).name() << "::sqr(" << A << ") mod " << n << " wrong result: " << r << std::endl;
            return false;
        }
        return true;
    }
    
    bool test_one_mul(Integer i, Integer j, Integer &n) const {
        Modulus m(n);
        Residue a(m), b(m);
        Integer r;
        m.set(a, i);
        m.set(b, j);
        m.mul(a, a, b);
        m.get(r, a);
        cxx_mpz A = (cxx_mpz) i, B = (cxx_mpz) j, N = (cxx_mpz) n;
        cxx_mpz R;
        mpz_mul(R, A, B);
        mpz_mod(R, R, N);
        if (R != (cxx_mpz) r) {
            std::cerr << typeid(T).name() << "::mul(" << A << ", " << B << ") mod " << N << " wrong result: " << r << std::endl;
            return false;
        }
        return true;
    }
    
    bool test_one_mul(Integer n) const {
        if (!T::valid(n))
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
    
    bool test_mul() const {
        Integer n;
        bool ok = true;
        ok &= test_one_mul(Integer(727));
        for (uint64_t i = 0; i < 10; i++) {
            ok &= test_one_mul(Integer(UINT64_MAX - i));
        }
        const uint64_t a[] = {1, 1};
        if (n.set(a, 2)) {
            ok &= test_one_mul(n);
            for (uint64_t i = 0; i < 10; i++) {
                ok &= test_one_mul(n + i);
            }
        }
        if (ok) {
            std::cout << "tests<" << typeid(T).name() << ">::tests_mul() passed" << std::endl;
        }
        return ok;
    }
    
    
    bool cmp_one_divn (const Residue &r, const Residue &a, const uint64_t b, const Modulus &m) const {
        /* Check that r * b == a */
        Residue _b(m), rb(m);
        m.set(_b, b);
        m.mul(rb, r, _b);
        if (!m.equal(a, rb)) {
            Integer ai, mi, ri;
            m.getmod(mi);
            m.get(ai, a);
            m.get(ri, r);
            std::cerr << typeid(T).name() << "::div" << b << "(" << ai << ") = " << ri << " wrong for modulus " << mi << std::endl;
            return false;
        }
        return true;
    }
    
    bool test_one_divn(const Residue &a, const Modulus &m) const {
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
    
    bool test_divn(const unsigned long iter) const {
        bool ok = true;
        for (unsigned long i = 0; i < iter; i++) {
            Modulus m = randomModulus();
            Integer ai = randomInteger<Integer>();
            Residue a(m);
            m.set(a, ai);
            ok &= test_one_divn(a, m);
        }
        if (ok) {
            std::cout << "tests<" << typeid(T).name() << ">::test_divn() passed" << std::endl;
        }
        return ok;
    }
    
    cxx_mpz compute_mpz_pow(const Integer &base, const Integer &exponent, const Integer &n) const {
        cxx_mpz B = (cxx_mpz) base, E = (cxx_mpz) exponent, N = (cxx_mpz) n, R;
        mpz_powm(R, B, E, N);
        return R;
    }

    cxx_mpz compute_mpz_pow(const Integer &base, const uint64_t *exponent, const size_t len, const Integer &n) const {
        cxx_mpz B = (cxx_mpz) base, E, R, N = (cxx_mpz) n;
        E.set(exponent, len);
        mpz_powm(R, B, E, N);
        return R;
    }

    bool test_one_pow(const Integer &base, const uint64_t *exponent, const size_t len, const Integer &n) const {
        if (!T::valid(n))
            return true;
        Modulus m(n);
        Residue p(m), b(m);
        Integer r, e;
        cxx_mpz R;
        R = compute_mpz_pow(base, exponent, len, n);
        
        const bool eFitsInt = e.set(exponent, len);
        
        m.set(b, base);
        m.pow(p, b, exponent, len);
        m.get(r, p);
        if (R != (cxx_mpz) r) {
            std::cerr << typeid(T).name() << "::pow(" << base << ", " << exponent << ", " << len << ") mod " << n << " wrong result: " << r << std::endl;
            return false;
        }
        if (eFitsInt) {
            m.pow(p, b, e);
            m.get(r, p);
            if (R != (cxx_mpz) r) {
                std::cerr << typeid(T).name() << "::pow(" << base << ", " << e << ") mod " << n << " wrong result: " << r << std::endl;
                return false;
            }
        }
        if (len == 1) {
            m.pow(p, b, exponent[0]);
            m.get(r, p);
            if (R != (cxx_mpz) r) {
                std::cerr << typeid(T).name() << "::pow(" << base << ", " << exponent << ") mod " << n << " wrong result: " << r << std::endl;
                return false;
            }
        }
        if (base == Integer(2)) {
            m.pow2(p, exponent, len);
            m.get(r, p);
            if (R != (cxx_mpz) r) {
                std::cerr << typeid(T).name() << "::pow2(" << exponent << ", " << len << ") mod " << n << " wrong result: " << r << std::endl;
                return false;
            }
            if (eFitsInt) {
                m.pow2(p, e);
                m.get(r, p);
                if (R != (cxx_mpz) r) {
                    std::cerr << typeid(T).name() << "::pow2(" << e << ") mod " << n << " wrong result: " << r << std::endl;
                    return false;
                }
            }
            if (len == 1) {
                m.pow2(p, exponent[0]);
                m.get(r, p);
                if (R != (cxx_mpz) r) {
                    std::cerr << typeid(T).name() << "::pow2(" << exponent << ") mod " << n << " wrong result: " << r << std::endl;
                    return false;
                }
            }
        }
        return true;
    }

    bool test_pow(const unsigned long iter) const {
        Integer n;
        uint64_t e = 3;
        bool ok = true;
        ok &= test_one_pow(Integer(2), &e, 1, Integer(727));
        for (e = 0; e < 20; e++) {
            for (uint64_t b = 2; b < 20; b++) {
                ok &= test_one_pow(Integer(b), &e, 1, Integer(727));
            }
        }
        uint64_t a1[] = {1,1};
        if (n.set(a1, 2)) {
            e = 3;
            ok &= test_one_pow(Integer(2), &e, 1, n);
        }

        for (unsigned long i = 0; i < iter; i++) {
            Modulus m = randomModulus();
            Integer b = randomInteger<Integer>();
            m.getmod(n);
            
            const size_t maxlen = 10;
            uint64_t e2[maxlen];
            size_t len = u64_random(state) % (maxlen - 1) + 1;
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
            std::cout << "tests<" << typeid(T).name() << ">::test_pow() passed" << std::endl;
        }        return ok;
    }
    
    bool test_one_sprp(const Modulus m, const bool isPrime) const {
        bool ok = true;
        const char *prime_str[2] = {"composite", "prime"};
        if (m.sprp2() != isPrime) {
            cxx_mpz N = modToMpz(m);
            std::cerr << N << " incorrectly declared " << prime_str[!isPrime] << " by " << typeid(T).name() << "::sprp2()" << std::endl;
            ok = false;
        }

        for (unsigned long b = 2; b < 10; b++) {
            Residue r(m);
            m.set(r, b);
            if (m.sprp(r) != isPrime) {
                cxx_mpz N = modToMpz(m);
                std::cerr << N << " incorrectly declared " << prime_str[!isPrime] << " by " << typeid(T).name() << "::sprp(" << b << ")" << std::endl;
                ok = false;
            }
        }
        return ok;
    }
    bool test_one_sprp(const uint64_t *a, const size_t l, const bool isPrime) const
    {
        Integer n;
        if (!n.set(a, l) || !T::valid(n)) {
            return true;
        }

        Modulus m(n);
        return test_one_sprp(m, isPrime);
    }

    bool test_sprp(const unsigned long iter MAYBE_UNUSED) const {
        bool ok = true;
        const uint64_t a1[1] = {727};
        const uint64_t a2[1] = {19*23};
        const uint64_t a3[2] = {13, 1};
        const uint64_t a4[2] = {27, 1};
        ok &= test_one_sprp(a1, sizeof(a1) / sizeof(a1[0]), true);
        ok &= test_one_sprp(a2, sizeof(a2) / sizeof(a2[0]), false);
        ok &= test_one_sprp(a3, sizeof(a3) / sizeof(a3[0]), true);
        ok &= test_one_sprp(a4, sizeof(a4) / sizeof(a4[0]), false);
        
        if (ok) {
            std::cout << "tests<" << typeid(T).name() << ">::tests_sprp() passed" << std::endl;
        }
        return ok;
    }
    
    bool test_one_isprime(const Modulus m, const bool isPrime) const {
        const char *prime_str[2] = {"composite", "prime"};
        if (m.isprime() != isPrime) {
            std::cerr << modToMpz(m) << " incorrectly declared " << prime_str[!isPrime] << " by " << typeid(T).name() << "::isprime()" << std::endl;
            return false;
        }
        return true;
    }
    
    bool test_isprime(const unsigned long iter) const {
        bool ok = true;
        for (unsigned long i = 0; i < iter; i++) {
            Modulus m = randomModulus();
            cxx_mpz M = modToMpz(m);
            bool gmpIsPrime = mpz_probab_prime_p(M, 10);
            ok &= test_one_isprime(m, gmpIsPrime);
        }
        if (ok) {
            std::cout << "tests<" << typeid(T).name() << ">::test_isprime() passed" << std::endl;
        }
        return ok;
    }

    bool test_one_gcd (const Integer &i, const Modulus m) const {
        bool ok = true;
        Integer g;
        Residue r(m);
        m.set(r, i);
        m.gcd(g, r);
        
        cxx_mpz I = (cxx_mpz) i, M = modToMpz(m), G;
        mpz_gcd(G, I, M);
        
        if (G != (cxx_mpz) g) {
            std::cerr << typeid(T).name() << "::gcd(" << I << ", " << M << ") = " << g << " wrong" << std::endl;
            ok = false;
        }
        
        return ok;
    }
    
    bool test_gcd(const unsigned long iter) const {
        bool ok = true;
        
        Modulus m = randomModulus();
        ok &= test_one_gcd(Integer(0), m);
        ok &= test_one_gcd(Integer(1), m);
        Integer a;
        m.getmod(a);
        ASSERT_ALWAYS(a > 0);
        ok &= test_one_gcd(a - 1, m);
        
        for (unsigned long i = 0; i < iter; i++) {
            Modulus m = randomModulus();
            Integer a = randomInteger<Integer>();
            ok &= test_one_gcd(a, m);
        }
        if (ok) {
            std::cout << "tests<" << typeid(T).name() << ">::test_gcd() passed" << std::endl;
        }
        return ok;
    }
    
    bool test_one_inv (const Integer &a, const Modulus m) const {
        bool ok = true;
        Integer i;
        Residue ar(m), ir(m);
        m.set(ar, a);
        bool invExists1 = m.inv(ir, ar);
        m.get(i, ir);
        
        cxx_mpz A = (cxx_mpz) a, M = modToMpz(m), I;
        int invExists2 = mpz_invert(I, A, M);

        if (invExists1 != (invExists2 != 0)) {
            std::cerr << typeid(T).name() << "::inv(" << A << ", " << M << ") wrongly thinks inverse " << (invExists1 ? "exists" : "does not exist") << std::endl;
            ok = false;
        } else if (invExists1 && I != (cxx_mpz) i) {
            std::cerr << typeid(T).name() << "::inv(" << A << ", " << M << ") = " << i << " wrong" << std::endl;
            ok = false;
        }
        
        return ok;
    }
    
    bool test_inv(const unsigned long iter) const {
        bool ok = true;
        
        Modulus m = randomModulus();
        ok &= test_one_inv(Integer(0), m);
        ok &= test_one_inv(Integer(1), m);
        Integer a;
        m.getmod(a);
        ASSERT_ALWAYS(a > 0);
        ok &= test_one_inv(a - 1, m);
        
        for (unsigned long i = 0; i < iter; i++) {
            Modulus m = randomModulus();
            Integer a = randomInteger<Integer>();
            ok &= test_one_inv(a, m);
        }
        if (ok) {
            std::cout << "tests<" << typeid(T).name() << ">::test_inv() passed" << std::endl;
        }
        return ok;
    }
    
    bool
    test_one_batchinv(const Residue *a, const size_t len, const Residue *c, const Modulus &m) const {
        ASSERT_ALWAYS(len <= 10);
        Residue r[10]{m, m, m, m, m, m, m, m, m, m};
        Residue t(m);
        const bool batchinv_valid = m.batchinv(r, a, len, c);
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
                std::cerr << typeid(T).name() << "::batchinv() wrongly thinks inverse exists" << std::endl;
                return false;
            }
        } else if (only_trivial_gcds) {
            std::cerr << typeid(T).name() << "::batchinv() wrongly thinks inverse does not exist" << std::endl;
            return false;
        } else {
            /* batchinv reported no inverse exists and there is, in fact, a
             * non-trivial gcd. All is well. */
            return true;
        }

        for (size_t i = 0; i < len; i++) {
            m.inv(t, a[i]);
            if (c != NULL) {
                m.mul(t, t, *c);
            }
            if (!m.equal(t, r[i])) {
                std::cerr << typeid(T).name() << "::batchinv() computed wrong inverse" << std::endl;
                return false;
            }
        }
        return true;
    }
    
    bool test_batchinv(const unsigned long iter) const {
        bool ok = true;
        
        for (unsigned long i_test = 0; i_test < iter; i_test++) {
            Modulus m = randomModulus();
            Residue a[10]{m, m, m, m, m, m, m, m, m, m}; /* Shoutout to C++ */
            Residue c(m);
            for (size_t i_residue = 0; i_residue < 10; i_residue++) {
                m.set(a[i_residue], randomInteger<Integer>());
            }
            m.set(c,  randomInteger<Integer>());
            for (size_t i_size = 0; i_size <= 10; i_size++) {
                ok &= test_one_batchinv(a, i_size, NULL, m);
                ok &= test_one_batchinv(a, i_size, &c, m);
            }
        }
        if (ok) {
            std::cout << "tests<" << typeid(T).name() << ">::test_batchinv() passed" << std::endl;
        }
        return ok;
    }
    
    bool test_one_jacobi (const Integer &a, const Modulus m) const {
        bool ok = true;
        Residue ar(m);
        m.set(ar, a);
        int jacobi1 = m.jacobi(ar);

        cxx_mpz A = (cxx_mpz) a, M = modToMpz(m);
        int jacobi2 = mpz_jacobi(A, M);

        if (jacobi1 != jacobi2) {
            std::cerr << typeid(T).name() << "::jacobi(" << A << ", " << M << ") = " << jacobi1 << " wrong" << std::endl;
            ok = false;
        }
        
        return ok;
    }
    
    bool test_jacobi(const unsigned long iter) const {
        bool ok = true;

        if (T::valid(Integer(1))) {
            Modulus m(Integer(1));
            ok &= test_one_jacobi(Integer(0), m);
        }
        
        Modulus m = randomModulus(true);
        ok &= test_one_jacobi(Integer(0), m);
        ok &= test_one_jacobi(Integer(1), m);
        Integer a;
        m.getmod(a);
        ASSERT_ALWAYS(a > 0);
        ok &= test_one_jacobi(a - 1, m);
        
        for (unsigned long i = 0; i < iter; i++) {
            Modulus m = randomModulus(true);
            Integer a = randomInteger<Integer>();
            ok &= test_one_jacobi(a, m);
        }
        if (ok) {
            std::cout << "tests<" << typeid(T).name() << ">::test_jacobi() passed" << std::endl;
        }
        return ok;
    }
    
    bool test_modop() const {
        Integer n;
        {
            Modulus m = randomModulus();
            ResidueOp r(m), s(m), t(m);
            n = randomInteger<Integer>();
        
            r = n;
            s = 2;
            t = r*s;
            if (t != r+r) {
                return false;
            }
            t += n;
            if (t != r+r+r) {
                return false;
            }
        }
        n = 727;
        if (Modulus::valid(n)) {
            Modulus m(n);
            ResidueOp r(m);

            r = 2;
            r = r.pow(n - 1);
            if ((Integer) r != 1)
                return false;
            r = 5;
            r = 5/r;
            if ((Integer) r != 1)
                return false;
        }
        std::cout << "tests<" << typeid(T).name() << ">::test_modop() passed" << std::endl;
        return true;
    }
    
    bool runTests(const unsigned long iter) const {
        bool ok = true;
        ok &= test_init();
        ok &= test_set(iter);
        ok &= test_neg();
        ok &= test_mul();
        ok &= test_divn(iter);
        ok &= test_pow(iter);
        ok &= test_sprp(iter);
        ok &= test_isprime(iter);
        ok &= test_gcd(iter);
        ok &= test_inv(iter);
        ok &= test_jacobi(iter);
        ok &= test_batchinv(iter);
        ok &= test_modop();
        return ok;
    }
};

template<>
Modulus64 Tests<Modulus64>::randomModulus(const bool odd) const {
    return Modulus64(u64_random(state) | (odd ? 1 : 0));
}

template<>
ModulusREDC64 Tests<ModulusREDC64>::randomModulus(const bool odd MAYBE_UNUSED) const {
    return ModulusREDC64(u64_random(state) | 1);
}

template<>
ModulusREDC126 Tests<ModulusREDC126>::randomModulus(const bool odd MAYBE_UNUSED) const {
    Integer128 m(u64_random(state) | 1, u64_random(state) & (UINT64_MAX >> 2));
    return ModulusREDC126(m);
}

template<>
ModulusMPZ Tests<ModulusMPZ>::randomModulus(const bool odd) const {
    ModulusMPZ::Integer i = randomInteger<ModulusMPZ::Integer>();
    if (odd && mpz_even_p(i))
        i += 1;
    ModulusMPZ m(i);
    return m;
}

int main(int argc, const char **argv) {
    unsigned long iter = 100;
    bool ok = true;
  
    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
    tests_common_get_iter(&iter);

    Tests<Modulus64> test1;
    ok &= test1.runTests(iter);
    
    Tests<ModulusREDC64> test2;
    ok &= test2.runTests(iter);
    
    Tests<ModulusREDC126> test3;
    ok &= test3.runTests(iter);
    
    Tests<ModulusMPZ> test4;
    ok &= test4.runTests(iter);
  
    tests_common_clear();
    exit(ok ? EXIT_SUCCESS : EXIT_FAILURE);
}
