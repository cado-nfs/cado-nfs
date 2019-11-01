#ifndef MODINT_HPP
#define MODINT_HPP

#include <cassert>
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include "macros.h"
#include "u64arith.h"

#ifndef ASSERT
#define ASSERT(x)	assert(x)
#endif

#ifndef ASSERT_ALWAYS
#define ASSERT_ALWAYS(x)	assert(x)
#endif

/* Even simple assertions are relatively expensive in very simple functions.
   If we want them anyway to hunt a bug, define WANT_ASSERT_EXPENSIVE */
#ifdef WANT_ASSERT_EXPENSIVE
#ifndef ASSERT_EXPENSIVE
#define ASSERT_EXPENSIVE(x) ASSERT_ALWAYS(x)
#endif
#else
#define ASSERT_EXPENSIVE(x)
#endif

/* Integers of 64 bits. We want additional conversion function from/to
   arrays of uint64_t and mpz_ts. Unfortunately, we can't inherit from
   standard integer types, so there's a lot of clutter here. */
class Integer64 {
    uint64_t v[1];
public:
    Integer64() : v{0} {}
    Integer64(const uint64_t a) : v{a} {}
    /* Use default copy and move constructors */
    
    static size_t maxsize() {return 1;}
    /* Set Integer64 from array of uint64_t of length n. Returns true if the
     * value fits in an Integer64 and false if not. If it does not fit, the
       contents of the Integer64 are undefined. */
    bool set(const uint64_t *s, const size_t n) {
        if(n > maxsize()) {
            return false;
        }
        if (n == 0)
            v[0] = 0;
        else
            v[0] = s[0];
        return true;
    }
    /* Set Integer64 from iterator over uint64_t. Returns true if the value
     * fits in an Integer64 and false if not. If it does not fit, the contents
       of the Integer64 are undefined.*/
    bool set(const uint64_t *begin, const uint64_t *end) {
        size_t i = 0;
        const uint64_t *iter;
        for (iter = begin; iter != end; iter++) {
            if (i >= maxsize())
                break;
            v[i++] = *iter;
        }
        for ( ; i < maxsize(); i++) {
            v[i] = 0;
        }
        return (iter == end);
    }

    /* FIXME: this method is terrible and needs to go */
    const uint64_t *get() const {return v;}
    /* Return the size in uint64_ts that is required in the output for 
     * get(uint64_t *, size_t) */
    size_t size() const {return 1;}
    void get (uint64_t *r, const size_t len) const {
        if (len > 0)
            r[0] = v[0];
        for (size_t i = 1; i < len; i++)
            r[i] = 0;
    }

    /* Typecast operators */
    explicit operator bool() const {return v[0] != 0;};
    explicit operator uint8_t() const {return (uint8_t) v[0];}
    explicit operator uint16_t() const {return (uint16_t) v[0];}
    explicit operator uint32_t() const {return (uint32_t) v[0];}
    explicit operator uint64_t() const {return v[0];}

    bool operator==(const Integer64 &a) const {return v[0] == a.v[0];}
    bool operator==(const uint64_t a) const {return v[0] == a;}
    bool operator!=(const Integer64 &a) const {return v[0] != a.v[0];}
    bool operator!=(const uint64_t a) const {return v[0] != a;}
    bool operator<(const Integer64 &a) const {return v[0] < a.v[0];}
    bool operator<(const uint64_t a) const {return v[0] < a;}
    bool operator>(const Integer64 &a) const {return v[0] > a.v[0];}
    bool operator>(const uint64_t a) const {return v[0] > a;}
    bool operator<=(const Integer64 &a) const {return v[0] <= a.v[0];}
    bool operator<=(const uint64_t a) const {return v[0] <= a;}
    bool operator>=(const Integer64 &a) const {return v[0] >= a.v[0];}
    bool operator>=(const uint64_t a) const {return v[0] >= a;}
    /* x.cmp(a) returns -1 for x<a, 0 for x==a, 1 for x>a.
       This is essentially the <=> operator from C++20. */
    int cmp(const Integer64 &a) const {return (*this < a) ? -1 : (*this == a) ? 0 : 1;}
    int cmp(const uint64_t a) const {return (*this < a) ? -1 : (*this == a) ? 0 : 1;}
    
    Integer64& operator=(const uint64_t a) {v[0] = a; return *this;}
    Integer64& operator++ () {++v[0]; return *this;}
    Integer64  operator++ (int) {Integer64 r = *this; v[0]++; return r;}
    Integer64  operator+ (const Integer64 &a) const {return Integer64(v[0] + a.v[0]);}
    Integer64  operator+ (const uint64_t a)   const {return Integer64(v[0] + a);}
    Integer64& operator+=(const Integer64 &a) {v[0] += a.v[0]; return *this;}
    Integer64& operator+=(const uint64_t a)   {v[0] += a; return *this;}
    Integer64& operator-- () {--v[0]; return *this;}
    Integer64  operator-- (int) {Integer64 r = *this; v[0]--; return r;}
    Integer64  operator- (const Integer64 &a) const {return Integer64(v[0] - a.v[0]);}
    Integer64  operator- (const uint64_t a)   const {return Integer64(v[0] - a);}
    Integer64& operator-=(const Integer64 &a) {v[0] -= a.v[0]; return *this;}
    Integer64& operator-=(const uint64_t a)   {v[0] -= a; return *this;}
    Integer64  operator>>(const int i) const {return Integer64(v[0] >> i);}
    Integer64& operator>>=(const int i)      {v[0] >>= i; return *this;}
    Integer64  operator<<(const int i) const {return Integer64(v[0] << i);}
    Integer64& operator<<=(const int i)      {v[0] <<= i; return *this;}
    Integer64  operator* (const Integer64 &a) const {return Integer64(v[0] * a.v[0]);}
    Integer64  operator* (const uint64_t a)   const {return Integer64(v[0] * a);}
    Integer64& operator*=(const Integer64 &a) {v[0] *= a.v[0]; return *this;}
    Integer64& operator*=(const uint64_t a)   {v[0] *= a; return *this;}
    Integer64  operator/ (const Integer64 &a) const {return Integer64(v[0] / a.v[0]);}
    Integer64  operator/ (const uint64_t a)   const {return Integer64(v[0] / a);}
    Integer64& operator/=(const Integer64 &a) {v[0] /= a.v[0]; return *this;}
    Integer64& operator/=(const uint64_t a)   {v[0] /= a; return *this;}
    Integer64  operator% (const Integer64 &a) const {return Integer64(v[0] % a.v[0]);}
    Integer64  operator% (const uint64_t a)   const {return Integer64(v[0] % a);}
    Integer64& operator%=(const Integer64 &a) {v[0] %= a.v[0]; return *this;}
    Integer64& operator%=(const uint64_t a)   {v[0] %= a; return *this;}
    Integer64  operator| (const Integer64 &a) const {return Integer64(v[0] | a.v[0]);}
    Integer64  operator| (const uint64_t a)   const {return Integer64(v[0] | a);}
    Integer64& operator|=(const Integer64 &a) {v[0] |= a.v[0]; return *this;}
    Integer64& operator|=(const uint64_t a)   {v[0] |= a; return *this;}
    Integer64  operator& (const Integer64 &a) const {return Integer64(v[0] & a.v[0]);}
    Integer64  operator& (const uint64_t a)   const {return Integer64(v[0] & a);}
    Integer64& operator&=(const Integer64 &a) {v[0] &= a.v[0]; return *this;}
    Integer64& operator&=(const uint64_t a)   {v[0] &= a; return *this;}
    Integer64  operator^ (const Integer64 &a) const {return Integer64(v[0] ^ a.v[0]);}
    Integer64  operator^ (const uint64_t a)   const {return Integer64(v[0] ^ a);}
    Integer64& operator^=(const Integer64 &a) {v[0] ^= a.v[0]; return *this;}
    Integer64& operator^=(const uint64_t a)   {v[0] ^= a; return *this;}
    Integer64  operator~ () const {return Integer64(~v[0]);}

    /* r = v/a. We require a|v. */
    Integer64 divexact(const Integer64 &a) const {ASSERT_EXPENSIVE(v[0] % a.v[0] == 0); return Integer64(v[0] / a.v[0]);}
    
    Integer64& operator=(const mpz_class &s) {
        ASSERT_ALWAYS(mpz_sizeinbase(s.get_mpz_t(), 2) <= 64);
        const int ORDER = -1;
        const size_t SIZE = sizeof(v[0]);
        const int ENDIAN = 0;
        const size_t NAILS = 0;
        mpz_export(&v, NULL, ORDER, SIZE, ENDIAN, NAILS, s.get_mpz_t());
        return *this;
    }
    
    void get(mpz_class &r) const {
        const size_t COUNT = 1;
        const int ORDER = -1;
        const size_t SIZE = sizeof(v[0]);
        const int ENDIAN = 0;
        const size_t NAILS = 0;
        mpz_import(r.get_mpz_t(), COUNT, ORDER, SIZE, ENDIAN, NAILS, &v[0]);
    }
    
    explicit operator mpz_class () const {
        mpz_class r;
        get(r);
        return r;
    }
    /* Returns the number of bits in a, that is, floor(log_2(n))+1. For n==0 returns 0. */
    size_t bits() const {
        if (v[0] == 0)
            return 0;
        return 64 - u64arith_clz (v[0]);
    }
    int ctz() const {
        if (v[0] == 0)
            return 64;
        return u64arith_ctz(v[0]);
    }
    friend std::ostream & operator << (std::ostream &out, const Integer64 &i) {
        out << i.v[0];
        return out;
    }
    
};


class Integer128 {
    uint64_t v[2]; /* Least significant word at index 0 */
public:
    Integer128() : v{0,0} {}
    Integer128(const uint64_t a) : v{a,0} {}
    Integer128(const uint64_t a0, const uint64_t a1) : v{a0,a1} {}
    Integer128(Integer64 &a) : v{(uint64_t) a, 0} {}
    
    static size_t maxsize() {return 2;}

    /* The two set/get() functions import/export Integer128 from/to an array of 
        uint64_t. For set, the size of the array is passed as a parameter n which
        must not exceed 2.
        For get(), the required array size can be determined via size(); if the
        Integer128 is zero, get() writes 0 to the first output uint64_t.
    */
    bool set(const uint64_t *s, const size_t n) {
        if (n > maxsize()) {
            return false;
        }
        v[0] = v[1] = 0;
        if (n > 0)
            v[0] = s[0];
        if (n > 1)
            v[1] = s[1];
        return true;
    }

    /* Set Integer128 from iterator over uint64_t. Returns true if the value
     * fits in an Integer128 and false if not. If it does not fit, the
     * contents of the Integer64 are undefined. */
    bool set(const uint64_t *begin, const uint64_t *end) {
        size_t i = 0;
        const uint64_t *iter;
        for (iter = begin; iter != end; iter++) {
            if (i >= maxsize())
                break;
            v[i++] = *iter;
        }
        for ( ; i < maxsize(); i++) {
            v[i] = 0;
        }
        return (iter == end);
    }

    /* FIXME: this method is terrible and needs to go */
    const uint64_t *get() const {return v;}
    /* Return the size in uint64_ts that is required in the output for get() */
    size_t size() const {return (v[1] != 0) ? 2 : 1;}
    void get (uint64_t *r, const size_t len) const {
        if (len > 0)
            r[0] = v[0];
        if (len > 1)
            r[1] = v[1];
        for (size_t i = 2; i < len; i++)
            r[i] = 0;
    }

    /* Typecast operators */
    explicit operator bool() const {return v[0] != 0 || v[1] != 0;}
    explicit operator uint8_t() const {return (uint8_t) v[0];}
    explicit operator uint16_t() const {return (uint16_t) v[0];}
    explicit operator uint32_t() const {return (uint32_t) v[0];}
    explicit operator uint64_t() const {return v[0];}

    bool operator==(const Integer128 &a) const {return (v[0] == a.v[0] && v[1] == a.v[1]);}
    bool operator==(const uint64_t a)   const {return (v[0] == a && v[1] == 0);}
    bool operator!=(const Integer128 &a) const {return !(*this == a);}
    bool operator!=(const uint64_t a)   const {return !(*this == a);}
    bool operator>=(const Integer128 &a) const {return u64arith_ge_2_2 (v[0], v[1], a.v[0], a.v[1]);}
    bool operator>=(const uint64_t a)   const {return u64arith_ge_2_2 (v[0], v[1], a, 0);}
    bool operator<=(const Integer128 &a) const {return u64arith_le_2_2 (v[0], v[1], a.v[0], a.v[1]);}
    bool operator<=(const uint64_t a)   const {return u64arith_le_2_2 (v[0], v[1], a, 0);}
    bool operator<(const Integer128 &a)  const {return u64arith_lt_2_2 (v[0], v[1], a.v[0], a.v[1]);}
    bool operator<(const uint64_t a)    const {return u64arith_lt_2_2 (v[0], v[1], a, 0);}
    bool operator>(const Integer128 &a)  const {return u64arith_gt_2_2 (v[0], v[1], a.v[0], a.v[1]);}
    bool operator>(const uint64_t a)    const {return u64arith_gt_2_2 (v[0], v[1], a, 0);}
    /* x.cmp(a) returns -1 for x<a, 0 for x==a, 1 for x>a */
    int cmp(const Integer128 &a) const {return (*this < a) ? -1 : (*this == a) ? 0 : 1;}
    int cmp(const uint64_t a) const {return (*this < a) ? -1 : (*this == a) ? 0 : 1;}
    
    Integer128& operator=(const uint64_t a) {v[0] = a; v[1] = 0; return *this;}
    Integer128& operator++ () {u64arith_add_2_2(&v[0], &v[1], 1, 0); return *this;}
    Integer128  operator++ (int) {Integer128 r = *this; u64arith_add_2_2(&v[0], &v[1], 1, 0); return r;}
    Integer128& operator+=(const Integer128 &a) {u64arith_add_2_2(&v[0], &v[1], a.v[0], a.v[1]); return *this;}
    Integer128& operator+=(const uint64_t a)   {u64arith_add_2_2(&v[0], &v[1], a, 0); return *this;}
    Integer128  operator+ (const Integer128 &a) const {Integer128 r = *this; r += a; return r;}
    Integer128  operator+ (const uint64_t a)   const {Integer128 r = *this; r += a; return r;}
    Integer128& operator-- () {u64arith_sub_2_2(&v[0], &v[1], 1, 0); return *this;}
    Integer128  operator-- (int) {Integer128 r = *this; u64arith_sub_2_2(&v[0], &v[1], 1, 0); return r;}
    Integer128& operator-=(const Integer128 &a) {u64arith_sub_2_2(&v[0], &v[1], a.v[0], a.v[1]); return *this;}
    Integer128& operator-=(const uint64_t a)   {u64arith_sub_2_2(&v[0], &v[1], a, 0); return *this;}
    Integer128  operator- (const Integer128 &a) const {Integer128 r = *this; r -= a; return r;}
    Integer128  operator- (const uint64_t a)   const {Integer128 r = *this; r -= a; return r;}
    Integer128& operator>>=(const int i) {
        ASSERT_EXPENSIVE(0 <= i && i < 128);
        if (i >= 64) {
            v[0] = v[1] >> (i - 64);
            v[1] = 0;
        } else {
            u64arith_shr_2(&v[0], &v[1], i);
        }
        return *this;
    }
    Integer128  operator>>(const int i) const {Integer128 r = *this; r >>= i; return r;}
    Integer128& operator<<=(const int i) {
        ASSERT_EXPENSIVE(0 <= i && i < 128);
        if (i >= 64) {
            v[1] = v[0] << (i - 64);
            v[0] = 0;
        } else {
            u64arith_shl_2(&v[0], &v[1], i);
        }
        return *this;
    }
    Integer128  operator<<(const int i) const {Integer128 r = *this; r <<= i; return r;}
    Integer128& operator*=(const Integer128 &a) {
        uint64_t r0, r1;
        u64arith_mul_1_1_2 (&r0, &r1, v[0], a.v[0]);
        r1 += v[0] * a.v[1] + v[1] * a.v[0];
        v[0] = r0; v[1] = r1;
        return *this;
    }
    Integer128& operator*=(const uint64_t a) {
        const uint64_t t = v[1];
        u64arith_mul_1_1_2 (&v[0], &v[1], v[0], a);
        v[1] += t * a;
        return *this;
    }
    Integer128  operator*(const Integer128 &a) const {Integer128 r = *this; r *= a; return r;}
    Integer128  operator*(const uint64_t a) const {Integer128 r = *this; r *= a; return r;}
    Integer128& operator/=(const uint64_t a) {
        uint64_t q, r;
        u64arith_divqr_2_1_1 (&q, &v[1], v[1], 0, a);
        u64arith_divqr_2_1_1 (&v[0], &r, v[0], v[1], a);
        v[1] = q;
        return *this;
    }
    Integer128& operator/=(const Integer128 &a) {
        if (a.v[1] == 0) {
            return *this /= a.v[0];
        } else {
            uint64_t q, r0, r1;
            u64arith_divqr_3_2_1(&q, &r0, &r1, v[0], v[1], 0, a.v[0], a.v[1]);
            v[0] = q;
            v[1] = 0;
        }
        return *this;
    }
    Integer128  operator/(const Integer128 &a) const {Integer128 r = *this; r /= a; return r;}
    Integer128  operator/(const uint64_t a) const {Integer128 r = *this; r /= a; return r;}
    
    Integer128& operator%=(const uint64_t a) {
        uint64_t q;
        u64arith_divqr_2_1_1 (&q, &v[1], v[1], 0, a);
        u64arith_divqr_2_1_1 (&q, &v[0], v[0], v[1], a);
        v[1] = 0;
        return *this;
    }
    Integer128& operator%=(const Integer128 &a) {
        if (a.v[1] == 0) {
            return *this %= a.v[0];
        } else {
            uint64_t q;
            u64arith_divqr_3_2_1(&q, &v[0], &v[1], v[0], v[1], 0, a.v[0], a.v[1]);
        }
        return *this;
    }
    Integer128 operator%(const Integer128 &a) const {Integer128 r = *this; r %= a; return r;}
    Integer128 operator%(const uint64_t a) const {Integer128 r = *this; r %= a; return r;}
    /* r = v/a. We require a|v. */
    Integer128 divexact(const Integer128 &a) const {
        Integer128 n1, d1, r;
        uint64_t invf, r0, k0, k1;
        int i;
        
        n1 = *this;
        d1 = a;

        /* Make d odd */
        if (d1.v[0] == 0) {
            ASSERT (n1.v[0] == 0); /* Implied by d|n */
            d1.v[0] = d1.v[1];
            d1.v[1] = 0;
            n1.v[0] = n1.v[1];
            n1.v[1] = 0;
        }
        ASSERT(d1.v[0] != 0);
        i = u64arith_ctz (d1.v[0]);
        u64arith_shr_2 (&d1.v[0], &d1.v[1], i);
        ASSERT((n1.v[0] & (((uint64_t)1 << i) - 1)) == 0); /* Implied by d|n */ 
        u64arith_shr_2 (&n1.v[0], &n1.v[1], i);
        
        invf = u64arith_invmod (d1.v[0]);
        r0 = invf * n1.v[0];
        u64arith_mul_1_1_2 (&k0, &k1, r0, d1.v[0]);
        u64arith_sub_2_2 (&n1.v[0], &n1.v[1], k0, k1);
        n1.v[1] -= r0 * d1.v[1];
        ASSERT (n1.v[0] == 0);
        r.v[0] = r0;
        r.v[1] = invf * n1.v[1];
#ifdef WANT_ASSERT_EXPENSIVE
        ASSERT_EXPENSIVE((*this) * a == r);
#endif
        return r;
    }
    Integer128 divexact(const uint64_t a) const {
        /* We could simplify the right-adjust and save a single-word mul here,
         * but is it worth it? */
        return divexact(Integer128(a));
    }
    
    Integer128  operator| (const Integer128 &a) const {return Integer128(v[0] | a.v[0], v[1] | a.v[1]);}
    Integer128  operator| (const uint64_t a)   const {return Integer128(v[0] | a, v[1]);}
    Integer128& operator|=(const Integer128 &a) {v[0] |= a.v[0]; v[1] |= a.v[1]; return *this;}
    Integer128& operator|=(const uint64_t a)   {v[0] |= a; return *this;}
    Integer128  operator& (const Integer128 &a) const {return Integer128(v[0] & a.v[0], v[1] & a.v[1]);}
    Integer128  operator& (const uint64_t a)   const {return Integer128(v[0] & a, 0);}
    Integer128& operator&=(const Integer128 &a) {v[0] &= a.v[0]; v[1] &= a.v[1]; return *this;}
    Integer128& operator&=(const uint64_t a)   {v[0] &= a; v[1] = 0; return *this;}
    Integer128  operator^ (const Integer128 &a) const {return Integer128(v[0] ^ a.v[0], v[1] ^ a.v[1]);}
    Integer128  operator^ (const uint64_t a)   const {return Integer128(v[0] ^ a, v[1]);}
    Integer128& operator^=(const Integer128 &a) {v[0] ^= a.v[0]; v[1] ^= a.v[1]; return *this;}
    Integer128& operator^=(const uint64_t a)   {v[0] ^= a; return *this;}
    Integer128  operator~ () const {return Integer128(~v[0], ~v[1]);}

    Integer128 operator=(const mpz_class s) {
        ASSERT_ALWAYS(mpz_sizeinbase(s.get_mpz_t(), 2) <= 128);
        const int ORDER = -1;
        const size_t SIZE = sizeof(v[0]);
        const int ENDIAN = 0;
        const size_t NAILS = 0;
        mpz_export(&v, NULL, ORDER, SIZE, ENDIAN, NAILS, s.get_mpz_t());
        return *this;
    }
    void get(mpz_class &r) const {
        const size_t COUNT = 2;
        const int ORDER = -1;
        const size_t SIZE = sizeof(v[0]);
        const int ENDIAN = 0;
        const size_t NAILS = 0;
        mpz_import(r.get_mpz_t(), COUNT, ORDER, SIZE, ENDIAN, NAILS, v);
    }
    explicit operator mpz_class () const {
        mpz_class r;
        get(r);
        return r;
    }
    
    /* Returns the number of bits in a, that is, floor(log_2(n))+1. For n==0 returns 0. */
    size_t bits() const {
        if (v[1] == 0) {
            if (v[0] == 0) {
                return 0;
            }
            return 64 - u64arith_clz (v[0]);
        }
        return 64 + 64 - u64arith_clz (v[1]);
    }
    int ctz() const {
        if (v[0] == 0) {
            if (v[1] == 0)
                return 128;
            return 64 + u64arith_ctz(v[1]);
        }
        return u64arith_ctz(v[0]);
    }
    friend std::ostream & operator << (std::ostream &out, const Integer128 &i) {
        out << i.v[0] << ":" << i.v[1];
        return out;
    }
    
};
#endif
