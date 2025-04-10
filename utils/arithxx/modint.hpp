#ifndef CADO_MODINT_HPP
#define CADO_MODINT_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <array>

#include <gmp.h>

#include "macros.h"
#include "u64arith.h"
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "gmp_auxx.hpp"
#include "misc.h"

template<size_t NN>
class Integer_base : public std::array<uint64_t, NN>
{
public:
    typedef std::array<uint64_t, NN> super;
    using super::begin;
    using super::end;
    using super::data;
    // https://brevzin.github.io/c++/2020/02/05/constexpr-array-size/
    static constexpr size_t word_bits = 64;
    static constexpr size_t max_bits = NN * word_bits;
    static constexpr size_t max_size_in_words = NN;

    typedef uint64_t value_type;
    /* value-initialization sets all trailing values to zero */
    Integer_base() : super {} {}
    Integer_base(const uint64_t a) : super {a} {}
    /* Use default copy and move constructors */
    
    /* Set Integer64 from array of uint64_t of length n. Returns true if the
     * value fits in an Integer64 and false if not. If it does not fit, the
       contents of the Integer64 are undefined. */
    bool set(const uint64_t *s, const size_t n) {
        if (n > NN)
            return false;
        std::copy_n(s, n, begin());
        std::fill_n(begin() + n, NN - n, 0);
        return true;
    }
    /* Set from iterator over uint64_t. Returns true if the value
     * fits and false if not. If it does not fit, the resulting contents
       are undefined.*/
    bool set(const uint64_t *begin, const uint64_t *end) {
        return set(begin, end - begin);
    }
    template<size_t N>
    bool set(std::array<uint64_t, N> const & s) {
        static_assert(N <= NN, "N is too large");
        std::copy_n(s.begin(), N, begin());
        std::fill_n(begin() + N, NN - N, 0);
        return true;
    }

    /** Return the size in uint64_ts that is required in the output for
     * get(uint64_t *, size_t) */
    size_t size_in_words() const {
        size_t s;
        for(s = NN ; s && (*this)[s-1] == 0 ; s--);
        return s;
    }
    /** Write the Integer to r. Exactly len words are written.
     * If len is less than the required size as given by size(), output is
     * truncated. If len is greater, output is padded with zeroes. */
    void get (uint64_t *r, const size_t len) const {
        std::copy_n(begin(), std::min(len, NN), r);
    }
    template<size_t N>
    void get (std::array<uint64_t, N> & r) const {
        std::copy_n(begin(), std::min(N, NN), r.begin());
    }
    super const & get() const { return *this; }
    /* same as operator[], but with implicit zero-extend */
    value_type getWord(const size_t i) const {return (i < NN) ? (*this)[i] : 0;}

    void get(cxx_mpz &r) const {
        const size_t COUNT = NN;
        const int ORDER = -1;
        const size_t SIZE = sizeof(value_type);
        const int ENDIAN = 0;
        const size_t NAILS = 0;
        mpz_import(r, COUNT, ORDER, SIZE, ENDIAN, NAILS, data());
    }
    
    explicit operator cxx_mpz () const {
        cxx_mpz r;
        get(r);
        return r;
    }
    friend std::ostream & operator << (std::ostream &out, const Integer_base &i) {
        out << i[0];
        return out;
    }
};

/* Integers of 64 bits. We want additional conversion function from/to
   arrays of uint64_t and mpz_ts. Unfortunately, we can't inherit from
   standard integer types, so there's a lot of clutter here. */
class Integer64 : public Integer_base<1>
{
public:
    Integer64() = default;
    Integer64(const uint64_t a) : Integer_base(a) {}
    /* Use default copy and move constructors */
    
    /* Typecast operators */
    explicit operator bool() const {return (*this)[0] != 0;};
    explicit operator uint8_t() const {return (uint8_t) (*this)[0];}
    explicit operator uint16_t() const {return (uint16_t) (*this)[0];}
    explicit operator uint32_t() const {return (uint32_t) (*this)[0];}
    explicit operator uint64_t() const {return (*this)[0];}

    bool operator==(const Integer64 &a) const {return (*this)[0] == a[0];}
    bool operator==(const uint64_t a) const {return (*this)[0] == a;}
    bool operator!=(const Integer64 &a) const {return (*this)[0] != a[0];}
    bool operator!=(const uint64_t a) const {return (*this)[0] != a;}
    bool operator<(const Integer64 &a) const {return (*this)[0] < a[0];}
    bool operator<(const uint64_t a) const {return (*this)[0] < a;}
    bool operator>(const Integer64 &a) const {return (*this)[0] > a[0];}
    bool operator>(const uint64_t a) const {return (*this)[0] > a;}
    bool operator<=(const Integer64 &a) const {return (*this)[0] <= a[0];}
    bool operator<=(const uint64_t a) const {return (*this)[0] <= a;}
    bool operator>=(const Integer64 &a) const {return (*this)[0] >= a[0];}
    bool operator>=(const uint64_t a) const {return (*this)[0] >= a;}
    /* x.cmp(a) returns -1 for x<a, 0 for x==a, 1 for x>a.
       This is essentially the <=> operator from C++20. */
    int cmp(const Integer64 &a) const { return (*this > a) - (*this < a); }
    int cmp(const uint64_t a) const { return (*this > a) - (*this < a); }
    
    Integer64& operator=(const uint64_t a) {(*this)[0] = a; return *this;}
    Integer64& operator++ () {++(*this)[0]; return *this;}
    Integer64  operator++ (int) {Integer64 r = *this; (*this)[0]++; return r;}
    Integer64  operator+ (const Integer64 &a) const {return { (*this)[0] + a[0] };}
    Integer64  operator+ (const uint64_t a)   const {return { (*this)[0] + a }; }
    Integer64& operator+=(const Integer64 &a) {(*this)[0] += a[0]; return *this;}
    Integer64& operator+=(const uint64_t a)   {(*this)[0] += a; return *this;}
    Integer64& operator-- ()                  {--(*this)[0]; return *this;}
    Integer64  operator-- (int) {Integer64 r = *this; (*this)[0]--; return r;}
    Integer64  operator- (const Integer64 &a) const {return { (*this)[0] - a[0] }; }
    Integer64  operator- (const uint64_t a)   const {return { (*this)[0] - a }; }
    Integer64& operator-=(const Integer64 &a) {(*this)[0] -= a[0]; return *this;}
    Integer64& operator-=(const uint64_t a)   {(*this)[0] -= a; return *this;}
    Integer64  operator>>(const int i) const  {return { (*this)[0] >> i }; }
    Integer64& operator>>=(const int i)       {(*this)[0] >>= i; return *this;}
    Integer64  operator<<(const int i) const  {return { (*this)[0] << i }; }
    Integer64& operator<<=(const int i)       {(*this)[0] <<= i; return *this;}
    Integer64  operator* (const Integer64 &a) const {return { (*this)[0] * a[0] }; }
    Integer64  operator* (const uint64_t a)   const {return { (*this)[0] * a }; }
    Integer64& operator*=(const Integer64 &a) {(*this)[0] *= a[0]; return *this;}
    Integer64& operator*=(const uint64_t a)   {(*this)[0] *= a; return *this;}
    Integer64  operator/ (const Integer64 &a) const {return { (*this)[0] / a[0] }; }
    Integer64  operator/ (const uint64_t a)   const {return { (*this)[0] / a }; }
    Integer64& operator/=(const Integer64 &a) {(*this)[0] /= a[0]; return *this;}
    Integer64& operator/=(const uint64_t a)   {(*this)[0] /= a; return *this;}
    Integer64  operator% (const Integer64 &a) const {return { (*this)[0] % a[0] }; }
    Integer64  operator% (const uint64_t a)   const {return { (*this)[0] % a }; }
    Integer64& operator%=(const Integer64 &a) {(*this)[0] %= a[0]; return *this;}
    Integer64& operator%=(const uint64_t a)   {(*this)[0] %= a; return *this;}
    Integer64  operator| (const Integer64 &a) const {return { (*this)[0] | a[0] }; }
    Integer64  operator| (const uint64_t a)   const {return { (*this)[0] | a }; }
    Integer64& operator|=(const Integer64 &a) {(*this)[0] |= a[0]; return *this;}
    Integer64& operator|=(const uint64_t a)   {(*this)[0] |= a; return *this;}
    Integer64  operator& (const Integer64 &a) const {return { (*this)[0] & a[0] }; }
    Integer64  operator& (const uint64_t a)   const {return { (*this)[0] & a }; }
    Integer64& operator&=(const Integer64 &a) {(*this)[0] &= a[0]; return *this;}
    Integer64& operator&=(const uint64_t a)   {(*this)[0] &= a; return *this;}
    Integer64  operator^ (const Integer64 &a) const {return { (*this)[0] ^ a[0] }; }
    Integer64  operator^ (const uint64_t a)   const {return { (*this)[0] ^ a }; }
    Integer64& operator^=(const Integer64 &a) {(*this)[0] ^= a[0]; return *this;}
    Integer64& operator^=(const uint64_t a)   {(*this)[0] ^= a; return *this;}
    Integer64  operator~ () const {return { ~(*this)[0] }; }

    /* r = (*this)/a. We require a|(*this). */
    Integer64 divexact(const Integer64 &a) const {ASSERT_EXPENSIVE((*this)[0] % a[0] == 0); return { (*this)[0] / a[0] }; }
    
    Integer64& operator=(const cxx_mpz &s) {
        (*this)[0] = mpz_get_uint64(s);
        return *this;
    }
    
    /* Returns the number of bits in a, that is, floor(log_2(n))+1. For n==0 returns 0. */
    size_t bits() const {
        if ((*this)[0] == 0)
            return 0;
        return 64 - u64arith_clz ((*this)[0]);
    }
    int ctz() const {
        if ((*this)[0] == 0)
            return 64;
        return u64arith_ctz((*this)[0]);
    }
    friend std::ostream & operator << (std::ostream &out, const Integer64 &i) {
        out << i[0];
        return out;
    }
    
};


class Integer128 : public Integer_base<2>
{
public:
    Integer128() = default;
    explicit Integer128(const uint64_t a) : Integer_base(a) {}
    Integer128(const uint64_t a0, const uint64_t a1) {
        (*this)[0] = a0;
        (*this)[1] = a1;
    }
    explicit Integer128(Integer64 &a) : Integer_base((uint64_t) a) {}
    
    /* Typecast operators */
    explicit operator bool() const {return (*this)[0] != 0 || (*this)[1] != 0;}
    explicit operator uint8_t() const {return (uint8_t) (*this)[0];}
    explicit operator uint16_t() const {return (uint16_t) (*this)[0];}
    explicit operator uint32_t() const {return (uint32_t) (*this)[0];}
    explicit operator uint64_t() const {return (*this)[0];}

    bool operator==(const Integer128 &a) const {return ((*this)[0] == a[0] && (*this)[1] == a[1]);}
    bool operator==(const uint64_t a)   const {return ((*this)[0] == a && (*this)[1] == 0);}
    bool operator!=(const Integer128 &a) const {return !(*this == a);}
    bool operator!=(const uint64_t a)   const {return !(*this == a);}
    bool operator>=(const Integer128 &a) const {return u64arith_ge_2_2 ((*this)[0], (*this)[1], a[0], a[1]);}
    bool operator>=(const uint64_t a)   const {return u64arith_ge_2_2 ((*this)[0], (*this)[1], a, 0);}
    bool operator<=(const Integer128 &a) const {return u64arith_le_2_2 ((*this)[0], (*this)[1], a[0], a[1]);}
    bool operator<=(const uint64_t a)   const {return u64arith_le_2_2 ((*this)[0], (*this)[1], a, 0);}
    bool operator<(const Integer128 &a)  const {return u64arith_lt_2_2 ((*this)[0], (*this)[1], a[0], a[1]);}
    bool operator<(const uint64_t a)    const {return u64arith_lt_2_2 ((*this)[0], (*this)[1], a, 0);}
    bool operator>(const Integer128 &a)  const {return u64arith_gt_2_2 ((*this)[0], (*this)[1], a[0], a[1]);}
    bool operator>(const uint64_t a)    const {return u64arith_gt_2_2 ((*this)[0], (*this)[1], a, 0);}
    /* x.cmp(a) returns -1 for x<a, 0 for x==a, 1 for x>a */
    int cmp(const Integer128 &a) const {return (*this > a) - (*this < a);}
    int cmp(const uint64_t a) const {return (*this > a) - (*this < a); }
    
    Integer128& operator=(const uint64_t a) {(*this)[0] = a; (*this)[1] = 0; return *this;}
    Integer128& operator++ () {u64arith_add_2_2(data(), data() + 1, 1, 0); return *this;}
    Integer128  operator++ (int) {Integer128 r = *this; u64arith_add_2_2(data(), data() + 1, 1, 0); return r;}
    Integer128& operator+=(const Integer128 &a) {u64arith_add_2_2(data(), data() + 1, a[0], a[1]); return *this;}
    Integer128& operator+=(const uint64_t a)   {u64arith_add_2_2(data(), data() + 1, a, 0); return *this;}
    Integer128  operator+ (const Integer128 &a) const {Integer128 r = *this; r += a; return r;}
    Integer128  operator+ (const uint64_t a)   const {Integer128 r = *this; r += a; return r;}
    Integer128& operator-- () {u64arith_sub_2_2(data(), data() + 1, 1, 0); return *this;}
    Integer128  operator-- (int) {Integer128 r = *this; u64arith_sub_2_2(data(), data() + 1, 1, 0); return r;}
    Integer128& operator-=(const Integer128 &a) {u64arith_sub_2_2(data(), data() + 1, a[0], a[1]); return *this;}
    Integer128& operator-=(const uint64_t a)   {u64arith_sub_2_2(data(), data() + 1, a, 0); return *this;}
    Integer128  operator- (const Integer128 &a) const {Integer128 r = *this; r -= a; return r;}
    Integer128  operator- (const uint64_t a)   const {Integer128 r = *this; r -= a; return r;}
    Integer128& operator>>=(const int i) {
        ASSERT_EXPENSIVE(0 <= i && i < 128);
        if (i >= 64) {
            (*this)[0] = (*this)[1] >> (i - 64);
            (*this)[1] = 0;
        } else {
            u64arith_shr_2(data(), data() + 1, i);
        }
        return *this;
    }
    Integer128  operator>>(const int i) const {Integer128 r = *this; r >>= i; return r;}
    Integer128& operator<<=(const int i) {
        ASSERT_EXPENSIVE(0 <= i && i < 128);
        if (i >= 64) {
            (*this)[1] = (*this)[0] << (i - 64);
            (*this)[0] = 0;
        } else {
            u64arith_shl_2(data(), data() + 1, i);
        }
        return *this;
    }
    Integer128  operator<<(const int i) const {Integer128 r = *this; r <<= i; return r;}
    Integer128& operator*=(const Integer128 &a) {
        uint64_t r0, r1;
        u64arith_mul_1_1_2 (&r0, &r1, (*this)[0], a[0]);
        r1 += (*this)[0] * a[1] + (*this)[1] * a[0];
        (*this)[0] = r0; (*this)[1] = r1;
        return *this;
    }
    Integer128& operator*=(const uint64_t a) {
        const uint64_t t = (*this)[1];
        u64arith_mul_1_1_2 (data(), data() + 1, (*this)[0], a);
        (*this)[1] += t * a;
        return *this;
    }
    Integer128  operator*(const Integer128 &a) const {Integer128 r = *this; r *= a; return r;}
    Integer128  operator*(const uint64_t a) const {Integer128 r = *this; r *= a; return r;}
    Integer128& operator/=(const uint64_t a) {
        uint64_t r;
        u64arith_divqr_2_1_1 (&(*this)[1], &r, (*this)[1], 0, a);
        u64arith_divqr_2_1_1 (data(), &r, (*this)[0], r, a);
        return *this;
    }
    Integer128& operator/=(const Integer128 &a) {
        if (a[1] == 0) {
            return *this /= a[0];
        } else {
            uint64_t q, r0, r1;
            u64arith_divqr_3_2_1(&q, &r0, &r1, (*this)[0], (*this)[1], 0, a[0], a[1]);
            (*this)[0] = q;
            (*this)[1] = 0;
        }
        return *this;
    }
    Integer128  operator/(const Integer128 &a) const {Integer128 r = *this; r /= a; return r;}
    Integer128  operator/(const uint64_t a) const {Integer128 r = *this; r /= a; return r;}
    
    Integer128& operator%=(const uint64_t a) {
        uint64_t q;
        u64arith_divqr_2_1_1 (&q, &(*this)[1], (*this)[1], 0, a);
        u64arith_divqr_2_1_1 (&q, data(), (*this)[0], (*this)[1], a);
        (*this)[1] = 0;
        return *this;
    }
    Integer128& operator%=(const Integer128 &a) {
        if (a[1] == 0) {
            return *this %= a[0];
        } else {
            uint64_t q;
            u64arith_divqr_3_2_1(&q, data(), data() + 1, (*this)[0], (*this)[1], 0, a[0], a[1]);
        }
        return *this;
    }
    Integer128 operator%(const Integer128 &a) const {Integer128 r = *this; r %= a; return r;}
    uint64_t operator%(const uint64_t a) const {Integer128 r = *this; r %= a; return (uint64_t) r;}
    /* r = (*this)/a. We require a|(*this). */
    Integer128 divexact(const Integer128 &a) const {
        Integer128 n1, d1, r;
        uint64_t invf, r0, k0, k1;
        int i;
        
        n1 = *this;
        d1 = a;

        /* Make d odd */
        if (d1[0] == 0) {
            ASSERT (n1[0] == 0); /* Implied by d|n */
            d1[0] = d1[1];
            d1[1] = 0;
            n1[0] = n1[1];
            n1[1] = 0;
        }
        ASSERT(d1[0] != 0);
        i = u64arith_ctz (d1[0]);
        u64arith_shr_2 (d1.data(), d1.data() + 1, i);
        ASSERT((n1[0] & (((uint64_t)1 << i) - 1)) == 0); /* Implied by d|n */ 
        u64arith_shr_2 (n1.data(), n1.data() + 1, i);
        
        invf = u64arith_invmod (d1[0]);
        r0 = invf * n1[0];
        u64arith_mul_1_1_2 (&k0, &k1, r0, d1[0]);
        u64arith_sub_2_2 (n1.data(), n1.data() + 1, k0, k1);
        n1[1] -= r0 * d1[1];
        ASSERT (n1[0] == 0);
        r[0] = r0;
        r[1] = invf * n1[1];
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
    
    Integer128  operator| (const Integer128 &a) const {return { (*this)[0] | a[0], (*this)[1] | a[1] }; }
    Integer128  operator| (const uint64_t a)   const {return { (*this)[0] | a, (*this)[1] }; }
    Integer128& operator|=(const Integer128 &a) {(*this)[0] |= a[0]; (*this)[1] |= a[1]; return *this;}
    Integer128& operator|=(const uint64_t a)   {(*this)[0] |= a; return *this;}
    Integer128  operator& (const Integer128 &a) const {return { (*this)[0] & a[0], (*this)[1] & a[1] }; }
    Integer128  operator& (const uint64_t a)   const {return { (*this)[0] & a, 0 }; }
    Integer128& operator&=(const Integer128 &a) {(*this)[0] &= a[0]; (*this)[1] &= a[1]; return *this;}
    Integer128& operator&=(const uint64_t a)   {(*this)[0] &= a; (*this)[1] = 0; return *this;}
    Integer128  operator^ (const Integer128 &a) const {return { (*this)[0] ^ a[0], (*this)[1] ^ a[1] }; }
    Integer128  operator^ (const uint64_t a)   const {return { (*this)[0] ^ a, (*this)[1] }; }
    Integer128& operator^=(const Integer128 &a) {(*this)[0] ^= a[0]; (*this)[1] ^= a[1]; return *this;}
    Integer128& operator^=(const uint64_t a)   {(*this)[0] ^= a; return *this;}
    Integer128  operator~ () const {return { ~(*this)[0], ~(*this)[1] }; }

    Integer128& operator=(const cxx_mpz & s) {
        s.get(data(), max_size_in_words);
        return *this;
    }
    
    /* Returns the number of bits in a, that is, floor(log_2(n))+1. For n==0 returns 0. */
    size_t bits() const {
        if ((*this)[1] == 0) {
            if ((*this)[0] == 0) {
                return 0;
            }
            return 64 - u64arith_clz ((*this)[0]);
        }
        return 64 + 64 - u64arith_clz ((*this)[1]);
    }
    int ctz() const {
        if ((*this)[0] == 0) {
            if ((*this)[1] == 0)
                return 128;
            return 64 + u64arith_ctz((*this)[1]);
        }
        return u64arith_ctz((*this)[0]);
    }
    friend std::ostream & operator << (std::ostream &out, const Integer128 &s) {

        if (s == 0) {
            return out << "0";
        } else if (s[1] == 0) {
            return out << s[0];
        }

        IoStreamFlagsRestorer dummy(out);

        constexpr uint64_t ten19 = UINT64_C(10000000000000000000);
        Integer128 t{s};
        const uint64_t lower19 = t % ten19;
        t -= lower19;
        t /= ten19;
        const uint64_t upper19 = t % ten19;
        t -= upper19;

        if (t != 0) {
            t /= ten19;
            ASSERT_ALWAYS(t < 4);
            out << t[0] << std::setfill('0') << std::setw(19) << upper19 << std::setfill('0') << std::setw(19) << lower19;
        } else {
            // This is apparently a bug in coverity. The
            // IoStreamFlagsRestorer *does* properly restore the stream
            // cflags. It's just odd that coverity doesn't seem to see
            // it.
            // coverity[format_changed]
            out << upper19 << std::setfill('0') << std::setw(19) << lower19;
        }
        return out;
    }
    
};

/* Make these types known to gmp_auxx::mpz_fits(). This also makes these types
   known to cxx_mpz::fits(). */
namespace gmp_auxx {

template <>
inline bool mpz_fits<Integer64> (mpz_srcptr v) {
    return mpz_fits_uint64_p(v);
}

template <>
inline bool mpz_fits<Integer128> (mpz_srcptr v) {
    return mpz_sizeinbase(v, 2) <= 128;
}

}

#endif
