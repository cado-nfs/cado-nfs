#ifndef CADO_MODINT_HPP
#define CADO_MODINT_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <iostream>
#include <algorithm>
#include <array>
#include <stdexcept>
#include <type_traits>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/ostream.h"

#include "macros.h"
#include "u64arith.h"
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "gmp_auxx.hpp"
#include "utils_cxx.hpp"

/* see also mpn_compile_time.hpp for a similar purpose (but a different
 * interface)
 */
template<typename T, size_t NN>
class Integer_base : public std::array<uint64_t, NN>
{
public:
    T& downcast() { return static_cast<T&>(*this); }
    T const & downcast() const { return static_cast<T const &>(*this); }

    friend T;

    typedef std::array<uint64_t, NN> super;
    using super::begin;
    using super::end;
    using super::rbegin;
    using super::rend;
    using super::data;
    // https://brevzin.github.io/c++/2020/02/05/constexpr-array-size/
    static constexpr size_t word_bits = 64;
    static constexpr size_t max_bits = NN * word_bits;
    static constexpr size_t max_size_in_words = NN;

    typedef uint64_t value_type;

    /* value-initialization sets all trailing values to zero */
private:
    Integer_base() : super {} {}
    explicit Integer_base(const uint64_t a) : super {a} {}
    template<size_t N, typename XX = typename std::enable_if<N <= NN>::type>
    explicit Integer_base(std::array<uint64_t, N> const & s)
#if GNUC_VERSION_ATMOST(8, 0, 0)
    /* It's probably a gcc bug. Definitely, the code below properly
     * initializes the std::array parent */
        : super()
#endif
    {
        std::copy_n(s.begin(), N, begin());
        std::fill_n(begin() + N, NN - N, 0);
    }
protected:
    /* Use default copy and move constructors */
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
    void get(cxx_mpz &r) const {
        const size_t COUNT = NN;
        const int ORDER = -1;
        const size_t SIZE = sizeof(value_type);
        const int ENDIAN = 0;
        const size_t NAILS = 0;
        mpz_import(r, COUNT, ORDER, SIZE, ENDIAN, NAILS, data());
    }
    
public:

    /** Return the size in uint64_ts that is required in the output for
     * get(uint64_t *, size_t) */
    size_t size_in_words() const {
        size_t s;
        for(s = NN ; s && (*this)[s-1] == 0 ; s--);
        return s;
    }
    bool fits_uint64_t() const { return size_in_words() <= 1; }

    super const & get() const { return *this; }
    /* same as operator[], but with implicit zero-extend */
    value_type getWord(const size_t i) const {return (i < NN) ? (*this)[i] : 0;}

    explicit operator cxx_mpz () const {
        cxx_mpz r;
        get(r);
        return r;
    }
    friend std::ostream & operator << (std::ostream &out, const Integer_base &i) {
        out << i[0];
        return out;
    }
    explicit operator double() const {
        return mpz_get_d(cxx_mpz(*this));
    }

    /* Typecast operators */
    explicit operator bool() const {return std::any_of(begin(), end(), convert_bool());}
    explicit operator uint8_t() const {return (uint8_t) (*this)[0];}
    explicit operator uint16_t() const {return (uint16_t) (*this)[0];}
    explicit operator uint32_t() const {return (uint32_t) (*this)[0];}
    explicit operator uint64_t() const {return (*this)[0];}

#if __cplusplus >= 202002L
    int operator<=>(const Integer_base & a) const { return std::lexicographical_compare_three_way(rbegin(), rend(), a.rbegin(), a.rend()); }
    int operator<=>(const uint64_t a) const {return (*this > a) - (*this < a); }
    int cmp(const Integer_base & a) const { return *this <=> a; }
    int cmp(const uint64_t a) const { return *this <=> a; }
#else
    int cmp(const Integer_base & a) const {return (*this > a) - (*this < a);}
    int cmp(const uint64_t a) const {return (*this > a) - (*this < a); }
#endif
    bool operator==(Integer_base const & a) const { return std::equal(begin(), end(), a.begin()); }
    bool operator< (Integer_base const & a) const { return std::lexicographical_compare(rbegin(), rend(), a.rbegin(), a.rend()); }
    bool operator!=(Integer_base const & a) const { return !operator==(a); }
    bool operator> (Integer_base const & a) const { return a < *this; }
    bool operator>=(Integer_base const & a) const { return !operator<(a); }
    bool operator<=(Integer_base const & a) const { return !(a < *this); }

    bool operator==(uint64_t const a) const { return (*this)[0] == a && std::none_of(begin() + 1, end(), convert_bool()); }
    bool operator< (uint64_t const a) const { return (*this)[0] < a && std::none_of(begin() + 1, end(), convert_bool()); }
    bool operator> (uint64_t const a) const { return std::any_of(begin() + 1, end(), convert_bool()) || (*this)[0] > a; }
    bool operator!=(uint64_t const a) const { return !operator==(a); }
    bool operator>=(uint64_t const a) const { return !operator<(a); }
    bool operator<=(uint64_t const a) const { return !operator<(a); }



    T& operator|=(const T &a) { auto c = a.begin(); for(auto & v: *this) v |= *c++; return downcast(); }
    T& operator&=(const T &a) { auto c = a.begin(); for(auto & v: *this) v &= *c++; return downcast(); }
    T& operator^=(const T &a) { auto c = a.begin(); for(auto & v: *this) v ^= *c++; return downcast(); }

    T& operator|=(uint64_t a) { (*this)[0] |= a; return downcast(); }
    T& operator&=(uint64_t a) { (*this)[0] &= a; std::fill_n(begin() + 1, NN - 1, 0); return downcast(); }
    T& operator^=(uint64_t a) { (*this)[0] ^= a; return downcast(); }

    T  operator| (const T &a) const { T r = downcast(); return r |= a; }
    T  operator& (const T &a) const { T r = downcast(); return r &= a; }
    T  operator^ (const T &a) const { T r = downcast(); return r ^= a; }

    T  operator| (uint64_t a) const { T r = downcast(); return r |= a; }
    T  operator& (uint64_t a) const { T r = downcast(); return r &= a; }
    T  operator^ (uint64_t a) const { T r = downcast(); return r ^= a; }

    T  operator~ () const { T res; for(auto& v: res) v=~v; return res; }

    size_t clz() const {
        size_t b = 0;
        for(auto c = super::rbegin() ; c != super::rend() ; ++c, b += word_bits)
            if (*c)
                return b + u64arith_clz(*c);
        return b;
    }
    /* Returns the number of bits in a, that is, floor(log_2(n))+1. For n==0 returns 0. */
    size_t bits() const {
        return max_bits - clz();
    }
    size_t ctz() const {
        size_t b = 0;
        for(auto c = super::begin() ; c != super::end() ; ++c, b += word_bits)
            if (*c)
                return b += u64arith_ctz(*c);
        return b;
    }

    /* Write all the element-creating operators as children of the
     * compound ones.  */
    T operator++ (int)                { T r = downcast(); ++r; return r;     }
    T operator-- (int)                { T r = downcast(); --r; return r;     }
    T operator+ (const T &a) const        { T r = downcast();  return r+=a;  }
    T operator- (const T &a) const        { T r = downcast();  return r-=a;  }
    T operator+ (const uint64_t a) const  { T r = downcast();  return r+=a;  }
    T operator- (const uint64_t a) const  { T r = downcast();  return r-=a;  }
    T operator>>(const int i) const       { T r = downcast();  return r>>=i; }
    T operator<<(const int i) const       { T r = downcast();  return r<<=i; }
    T operator* (const T &a) const        { T r = downcast();  return r*=a;  }
    T operator/ (const T &a) const        { T r = downcast();  return r/=a;  }
    T operator% (const T &a) const        { T r = downcast();  return r%=a;  }
    T operator* (const uint64_t a) const  { T r = downcast();  return r*=a;  }
    T operator/ (const uint64_t a) const  { T r = downcast();  return r/=a;  }
    uint64_t operator% (const uint64_t a) const  { T r = downcast();  return uint64_t(r%=a);  }

    friend std::ostream & operator << (std::ostream &out, const T &s) {
        return out << cxx_mpz(s);
    }
};

namespace fmt {
    template <typename T, size_t NN> struct formatter<Integer_base<T, NN>>: ostream_formatter {};
}

/* Integers of 64 bits. We want additional conversion function from/to
   arrays of uint64_t and mpz_ts. Unfortunately, we can't inherit from
   standard integer types, so there's a lot of clutter here. */
class Integer64 : public Integer_base<Integer64, 1>
{
public:
    typedef Integer_base<Integer64, 1> super;

    /* This blob must carry over unchanged for all instantiations */
    Integer64() = default;
    explicit Integer64(const uint64_t a) : super(a) {}
    explicit Integer64(cxx_mpz const & s) {
        if (s < 0 || mpz_sizeinbase(s, 2) > max_bits)
            throw std::invalid_argument("input does not fit");
        s.get(data(), max_size_in_words);
    }
    Integer64(const uint64_t * begin, const uint64_t * end) {
        if (!super::set(begin, end))
            throw std::invalid_argument("input does not fit");
    }
    Integer64(const uint64_t * begin, const size_t n) {
        if (!super::set(begin, n))
            throw std::invalid_argument("input does not fit");
    }

    template<size_t N, typename XX = typename std::enable_if<N <= super::max_size_in_words>::type>
    explicit Integer64(std::array<uint64_t, N> const & s) : super(s) {}
    template<size_t N, typename XX = typename std::enable_if<N <= super::max_size_in_words>::type>
    Integer64& operator=(std::array<uint64_t, N> const & s) { set(s); return *this; }
    Integer64& operator=(const cxx_mpz & s) { s.get(data(), max_size_in_words); return *this; }
    Integer64& operator=(const uint64_t a) { set(&a, 1); return *this; }

    
    Integer64& operator++ () {++(*this)[0]; return *this;}
    Integer64& operator-- () {--(*this)[0]; return *this;}
    Integer64& operator+=(const Integer64 &a) {(*this)[0] += a[0]; return *this;}
    Integer64& operator-=(const Integer64 &a) {(*this)[0] -= a[0]; return *this;}
    Integer64& operator+=(const uint64_t a)   {(*this)[0] += a; return *this;}
    Integer64& operator-=(const uint64_t a)   {(*this)[0] -= a; return *this;}
    Integer64& operator>>=(const int i)       {(*this)[0] >>= i; return *this;}
    Integer64& operator<<=(const int i)       {(*this)[0] <<= i; return *this;}
    Integer64& operator*=(const Integer64 &a) {(*this)[0] *= a[0]; return *this;}
    Integer64& operator/=(const Integer64 &a) {(*this)[0] /= a[0]; return *this;}
    Integer64& operator%=(const Integer64 &a) {(*this)[0] %= a[0]; return *this;}
    Integer64& operator*=(const uint64_t a)   {(*this)[0] *= a; return *this;}
    Integer64& operator/=(const uint64_t a)   {(*this)[0] /= a; return *this;}
    Integer64& operator%=(const uint64_t a)   {(*this)[0] %= a; return *this;}


    /* r = (*this)/a. We require a|(*this). */
    Integer64 divexact(const Integer64 &a) const {ASSERT_EXPENSIVE((*this)[0] % a[0] == 0); return Integer64((*this)[0] / a[0]); }
    
};

class Integer128 : public Integer_base<Integer128, 2>
{
public:
    typedef Integer_base<Integer128, 2> super;

    /* This blob must carry over unchanged for all instantiations */
    Integer128() = default;
    explicit Integer128(const uint64_t a) : super(a) {}
    explicit Integer128(cxx_mpz const & s) {
        if (s < 0 || mpz_sizeinbase(s, 2) > max_bits)
            throw std::invalid_argument("input does not fit");
        s.get(data(), max_size_in_words);
    }
    Integer128(const uint64_t * begin, const uint64_t * end) {
        if (!super::set(begin, end))
            throw std::invalid_argument("input does not fit");
    }
    Integer128(const uint64_t * begin, const size_t n) {
        if (!super::set(begin, n))
            throw std::invalid_argument("input does not fit");
    }
    template<size_t N, typename XX = typename std::enable_if<N <= super::max_size_in_words>::type>
    explicit Integer128(std::array<uint64_t, N> const & s) : super(s) {}
    template<size_t N, typename XX = typename std::enable_if<N <= super::max_size_in_words>::type>
    Integer128& operator=(std::array<uint64_t, N> const & s) { set(s); return *this; }
    Integer128& operator=(const cxx_mpz & s) { s.get(data(), max_size_in_words); return *this; }
    Integer128& operator=(const uint64_t a) { set(&a, 1); return *this; }

    /* for convenience only */
    Integer128(const uint64_t a0, const uint64_t a1) : super { super::super { a0, a1 } } {}
    
    Integer128& operator++ () {u64arith_add_2_2(data(), data() + 1, 1, 0); return *this;}
    Integer128& operator+=(const Integer128 &a) {u64arith_add_2_2(data(), data() + 1, a[0], a[1]); return *this;}
    Integer128& operator+=(const uint64_t a)   {u64arith_add_2_2(data(), data() + 1, a, 0); return *this;}
    Integer128& operator-- () {u64arith_sub_2_2(data(), data() + 1, 1, 0); return *this;}
    Integer128& operator-=(const Integer128 &a) {u64arith_sub_2_2(data(), data() + 1, a[0], a[1]); return *this;}
    Integer128& operator-=(const uint64_t a)   {u64arith_sub_2_2(data(), data() + 1, a, 0); return *this;}
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
    
};

namespace fmt {
    template <> struct formatter<Integer64>: ostream_formatter {};
    template <> struct formatter<Integer128>: ostream_formatter {};
}

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
