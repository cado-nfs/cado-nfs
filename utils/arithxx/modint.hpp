#ifndef CADO_MODINT_HPP
#define CADO_MODINT_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <climits>
#include <algorithm>
#include <array>
#include <stdexcept>
#include <compare>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/ostream.h"

#include "macros.h"
#include "u64arith.h"
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "gmp_auxx.hpp"
#include "utils_cxx.hpp"
#include "cado_compile_time_hacks.hpp"

/* see also mpn_compile_time.hpp for a similar purpose (but a different
 * interface)
 */

namespace Integer_details {

    template<typename Integer, int n, bool has_carry=true>
        struct reduce_multiple_impl;

    template<typename Integer, int n, size_t k = Integer::max_size_in_words>
        struct mod_n_impl;

template<typename T, size_t NN>
class Integer_base : public std::array<uint64_t, NN>
{
    static_assert(0 < NN && NN < (size_t) INT_MAX, "absurd size for integers");
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
    template<size_t N>
    explicit Integer_base(std::array<uint64_t, N> const & s)
    requires(N <= NN)
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

    std::strong_ordering operator<=>(const Integer_base & a) const
    {
        return std::lexicographical_compare_three_way(rbegin(), rend(), a.rbegin(), a.rend());
    }
    std::strong_ordering operator<=>(const uint64_t a) const {
        if (std::any_of(begin() + 1, end(), convert_bool()))
            return std::strong_ordering::greater;
        return (*this)[0] <=> a;
    }
    // int cmp(const Integer_base & a) const { return *this <=> a; }
    // int cmp(const uint64_t a) const { return *this <=> a; }

    bool operator==(Integer_base const & a) const { return std::ranges::equal(*this, a); }
    bool operator==(uint64_t const a) const { return (*this)[0] == a && std::none_of(begin() + 1, end(), convert_bool()); }

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
    T operator+ (const uint64_t a) const  { T r = downcast();  return r+=a;  }
    T operator- (const T &a) const        { T r = downcast();  return r-=a;  }
    T operator- (const uint64_t a) const  { T r = downcast();  return r-=a;  }
    /* we use another name so that we don't get overload conflicts on
     * operator- */
    T operator- () const  { return downcast().unary_minus(); }
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

    template<int n, bool has_carry=true, typename chooser>
        static T reduce_multiple(T & t, chooser const & c = {})
        {
            return Integer_details::reduce_multiple_impl<T, n, has_carry>::reduce(t, c);
        }

    template<int n>
        uint64_t mod_n() const
        {
            return Integer_details::mod_n_impl<T, n>::value(downcast());
        }

    uint64_t high_word() const { return super::back(); }
};
}

namespace fmt {
    template <typename T, size_t NN> struct formatter<Integer_details::Integer_base<T, NN>>: ostream_formatter {};
}

/* Integers of 64 bits. We want additional conversion function from/to
   arrays of uint64_t and mpz_ts. Unfortunately, we can't inherit from
   standard integer types, so there's a lot of clutter here. */
class Integer64 : public Integer_details::Integer_base<Integer64, 1>
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

    template<size_t N>
    explicit Integer64(std::array<uint64_t, N> const & s)
    requires(N <= super::max_size_in_words)
    : super(s) {}
    template<size_t N>
    Integer64& operator=(std::array<uint64_t, N> const & s)
        requires(N <= super::max_size_in_words)
    { set(s); return *this; }
    Integer64& operator=(const cxx_mpz & s) { s.get(data(), max_size_in_words); return *this; }
    Integer64& operator=(const uint64_t a) { set(&a, 1); return *this; }

    
    Integer64& operator++ () {++(*this)[0]; return *this;}
    Integer64& operator-- () {--(*this)[0]; return *this;}
    Integer64& operator+=(const Integer64 &a) {(*this)[0] += a[0]; return *this;}
    Integer64& operator-=(const Integer64 &a) {(*this)[0] -= a[0]; return *this;}
    Integer64& operator+=(const uint64_t a)   {(*this)[0] += a; return *this;}
    Integer64& operator-=(const uint64_t a)   {(*this)[0] -= a; return *this;}
    Integer64& operator>>=(const int i)       {(*this)[0] >>= i; return *this;}
    Integer64& signed_shift_right(const int i) {
        ASSERT_EXPENSIVE(0 <= i && i < 64);
        auto hi = int64_t((*this)[0]);
        (*this)[0] = uint64_t(hi >> i);
        return *this;
    }
    bool sign_bit() const { return ((int64_t)(*this)[0]) < 0; }

    Integer64& operator<<=(const int i)       {(*this)[0] <<= i; return *this;}
    Integer64& operator*=(const Integer64 &a) {(*this)[0] *= a[0]; return *this;}
    Integer64& operator/=(const Integer64 &a) {(*this)[0] /= a[0]; return *this;}
    Integer64& operator%=(const Integer64 &a) {(*this)[0] %= a[0]; return *this;}
    Integer64& operator*=(const uint64_t a)   {(*this)[0] *= a; return *this;}
    Integer64& operator/=(const uint64_t a)   {(*this)[0] /= a; return *this;}
    Integer64& operator%=(const uint64_t a)   {(*this)[0] %= a; return *this;}
    Integer64 unary_minus() const  { return Integer64 { -(*this)[0] }; }


    /* r = (*this)/a. We require a|(*this). */
    Integer64 divexact(const Integer64 &a) const {ASSERT_EXPENSIVE((*this)[0] % a[0] == 0); return Integer64((*this)[0] / a[0]); }
    
};

class Integer128 : public Integer_details::Integer_base<Integer128, 2>
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
    template<size_t N>
    explicit Integer128(std::array<uint64_t, N> const & s)
        requires(N <= super::max_size_in_words)
        : super(s) {}
    template<size_t N>
    Integer128& operator=(std::array<uint64_t, N> const & s)
        requires(N <= super::max_size_in_words)
    { set(s); return *this; }
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
    Integer128& signed_shift_right(const int i) {
        ASSERT_EXPENSIVE(0 <= i && i < 128);
        auto hi = int64_t((*this)[1]);
        if (i >= 64) {
            (*this)[0] = uint64_t(hi >> (i - 64));
            (*this)[1] = hi < 0 ? uint64_t(-1) : 0;
        } else {
            u64arith_shrd(data(), (*this)[1], (*this)[0], i);
            (*this)[1] = uint64_t(hi >> i);
        }
        return *this;
    }
    bool sign_bit() const { return ((int64_t)(*this)[1]) < 0; }

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

    Integer128 unary_minus() const  {
        auto a = *this;
        a[1] = -a[1];
        if (a[0] != 0)
            a[1]--;
        a[0] = -a[0];
        return a;
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

} /* namespace gmp_auxx */

namespace Integer_details {
    /* When we divide by small constants, there's a point where we create an
     * Integer type that is in the congruence class of the input, and that is
     * a multiple of the dividend. Or maybe just *congruent* to such a
     * multiple, in case the actual multiple involves a carry limb that was
     * discarded by wraparound.
     *
     * It's typically never a problem for Integer64, but it becomes
     * awkward for Integer128 and above.
     *
     * To cater for this situation, we define two distinct algorithms in
     * order to find the quotient. The first "simple" one is when there's no
     * carry. A truncating division + a multiplication are enough.
     *
     */
    template<int n>
        struct reduce_multiple_impl<Integer128, n, false> {
            template<typename ignore>
                static Integer128 reduce(Integer128 const & t, ignore const &) {
                    using cado_math_aux::invmod;
                    constexpr uint64_t c = invmod<n, uint64_t>::value;

                    /* a = a1 * 2^w + a0, n|a
                     * Let a = a' * n * 2^w + a'', a'' < n * 2^w.
                     * n | a'', a'' / n < 2^w
                     * So a / n = a' * w + a'' / n
                     * a' = trunc(a1 / n)
                     * a'' = a0 * n^{-1} (mod 2^w)
                     * Hence we get the correct result with one one-word
                     * multiplication and one one-word truncating
                     * division by a small constant.
                     *
                     * note that the truncating division by n is
                     * optimized by the compiler as a multiplication by c
                     * plus whatever is needed to get the quotient right.
                     */

                    return { t[0] * c, t[1] / n };
                }
        };

    /* If we do encounter a carry, then we're in for three
     * multiplications (one of them by n which is a small constant) */
    template<int n>
        struct reduce_multiple_impl<Integer128, n, true> {
            template<typename chooser_mul>
                static Integer128 reduce(Integer128 & t, chooser_mul const & cm) {
                    using cado_math_aux::invmod;
                    constexpr uint64_t c = invmod<n, uint64_t>::value;

                    Integer128 r, t2;

                    r[0] = t[0] * c;

                    /* r0 == (a+km)/n (mod w)
                       (r1*w + r0) * n = (a+km)
                       (r1*w + r0) * n == t (mod w^2)
                       r1*w*n == t - n*r0 (mod w^2)
                       t - n*r0 == 0 (mod w), thus
                       r1*n == (t - n*r0)/w (mod w) */

                    mul_c<n>(t2, r[0], cm);;

                    ASSERT_EXPENSIVE(t2[1] < n);
                    t -= t2;

                    ASSERT_EXPENSIVE(t[0] == 0);
                    r[1] = t[1] * c;

                    return r;
                }
        };

    /* implementations for Integer64 are easy */
    template<int n, bool b>
        struct reduce_multiple_impl<Integer64, n, b> {
            template<typename ignore>
                static Integer64 reduce(Integer64 & t, ignore const &) {
                    using cado_math_aux::invmod;
                    constexpr uint64_t c = invmod<n, uint64_t>::value;
                    return Integer64 { t[0] * c };
                }
        };

    /* This works only if n * NN does not overflow
     */
    template<typename Integer, int n, size_t k>
        struct mod_n_impl {
            static uint64_t value(Integer const & r) {
                using cado_math_aux::pow2_mod;
                constexpr uint64_t w_mod_n = pow2_mod<64, n>::value;
                return (mod_n_impl<Integer, n, k-1>::value(r) * w_mod_n + r[Integer::max_size_in_words-k] % n) % n;
            }
        };

    /* this specialization is probably not necessary, but I don't want
     * the compiler to do silly things. */
    template<typename Integer, int n>
        struct mod_n_impl<Integer, n, 1> {
            static uint64_t value(Integer const & r) {
                return r[Integer::max_size_in_words-1] % n;
            }
        };
    template<typename Integer, int n>
        struct mod_n_impl<Integer, n, 0> {
            static uint64_t value(Integer const &) {
                return 0;
            }
        };
} /* namespace Integer_details */

#endif
