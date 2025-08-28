#ifndef CADO_CXX_MPZ_HPP
#define CADO_CXX_MPZ_HPP

#include <cstdint>
#include <cstdlib>

#include <istream>
#include <ostream>
#include <type_traits>
#include <memory>

#include <gmp.h>
#include "fmt/ostream.h"
#include "fmt/base.h"

#include "is_non_narrowing_conversion.hpp"
#include "gmp_aux.h"
#include "gmp_auxx.hpp"
#include "macros.h"

struct cxx_mpz {
public:
    typedef mp_limb_t WordType;
    mpz_t x;
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
    cxx_mpz() { mpz_init(x); }

    template<typename T, typename U>
    static constexpr bool converts_via =
        integral_fits_v<T, U> &&
        cado_math_aux::is_non_narrowing_conversion_v<T, U>;

    template <typename T>
        // NOLINTNEXTLINE(hicpp-explicit-conversions)
        cxx_mpz (const T & rhs)
        requires converts_via<T, int64_t>
        {
            gmp_auxx::mpz_init_set(x, int64_t(rhs));
        }
    template <typename T>
        cxx_mpz & operator=(const T a)
        requires converts_via<T, int64_t>
        {
            gmp_auxx::mpz_set(x, int64_t(a));
            return *this;
        }
    template <typename T>
        // NOLINTNEXTLINE(hicpp-explicit-conversions)
        cxx_mpz (const T & rhs)
        requires converts_via<T, uint64_t>
        {
            gmp_auxx::mpz_init_set(x, uint64_t(rhs));
        }
    template <typename T>
        cxx_mpz & operator=(const T a)
        requires converts_via<T, uint64_t>
        {
            gmp_auxx::mpz_set(x, uint64_t(a));
            return *this;
        }

    ~cxx_mpz() { mpz_clear(x); }
    cxx_mpz(cxx_mpz const & o) {
        mpz_init_set(x, o.x);
    }
    // NOLINTNEXTLINE(hicpp-explicit-conversions)
    cxx_mpz(mpz_srcptr a) {
        mpz_init_set(x, a);
    }
    cxx_mpz & operator=(cxx_mpz const & o) {
        if (&o != this)
            mpz_set(x, o.x);
        return *this;
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init,hicpp-member-init)

    cxx_mpz(cxx_mpz && o) noexcept
        : cxx_mpz()
    {
        mpz_swap(x, o.x);
    }
    cxx_mpz& operator=(cxx_mpz && o) noexcept {
        if (&o != this)
            mpz_swap(x, o.x);
        return *this;
    }
    // NOLINTBEGIN(hicpp-explicit-conversions)
    operator mpz_ptr() { return x; }
    operator mpz_srcptr() const { return x; }
    /* it is very impotant to have the conversion to bool, otherwise the
     * implicit conversion to mpz_ptr wins!
     */
    explicit operator bool() { return mpz_size(x) != 0; }
    explicit operator bool() const { return mpz_size(x) != 0; }
    // NOLINTEND(hicpp-explicit-conversions)
    mpz_ptr operator->() { return x; }
    mpz_srcptr operator->() const { return x; }
    explicit operator uint64_t() const {return mpz_get_uint64(x);}
    
    /** Set the value of the cxx_mpz to that of the uint64_t array s,
     * least significant word first. 
     */
private:
    bool set(const uint64_t *s, const size_t len) {
        mpz_import(x, len, -1, sizeof(uint64_t), 0, 0, s);
        return true;
    }
public:
    cxx_mpz(const uint64_t *s, const size_t len)
        : cxx_mpz()
    {
        set(s, len);
    }

    /** Return the size in uint64_ts that is required in the output for
     * get(uint64_t *, size_t) */
    size_t size() const {return iceildiv(mpz_sizeinbase(x, 2), 64);}

    /** Write the absolute value of the cxx_mpz to r.
     * The least significant word is written first. Exactly len words are
     * written. If len is less than the required size as given by size(),
     * output is truncated. If len is greater, output is padded with zeroes.
     */
    void get(uint64_t *r, const size_t len) const {
        size_t written;
        if (len < size()) {
            auto t = std::unique_ptr<uint64_t[]>(static_cast<uint64_t*>(mpz_export(nullptr, &written, -1, sizeof(uint64_t), 0, 0, x)));
            /* Here, len < written. Write only the len least significant words
             * to r */
            for (size_t i = 0; i < len; i++)
                r[i] = t[i];
        } else {
            mpz_export(r, &written, -1, sizeof(uint64_t), 0, 0, x);
            ASSERT_ALWAYS(written <= len);
            for (size_t i = written; i < len; i++)
                r[i] = 0;
        }
    }

    /* Should use a C++ iterator instead? Would that be slower? */
    static constexpr size_t word_bits = GMP_NUMB_BITS;
    size_t size_in_words() const {return mpz_size(x);}
    const mp_limb_t * data() const { return mpz_limbs_read(x); }
    mp_bitcnt_t ctz() const { return mpz_scan1(x, 0); }
    WordType getWord(const size_t i) const {return mpz_getlimbn(x, i);}

    template <typename T>
    bool fits() const {
        return gmp_auxx::mpz_fits<T>(x);
    }
    size_t bits() const { return mpz_sizeinbase(x, 2); }

    /* Honestly I don't know. Having a few helper functions for
     * convenience is nice, but it's not clear at all where we should
     * stop, so there's a point in not starting at all. On the other
     * hand, we do have some operator overloads already, so adding helper
     * functions in the class doesn't look like an heresy to me.
     */
    cxx_mpz invmod(cxx_mpz const & p) const {
        cxx_mpz ri;
        mpz_invert (ri, x, p);
        return ri;
    }
};

template <>
inline bool cxx_mpz::fits<cxx_mpz>() const {
    return true;
}

#if GNUC_VERSION_ATLEAST(4,3,0)
extern void mpz_init(cxx_mpz & pl) __attribute__((error("mpz_init must not be called on a mpz reference -- it is the caller's business (via a ctor)")));
extern void mpz_clear(cxx_mpz & pl) __attribute__((error("mpz_clear must not be called on a mpz reference -- it is the caller's business (via a dtor)")));
#endif

struct cxx_mpq{
    mpq_t x;
    cxx_mpq() {mpq_init(x);}
    ~cxx_mpq() {mpq_clear(x);}
    cxx_mpq(unsigned long a, unsigned long b = 1) { mpq_init(x); mpq_set_ui(x, a,b); }
    cxx_mpq(cxx_mpq const & o) {
        mpq_init(x);
        mpq_set(x, o.x);
    }
    cxx_mpq & operator=(cxx_mpq const & o) {
        mpq_set(x, o.x);
        return *this;
    }
    cxx_mpq(cxx_mpq && o) {
        mpq_init(x);
        mpq_swap(x, o.x);
    }
    cxx_mpq& operator=(cxx_mpq && o) {
        mpq_swap(x, o.x);
        return *this;
    }
    template<typename T, typename U>
    static constexpr bool converts_via =
        integral_fits_v<T, U> &&
        cado_math_aux::is_non_narrowing_conversion_v<T, U>;

    template <typename T>
        // NOLINTNEXTLINE(hicpp-explicit-conversions)
        cxx_mpq (const T & rhs)
        requires converts_via<T, int64_t>
        {
            mpq_init(x);
            gmp_auxx::mpz_init_set(mpq_numref(x), int64_t(rhs));
            mpz_set_ui(mpq_denref(x), 1);
            mpq_canonicalize(x);
        }
    template <typename T>
        cxx_mpq & operator=(const T a)
        requires converts_via<T, int64_t>
        {
            gmp_auxx::mpz_set(mpq_numref(x), int64_t(a));
            mpz_set_ui(mpq_denref(x), 1);
            mpq_canonicalize(x);
            return *this;
        }
    template <typename T>
        // NOLINTNEXTLINE(hicpp-explicit-conversions)
        cxx_mpq (const T & rhs)
        requires converts_via<T, uint64_t>
        {
            mpq_init(x);
            gmp_auxx::mpz_set(mpq_numref(x), uint64_t(rhs));
            mpz_set_ui(mpq_denref(x), 1);
            mpq_canonicalize(x);
        }
    template <typename T>
        cxx_mpq & operator=(const T a)
        requires converts_via<T, uint64_t>
        {
            gmp_auxx::mpz_set(mpq_numref(x), uint64_t(a));
            mpz_set_ui(mpq_denref(x), 1);
            mpq_canonicalize(x);
            return *this;
        }
    // NOLINTNEXTLINE(hicpp-explicit-conversions)
    cxx_mpq(mpq_srcptr a) {
        mpq_init(x);
        mpq_set(x, a);
    }
    operator mpq_ptr() { return x; }
    operator mpq_srcptr() const { return x; }
    mpq_ptr operator->() { return x; }
    mpq_srcptr operator->() const { return x; }
};
#if GNUC_VERSION_ATLEAST(4,3,0)
extern void mpq_init(cxx_mpq & pl) __attribute__((error("mpq_init must not be called on a mpq reference -- it is the caller's business (via a ctor)")));
extern void mpq_clear(cxx_mpq & pl) __attribute__((error("mpq_clear must not be called on a mpq reference -- it is the caller's business (via a dtor)")));
#endif

#define CXX_MPZ_DEFINE_CMP(OP) \
inline bool operator OP(cxx_mpz const & a, cxx_mpz const & b) { return mpz_cmp(a, b) OP 0; } \
inline bool operator OP(mpz_srcptr a, cxx_mpz const & b) { return mpz_cmp(a, b) OP 0; } \
inline bool operator OP(cxx_mpz const & a, mpz_srcptr b) { return mpz_cmp(a, b) OP 0; } \
template <typename T> inline bool operator OP(cxx_mpz const & a, const T b) requires std::is_integral_v<T> { return gmp_auxx::mpz_cmp(a, b) OP 0; } \
template <typename T> inline bool operator OP(const T a, cxx_mpz const & b) requires std::is_integral_v<T> { return 0 OP gmp_auxx::mpz_cmp(b, a); }

CXX_MPZ_DEFINE_CMP(==)
CXX_MPZ_DEFINE_CMP(!=)
CXX_MPZ_DEFINE_CMP(<)
CXX_MPZ_DEFINE_CMP(>)
CXX_MPZ_DEFINE_CMP(<=)
CXX_MPZ_DEFINE_CMP(>=)

inline bool operator<=>(cxx_mpz const & a, cxx_mpz const & b) { return gmp_auxx::mpz_cmp(a, b); }

inline cxx_mpz operator+(cxx_mpz const & a, cxx_mpz const & b) { cxx_mpz r; mpz_add(r, a, b); return r; }
template <typename T>
inline cxx_mpz operator+(cxx_mpz const & a, const T b)
    requires std::is_integral_v<T>
{ cxx_mpz r; gmp_auxx::mpz_add(r, a, b); return r; }

template <typename T> inline cxx_mpz operator+(const T a, cxx_mpz const & b)
requires std::is_integral_v<T>
 { cxx_mpz r; gmp_auxx::mpz_add(r, b, a); return r; }

inline cxx_mpz & operator+=(cxx_mpz & a, cxx_mpz const & b) { mpz_add(a, a, b); return a; }
template <typename T> inline cxx_mpz & operator+=(cxx_mpz & a, const T b)
requires std::is_integral_v<T>
 { gmp_auxx::mpz_add(a, a, b); return a; }

inline cxx_mpz operator-(cxx_mpz const & a) { cxx_mpz r; mpz_neg(r, a); return r; }
inline cxx_mpz operator-(cxx_mpz const & a, cxx_mpz const & b) { cxx_mpz r; mpz_sub(r, a, b); return r; }
template <typename T> inline cxx_mpz operator-(cxx_mpz const & a, const T b)
requires std::is_integral_v<T>
 { cxx_mpz r; gmp_auxx::mpz_sub(r, a, b); return r; }
template <typename T> inline cxx_mpz operator-(const T a, cxx_mpz const & b)
requires std::is_integral_v<T>
 { cxx_mpz r; gmp_auxx::mpz_sub(r, a, b); return r; }

inline cxx_mpz & operator-=(cxx_mpz & a, cxx_mpz const & b) { mpz_sub(a, a, b); return a; }
template <typename T> inline cxx_mpz & operator-=(cxx_mpz & a, const T b) 
requires std::is_integral_v<T>
 { gmp_auxx::mpz_sub(a, a, b); return a; }

inline cxx_mpz operator*(cxx_mpz const & a, cxx_mpz const & b) { cxx_mpz r; mpz_mul(r, a, b); return r; }
template <typename T> inline cxx_mpz operator*(cxx_mpz const & a, const T b) 
requires std::is_integral_v<T>
 { cxx_mpz r; gmp_auxx::mpz_mul(r, a, b); return r; }
template <typename T> inline cxx_mpz operator*(const T a, cxx_mpz const & b) 
requires std::is_integral_v<T>
 { cxx_mpz r; gmp_auxx::mpz_mul(r, b, a); return r; }

inline cxx_mpz & operator*=(cxx_mpz & a, cxx_mpz const & b) { mpz_mul(a, a, b); return a; }
template <typename T> inline cxx_mpz & operator*=(cxx_mpz & a, const T b) 
requires std::is_integral_v<T>
 { gmp_auxx::mpz_mul(a, a, b); return a; }

inline cxx_mpz operator/(cxx_mpz const & a, cxx_mpz const & b) { cxx_mpz r; mpz_tdiv_q(r, a, b); return r; }
template <typename T> inline cxx_mpz operator/(cxx_mpz const & a, const T b) 
requires(std::is_integral_v<T> && std::is_unsigned_v<T>)
 { cxx_mpz r; mpz_tdiv_q_uint64(r, a, b); return r; }

inline cxx_mpz & operator/=(cxx_mpz & a, cxx_mpz const & b) { mpz_tdiv_q(a, a, b); return a; }
template <typename T> inline cxx_mpz & operator/=(cxx_mpz & a, const T b) 
requires(std::is_integral_v<T> && std::is_unsigned_v<T>)
 { mpz_tdiv_q_uint64(a, a, b); return a; }

inline cxx_mpz operator%(cxx_mpz const & a, cxx_mpz const & b)  { cxx_mpz r; mpz_tdiv_r(r, a, b); return r; }
template <typename T> inline cxx_mpz operator%(cxx_mpz const & a, const T b) 
requires(std::is_integral_v<T> && std::is_unsigned_v<T>)
 { cxx_mpz r; mpz_tdiv_r_uint64(r, a, b); return r; }

inline cxx_mpz & operator%=(cxx_mpz & a, cxx_mpz const & b)  { mpz_tdiv_r(a, a, b); return a; }
template <typename T> inline cxx_mpz & operator%=(cxx_mpz & a, const T b)  
requires(std::is_integral_v<T> && std::is_unsigned_v<T>)
 { mpz_tdiv_r_uint64(a, a, b); return a; }

inline cxx_mpz & operator<<=(cxx_mpz & a, const mp_bitcnt_t s)  { mpz_mul_2exp(a, a, s); return a; }
inline cxx_mpz operator<<(cxx_mpz const & a, const mp_bitcnt_t s)  { cxx_mpz r{a}; mpz_mul_2exp(r, r, s); return r; }

inline cxx_mpz & operator>>=(cxx_mpz & a, const mp_bitcnt_t s)  { mpz_tdiv_q_2exp(a, a, s); return a; }
inline cxx_mpz operator>>(cxx_mpz const & a, const mp_bitcnt_t s)  { cxx_mpz r{a}; mpz_tdiv_q_2exp(r, r, s); return r; }


inline cxx_mpz operator~(cxx_mpz const & a) { cxx_mpz r; mpz_com(r, a); return r; }
inline cxx_mpz operator|(cxx_mpz const & a, cxx_mpz const & b)  { cxx_mpz r; mpz_ior(r, a, b); return r; }
inline cxx_mpz operator|(cxx_mpz const & a, const unsigned long b)  { cxx_mpz r; mpz_ior(r, a, cxx_mpz(b)); return r; }
inline cxx_mpz & operator|=(cxx_mpz & a, cxx_mpz const & b)  { mpz_ior(a, a, b); return a; }
inline cxx_mpz & operator|=(cxx_mpz & a, const unsigned long b)  { mpz_ior(a, a, cxx_mpz(b)); return a; }

inline cxx_mpz operator^(cxx_mpz const & a, cxx_mpz const & b)  { cxx_mpz r; mpz_xor(r, a, b); return r; }
inline cxx_mpz operator^(cxx_mpz const & a, const unsigned long b)  { cxx_mpz r; mpz_xor(r, a, cxx_mpz(b)); return r; }
inline cxx_mpz & operator^=(cxx_mpz & a, cxx_mpz const & b)  { mpz_xor(a, a, b); return a; }
inline cxx_mpz & operator^=(cxx_mpz & a, const unsigned long b)  { mpz_xor(a, a, cxx_mpz(b)); return a; }

inline cxx_mpz operator&(cxx_mpz const & a, cxx_mpz const & b)  { cxx_mpz r; mpz_and(r, a, b); return r; }
inline cxx_mpz operator&(cxx_mpz const & a, const unsigned long b)  { cxx_mpz r; mpz_and(r, a, cxx_mpz(b)); return r; }
inline cxx_mpz & operator&=(cxx_mpz & a, cxx_mpz const & b)  { mpz_and(a, a, b); return a; }
inline cxx_mpz & operator&=(cxx_mpz & a, const unsigned long b)  { mpz_and(a, a, cxx_mpz(b)); return a; }

inline bool operator==(cxx_mpq const & a, cxx_mpq const & b) { return mpq_cmp(a, b) == 0; }
inline bool operator!=(cxx_mpq const & a, cxx_mpq const & b) { return mpq_cmp(a, b) != 0; }
inline bool operator<(cxx_mpq const & a, cxx_mpq const & b) { return mpq_cmp(a, b) < 0; }
inline bool operator>(cxx_mpq const & a, cxx_mpq const & b) { return mpq_cmp(a, b) > 0; }
inline std::ostream& operator<<(std::ostream& os, cxx_mpz const& x) { return os << (mpz_srcptr) x; }
inline std::ostream& operator<<(std::ostream& os, cxx_mpq const& x) { return os << (mpq_srcptr) x; }
inline std::istream& operator>>(std::istream& is, cxx_mpz & x) { return is >> (mpz_ptr) x; }
inline std::istream& operator>>(std::istream& is, cxx_mpq & x) { return is >> (mpq_ptr) x; }

namespace fmt {
    template <> struct formatter<cxx_mpz>: ostream_formatter {};
    template <> struct formatter<cxx_mpq>: ostream_formatter {};
}

/* a shorthand so that we can use user-defined literals */
static inline cxx_mpz operator"" _mpz(char const * str, size_t)
{
    cxx_mpz res;
    mpz_set_str(res, str, 0);
    return res;
}

#endif	/* CADO_CXX_MPZ_HPP */
