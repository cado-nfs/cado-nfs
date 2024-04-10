#ifndef CXX_MPZ_HPP_
#define CXX_MPZ_HPP_

#include "macros.h"

#include <gmp.h>
#include <istream>
#include <ostream>
#include <limits>
#include <type_traits>
#include <stdlib.h>
#include <fmt/core.h>
#include <sstream>
#include "gmp_aux.h"
#include "gmp_auxx.hpp"

struct cxx_mpz {
public:
    typedef mp_limb_t WordType;
    mpz_t x;
    cxx_mpz() { mpz_init(x); }
    template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
    cxx_mpz (const T & rhs) {
        gmp_auxx::mpz_init_set(x, rhs);
    }

    ~cxx_mpz() { mpz_clear(x); }
    cxx_mpz(cxx_mpz const & o) {
        mpz_init_set(x, o.x);
    }
    cxx_mpz(mpz_srcptr a) {
        mpz_init_set(x, a);
    }
    cxx_mpz & operator=(cxx_mpz const & o) {
        mpz_set(x, o.x);
        return *this;
    }
    template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
    cxx_mpz & operator=(const T a) {
        gmp_auxx::mpz_set(x, a);
        return *this;
    }

#if __cplusplus >= 201103L
    cxx_mpz(cxx_mpz && o) {
        mpz_init(x);
        mpz_swap(x, o.x);
    }
    cxx_mpz& operator=(cxx_mpz && o) {
        mpz_swap(x, o.x);
        return *this;
    }
#endif
    operator mpz_ptr() { return x; }
    operator mpz_srcptr() const { return x; }
    mpz_ptr operator->() { return x; }
    mpz_srcptr operator->() const { return x; }
    explicit operator uint64_t() const {return mpz_get_uint64(x);}
    
    /** Set the value of the cxx_mpz to that of the uint64_t array s,
     * least significant word first. 
     */
    bool set(const uint64_t *s, const size_t len) {
        mpz_import(x, len, -1, sizeof(uint64_t), 0, 0, s);
        return true;
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
        const bool useTemp = len < size();
        size_t written;
        uint64_t *t = (uint64_t *) mpz_export(useTemp ? NULL : r, &written,
                                              -1, sizeof(uint64_t), 0, 0, x);
        if (useTemp) {
            /* Here, len < written. Write only the len least significant words
             * to r */
            for (size_t i = 0; i < len; i++)
                r[i] = t[i];
            free(t);
        } else {
            ASSERT_ALWAYS(written <= len);
            for (size_t i = written; i < len; i++)
                r[i] = 0;
        }
    }

    /* Should use a C++ iterator instead? Would that be slower? */
    static int getWordSize() {return GMP_NUMB_BITS;}
    size_t getWordCount() const {return mpz_size(x);}
    WordType getWord(const size_t i) const {return mpz_getlimbn(x, i);}

    template <typename T>
    bool fits() const {
        return gmp_auxx::mpz_fits<T>(x);
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
#if __cplusplus >= 201103L
    cxx_mpq(cxx_mpq && o) {
        mpq_init(x);
        mpq_swap(x, o.x);
    }
    cxx_mpq& operator=(cxx_mpq && o) {
        mpq_swap(x, o.x);
        return *this;
    }
#endif
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
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >   \
inline bool operator OP(cxx_mpz const & a, const T b) { return gmp_auxx::mpz_cmp(a, b) OP 0; } \
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >   \
inline bool operator OP(const T a, cxx_mpz const & b) { return 0 OP gmp_auxx::mpz_cmp(b, a); }

CXX_MPZ_DEFINE_CMP(==)
CXX_MPZ_DEFINE_CMP(!=)
CXX_MPZ_DEFINE_CMP(<)
CXX_MPZ_DEFINE_CMP(>)
CXX_MPZ_DEFINE_CMP(<=)
CXX_MPZ_DEFINE_CMP(>=)

inline cxx_mpz operator+(cxx_mpz const & a, cxx_mpz const & b) { cxx_mpz r; mpz_add(r, a, b); return r; }
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
inline cxx_mpz operator+(cxx_mpz const & a, const T b) { cxx_mpz r; gmp_auxx::mpz_add(r, a, b); return r; }
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
inline cxx_mpz operator+(const T a, cxx_mpz const & b) { cxx_mpz r; gmp_auxx::mpz_add(r, b, a); return r; }

inline cxx_mpz & operator+=(cxx_mpz & a, cxx_mpz const & b) { mpz_add(a, a, b); return a; }
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
inline cxx_mpz & operator+=(cxx_mpz & a, const T b) { gmp_auxx::mpz_add(a, a, b); return a; }

inline cxx_mpz operator-(cxx_mpz const & a, cxx_mpz const & b) { cxx_mpz r; mpz_sub(r, a, b); return r; }
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
inline cxx_mpz operator-(cxx_mpz const & a, const T b) { cxx_mpz r; gmp_auxx::mpz_sub(r, a, b); return r; }
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
inline cxx_mpz operator-(const T a, cxx_mpz const & b) { cxx_mpz r; gmp_auxx::mpz_sub(r, a, b); return r; }

inline cxx_mpz & operator-=(cxx_mpz & a, cxx_mpz const & b) { mpz_sub(a, a, b); return a; }
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
inline cxx_mpz & operator-=(cxx_mpz & a, const T b)  { gmp_auxx::mpz_sub(a, a, b); return a; }

inline cxx_mpz operator*(cxx_mpz const & a, cxx_mpz const & b) { cxx_mpz r; mpz_mul(r, a, b); return r; }
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
inline cxx_mpz operator*(cxx_mpz const & a, const T b)  { cxx_mpz r; gmp_auxx::mpz_mul(r, a, b); return r; }
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
inline cxx_mpz operator*(const T a, cxx_mpz const & b)  { cxx_mpz r; gmp_auxx::mpz_mul(r, b, a); return r; }

inline cxx_mpz & operator*=(cxx_mpz & a, cxx_mpz const & b) { mpz_mul(a, a, b); return a; }
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
inline cxx_mpz & operator*=(cxx_mpz & a, const T b)  { gmp_auxx::mpz_mul(a, a, b); return a; }

inline cxx_mpz operator/(cxx_mpz const & a, cxx_mpz const & b) { cxx_mpz r; mpz_tdiv_q(r, a, b); return r; }
template <typename T, std::enable_if_t<std::is_integral<T>::value && std::is_unsigned<T>::value, int> = 0 >
inline cxx_mpz operator/(cxx_mpz const & a, const T b)  { cxx_mpz r; mpz_tdiv_q_uint64(r, a, b); return r; }

inline cxx_mpz & operator/=(cxx_mpz & a, cxx_mpz const & b) { mpz_tdiv_q(a, a, b); return a; }
template <typename T, std::enable_if_t<std::is_integral<T>::value && std::is_unsigned<T>::value, int> = 0 >
inline cxx_mpz & operator/=(cxx_mpz & a, const T b)  { mpz_tdiv_q_uint64(a, a, b); return a; }

inline cxx_mpz operator%(cxx_mpz const & a, cxx_mpz const & b)  { cxx_mpz r; mpz_tdiv_r(r, a, b); return r; }
template <typename T, std::enable_if_t<std::is_integral<T>::value && std::is_unsigned<T>::value, int> = 0 >
inline cxx_mpz operator%(cxx_mpz const & a, const T b)  { cxx_mpz r; mpz_tdiv_r_uint64(r, a, b); return r; }

inline cxx_mpz & operator%=(cxx_mpz & a, cxx_mpz const & b)  { mpz_tdiv_r(a, a, b); return a; }
template <typename T, std::enable_if_t<std::is_integral<T>::value && std::is_unsigned<T>::value, int> = 0 >
inline cxx_mpz & operator%=(cxx_mpz & a, const T b)   { mpz_tdiv_r_uint64(a, a, b); return a; }

inline cxx_mpz & operator<<=(cxx_mpz & a, const mp_bitcnt_t s)  { mpz_mul_2exp(a, a, s); return a; }
inline cxx_mpz operator<<(cxx_mpz & a, const mp_bitcnt_t s)  { cxx_mpz r{a}; mpz_mul_2exp(r, r, s); return r; }

inline cxx_mpz & operator>>=(cxx_mpz & a, const mp_bitcnt_t s)  { mpz_tdiv_q_2exp(a, a, s); return a; }
inline cxx_mpz operator>>(cxx_mpz & a, const mp_bitcnt_t s)  { cxx_mpz r{a}; mpz_tdiv_q_2exp(r, r, s); return r; }


#if 0
inline cxx_mpz operator|(cxx_mpz const & a, cxx_mpz const & b)  { cxx_mpz r; mpz_ior(r, a, b); return r; }
inline cxx_mpz operator|(cxx_mpz const & a, const unsigned long b)  { cxx_mpz r; mpz_ior(r, a, cxx_mpz(b)); return r; }
inline cxx_mpz & operator|=(cxx_mpz & a, cxx_mpz const & b)  { mpz_ior(a, a, b); return a; }
inline cxx_mpz & operator|=(cxx_mpz & a, const unsigned long b)  { mpz_ior(a, a, cxx_mpz(b)); return a; }
#endif

inline bool operator==(cxx_mpq const & a, cxx_mpq const & b) { return mpq_cmp(a, b) == 0; }
inline bool operator!=(cxx_mpq const & a, cxx_mpq const & b) { return mpq_cmp(a, b) != 0; }
inline bool operator<(cxx_mpq const & a, cxx_mpq const & b) { return mpq_cmp(a, b) < 0; }
inline bool operator>(cxx_mpq const & a, cxx_mpq const & b) { return mpq_cmp(a, b) > 0; }
inline std::ostream& operator<<(std::ostream& os, cxx_mpz const& x) { return os << (mpz_srcptr) x; }
inline std::ostream& operator<<(std::ostream& os, cxx_mpq const& x) { return os << (mpq_srcptr) x; }
inline std::istream& operator>>(std::istream& is, cxx_mpz & x) { return is >> (mpz_ptr) x; }
inline std::istream& operator>>(std::istream& is, cxx_mpq & x) { return is >> (mpq_ptr) x; }

namespace fmt {
    template <> struct /* fmt:: */ formatter<cxx_mpz>: formatter<string_view> {
    // only allow {} for formatting. No :, no :x, etc. It could be nice
    // to allow them, though. Note that this should be constexpr with
    // c++-14 or later
    // auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) { return ctx.begin(); }
    template <typename FormatContext>
        auto format(cxx_mpz const & c, FormatContext& ctx) -> decltype(ctx.out())
        {
            std::ostringstream os;
            os << c;
            return formatter<string_view>::format( string_view(os.str()), ctx);
        }
};
}

#endif	/* CXX_MPZ_HPP_ */
