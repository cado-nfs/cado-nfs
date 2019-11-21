#ifndef CXX_MPZ_HPP_
#define CXX_MPZ_HPP_
#include "macros.h"

#include <gmp.h>
#include <istream>
#include <ostream>
#include <limits>

#include "gmp_aux.h"

struct cxx_mpz {
    mpz_t x;
    cxx_mpz() { mpz_init(x); }
    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0 >
    cxx_mpz (const T & rhs) {
        static_assert(std::numeric_limits<T>::max() <= std::numeric_limits<uint64_t>::max(), 
                      "Integer type with more than 64 bits currently not supported");
        if (std::is_signed<T>::value) {
            mpz_init_set_int64(x, rhs);
        } else {
            mpz_init_set_uint64(x, (uint64_t) rhs);
        }
    }

    ~cxx_mpz() { mpz_clear(x); }
    cxx_mpz(cxx_mpz const & o) {
        mpz_init_set(x, o.x);
    }
    cxx_mpz & operator=(cxx_mpz const & o) {
        mpz_set(x, o.x);
        return *this;
    }
    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0 >
    cxx_mpz & operator=(const T a) {
        if (std::is_signed<T>::value)
            mpz_set_int64(x, a);
        else
            mpz_set_uint64(x, a);
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
};
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
inline bool operator OP(cxx_mpz const & a, cxx_mpz const & b) { return mpz_cmp(a, b) OP 0; }        \
inline bool operator OP(cxx_mpz const & a, const uint64_t b)  { return mpz_cmp_uint64(a, b) OP 0; } \
inline bool operator OP(const uint64_t a, cxx_mpz const & b)  { return 0 OP mpz_cmp_uint64(b, a); } \
inline bool operator OP(cxx_mpz const & a, const int64_t b)   { return mpz_cmp_int64(a, b) OP 0; }  \
inline bool operator OP(const int64_t a, cxx_mpz const & b)   { return 0 OP mpz_cmp_int64(b, a); }

CXX_MPZ_DEFINE_CMP(==)
CXX_MPZ_DEFINE_CMP(!=)
CXX_MPZ_DEFINE_CMP(<)
CXX_MPZ_DEFINE_CMP(>)
CXX_MPZ_DEFINE_CMP(<=)
CXX_MPZ_DEFINE_CMP(>=)

inline cxx_mpz operator+(cxx_mpz const & a, cxx_mpz const & b) { cxx_mpz r; mpz_add(r, a, b); return r; }
inline cxx_mpz operator+(cxx_mpz const & a, const uint64_t b)  { cxx_mpz r; mpz_add_uint64(r, a, b); return r; }
inline cxx_mpz operator+(const uint64_t a, cxx_mpz const & b)  { cxx_mpz r; mpz_add_uint64(r, b, a); return r; }
inline cxx_mpz operator+(cxx_mpz const & a, const int64_t b)   { cxx_mpz r; mpz_add_int64(r, a, b); return r; }
inline cxx_mpz operator+(const int64_t a, cxx_mpz const & b)   { cxx_mpz r; mpz_add_int64(r, b, a); return r; }

inline cxx_mpz & operator+=(cxx_mpz & a, cxx_mpz const & b) { mpz_add(a, a, b); return a; }
inline cxx_mpz & operator+=(cxx_mpz & a, const uint64_t b)  { mpz_add_uint64(a, a, b); return a; }
inline cxx_mpz & operator+=(cxx_mpz & a, const int64_t b)   { mpz_add_int64(a, a, b); return a; }

inline cxx_mpz operator-(cxx_mpz const & a, cxx_mpz const & b) { cxx_mpz r; mpz_sub(r, a, b); return r; }
inline cxx_mpz operator-(cxx_mpz const & a, const uint64_t b)  { cxx_mpz r; mpz_sub_uint64(r, a, b); return r; }
inline cxx_mpz operator-(const uint64_t a, cxx_mpz const & b)  { cxx_mpz r; mpz_uint64_sub(r, a, b); return r; }
inline cxx_mpz operator-(cxx_mpz const & a, const int64_t b)   { cxx_mpz r; mpz_sub_int64(r, a, b); return r; }
inline cxx_mpz operator-(const int64_t a, cxx_mpz const & b)   { cxx_mpz r; mpz_sub_int64(r, b, a); mpz_neg(r, r); return r; }

inline cxx_mpz & operator-=(cxx_mpz & a, cxx_mpz const & b) { mpz_sub(a, a, b); return a; }
inline cxx_mpz & operator-=(cxx_mpz & a, const uint64_t b)  { mpz_sub_uint64(a, a, b); return a; }
inline cxx_mpz & operator-=(cxx_mpz & a, const int64_t b)   { mpz_sub_int64(a, a, b); return a; }

inline cxx_mpz operator*(cxx_mpz const & a, cxx_mpz const & b) { cxx_mpz r; mpz_mul(r, a, b); return r; }
inline cxx_mpz operator*(cxx_mpz const & a, const uint64_t b)  { cxx_mpz r; mpz_mul_uint64(r, a, b); return r; }
inline cxx_mpz operator*(const uint64_t a, cxx_mpz const & b)  { cxx_mpz r; mpz_mul_uint64(r, b, a); return r; }
inline cxx_mpz operator*(cxx_mpz const & a, const int64_t b)   { cxx_mpz r; mpz_mul_int64(r, a, b); return r; }
inline cxx_mpz operator*(const int64_t a, cxx_mpz const & b)   { cxx_mpz r; mpz_mul_int64(r, b, a); return r; }

inline cxx_mpz & operator*=(cxx_mpz & a, cxx_mpz const & b) { mpz_mul(a, a, b); return a; }
inline cxx_mpz & operator*=(cxx_mpz & a, const uint64_t b)  { mpz_mul_uint64(a, a, b); return a; }
inline cxx_mpz & operator*=(cxx_mpz & a, const int64_t b)   { mpz_mul_int64(a, a, b); return a; }

inline cxx_mpz operator/(cxx_mpz const & a, cxx_mpz const & b) { cxx_mpz r; mpz_tdiv_q(r, a, b); return r; }
inline cxx_mpz operator/(cxx_mpz const & a, const uint64_t b)  { cxx_mpz r; mpz_tdiv_q_uint64(r, a, b); return r; }

inline cxx_mpz & operator/=(cxx_mpz & a, cxx_mpz const & b) { mpz_tdiv_q(a, a, b); return a; }
inline cxx_mpz & operator/=(cxx_mpz & a, const uint64_t b)  { mpz_tdiv_q_uint64(a, a, b); return a; }

inline cxx_mpz operator%(cxx_mpz const & a, cxx_mpz const & b)  { cxx_mpz r; mpz_tdiv_r(r, a, b); return r; }
inline cxx_mpz operator%(cxx_mpz const & a, const uint64_t b)  { cxx_mpz r; mpz_tdiv_r_uint64(r, a, b); return r; }

inline cxx_mpz & operator%=(cxx_mpz & a, cxx_mpz const & b)  { mpz_tdiv_r(a, a, b); return a; }
inline cxx_mpz & operator%=(cxx_mpz & a, const uint64_t b)   { mpz_tdiv_r_uint64(a, a, b); return a; }

inline cxx_mpz operator|(cxx_mpz const & a, cxx_mpz const & b)  { cxx_mpz r; mpz_ior(r, a, b); return r; }
inline cxx_mpz operator|(cxx_mpz const & a, const unsigned long b)  { cxx_mpz r; mpz_ior(r, a, cxx_mpz(b)); return r; }
inline cxx_mpz & operator|=(cxx_mpz & a, cxx_mpz const & b)  { mpz_ior(a, a, b); return a; }
inline cxx_mpz & operator|=(cxx_mpz & a, const unsigned long b)  { mpz_ior(a, a, cxx_mpz(b)); return a; }

inline bool operator==(cxx_mpq const & a, cxx_mpq const & b) { return mpq_cmp(a, b) == 0; }
inline bool operator!=(cxx_mpq const & a, cxx_mpq const & b) { return mpq_cmp(a, b) != 0; }
inline bool operator<(cxx_mpq const & a, cxx_mpq const & b) { return mpq_cmp(a, b) < 0; }
inline bool operator>(cxx_mpq const & a, cxx_mpq const & b) { return mpq_cmp(a, b) > 0; }
inline std::ostream& operator<<(std::ostream& os, cxx_mpz const& x) { return os << (mpz_srcptr) x; }
inline std::ostream& operator<<(std::ostream& os, cxx_mpq const& x) { return os << (mpq_srcptr) x; }
inline std::istream& operator>>(std::istream& is, cxx_mpz & x) { return is >> (mpz_ptr) x; }
inline std::istream& operator>>(std::istream& is, cxx_mpq & x) { return is >> (mpq_ptr) x; }
#endif	/* CXX_MPZ_HPP_ */
