#ifndef GMP_AUXX_HPP_
#define GMP_AUXX_HPP_

#include "cxx_misc.hpp"
#include "gmp_aux.h"

/* A C++ wrapper around functions in gmp_aux.h. The functions here delegate
   to the proper C function depending on the argument type. We keep them
   out of global name space to avoid obscuring errors */

namespace gmp_auxx {
static inline void mpz_set (mpz_ptr a, mpz_srcptr b) {
    ::mpz_set(a, b);
}

template <typename T, integral_fits_t<T, long> = 0>
static inline void mpz_set (mpz_ptr a, const T b) {
    mpz_set_si(a, b);
}

template <typename T, integral_fits_t<T, unsigned long> = 0 >
static inline void mpz_set (mpz_ptr a, const T b) {
    mpz_set_ui(a, b);
}

inline void mpz_set (mpz_ptr a, const int64_t b) {
    return mpz_set_int64(a, b);
}

inline void mpz_set (mpz_ptr a, const uint64_t b) {
    mpz_set_uint64(a, b);
}

static inline void mpz_init_set (mpz_ptr a, mpz_srcptr b) {
    ::mpz_init_set(a, b);
}

template <typename T, integral_fits_t<T, long> = 0 >
static inline void mpz_init_set (mpz_ptr a, const T b) {
    mpz_init_set_si(a, b);
}

template <typename T, integral_fits_t<T, unsigned long> = 0 >
static inline void mpz_init_set (mpz_ptr a, const T b) {
    mpz_init_set_ui(a, b);
}

inline void mpz_init_set (mpz_ptr a, const int64_t b) {
    mpz_init_set_int64(a, b);
}

inline void mpz_init_set (mpz_ptr a, const uint64_t b) {
    mpz_init_set_uint64(a, b);
}

static inline int mpz_cmp (mpz_srcptr a, mpz_srcptr b) {
    return ::mpz_cmp(a, b);
}

template <typename T, integral_fits_t<T, long> = 0 >
static inline int mpz_cmp (mpz_srcptr a, const T b) {
    return mpz_cmp_si(a, b);
}

template <typename T, integral_fits_t<T, unsigned long> = 0 >
static inline int mpz_cmp (mpz_srcptr a, const T b) {
    return mpz_cmp_ui(a, b);
}

inline int mpz_cmp (mpz_srcptr a, const int64_t b) {
    return mpz_cmp_int64(a, b);
}

inline int mpz_cmp (mpz_srcptr a, const uint64_t b) {
    return mpz_cmp_uint64(a, b);
}

#define GMP_AUXX_DEFINE_FUNC3_U(OP) \
static inline void mpz_##OP (mpz_ptr a, mpz_srcptr b, mpz_srcptr c) { ::mpz_##OP(a, b, c); } \
template <typename T, integral_fits_t<T, unsigned long> = 0 > \
static inline void mpz_##OP (mpz_ptr a, mpz_srcptr b, const T c) { mpz_##OP##_ui(a, b, c); } \
inline void mpz_##OP (mpz_ptr a, mpz_srcptr b, const uint64_t c) { mpz_##OP##_uint64(a, b, c); }
#define GMP_AUXX_DEFINE_FUNC3_S(OP) \
template <typename T, integral_fits_t<T, long> = 0 > \
static inline void mpz_##OP (mpz_ptr a, mpz_srcptr b, const T c) { mpz_##OP##_si(a, b, c); } \
inline void mpz_##OP (mpz_ptr a, mpz_srcptr b, const int64_t c) { mpz_##OP##_int64(a, b, c); } \

GMP_AUXX_DEFINE_FUNC3_U(add)
GMP_AUXX_DEFINE_FUNC3_S(add)
GMP_AUXX_DEFINE_FUNC3_U(sub)
GMP_AUXX_DEFINE_FUNC3_S(sub)
GMP_AUXX_DEFINE_FUNC3_U(mul)
GMP_AUXX_DEFINE_FUNC3_S(mul)
GMP_AUXX_DEFINE_FUNC3_U(addmul)
GMP_AUXX_DEFINE_FUNC3_S(addmul)
GMP_AUXX_DEFINE_FUNC3_U(submul)
GMP_AUXX_DEFINE_FUNC3_S(submul)
GMP_AUXX_DEFINE_FUNC3_U(divexact)

static inline int mpz_divisible_p (mpz_ptr a, mpz_srcptr c) {
    return ::mpz_divisible_p(a, c);
}

template <typename T, integral_fits_t<T, unsigned long> = 0 >
static inline int mpz_divisible_p (mpz_ptr a, const T c) {
    return mpz_divisible_ui_p(a, c);
}

inline int mpz_divisible_p (mpz_ptr a, const uint64_t c) {
    return mpz_divisible_uint64_p(a, c);
}

template <typename T, integral_fits_t<T, unsigned long> = 0 >
static inline T
mpz_tdiv_qr (mpz_ptr Q, mpz_ptr R, mpz_srcptr N, T d) {
    return (T) mpz_tdiv_qr_ui(Q, R, N, d);
}

inline uint64_t
mpz_tdiv_qr (mpz_ptr Q, mpz_ptr R, mpz_srcptr N, uint64_t d) {
    return mpz_tdiv_qr_uint64(Q, R, N, d);
}

static inline void
mpz_tdiv_q (mpz_ptr Q, mpz_srcptr N, mpz_srcptr D) {
    return ::mpz_tdiv_q(Q, N, D);
}

template <typename T, integral_fits_t<T, unsigned long> = 0 >
static inline T
mpz_tdiv_q (mpz_ptr Q, mpz_srcptr N, T d) {
    return (T) mpz_tdiv_q_ui(Q, N, d);
}

inline uint64_t
mpz_tdiv_q (mpz_ptr Q, mpz_srcptr N, uint64_t d) {
    return mpz_tdiv_q_uint64(Q, N, d);
}

static inline void
mpz_tdiv_r (mpz_ptr R, mpz_srcptr N, mpz_srcptr D) {
    return ::mpz_tdiv_r(R, N, D);
}

template <typename T, integral_fits_t<T, unsigned long> = 0 >
static inline T
mpz_tdiv_r (mpz_ptr R, mpz_srcptr N, T d) {
    return (T) mpz_tdiv_r_ui(R, N, d);
}

inline uint64_t
mpz_tdiv_r (mpz_ptr R, mpz_srcptr N, uint64_t d) {
    return mpz_tdiv_r_uint64(R, N, d);
}

template <typename T, integral_fits_t<T, unsigned long> = 0 >
static inline T
mpz_tdiv (mpz_srcptr N, T d) {
    return (T) mpz_tdiv_ui(N, d);
}

inline uint64_t
mpz_tdiv (mpz_srcptr N, uint64_t d) {
    return mpz_tdiv_uint64(N, d);
}

}

#endif
