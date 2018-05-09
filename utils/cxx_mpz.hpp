#ifndef CXX_MPZ_HPP_
#define CXX_MPZ_HPP_
#include "macros.h"

#include <gmp.h>
#include <istream>
#include <ostream>

struct cxx_mpz {
    mpz_t x;
    cxx_mpz() { mpz_init(x); }
    cxx_mpz(unsigned long p) { mpz_init_set_ui(x, p); }
    ~cxx_mpz() { mpz_clear(x); }
    cxx_mpz(cxx_mpz const & o) {
        mpz_init(x);
        mpz_set(x, o.x);
    }
    cxx_mpz & operator=(cxx_mpz const & o) {
        mpz_set(x, o.x);
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

inline std::ostream& operator<<(std::ostream& os, cxx_mpz const& x) { return os << (mpz_srcptr) x; }
inline std::ostream& operator<<(std::ostream& os, cxx_mpq const& x) { return os << (mpq_srcptr) x; }
inline std::istream& operator>>(std::istream& is, cxx_mpz & x) { return is >> (mpz_ptr) x; }
inline std::istream& operator>>(std::istream& is, cxx_mpq & x) { return is >> (mpq_ptr) x; }
#endif	/* CXX_MPZ_HPP_ */
