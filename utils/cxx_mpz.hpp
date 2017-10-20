#ifndef CXX_MPZ_HPP_
#define CXX_MPZ_HPP_

#include <gmp.h>

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
#endif	/* CXX_MPZ_HPP_ */
