#ifndef CADO_MPZ_POLYMODF_H
#define CADO_MPZ_POLYMODF_H

#include "mpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Let F(x) be a (non-monic) polynomial of degree d:
   F(x) = f_d x^d + f_{d-1} x^{d-1} + .... + f_1 x + f_0
   The following type represents a polynomial modulo F(x):
   If P is such an element, it means: P = P.p / f_d^P.v */
struct mpz_polymodF_s {
  mpz_poly p;
  unsigned long v;
};

typedef struct mpz_polymodF_s mpz_polymodF[1];
typedef struct mpz_polymodF_s * mpz_polymodF_ptr;
typedef const struct mpz_polymodF_s * mpz_polymodF_srcptr;

void mpz_polymodF_init(mpz_polymodF_ptr P, int x);
void mpz_polymodF_clear(mpz_polymodF_ptr P);
void mpz_polymodF_swap(mpz_polymodF_ptr P, mpz_polymodF_ptr Q);
void mpz_polymodF_set_ui(mpz_polymodF_ptr P, unsigned long x);
void mpz_polymodF_set(mpz_polymodF_ptr P, mpz_polymodF_srcptr Q);
void mpz_polymodF_set_from_ab(mpz_polymodF_ptr P, mpz_srcptr a, mpz_srcptr b);
void mpz_poly_reducemodF(mpz_polymodF_ptr P, mpz_poly_ptr p, mpz_poly_srcptr F);
void mpz_polymodF_mul(mpz_polymodF_ptr Q, mpz_polymodF_srcptr P1, mpz_polymodF_srcptr P2, mpz_poly_srcptr F);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/* This is sort of a generic way to write a c++ equivalent to the C type.
 * The first-class citizen in the cado-nfs code is (still) the C type, so
 * we're definitely bound to have a few infelicities here:
 *  - the type name can't be the same because of the size-1 array trick
 *    in C.
 *  - the C type is embedded as a member x for the same reason.
 *  - most operations on the C type should go through the member x
 *    (however, the conversions we have to _ptr and _srcptr can ease
 *    things a bit).
 */
struct cxx_mpz_polymodF {
    mpz_polymodF x;
    cxx_mpz_polymodF() { mpz_polymodF_init(x, -1); }
    cxx_mpz_polymodF(int deg) { mpz_polymodF_init(x, deg); }
    inline int degree() const { return x->p->deg; } /* handy */
    cxx_mpz_polymodF(mpz_poly_srcptr f, unsigned long scale=0) {
        mpz_polymodF_init(x, -1);
        mpz_poly_set(x->p, f);
        x->v = scale;
    }
    cxx_mpz_polymodF(mpz_polymodF_srcptr f) { mpz_polymodF_init(x, -1); mpz_polymodF_set(x, f); }
    ~cxx_mpz_polymodF() { mpz_polymodF_clear(x); }
    cxx_mpz_polymodF(cxx_mpz_polymodF const & o) {
        mpz_polymodF_init(x, -1);
        mpz_polymodF_set(x, o.x);
    }
    cxx_mpz_polymodF & operator=(cxx_mpz_polymodF const & o) {
        mpz_polymodF_set(x, o.x);
        return *this;
    }
    cxx_mpz_polymodF(cxx_mpz_polymodF && o) {
        mpz_polymodF_init(x, -1);
        mpz_polymodF_swap(x, o.x);
    }
    cxx_mpz_polymodF& operator=(cxx_mpz_polymodF && o) {
        mpz_polymodF_swap(x, o.x);
        return *this;
    }
    operator mpz_polymodF_ptr() { return x; }
    operator mpz_polymodF_srcptr() const { return x; }
    mpz_polymodF_ptr operator->() { return x; }
    mpz_polymodF_srcptr operator->() const { return x; }
    std::string print_poly(std::string const& var) const;
};

std::ostream& operator<<(std::ostream& o, cxx_mpz_polymodF const & f);
std::istream& operator>>(std::istream& in, cxx_mpz_polymodF & f);

#if GNUC_VERSION_ATLEAST(4,3,0)
extern void mpz_poly_init(cxx_mpz_poly & pl, int) __attribute__((error("mpz_poly_init must not be called on a mpz_poly reference -- it is the caller's business (via a ctor)")));
extern void mpz_poly_clear(cxx_mpz_poly & pl) __attribute__((error("mpz_poly_clear must not be called on a mpz_poly reference -- it is the caller's business (via a dtor)")));
#endif

#endif


/* implementation is in mpz_poly.cpp */



#endif	/* CADO_MPZ_POLYMODF_H */
