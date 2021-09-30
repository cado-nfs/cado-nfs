#ifndef POLYSELECT_POLY_HEADER_H_
#define POLYSELECT_POLY_HEADER_H_

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

/* structure to store information on N, d, ad, etc... */
struct polyselect_poly_header_s {
  mpz_srcptr N;         // a copy (const pointer), for convenience
  unsigned long d;
  mpz_srcptr ad;        // a copy (const pointer), for convenience
  mpz_t Ntilde;
  mpz_t m0;
};
typedef struct polyselect_poly_header_s polyselect_poly_header_t[1];
typedef struct polyselect_poly_header_s * polyselect_poly_header_ptr;
typedef const struct polyselect_poly_header_s * polyselect_poly_header_srcptr;

void polyselect_poly_header_init (polyselect_poly_header_ptr, mpz_ptr, unsigned long, mpz_ptr);
void polyselect_poly_header_clear (polyselect_poly_header_ptr);
int polyselect_poly_header_skip (polyselect_poly_header_srcptr, unsigned long);


#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_POLY_HEADER_H_ */
