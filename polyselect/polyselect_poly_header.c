#include "cado.h"
#include <gmp.h>
#include "polyselect_poly_header.h"


/* init the header struct */
void
polyselect_poly_header_init (polyselect_poly_header_ptr header)
{
  /* compute Ntilde, m0 */
  mpz_init(header->N);
  mpz_init (header->Ntilde);
  mpz_init (header->m0);
  header->d = -1;
  mpz_init(header->ad);
}

void
polyselect_poly_header_set_Nd (polyselect_poly_header_ptr header, mpz_srcptr N, int d)
{
  mpz_set(header->N, N);
  header->d = d;
}

void
polyselect_poly_header_set_ad (polyselect_poly_header_ptr header,
              mpz_srcptr ad )
{
  /* compute Ntilde, m0 */
  mpz_set(header->ad, ad);

  /* compute Ntilde, ... from N, ... */
  mpz_set (header->Ntilde, header->ad);
  mpz_mul_ui (header->Ntilde, header->Ntilde, header->d);
  mpz_pow_ui (header->Ntilde, header->Ntilde, header->d - 1);
  mpz_mul_ui (header->Ntilde, header->Ntilde, header->d);
  mpz_mul (header->Ntilde, header->Ntilde, header->N); /* d^d * ad^(d-1) * N */
  mpz_root (header->m0, header->Ntilde, header->d);
}


void polyselect_poly_header_set (polyselect_poly_header_ptr to, polyselect_poly_header_srcptr from)
{
    mpz_set(to->N, from->N);
    mpz_set(to->Ntilde, from->Ntilde);
    mpz_set(to->m0, from->m0);
    mpz_set(to->ad, from->ad);
    to->d = from->d;
}

/* clear header struct */
void
polyselect_poly_header_clear (polyselect_poly_header_ptr header )
{
  mpz_clear(header->N);
  mpz_clear (header->m0);
  mpz_clear (header->Ntilde);
  mpz_clear (header->ad);
}

int
polyselect_poly_header_skip (polyselect_poly_header_srcptr header, unsigned long p)
{
  return header->d % p == 0 || mpz_divisible_ui_p (header->ad, p);
}

