#include "cado.h"
#include <gmp.h>
#include "polyselect_poly_header.h"


/* init the header struct */
void
polyselect_poly_header_init (polyselect_poly_header_ptr header,
              mpz_ptr N,
              unsigned long d,
              mpz_ptr ad )
{
  /* compute Ntilde, m0 */
  mpz_init_set (header->N, N);
  mpz_init (header->Ntilde);
  mpz_init (header->m0);
  header->d = d;
  mpz_init_set (header->ad, ad);

  /* compute Ntilde, ... from N, ... */
  mpz_set (header->Ntilde, header->ad);
  mpz_mul_ui (header->Ntilde, header->Ntilde, header->d);
  mpz_pow_ui (header->Ntilde, header->Ntilde, header->d - 1);
  mpz_mul_ui (header->Ntilde, header->Ntilde, header->d);
  mpz_mul (header->Ntilde, header->Ntilde, header->N); /* d^d * ad^(d-1) * N */
  mpz_root (header->m0, header->Ntilde, header->d);
}


/* clear header struct */
void
polyselect_poly_header_clear (polyselect_poly_header_ptr header )
{
  mpz_clear (header->m0);
  mpz_clear (header->Ntilde);
  mpz_clear (header->N);
  mpz_clear (header->ad);
}

int
polyselect_poly_header_skip (polyselect_poly_header_srcptr header, unsigned long p)
{
  return header->d % p == 0 || mpz_divisible_ui_p (header->ad, p);
}

