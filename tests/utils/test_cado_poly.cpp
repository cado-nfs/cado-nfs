#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <gmp.h>
#include <sstream>
#include <istream>
#include <fstream>
#include "tests_common.h"
#include "cado_poly.h"
#include "mpz_poly.h"
#include "macros.h"

void
test_cado_poly_set ()
{
  cado_poly p, q;

  cado_poly_init (p);
  cado_poly_init (q);

  const double s = 3.1415;

  mpz_set_ui (q->n, 1000000007);
  q->skew = s;
  cado_poly_provision_new_poly(q);
  cado_poly_provision_new_poly(q);
  q->pols[0]->deg = 1;
  mpz_poly_setcoeff_si(q->pols[0], 0, -123128869);
  mpz_poly_setcoeff_si(q->pols[0], 1, 1000000008);
  q->pols[1]->deg = 2;
  mpz_poly_setcoeff_ui(q->pols[1], 0, 228868283);
  mpz_poly_setcoeff_ui(q->pols[1], 1, 887036294);
  mpz_poly_setcoeff_ui(q->pols[1], 2, 429156742);

  cado_poly_set (p, q);

  ASSERT_ALWAYS (mpz_cmp_ui (p->n, 1000000007) == 0);
  ASSERT_ALWAYS (p->skew == s);
  ASSERT_ALWAYS (mpz_cmp_si (p->pols[0]->coeff[0], -123128869) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (p->pols[0]->coeff[1], 1000000008) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (p->pols[1]->coeff[0], 228868283) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (p->pols[1]->coeff[1], 887036294) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (p->pols[1]->coeff[2], 429156742) == 0);
  mpz_t m;
  mpz_init(m);
  cado_poly_getm(m, p, p->n);
  ASSERT_ALWAYS (mpz_cmp_ui (m, 123128869) == 0);
  mpz_clear(m);
  cado_poly_clear (p);
  cado_poly_clear (q);
}

// this could go to the public cado_poly interface, but that interface is
// C only at the moment...
//
// returns 0 on failure, 1 on success.
extern "C" int cado_poly_set_plist(cado_poly_ptr cpoly, param_list_ptr pl);

int cado_poly_read (cxx_cado_poly & poly, std::istream& is)
{
  cxx_param_list pl;
  param_list_read(pl, is, 0);
  int r = cado_poly_set_plist (poly, pl);
  return r;
}


void test_cado_poly_sanitycheck_stream(std::istream & is)
{
    cxx_cado_poly cpoly;
    cado_poly_read(cpoly, is);
    cxx_mpz_poly G;
    int ret = cado_poly_check_mapping(G, cpoly, cpoly->n);
    ASSERT_ALWAYS(ret != 0);
}

void test_cado_poly_sanitychecks_static()
{
    std::string s =
        "n: 54022122323205311359700529131254845253584832080092810873601245077747279904751944559089001546838958178759103\n"
        "skew: 1\n"
        "c0: -3812358699277286779054168\n"
        "c1: -565926392227485137902\n"
        "c2: 182765314800266891\n"
        "c3: 6921015189226\n"
        "c4: -494461352\n"
        "c5: 12480\n"
        "Y0: -337674472307257214145\n"
        "Y1: 120135550263449\n"
        ;

    std::istringstream is(s);
    test_cado_poly_sanitycheck_stream(is);
}

void test_cado_poly_sanitycheck_file(const char * file)
{
    std::ifstream is(file);
    ASSERT_ALWAYS(is);
    test_cado_poly_sanitycheck_stream(is);
}

// coverity[root_function]
int
main (int argc, const char *argv[])
{
  tests_common_cmdline(&argc, &argv, 0);
  test_cado_poly_set ();

  test_cado_poly_sanitychecks_static();

  for( ; argc > 1 ; argc--,argv++) {
      /* Test another poly file */
      test_cado_poly_sanitycheck_file(argv[1]);
  }

  tests_common_clear();
  exit (EXIT_SUCCESS);
}
