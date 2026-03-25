#include "cado.h" // IWYU pragma: keep

#include <cstdlib>

#include <sstream>
#include <istream>
#include <fstream>
#include <string>

#include <gmp.h>

#include "params.h"
#include "tests_common.h"
#include "cado_poly.hpp"
#include "mpz_poly.h"
#include "macros.h"

static void
test_cado_poly_set ()
{
  cxx_cado_poly p, q;

  const double s = 3.1415;

  mpz_set_ui (q.n, 1000000007);
  q.skew = s;
  q.provision_new_poly();
  q.provision_new_poly();
  q[0]->deg = 1;
  mpz_poly_setcoeff_si(q[0], 0, -123128869);
  mpz_poly_setcoeff_si(q[0], 1, 1000000008);
  q[1]->deg = 2;
  mpz_poly_setcoeff_ui(q[1], 0, 228868283);
  mpz_poly_setcoeff_ui(q[1], 1, 887036294);
  mpz_poly_setcoeff_ui(q[1], 2, 429156742);

  p = q;

  ASSERT_ALWAYS (mpz_cmp_ui (p.n, 1000000007) == 0);
  ASSERT_ALWAYS (p.skew == s);
  ASSERT_ALWAYS (mpz_cmp_si (mpz_poly_coeff_const(p[0], 0), -123128869) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (mpz_poly_coeff_const(p[0], 1), 1000000008) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (mpz_poly_coeff_const(p[1], 0), 228868283) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (mpz_poly_coeff_const(p[1], 1), 887036294) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (mpz_poly_coeff_const(p[1], 2), 429156742) == 0);
  cxx_mpz m;
  p.getm(m, p.n);
  ASSERT_ALWAYS (mpz_cmp_ui (m, 123128869) == 0);
}

static int cado_poly_read (cxx_cado_poly & poly, std::istream& is)
{
  cxx_param_list pl;
  param_list_read(pl, is, false);
  poly = cxx_cado_poly(cxx_cado_poly::plist {}, pl);
  int const r = 1;
  return r;
}


static void test_cado_poly_sanitycheck_stream(std::istream & is)
{
    cxx_cado_poly cpoly;
    cado_poly_read(cpoly, is);
    cxx_mpz_poly G;
    int const ret = cpoly.check_mapping(G, cpoly.n);
    ASSERT_ALWAYS(ret != 0);
}

static void test_cado_poly_sanitychecks_static()
{
    std::string const s =
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

static void test_cado_poly_sanitycheck_file(const char * file)
{
    std::ifstream is(file);
    ASSERT_ALWAYS(is);
    test_cado_poly_sanitycheck_stream(is);
}

// coverity[root_function]
int main(int argc, char const * argv[])
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
