#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "macros.h"
#include "tests_common.h"
#include "usp.h"

#define MAX_DEGREE 100

void
test_usp ()
{
  mpz_t p[MAX_DEGREE], u;
  int n, i, d;
  usp_root_data R[MAX_DEGREE];
  double root;

  for (i = 0; i < MAX_DEGREE; i++)
    {
      mpz_init (p[i]);
      mpz_init (R[i].a);
      mpz_init (R[i].b);
    }

  /* polynomial x */
  mpz_set_ui (p[0], 0);
  mpz_set_ui (p[1], 1);
  n = numberOfRealRoots (p, 1, 2, 1, R);
  ASSERT_ALWAYS (n == 1);
  /* check [a/2^ka, b/2^kb] contains root 0 */
  ASSERT_ALWAYS (mpz_cmp_ui (R[0].a, 0) <= 0);
  ASSERT_ALWAYS (0 <= mpz_cmp_ui (R[0].b, 0));

  /* polynomial 2*x */
  mpz_set_ui (p[0], 0);
  mpz_set_ui (p[1], 2);
  n = numberOfRealRoots (p, 1, 2, 1, R);
  ASSERT_ALWAYS (n == 1);
  root = rootRefine (R, p, 1, 1e-9);
  ASSERT_ALWAYS(-0.000001 <= root && root <= 0.000001);

  /* polynomial x+1 */
  mpz_set_ui (p[0], 1);
  mpz_set_ui (p[1], 1);
  n = numberOfRealRoots (p, 1, 0, 0, NULL);
  ASSERT_ALWAYS (n == 1);

  /* polynomial x-1 */
  mpz_set_si (p[0], -1);
  mpz_set_ui (p[1], 1);
  n = numberOfRealRoots (p, 1, 0, 0, NULL);
  ASSERT_ALWAYS (n == 1);

  /* polynomial x^2+1 */
  mpz_set_si (p[0], 1);
  mpz_set_ui (p[1], 0);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 0);

  /* polynomial x^2+2 */
  mpz_set_si (p[0], 2);
  mpz_set_ui (p[1], 0);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 0);

  /* polynomial x^2-1 */
  mpz_set_si (p[0], -1);
  mpz_set_ui (p[1], 0);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 2);

  /* polynomial x^2-3*x+2 */
  mpz_set_si (p[0], 2);
  mpz_set_si (p[1], -3);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 2);

  /* polynomial x^2-2*x+2 */
  mpz_set_si (p[0], 2);
  mpz_set_si (p[1], -2);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 0);

  /* polynomial x^2+1000*x+2 */
  mpz_set_si (p[0], 2);
  mpz_set_si (p[1], 1000);
  mpz_set_ui (p[2], 1);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 2);

  /* polynomial (x-1)*(x-2)*(x-3) */
  mpz_set_si (p[0], -6);
  mpz_set_si (p[1], 11);
  mpz_set_si (p[2], -6);
  mpz_set_si (p[3], 1);
  n = numberOfRealRoots (p, 3, 0, 0, NULL);
  ASSERT_ALWAYS (n == 3);

  /* polynomial (x-1)*(x-2)*(x-3)*(x-4) */
  mpz_set_si (p[0], 24);
  mpz_set_si (p[1], -50);
  mpz_set_si (p[2], 35);
  mpz_set_si (p[3], -10);
  mpz_set_si (p[4], 1);
  n = numberOfRealRoots (p, 4, 0, 0, NULL);
  ASSERT_ALWAYS (n == 4);

  /* polynomial (x-1)*(x-2)*(x-3)*(x-4)*(x-5) */
  mpz_set_si (p[0], -120);
  mpz_set_si (p[1], 274);
  mpz_set_si (p[2], -225);
  mpz_set_si (p[3], 85);
  mpz_set_si (p[4], -15);
  mpz_set_si (p[5], 1);
  n = numberOfRealRoots (p, 5, 0, 0, NULL);
  ASSERT_ALWAYS (n == 5);

  /* polynomial (x-1)*(x-2)*(x-3)*(x-4)*(x-5)*(x-6) */
  mpz_set_si (p[0], 720);
  mpz_set_si (p[1], -1764);
  mpz_set_si (p[2], 1624);
  mpz_set_si (p[3], -735);
  mpz_set_si (p[4], 175);
  mpz_set_si (p[5], -21);
  mpz_set_si (p[6], 1);
  n = numberOfRealRoots (p, 6, 0, 0, NULL);
  ASSERT_ALWAYS (n == 6);

  /* polynomial 2*(x-1)*(x-2)*(x-3/2) */
  mpz_set_si (p[0], -6);
  mpz_set_si (p[1], 13);
  mpz_set_si (p[2], -9);
  mpz_set_si (p[3], 2);
  n = numberOfRealRoots (p, 3, 0, 0, NULL);
  ASSERT_ALWAYS (n == 3);

  /* polynomial 8*(x-1)*(x-2)*(x-3/2)*(x-5/4) */
  mpz_set_si (p[0], 30);
  mpz_set_si (p[1], -89);
  mpz_set_si (p[2], 97);
  mpz_set_si (p[3], -46);
  mpz_set_si (p[4], 8);
  n = numberOfRealRoots (p, 4, 0, 0, NULL);
  ASSERT_ALWAYS (n == 4);

  /* polynomial 2*(x-1)*(x-3/2)*(x-2)*(x-3)*(x-4) */
  mpz_set_si (p[0], -72);
  mpz_set_si (p[1], 198);
  mpz_set_si (p[2], -205);
  mpz_set_si (p[3], 100);
  mpz_set_si (p[4], -23);
  mpz_set_si (p[5], 2);
  n = numberOfRealRoots (p, 5, 0, 0, NULL);
  ASSERT_ALWAYS (n == 5);

  mpz_set_si (p[0], -3);
  mpz_set_si (p[1], 2);
  n = numberOfRealRoots (p, 1, 0, 0, NULL);
  ASSERT_ALWAYS (n == 1);

  /* 2*(x-17)*(x-17/5) */
  mpz_set_si (p[0], 595);
  mpz_set_si (p[1], -69);
  mpz_set_si (p[2], 2);
  n = numberOfRealRoots (p, 2, 0, 0, NULL);
  ASSERT_ALWAYS (n == 2);

  /* 2000*x^3 - 35000*x^2 + 2*x - 35 */
  mpz_set_si (p[0], -35);
  mpz_set_si (p[1], 2);
  mpz_set_si (p[2], -35000);
  mpz_set_si (p[3], 2000);
  n = numberOfRealRoots (p, 3, 0, 0, NULL);
  ASSERT_ALWAYS (n == 1);

  /* 33*x^4-124*x^3+137*x^2-32*x-11 */
  mpz_set_si (p[0], -11);
  mpz_set_si (p[1], -32);
  mpz_set_si (p[2], 137);
  mpz_set_si (p[3], -124);
  mpz_set_si (p[4], 33);
  n = numberOfRealRoots (p, 4, 0, 0, NULL);
  ASSERT_ALWAYS (n == 4);

  /* 2*x^4 + 4*x^3 - 4*x^2 + 2*x + 3 */
  mpz_set_si (p[0], 3);
  mpz_set_si (p[1], 2);
  mpz_set_si (p[2], -4);
  mpz_set_si (p[3], 4);
  mpz_set_si (p[4], 2);
  n = numberOfRealRoots (p, 4, 0, 0, NULL);
  ASSERT_ALWAYS (n == 2);

  /* 2*x^4 + x^3 - 3*x^2 + 2*x + 3 */
  mpz_set_si (p[0], 3);
  mpz_set_si (p[1], 2);
  mpz_set_si (p[2], -3);
  mpz_set_si (p[3], 1);
  mpz_set_si (p[4], 2);
  n = numberOfRealRoots (p, 4, 0, 0, NULL);
  ASSERT_ALWAYS (n == 2);

  /* infinite loop that appeared on 32-bit computers, 2nd April 2014 */
  mpz_set_str (p[0], "202156081496031128767910903808", 10);
  mpz_set_str (p[1], "-1014774808369763543554925264896", 10);
  mpz_set_str (p[2], "305592973565950609559904059392", 10);
  mpz_set_str (p[3], "532741928198739321928042938368", 10);
  mpz_set_str (p[4], "-265275048860303928171627020288", 10);
  mpz_set_str (p[5], "152821965764546794524860481536", 10);
  mpz_set_str (p[6], "-40951570592501524318484692992", 10);
  mpz_set_ui (R[0].a, 2);
  mpz_set_ui (R[0].b, 4);
  R[0].ka = R[0].kb = 0;
  root = rootRefine (R, p, 6, 1e-9); /* root is near 3.00763029864372 */
  ASSERT_ALWAYS(3.00763 <= root && root <= 3.00764);

  mpz_init (u);
  mpz_set_ui (u, 1);
  mpz_mul_2exp (u, u, 127);
  for (d = 1; d < MAX_DEGREE; d++)
    {
      for (i = 0; i <= d; i++)
        {
          mpz_urandomb (p[i], state, 128);
          mpz_sub (p[i], p[i], u);
        }
      if (mpz_cmp_ui (p[d], 0) == 0)
        mpz_set_ui (p[d], 1);
      n = numberOfRealRoots (p, d, 0, 0, NULL);
      ASSERT_ALWAYS (0 <= n && n <= d);
    }

  /* check with large coefficients */
  mpz_urandomb (p[0], state, 2048);
  mpz_urandomb (p[1], state, 2048);
  mpz_urandomb (p[2], state, 2048);
  mpz_urandomb (p[3], state, 2048);
  n = numberOfRealRoots (p, 3, 0, 0, NULL);
  ASSERT_ALWAYS (0 <= n && n <= 3);

  /* bug #18879 */
  mpz_set_str (p[0], "-399370140104138021004049555602050048725505629855192896740786176", 10);
  mpz_set_str (p[1], "6152454519896677691538310427796089495118818146087045879381884928", 10);
  mpz_set_str (p[2], "-44005628564777919654831813189277128400998342025364587400350662656", 10);
  mpz_set_str (p[3], "193693587647452249746878956308340478214944212045235974402799042560", 10);
  mpz_set_str (p[4], "-586131524042923039573328295315283165577996285519208621852726394880", 10);
  mpz_set_str (p[5], "1289948018788260792668169964588439385725345410925233394477273448448", 10);
  mpz_set_str (p[6], "-2129171214885956678206724080558522465836045133020822627924230275072", 10);
  mpz_set_str (p[7], "2677624344970803918386735215388265920949327908075617878580737867776", 10);
  mpz_set_str (p[8], "-2578130043460498861049699336930493887253985498977331087052720046080", 10);
  mpz_set_str (p[9], "1891301125005133082116243679283866801297483539999500792756728496128", 10);
  mpz_set_str (p[10], "-1040585580532476627932685878679921764622233387744782090307763175424", 10);
  mpz_set_str (p[11], "416382268948646574906631753282253840199545851564089327142421135360", 10);
  mpz_set_str (p[12], "-114545848993024762527402482645186791330188904180948233736106803200", 10);
  mpz_set_str (p[13], "19391576494070251913587093980476025365869714617333619797109243904", 10);
  mpz_set_str (p[14], "-1524165754529617102335904328302294916528952457696967340985417728", 10);
  n = numberOfRealRoots (p, 14, 0, 0, R);
  ASSERT_ALWAYS (n == 2);
  root = rootRefine (R, p, 14, 6.103516e-05);
  root = rootRefine (R + 1, p, 14, 6.103516e-05);

  mpz_set_str (p[0], "-1907753428620166438197776483820542501644580146356003274752", 10);
  mpz_set_str (p[1], "15452436600339983267007229667323593777349076211560511176704", 10);
  mpz_set_str (p[2], "-54758274613009131316199852690823331935974939002324115259392", 10);
  mpz_set_str (p[3], "110882878534980127672983781446528375367066401782722337964032", 10);
  mpz_set_str (p[4], "-140332817723613361637526395534759643288947238768788031668224", 10);
  mpz_set_str (p[5], "113666888827955590415811847475276525451763140364809078308864", 10);
  mpz_set_str (p[6], "-57542498902827796543386171783637664156309699559757732380672", 10);
  mpz_set_str (p[7], "16645828445978446377629358188283303394024754747896961695744", 10);
  mpz_set_str (p[8], "-2106687741183691798122020284370616107326915629763006562304", 10);
  mpz_set_ui (R[0].a, 63);
  mpz_set_ui (R[0].b, 64);
  R[0].ka = R[0].kb = 6;
  root = rootRefine (R, p, 8, 0.000244);

  for (i = 0; i < MAX_DEGREE; i++)
    {
      mpz_clear (p[i]);
      mpz_clear (R[i].a);
      mpz_clear (R[i].b);
    }
  mpz_clear (u);
}

int
main (int argc, const char *argv[])
{
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  test_usp ();
  tests_common_clear();
  exit (EXIT_SUCCESS);
}
