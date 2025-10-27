#include "cado.h" // IWYU pragma: keep

#include <stdlib.h>

#include <gmp.h>

#include "macros.h"
#include "mpz_mat.h"

void test_mpz_mat_LLL()
{
  mpz_mat Mroot;
  mpz_mat M;
  mpz_mat U;
  mpz_mat Mtmp;
  mpz_t det;
  mpz_t a;
  mpz_t b;

  mpz_mat_init(Mroot, 5, 5);
  mpz_mat_init(M, 5, 5);
  mpz_mat_init(U, 5, 5);
  mpz_mat_init(Mtmp, 5, 5);
  mpz_init(det);
  mpz_init(a);
  mpz_init(b);

  mpz_set_str(mpz_mat_entry(M, 0, 0), "303263745070282", 10);
  mpz_set_str(mpz_mat_entry(M, 0, 1), "534412315049442", 10);
  mpz_set_str(mpz_mat_entry(M, 0, 2), "420777937938826", 10);
  mpz_set_str(mpz_mat_entry(M, 0, 3), "1072079634362711", 10);
  mpz_set_str(mpz_mat_entry(M, 0, 4), "292981222244052", 10);
  mpz_set_str(mpz_mat_entry(M, 1, 0), "-276173394456956", 10);
  mpz_set_str(mpz_mat_entry(M, 1, 1), "-519451951511957", 10);
  mpz_set_str(mpz_mat_entry(M, 1, 2), "-601613904582636", 10);
  mpz_set_str(mpz_mat_entry(M, 1, 3), "143493794871706", 10);
  mpz_set_str(mpz_mat_entry(M, 1, 4), "1026130037924106", 10);
  mpz_set_str(mpz_mat_entry(M, 2, 0), "-204300942856958", 10);
  mpz_set_str(mpz_mat_entry(M, 2, 1), "-257594425461071", 10);
  mpz_set_str(mpz_mat_entry(M, 2, 2), "-84231523408738", 10);
  mpz_set_str(mpz_mat_entry(M, 2, 3), "-755649773432484", 10);
  mpz_set_str(mpz_mat_entry(M, 2, 4), "-845162459661229", 10);
  mpz_set_str(mpz_mat_entry(M, 3, 0), "497878367605355", 10);
  mpz_set_str(mpz_mat_entry(M, 3, 1), "-1088691215483927", 10);
  mpz_set_str(mpz_mat_entry(M, 3, 2), "-578088569966458", 10);
  mpz_set_str(mpz_mat_entry(M, 3, 3), "896932555709611", 10);
  mpz_set_str(mpz_mat_entry(M, 3, 4), "738312266252111", 10);
  mpz_set_str(mpz_mat_entry(M, 4, 0), "-375049643008522", 10);
  mpz_set_str(mpz_mat_entry(M, 4, 1), "-104542010743550", 10);
  mpz_set_str(mpz_mat_entry(M, 4, 2), "594090013097381", 10);
  mpz_set_str(mpz_mat_entry(M, 4, 3), "-777686688213268", 10);
  mpz_set_str(mpz_mat_entry(M, 4, 4), "-602810852352675", 10);
  mpz_mat_set(Mroot, M);
  mpz_set_ui(a, 99);
  mpz_set_ui(b, 100);

  mpz_mat_LLL(det, M, U, a, b);

  mpz_set_str(a,
      "256607830916257875370010914282002409827349022739443960549533444356139981381",
      10);
  mpz_mul(a, a, a);
  ASSERT_ALWAYS(mpz_cmpabs(a, det) == 0);

  mpz_mat_mul(Mtmp, U, Mroot);
  ASSERT_ALWAYS(mpz_mat_cmp(M, Mtmp) == 0);

  //Maybe a little bit strict.
  /*mpz_set_str(a, "-177210592243632", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 0, 0), a) == 0);*/
  /*mpz_set_str(a, "-242634061923586", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 0, 1), a) == 0);*/
  /*mpz_set_str(a, "-265067490052548", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 0, 2), a) == 0);*/
  /*mpz_set_str(a, "459923655801933", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 0, 3), a) == 0);*/
  /*mpz_set_str(a, "473948800506929", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 0, 4), a) == 0);*/
  /*mpz_set_str(a, "98962802213324", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 1, 0), a) == 0);*/
  /*mpz_set_str(a, "276817889588371", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 1, 1), a) == 0);*/
  /*mpz_set_str(a, "336546414530088", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 1, 2), a) == 0);*/
  /*mpz_set_str(a, "316429860930227", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 1, 3), a) == 0);*/
  /*mpz_set_str(a, "-552181237417177", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 1, 4), a) == 0);*/
  /*mpz_set_str(a, "-170748700151564", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 2, 0), a) == 0);*/
  /*mpz_set_str(a, "153052414717521", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 2, 1), a) == 0);*/
  /*mpz_set_str(a, "678321536506119", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 2, 2), a) == 0);*/
  /*mpz_set_str(a, "-22036914780784", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 2, 3), a) == 0);*/
  /*mpz_set_str(a, "242351607308554", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 2, 4), a) == 0);*/
  /*mpz_set_str(a, "552260235252154", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 3, 0), a) == 0);*/
  /*mpz_set_str(a, "347176072667136", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 3, 1), a) == 0);*/
  /*mpz_set_str(a, "-329022523044833", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 3, 2), a) == 0);*/
  /*mpz_set_str(a, "317763032411335", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 3, 3), a) == 0);*/
  /*mpz_set_str(a, "128862051845746", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 3, 4), a) == 0);*/
  /*mpz_set_str(a, "681550851941055", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 4, 0), a) == 0);*/
  /*mpz_set_str(a, "-450370676919234", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 4, 1), a) == 0);*/
  /*mpz_set_str(a, "630367946644757", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 4, 2), a) == 0);*/
  /*mpz_set_str(a, "-44951670675039", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 4, 3), a) == 0);*/
  /*mpz_set_str(a, "32766272546807", 10);*/
  /*ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry_const(M, 4, 4), a) == 0);*/

  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 0, 0), 1));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 0, 1), 1));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 0, 2), 1));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 0, 3), 0));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 0, 4), 0));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 1, 0), 1));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 1, 1), 0));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 1, 2), 1));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 1, 3), 0));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 1, 4), 0));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 2, 0), 0));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 2, 1), 0));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 2, 2), -1));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 2, 3), 0));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 2, 4), 1));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 3, 0), -1));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 3, 1), -1));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 3, 2), -1));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 3, 3), 0));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 3, 4), -1));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 4, 0), -2));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 4, 1), -2));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 4, 2), -3));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 4, 3), 1));*/
  /*ASSERT_ALWAYS(mpz_cmp_si(mpz_mat_entry_const(M, 4, 4), 1));*/

  mpz_clear(det);
  mpz_clear(a);
  mpz_clear(b);
  mpz_mat_clear(M);
  mpz_mat_clear(Mroot);
  mpz_mat_clear(U);
  mpz_mat_clear(Mtmp);
}

int main()
{
  test_mpz_mat_LLL();
  return EXIT_SUCCESS;
}
