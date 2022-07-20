#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include "bit_vector.h"
#include "tests_common.h"
#include "macros.h"
#include "portability.h" // IWYU pragma: keep

void
test_bit_vector (unsigned long iter)
{
  bit_vector b, c;
  size_t n, i, w;
  int ret;

  while (iter--)
    {
      n = 1 + gmp_urandomm_ui(state, 3 * BV_BITS);
      w = (n + BV_BITS - 1) / BV_BITS;
      bit_vector_init (b, n);
      ASSERT_ALWAYS(bit_vector_memory_footprint (b) == w * sizeof (bv_t));
      bit_vector_init (c, n);
      bit_vector_set (b, 0);
      bit_vector_neg (c, b);
      for (i = 0; i < n; i++)
        {
          ASSERT_ALWAYS(bit_vector_getbit (b, i) == 0);
          ASSERT_ALWAYS(bit_vector_getbit (c, i) == 1);
        }
      for (i = n; i < w * BV_BITS; i++)
        {
          ASSERT_ALWAYS(bit_vector_getbit (b, i) == 0);
          ASSERT_ALWAYS(bit_vector_getbit (c, i) == 0);
        }
      bit_vector_clear (b);
      bit_vector_init_set (b, n, 1);
      bit_vector_neg (c, b);
      for (i = 0; i < n; i++)
        {
          ASSERT_ALWAYS(bit_vector_getbit (b, i) == 1);
          ASSERT_ALWAYS(bit_vector_getbit (c, i) == 0);
        }
      for (i = n; i < w * BV_BITS; i++)
        {
          ASSERT_ALWAYS(bit_vector_getbit (b, i) == 0);
          ASSERT_ALWAYS(bit_vector_getbit (c, i) == 0);
        }
      i = gmp_urandomm_ui(state, n);
      bit_vector_setbit (b, i);
      ASSERT_ALWAYS(bit_vector_getbit (b, i) == 1);
      ret = bit_vector_clearbit (b, i);
      ASSERT_ALWAYS(bit_vector_getbit (b, i) == 0);
      ASSERT_ALWAYS(ret == 1);
      ret = bit_vector_flipbit (b, i);
      ASSERT_ALWAYS(bit_vector_getbit (b, i) == 1);
      ASSERT_ALWAYS(ret == 0);
      for (i = 0; i < n; i++)
        bit_vector_flipbit (b, gmp_urandomm_ui(state, n));
      for (w = 0, i = 0; i < n; i++)
        w += bit_vector_getbit (b, i);
      ASSERT_ALWAYS(w == bit_vector_popcount (b));
      bit_vector_clear (b);
      bit_vector_clear (c);
    }
}

void
test_bit_vector_read_from_file (void)
{
  bit_vector b, c;
  size_t n, i;

  /* create a random bit vector */
  n = 1 + gmp_urandomm_ui(state, 3 * BV_BITS);
  bit_vector_init_set (b, n, 0);
  for (i = 0; i < n; i++)
    if (gmp_urandomb_ui(state, 1))
      bit_vector_flipbit (b, i);


  // coverity complains about insecure temp files. For tests, I don't
  // think it's a problem, really.
  // coverity[secure_temp]
  FILE * f = tmpfile();
  bit_vector_write_to_stream (b, f);
  fseek(f, 0L, SEEK_SET);
  bit_vector_init (c, n);
  bit_vector_read_from_stream (c, f);
  fclose(f);

  /* check both vectors are equal */
  for (i = 0; i < n; i++)
    ASSERT_ALWAYS(bit_vector_getbit (b, i) == bit_vector_getbit (c, i));

  bit_vector_clear (b);
  bit_vector_clear (c);
}

int
main (int argc, const char *argv[])
{
  unsigned long iter = 20000;
  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter (&iter);
  test_bit_vector (iter);
  test_bit_vector_read_from_file ();
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}
