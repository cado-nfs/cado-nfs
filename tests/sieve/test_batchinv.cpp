/* test batchinvredc_u32 */

#include "cado.h" // IWYU pragma: keep
#include <cinttypes>               // for PRId64, PRIu64
#include <cstdint>                 // for uint32_t, uint64_t
#include <cstdio>
#include "fb-types.h"               // for fbprime_t, FBPRIME_FORMAT
#include "las-arith.hpp"            // for batchinvredc_u32
#include "macros.h"
#include "gcd.h"
#include "tests_common.h"

static inline uint32_t
mulmod(const uint32_t a, const uint32_t b, const uint32_t p) {
    return ((uint64_t) a * (uint64_t) b) % p;
}

bool
test_batchinvredc_u32(const uint32_t *a, const uint32_t p, const uint32_t invp, const size_t n, const bool should_work)
{
    const uint32_t two32 = UINT32_MAX % p + 1; /* 2^32 % p */
    uint32_t r[n];
    bool worked;
    
    if (redc_no_carry(p))
        worked = batchinvredc_u32<false> (r, a, n, p, invp);
    else
        worked = batchinvredc_u32<true> (r, a, n, p, invp);

    if (worked != should_work) {
        fprintf (stderr, "Inverse should %shave worked but %s\n",
                 should_work ? "" : "not ", worked ? "did" : "didn't");
        return false;
    } else if (worked == true) {
        for (size_t i = 0; i < n; i++) {
            if (mulmod(a[i], r[i], p) != two32) {
                fprintf (stderr, "2^32 * %" PRIu32 "^-1 (mod %" PRIu32 ") wrong: %" PRIu32 "\n",
                    a[i], p, r[i]);
                return false;
            }
        }
        return true;
    } else {
        /* Should not work and didn't: all is well */
        return true;
    }
}

/* The usual tests command line parameters "-seed" and "-iter" are accepted. */

int
main (int argc, const char *argv[])
{
  unsigned long N = 100;
  bool all_ok = true;

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&N);
  
  for (size_t n = 0; all_ok && n < N; n++) {
      const uint32_t p = (uint64_t) u64_random (state) | 1;
      const uint32_t invp = -invmod_po2 (p);
      bool should_work = true;
      
      uint32_t a[n];
      /* Try one that may or may not work, depending on which residues RNG
       * gives us */
      for (size_t i = 0; i < n; i++) {
          a[i] = (uint64_t) u64_random(state);
          should_work &= gcd_ul(a[i], p) == 1;
      }
      
      all_ok &= test_batchinvredc_u32 (a, p, invp, n, should_work);
      
      /* Make all residues coprime to p so that we can test one that should
       * indeed work */
      for (size_t i = 0; i < n; i++) {
          while (gcd_ul(a[i], p) != 1) {
              a[i]++;
          }
      }

      all_ok &= test_batchinvredc_u32 (a, p, invp, n, true);
  }

  exit (all_ok ? EXIT_SUCCESS : EXIT_FAILURE);
}
