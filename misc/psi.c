#include "cado.h"       // IWYU pragma: keep

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include "getprime.h"

/*
  Results
  -k 0 100000000000000000 1000 10000: 7991165172719
  -k 0 1000000000000000000 1000 10000: 30538908889103
  -k 0 10000000000000000000 1000 10000: 113981832007394
  -k 0 10000000000000000 10000 10000: 57822909092253
*/

uint64_t Psi(uint64_t x, uint64_t y);

uint64_t *primes;
size_t nr_primes;

static inline uint64_t
min_u64(uint64_t a, uint64_t b)
{
  return (a < b) ? a : b;
}

static inline uint64_t
max_u64(uint64_t a, uint64_t b)
{
  return (a > b) ? a : b;
}

/* Returns floor(log_2(x)). Requires x > 0.
   Assumes unsigned long long has 64 bits. */
static uint64_t
intlog2(const uint64_t x)
{
  return 63 - __builtin_clzll(x);
}

/* A look-up table for Psi(x, y) with small x and y.
   It serves two purposes: speeding up the recursive computation of Psi(x, y),
   and testing that Psi(x, y) is correct for small x, y.
   The LUT is computed by factoring each x to determine its larges prime factor.
   Each entry of the LUT so obtained is compared to the value produced by the
   Psi(x, y) function (which for this purposed is computed without referencing
   the LUT) */

/* To use LUT_X > 256 or LUT_Y > 256, the data type of a table entry must
   be changed to short in order to avoid overflow  */
#define LUT_X 256
#define LUT_Y 128
static const int use_LUT = 1; /* enable or disable LUT for speed comparison */
static int in_check_LUT = 0; /* disable use of LUT in Psi() while checking correctness */
static unsigned short Psi_LUT[LUT_X][LUT_Y];

/* Return the largest prime factor of x.
   For x == 0 or x == 1, returns 0.
   Does not need to be fast, and isn't. */
static inline uint64_t find_lpf(const uint64_t x)
{
  uint64_t p, c = x, lpf = 0;
  for (p = 2; c > 1; p++) {
    while(c % p == 0) {
      c /= p;
      lpf = p;
    }
  }
  /* printf("find_lpf(%" PRIu64 ") = %" PRIu64 "\n", x, lpf);
  fflush(stdout); */
  return lpf;
}

static void build_LUT()
{
  uint64_t lpf[LUT_X];
  for (size_t x = 0; x < LUT_X; x++)
    lpf[x] = find_lpf(x);

  for (size_t y = 0; y < LUT_Y; y++) {
    Psi_LUT[0][y] = 0; /* There are 0 integers in [1, 0] */
    for (size_t x = 1; x < LUT_X; x++) {
      Psi_LUT[x][y] = Psi_LUT[x - 1][y] + (lpf[x] <= y);
    }
  }
}

static void check_LUT()
{
  in_check_LUT = 1;
  for (size_t y = 0; y < LUT_Y; y++) {
    for (size_t x = 1; x < LUT_X; x++) {
      if (Psi_LUT[x][y] != Psi(x, y)) {
        fprintf (stderr,
                 "Error, Psi_LUT[%zu][%zu] = %u, but Psi(x, y) = %" PRIu64 "\n",
                 x, y, (unsigned int) Psi_LUT[x][y], Psi(x, y));
        abort();
      }
    }
  }
  in_check_LUT = 0;
}


/* Count integers in [1, x] that are y-smooth.
   An integer is y-smooth if no prime > y divides it. */
uint64_t
Psi(const uint64_t x, const uint64_t y)
{
  uint64_t r;
  if (x == 0) {
    r = 0; /* [1, 0] is the empty interval */
  } else if (x <= y) {
    r = x; /* All integers in [1, x] are smooth */
  } else if (y <= 1) {
    r = 1; /* Only 1 is smooth */
  } else if (use_LUT && !in_check_LUT && x < LUT_X && y < LUT_Y) {
    r = Psi_LUT[x][y];
  } else {
    r = 1;
    size_t i = 0;

    if (y >= 2) {
      /* Handle p = primes[0] = 2 without recursion */
      r += intlog2(x);
      i++;
    }

    if (y >= 3) {
      /* Handle p = primes[1] = 3 without recursion */
      for (uint64_t t = x / 3; t > 0; t /= 3)
        r += intlog2(t) + 1;
      i++;
    }

    const uint64_t m_max = 1;
    uint64_t l_max, m;

    if (m_max == 1) {
      m = 1;
      l_max = y;
    } else {
      m = min_u64(m_max, y);

      /* At first, we process all p such that floor(x/p) >= m
         which is true iff x/p >= m */

      l_max = min_u64(x / m, y);
    }

    /* We partition the set of y-smooth integers up to x according to their
       largest prime factor, p_1. For each such p_1 <= y, we count the
       elements in it partition, which is Psi(x / p_1, p_1) */
    for ( ; primes[i] <= l_max; i++) {
      /* This integer division is where we spend most of the time */
      const uint64_t t = x / primes[i];
      ASSERT(t >= m);
      r += Psi(t, primes[i]);
    }

    for (m-- ; m > 0; m--) {
      /* Count Psi(x/p, p) where floor(x/p) = m.
         Could use a primepi() function here that does bisect search,
         but this part of the code is currently not the bottleneck */
      l_max = min_u64(x / m, y);
      for ( ; primes[i] <= l_max; i++) {
        ASSERT (x / primes[i] == m);
        r += m; /* With y >= l_max, Psi(m, y) = m */
      }
    }
  }

  if (0 && !in_check_LUT)
    printf("Psi(%" PRIu64 ", %" PRIu64 ") = %" PRIu64 "\n", x, y, r);
  return r;
}

/* Count integers in [1, x] that have exactly k prime factors
   (counting multiplicity) in ]y, z], and all other prime factors <= y. */

uint64_t
Psi_k(const uint64_t x, const uint64_t y, const uint64_t z,
      const unsigned long k)
{
  uint64_t s = 0;
  size_t i;

  if (k == 0) {
    return Psi(x, y);
  }

  /* Find first prime > y. TODO: bisect */
  for (i = 0; primes[i] <= y; i++);

  if (k == 1) {
    for ( ; primes[i] <= z; i++)
      s += Psi(x / primes[i], y);
  } else {
    for ( ; primes[i] <= z; i++)
      s += Psi_k(x / primes[i], y, primes[i], k-1);
  }

  return s;
}

static uint64_t factors[256];

static void
print_factors(const size_t nr_factors)
{
  if (nr_factors == 0)
    printf("1");
  for (size_t i = 0; i < nr_factors; i++)
    printf("%s%" PRIu64, (i == 0) ? "" : " ", factors[i]);
  printf ("\n");
}

/* Counts smooth integers and also prints each one as a list of primes */

uint64_t
Psi_print(uint64_t x, uint64_t y, const size_t nr_factors)
{
  if (x == 0) return 0;
  if (y == 1) {
    print_factors(nr_factors);
    return 1;
  };

  uint64_t r = Psi_print(x, 1, nr_factors);

  for (size_t i = 0; primes[i] <= y; i++) {
    factors[nr_factors] = primes[i];
    r += Psi_print(x / primes[i], primes[i], nr_factors + 1);
  }

  return r;
}

uint64_t
Psi_k_print(const uint64_t x, const uint64_t y, const uint64_t z,
            const unsigned long k, const size_t nr_factors)
{
  uint64_t s = 0;
  size_t i;

  if (k == 0) {
    return Psi_print(x, y, nr_factors);
  }

  /* Find first prime > y. Could use bisect */
  for (i = 0; primes[i] <= y; i++);

  for ( ; primes[i] <= z; i++) {
    factors[nr_factors] = primes[i];
    s += Psi_k_print(x / primes[i], y, primes[i], k-1, nr_factors + 1);
  }

  return s;
}



int main(int argc, char * argv[])  // not 'char const' because of strtok
{
  uint64_t x, y, z, result = 0, sum = 0;
  int print = 0;
  char k_default[4] = {'0',',','1',0}; /* Need this non-const, thus no literal string */
  char * k_values = k_default;

  while (argc > 1) {
    if (strcmp(argv[1], "-k") == 0) {
      k_values = argv[2];
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-print") == 0) {
      print = 1;
      argc--;
      argv++;
    } else
      break;
  }

  if (argc < 4){
    printf("Usage: psi [-k 0,1,2,5] [-print] x y z\n");
    printf("Computes Psi_k(x, y, z), the number of intgers in [1, x] with exactly k prime\n"
           "factors in ]y, z], and all others in [1, y]\n");
    printf("Default for k is 0,1; i.e., integers with at most one large prime.\n");
    printf("The -k parameter takes a comma-separated list of k-values to process.\n");
    printf("With -print, prints each smooth integer as a list of primes\n");
    exit(EXIT_FAILURE);
  }
  
  x = strtoull(argv[1], NULL, 10);
  y = strtoull(argv[2], NULL, 10);
  z = strtoull(argv[3], NULL, 10);

  /* With small z, we need to init the primes table up to LUT_Y,
     or build_LUT() crashes */
  uint64_t prime_table_max = max_u64(z, LUT_Y) + 1;
  size_t nr_primes = 1; /* 1 extra at the end */
  prime_info pi;
  prime_info_init (pi);
  for (uint64_t p = 2; p <= prime_table_max; p = getprime_mt (pi))
    nr_primes++;
  
  primes = malloc(nr_primes * sizeof(uint64_t));
  prime_info_clear (pi);
  primes[0] = 2;
  prime_info_init (pi);
  for (size_t i = 1; i < nr_primes; i++)
    primes[i] = getprime_mt (pi);
  prime_info_clear (pi);

  build_LUT();
  check_LUT();

  char * saveptr = NULL;

  for (char * k_str = strtok_r(k_values, ",", &saveptr); k_str != NULL; k_str = strtok_r(NULL, ",", &saveptr)) {
    uint64_t k = strtoul(k_str, NULL, 10);
    if (print) {
      result = Psi_k_print(x, y, z, k, 0);
    } else {
      result = Psi_k(x, y, z, k);
    }
    sum += result;
    printf("There are %" PRIu64 " integers in [1, %" PRIu64 "] with exactly %" PRIu64
           " prime factor%s in ]%" PRIu64 ", %" PRIu64"] and all other prime factors <= %" PRIu64 "\n",
           result, x, k, k == 1 ? "" : "s", y, z, y);
  }
  if (sum != result) {
    printf("Sum over all processed k-values: %" PRIu64 "\n", sum);
  }
  free(primes);
  return 0;
}
