#include "cado.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>
#ifdef HAVE_SSSE3
#include "tmmintrin.h"
#endif
#ifdef HAVE_AVX2
#include "immintrin.h"
#endif
#ifdef HAVE_ARM_NEON
#include <arm_neon.h>
#endif
//#define DO_TIMING 1
#ifdef DO_TIMING
// The Jevents library is part of the PMU tools
// https://github.com/andikleen/pmu-tools
#ifdef HAVE_JEVENTS
#define USE_JEVENTS 1
#include "rdtsc.h"
#endif
#endif
#include "utils.h"
#include "tests_common.h"
#include "las-sieve2357.hpp"
#include "las-debug.hpp"

template<typename T>
class gettypename {
public:
    static constexpr const char * name = "unknown";
};

template<>
class gettypename<unsigned long> {
public:
    static constexpr const char * name = "unsigned long";
};

template<>
class gettypename<unsigned char> {
public:
    static constexpr const char * name = "unsigned char";
};

#ifdef HAVE_SSSE3
template<>
class gettypename<__m128i> {
public:
    static constexpr const char * name = "__m128i";
};
#endif

#ifdef HAVE_AVX2
template<>
class gettypename<__m256i> {
public:
    static constexpr const char * name = "__m256i";
};
#endif


template<typename T>
T * tolerant_malloc_aligned(size_t size)
{
    size_t s = sizeof(T);
    for( ; s % 8 ; s<<=1);
    return (T*) malloc_aligned(size, s);
}


template <typename SIMDTYPE, typename ELEMTYPE>
bool test(const unsigned long iter, const size_t arraysize,
  const sieve2357base::prime_t *all_primes, const size_t nr_all_primes)
{
  const size_t N = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
  bool ok = true;

  if(arraysize % N != 0) {abort();}
#define ARRAY_ON_HEAP 1
#if ARRAY_ON_HEAP
  SIMDTYPE *sievearray = tolerant_malloc_aligned<SIMDTYPE>(arraysize);
  ELEMTYPE *sievearray2 = tolerant_malloc_aligned<ELEMTYPE>(arraysize);
#else
#ifdef HAVE_ALIGNAS
  alignas(sizeof(SIMDTYPE))
#endif
  SIMDTYPE sievearray[arraysize / N];
  ELEMTYPE sievearray2[arraysize];
#endif
  sieve2357base::prime_t use_primes[nr_all_primes + 1]; /* +1 for end marker */

  /* Skip those prime (powers) this SIMD type can't sieve */
  size_t j = 0;
  for (size_t i = 0; i < nr_all_primes; i++) {
    if (sieve2357<SIMDTYPE, ELEMTYPE>::can_sieve(all_primes[i].q)) {
      use_primes[j++] = all_primes[i];
    }
  }
  use_primes[j] = sieve2357base::prime_t{0,0,0};

  where_am_I w;

#ifdef DO_TIMING
  sieve2357<SIMDTYPE, ELEMTYPE>:sieve(sievearray, arraysize, use_primes, false, sieve2357base::update_set, w);
  sieve2357<SIMDTYPE, ELEMTYPE>::sieve(sievearray, arraysize, use_primes, false, sieve2357base::update_set, w);
  start_timing();
#endif
  for (unsigned long i = 0; i < iter; i++) {
    sieve2357<SIMDTYPE, ELEMTYPE>::sieve(sievearray, arraysize, use_primes, false, sieve2357base::update_set, w);
  }
#ifdef DO_TIMING
  end_timing();
  if (iter > 1) {
    printf("%lu calls of sieve2357::sieve<%s, %s>(%zu) took %lu cycles, %lu per call, %.2f per %s\n",
           iter, gettypename<SIMDTYPE>::name, gettypename<ELEMTYPE>::name, arraysize,
           (unsigned long) get_diff_timing(), (unsigned long) get_diff_timing() / iter,
           (float)get_diff_timing() / iter / (arraysize / sizeof(SIMDTYPE)),
           gettypename<SIMDTYPE>::name);
  }
#endif

  /* Do the same thing again, using simple sieve code */
  memset(sievearray2, 0, arraysize);
  for (sieve2357base::prime_t *p = use_primes; p->q != 0; p++) {
    for (size_t i = p->idx; i < arraysize; i += p->q) {
      sievearray2[i] += p->logp;
    }
  }

  /* Compare the two arrays */
  unsigned long nr_errors = 0;
  for (size_t i = 0; i < arraysize && nr_errors < 10; i++) {
    if (((ELEMTYPE *)&sievearray[0])[i] != (sievearray2)[i]) {
      printf("Mismatch, sievearray<%s,%s>[%zu] = %hhu, sievearray2[%zu] = %hhu\n",
             gettypename<SIMDTYPE>::name, gettypename<ELEMTYPE>::name, i,
             ((ELEMTYPE *)&sievearray[0])[i], i, sievearray2[i]);
      nr_errors++;
    }
  }
  ok &= (nr_errors == 0);
#if ARRAY_ON_HEAP
  free(sievearray);
  free(sievearray2);
#endif
  return ok;
}

bool test_all(const unsigned long iter, const size_t arraysize,
  const sieve2357base::prime_t *all_primes, const size_t nr_all_primes)
{
  bool ok = true;
  ok &= test<uint32_t, unsigned char>(iter, arraysize, all_primes, nr_all_primes);
  ok &= test<uint64_t, unsigned char>(iter, arraysize, all_primes, nr_all_primes);
#ifdef HAVE_SSSE3
  ok &= test<__m128i, unsigned char>(iter, arraysize, all_primes, nr_all_primes);
#endif
#ifdef HAVE_AVX2
  ok &= test<__m256i, unsigned char>(iter, arraysize, all_primes, nr_all_primes);
#endif
#ifdef HAVE_ARM_NEON
  ok &= test<uint8x16_t, unsigned char>(iter, arraysize, all_primes, nr_all_primes);
#endif
  return ok;
}

int main(int argc, const char **argv)
{
  unsigned long iter = 1;
  size_t arraysize = 65536;
  bool ok = true;
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);

#ifdef DO_TIMING
  init_timing();
#endif

  {
    sieve2357base::prime_t all_primes[] = {
      /* q idx logp */
      {1, 0, 2},
      {2, 1, 1},
      {4, 3, 2},
      {8, 0, 6},
      {16, 14, 12},
      {32, 21, 24},
      {3, 2, 3},
      {3, 2, 3},
      {5, 4, 4},
      {7, 5, 5},
    };

    ok &= test_all(iter, arraysize, &all_primes[0], sizeof(all_primes) / sizeof(all_primes[0]));
  }

  for (int k = 0; k <= 6; k++) {
    const fbprime_t p = 1 << k;
    for(fbroot_t i = 0; i < p; i++) {
      sieve2357base::prime_t cur_prime[1] = {{p, i, 1}};
      ok &= test_all(1, arraysize, cur_prime, 1);
    }
  }

  for (fbprime_t p = 3; p <= 7; p += 2) {
    for(fbroot_t i = 0; i < p; i++) {
      sieve2357base::prime_t cur_prime[1] = {{p, i, 1}};
      ok &= test_all(1, arraysize, cur_prime, 1);
    }
  }

#ifdef DO_TIMING
  clear_timing();
#endif
  tests_common_clear();
  exit (ok ? EXIT_SUCCESS : EXIT_FAILURE);
}
