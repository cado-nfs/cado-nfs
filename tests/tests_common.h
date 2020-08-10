#ifndef TESTS_COMMON_H
#define TESTS_COMMON_H

#include <stdint.h>
#include <gmp.h>

#define PARSE_SEED 1
#define PARSE_ITER 2
#define PARSE_VERBOSE 4
#define PARSE_QUIET 8

extern gmp_randstate_t state;

#ifdef __cplusplus
extern "C" {
#endif

int cmp_double(double, double, double);
int64_t random_int64 ();
uint64_t random_uint64 ();
/** Generate a uniformly distributed random integer in the range 0 to
 *  2^N-1, inclusive. (Copied from GMP info page) */
void tests_common_urandomb (mpz_t R, mp_bitcnt_t N);
/** Generate a uniform random integer in the range 0 to N-1,
 *  inclusive. (Copied from GMP info page) */
void tests_common_urandomm (mpz_t R, mpz_t N);
/** Generate a random integer with long strings of zeros and ones in
 *  the binary representation. (Copied from GMP info page) */
void tests_common_rrandomb (mpz_t R, mp_bitcnt_t N);
/** If the "-iter" command line parameter was parsed by tests_common_cmdline(),
 *  then write the iter value to *output, otherwise do nothing. */
void tests_common_get_iter(unsigned long *output);
int tests_common_get_verbose();
int tests_common_get_quiet();
void tests_common_cmdline(int *, const char ***, uint64_t);
void tests_common_clear();

#ifdef __cplusplus
}
#endif

#endif
