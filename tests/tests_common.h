#ifndef TESTS_COMMON_H
#define TESTS_COMMON_H

#include <stdint.h>
#include <gmp.h>

#define PARSE_SEED 1
#define PARSE_ITER 2
#define PARSE_VERBOSE 4
#define PARSE_QUIET 8
#define PARSE_CHECK 16
#define PARSE_TIME 32

extern gmp_randstate_t state;

#ifdef __cplusplus
extern "C" {
#endif

int cmp_double(double, double, double);
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
/** If the "-v" command line parameter was parsed by tests_common_cmdline(),
 * then return non-zero, otherwise return zero. */
int tests_common_get_verbose();
/** If the "-q" command line parameter was parsed by tests_common_cmdline(),
 * then return non-zero, otherwise return zero. */
int tests_common_get_quiet();
/** If the "-check" command line parameter was parsed by tests_common_cmdline(),
 * then return non-zero, otherwise return zero. */
int tests_common_get_check();
/** If the "-time" command line parameter was parsed by tests_common_cmdline(),
 * then return non-zero, otherwise return zero. */
int tests_common_get_time();
/* If neither "-check" nor "-time" command line parameters were given, this
   function does nothing. Otherwise, do_check and do_time are set according
   to whether "-check" and "-time" was given, respectively. */
void tests_common_get_check_and_time(int *do_check, int *do_time);
void tests_common_cmdline(int *, const char ***, uint64_t);
void tests_common_clear();

#ifdef __cplusplus
}
#endif

#endif
