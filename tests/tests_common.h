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
void tests_common_get_iter(unsigned long *);
int tests_common_get_verbose();
int tests_common_get_quiet();
void tests_common_cmdline(int *, const char ***, uint64_t);
void tests_common_clear();

#ifdef __cplusplus
}
#endif

#endif
