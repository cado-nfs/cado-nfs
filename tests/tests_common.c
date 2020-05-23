#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include "tests_common.h"
#include "macros.h"
#include "portability.h"

static int rng_state_inited = 0;
gmp_randstate_t state;
static int parsed_iter = 0;
unsigned long iter;
static int verbose = 0, quiet = 0;

/* Return non-zero iff |d2| is in the interval |d1| * (1 +- err_margin) */
int
cmp_double(const double d1, const double d2, const double err_margin)
{
  return fabs(d1) * (1. - err_margin) <= fabs(d2) && fabs(d2) <= fabs(d1) * (1. + err_margin);
}

uint64_t
random_uint64 ()
{
#ifdef HAVE_LRAND48
  return ((uint64_t) lrand48 () << 33) + ((uint64_t) lrand48 () << 2) + ((uint64_t) lrand48 () & 3);
#else
  ASSERT_ALWAYS(RAND_MAX < UINT64_MAX - 1);
  uint64_t r = rand(), q = UINT64_MAX, b = (uint64_t) RAND_MAX + 1;
  /* +1 because range of rand() is [0, RAND_MAX] */
  /* We want ceil(log_b(UINT64_MAX + 1)) iterations
     = ceil(log_b(UINT64_MAX)) unless log_b(UINT64_MAX) is an integer, which it
     is only if b = UINT64_MAX, but the return type of rand() is int. */
  while (q > b) {
    r = r * b + (uint64_t) rand();
    q = (q - 1) / b + 1;
  }
  return r;
#endif
}


int64_t
random_int64 ()
{
#ifdef HAVE_LRAND48
  return ((int64_t) mrand48 () << 32) + (int64_t) mrand48 ();
#else
  return random_uint64() + INT64_MIN;
#endif
}

/* Set *output to the value from -iter, if -iter was given,
   otherwise do nothing */
void
tests_common_get_iter(unsigned long *output)
{
  if (parsed_iter)
    *output = iter;
}

int
tests_common_get_verbose()
{
  return verbose;
}

int
tests_common_get_quiet()
{
  return quiet;
}

static int
parse_l(long *result, const char *s)
{
  char *endptr;
  long r = strtol(s, &endptr, 10);
  if (s[0] == '\0' || endptr[0] != '\0' || errno == ERANGE)
    return 0;
  *result = r;
  return 1;
}

static int
parse_ul(unsigned long *result, const char *s)
{
  char *endptr;
  unsigned long r = strtoul(s, &endptr, 10);
  if (s[0] == '\0' || endptr[0] != '\0' || errno == ERANGE)
    return 0;
  *result = r;
  return 1;
}

void
tests_common_cmdline(int *argc, const char ***argv, const uint64_t flags)
{
  long seed = 0;

  if ((flags & PARSE_SEED) != 0)
    seed = time (NULL);

  while (1) {
    const char *name = "-seed";
    if ((flags & PARSE_SEED) != 0 && (*argc) > 1 && 
        strcmp(name, (*argv)[1]) == 0) {
      if ((*argc) <= 2) {
        fprintf (stderr, "No value given for %s parameter\n", name);
        exit(EXIT_FAILURE);
      }
      if (!parse_l(&seed, (*argv)[2])) {
        fprintf (stderr, "Invalid value \"%s\" given for %s parameter\n",
                 (*argv)[2], name);
        exit(EXIT_FAILURE);
      }
      *argc -= 2;
      *argv += 2;
      continue;
    } 

    name = "-iter";
    if ((flags & PARSE_ITER) != 0 && (*argc) > 1 && 
        strcmp(name, (*argv)[1]) == 0) {
      if ((*argc) <= 2) {
        fprintf (stderr, "No value given for %s parameter\n", name);
        exit(EXIT_FAILURE);
      }
      if (!parse_ul(&iter, (*argv)[2])) {
        fprintf (stderr, "Invalid value \"%s\" given for %s parameter\n",
                 (*argv)[2], name);
        exit(EXIT_FAILURE);
      }
      *argc -= 2;
      *argv += 2;
      parsed_iter = 1;
      printf ("Using %lu iterations\n", iter);
      continue;
    }

    name = "-v";
    if ((flags & PARSE_VERBOSE) != 0 && (*argc) > 1 && 
        strcmp(name, (*argv)[1]) == 0) {
      verbose = 1;
      printf ("Using verbose output\n");
      *argc -= 1;
      *argv += 1;
      continue;
    }

    name = "-q";
    if ((flags & PARSE_QUIET) != 0 && (*argc) > 1 && 
        strcmp(name, (*argv)[1]) == 0) {
      quiet = 1;
      *argc -= 1;
      *argv += 1;
      continue;
    }
    break;
  }

  if ((flags & PARSE_SEED) != 0) {
    printf ("Using random seed=%ld\n", seed);
#ifdef HAVE_LRAND48
    srand48 (seed);
#else
    unsigned int s = labs(seed);
    if (seed < 0 || s != (unsigned long) seed)
      printf ("Warning, seed truncated to %u\n", s);
    srand ((unsigned int) labs(seed));
#endif
    fflush (stdout);
    gmp_randinit_default (state);
    gmp_randseed_ui (state, (unsigned long) labs(seed));
    rng_state_inited = 1;
  }
}

/* Clean up, free any memory that may have been allocated */
void
tests_common_clear()
{
  if (rng_state_inited) {
    gmp_randclear (state);
    rng_state_inited = 0;
  }
}
