/* check_rels: check the factorization of the norm on both alg and rat side.
   Can also, in option, check the primality of ideal and correct uncomplete
   relation or relation with non prime ideal */

#include "cado.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <limits.h> /* for CHAR_BIT */
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>   /* for _O_BINARY */

#include "portability.h"
#include "macros.h"
#include "utils_with_io.h"

#define FACTOR_DO_NOT_DIVIDE 1UL
#define FACTOR_NOT_PRIME 2UL
#define FACTORIZATION_NOT_COMPLETE 4UL
#define FACTOR_ABOVE_LPB 8UL
#define REL_FULLY_FIXED 16UL

uint64_t nrels_read = 0, nrels_ok = 0, nrels_err = 0, nrels_fullyfixed = 0,
         nrels_donotdivide = 0, nrels_notprime = 0, nrels_notcomplete = 0,
         nrels_abovelpb = 0, nrels_fullycompleted = 0, nrels_fixednotprime = 0;
cado_poly cpoly;
unsigned long lpb[NB_POLYS_MAX] = {0, 0, 0, 0, 0, 0, 0, 0};
unsigned long lpb_max = 0;
int verbose = 0;
int abhexa = 0;
int check_primality = 0; /* By default no primality check */
int fix_it = 0; /* By default, we just check the rels */

/* used for counting time in different processes */
timingstats_dict_t stats;

void
rel_add_prime (earlyparsed_relation_ptr rel, unsigned int side, p_r_values_t p,
               exponent_t e)
{
  for(weight_t i = 0; i < rel->nb ; i++)
  {
    if (rel->primes[i].p == p && rel->primes[i].side == side)
    {
      rel->primes[i].e += e;
      return;
    }
  }

  if (rel->nb == rel->nb_alloc)
    realloc_buffer_primes_c (rel);
  rel->primes[rel->nb] = (prime_t) {.h = 0, .p = p, .e = e, .side = side};
  rel->nb++;
}

static inline void
print_error_line (prime_t prime, mpz_t norm[],
                  unsigned long err_type, int will_be_fixed)
{
  unsigned int side = prime.side;
  p_r_values_t p = prime.p;
  exponent_t e = prime.e;
  char *str = (will_be_fixed) ? "Warning" : "Error";

  if (err_type == FACTOR_DO_NOT_DIVIDE)
  {
    fprintf (stderr, "#   Error, given factor %" PRpr " with exponent %u does "
                     "not divide the norm on side %u\n", p, e, side);
  }
  else if (err_type == FACTOR_NOT_PRIME)
  {
    fprintf (stderr, "#   %s, given factor %" PRpr " is not prime on side "
                     "%u\n", str, p, side);
  }
  else if (err_type == FACTOR_ABOVE_LPB)
  {
    fprintf (stderr, "#   %s, given factor %" PRpr " is greater than lpb "
                     "(=%lu) on side %u\n", str, p, lpb[side], side);
  }
  else if (err_type == FACTORIZATION_NOT_COMPLETE)
  {
    gmp_fprintf (stderr, "#   %s, factorization of the norm on side %u is not "
                         "complete (%Zu is missing)\n", str, side, norm[side]);
  }
}

static inline void
factor_nonprime_ideal (earlyparsed_relation_ptr rel, weight_t i)
{
  exponent_t e = rel->primes[i].e;
  p_r_values_t p = rel->primes[i].p;
  unsigned int side = rel->primes[i].side, pr = 2;
  modulusul_t m;
  prime_info pi;
  prime_info_init (pi);
  do
  {
    exponent_t e_pr_in_p = 0;
    while (p % pr == 0)
    {
      p = p / pr;
      e_pr_in_p++;
    }

    if (p == 1)
    {
      rel->primes[i] = (prime_t) {.h= 0, .p= pr, .e= e * e_pr_in_p, .side = side};
      break;
    }
    else if (e_pr_in_p != 0)
      rel_add_prime (rel, side, pr, e * e_pr_in_p);

    pr = getprime_mt (pi);
    modul_initmod_ul (m, p);
  } while (!modul_isprime(m));
  prime_info_clear (pi);
  modul_clearmod (m);
  if (p != 1) /* means remaining p is prime */
    rel->primes[i] = (prime_t) {.h = 0, .p = p, .e = e, .side = side};
}

static int
more_job (mpz_t norm[], unsigned int nb_poly)
{
  for(unsigned int side = 0; side < nb_poly; side++)
    if(mpz_cmp_ui (norm[side], 1) != 0)
      return 1;
  return 0;
}

/* return 0 if everything is ok (factorization, primality, and complete)
 * else the ith bit of the return value is set if
 *  i = 1: a factor does not divide the norms (the check is then stopped)
 *  i = 2: a factor is not prime
 *  i = 3: the factorization of a norm is not complete
 *  i = 4: a factor is above a lpb
 *  i = 5: (only with fix_it != 0) non-prime factors were successfully factored
/bin/bash: q: command not found
 */
unsigned long
process_one_relation (earlyparsed_relation_ptr rel)
{
  mpz_t norm[NB_POLYS_MAX];
  char used[NB_POLYS_MAX];
  memset(used, 0, NB_POLYS_MAX); /* which sides are in use */
  for(int side = 0 ; side < cpoly->nb_polys ; side++)
    mpz_init (norm[side]);
  unsigned long err = 0;

  if (verbose)
    fprintf (stderr, "# relation %" PRIu64 " with (a,b) = (%" PRId64 ","
                     "%" PRIu64 "):\n", rel->num, rel->a, rel->b);

  /* Look for which sides are in use */
  for(weight_t i = 0; i < rel->nb ; i++)
    used[rel->primes[i].side] = 1;

  /* compute the norm on alg and rat sides */
  for(int side = 0 ; side < cpoly->nb_polys ; side++)
  {
    if (used[side])
    {
      mpz_poly_ptr ps = cpoly->pols[side];
      mpz_poly_homogeneous_eval_siui (norm[side], ps, rel->a, rel->b);
      if (verbose)
        gmp_fprintf (stderr, "#   norm on side %d = %Zu\n", side, norm[side]);
    }
    else
      mpz_set_ui (norm[side], 1);
  }

  /* check for correctness of the factorization of the norms */
  for(weight_t i = 0; i < rel->nb ; i++)
  {
    unsigned int side = rel->primes[i].side;
    p_r_values_t p = rel->primes[i].p;
    exponent_t e = rel->primes[i].e;
    ASSERT_ALWAYS(p != 0); /* could reveal a problem in parsing */
    ASSERT_ALWAYS(e > 0); /* non positive exponent is not possible */
    for (int j = 0; j < e; j++)
    {
      if (!mpz_divisible_ui_p (norm[side], p))
      {
        err |= FACTOR_DO_NOT_DIVIDE;
        if (verbose != 0)
          print_error_line (rel->primes[i], norm, FACTOR_DO_NOT_DIVIDE, fix_it);
      }
      else
        mpz_divexact_ui (norm[side], norm[side], p);
    }
  }

  /* With an error in the factorization of the norm, no need to continue */
  if (!err)
  {
    /* check primality of all ideals appearing in the relations */
    if (check_primality != 0)
    {
      weight_t len = rel->nb; /* no need to check the prime that can be added */
      for(weight_t i = 0; i < len ; i++)
      {
        modulusul_t p;
        modul_initmod_ul (p, rel->primes[i].p);
        if (!modul_isprime(p))
        {
          err |= FACTOR_NOT_PRIME;
          if (verbose != 0)
            print_error_line (rel->primes[i], norm, FACTOR_NOT_PRIME, fix_it);
          if (fix_it) /* if fix_it = 1, we factor it */
            factor_nonprime_ideal(rel, i);
        }
        modul_clearmod (p);
      }
    }

    /* Check that the product of the factors is equal to the norm. */
    for(int side = 0 ; side < cpoly->nb_polys ; side++)
    {
      if (mpz_cmp_ui (norm[side], 1) != 0)
      {
        err |= FACTORIZATION_NOT_COMPLETE;
        if (verbose != 0)
        {
          prime_t fake = {.h = 0, .p = 0, .e = 0, .side = side};
          print_error_line (fake, norm, FACTORIZATION_NOT_COMPLETE, fix_it);
        }
      }
    }

    /* complete relation if it is asked and necessary */
    if (fix_it)
    {
      /* complete at least for primes up to 10000 (because of GGNFS and Msieve
       * that skip these primes) */
      unsigned long max_p = MAX (lpb_max, 10000);
      prime_info pi;
      prime_info_init (pi);
      for (unsigned long p = 2; more_job (norm, cpoly->nb_polys) && p < max_p ;
           p = getprime_mt (pi))
      {
        for(int side = 0 ; side < cpoly->nb_polys ; side++)
        {
          exponent_t e = 0;
          while (mpz_divisible_ui_p (norm[side], p))
          {
            e++;
            mpz_divexact_ui (norm[side], norm[side], p);
          }
          if (e != 0)
          {
            rel_add_prime (rel, side, p, e);
            if (verbose != 0)
            {
              gmp_fprintf (stderr, "#   factorization of the norm on side %u is "
                                   "not complete, need to add %" PRpr "^%u\n",
                                   side, p, e);
            }
          }
        }
      }
      prime_info_clear (pi);

      if (more_job (norm, cpoly->nb_polys))
        gmp_fprintf (stderr, "#   factorization of the norm is still not "
                             "complete on at least one side\n");
    }

    /* check that ideals appearing in the relations are below the lpb. We
     * skipped this check if all the primes dividing the norms are not known
     * (because we do not know if the missing primes are below or above the
     * lpbs). */
    if (!more_job (norm, cpoly->nb_polys))
    {
      for(weight_t i = 0; i < rel->nb ; i++)
      {
        p_r_values_t p = rel->primes[i].p;
        unsigned int side = rel->primes[i].side;
        if (p > lpb[side])
        {
          err |= FACTOR_ABOVE_LPB;
          if (verbose != 0)
            print_error_line (rel->primes[i], norm, FACTOR_ABOVE_LPB, fix_it);
        }
      }
    }

    if (fix_it && (err & FACTOR_NOT_PRIME || err & FACTORIZATION_NOT_COMPLETE))
      if (!more_job (norm, cpoly->nb_polys))
        err |= REL_FULLY_FIXED;
  }

  for(int side = 0 ; side < cpoly->nb_polys ; side++)
    mpz_clear(norm[side]);

  if (verbose)
    fprintf (stderr, "#   end with error status %lx\n", err);
  return err;
}

static inline void
print_relation (FILE *outfile, earlyparsed_relation_ptr rel)
{
  char buf[1 << 12], *p, *op;
  size_t t;
  unsigned int i, j;

  if (!abhexa)
  {
    p = d64toa10(buf, rel->a);
    *p++ = ',';
    p = u64toa10(p, rel->b);
    *p++ = ':';
  }
  else
  {
    p = d64toa16(buf, rel->a);
    *p++ = ',';
    p = u64toa16(p, rel->b);
    *p++ = ':';
  }

  for(unsigned int side = 0 ; side < 2 ; side++)
  {
    *(--p) = ':';
    p++;
    for (i = 0; i < rel->nb; i++)
    {
      if (rel->primes[i].side == side && rel->primes[i].e > 0)
      {
        op = p;
        p = u64toa16(p, (uint64_t) rel->primes[i].p);
        *p++ = ',';
        t = p - op;
        for (j = (unsigned int) ((rel->primes[i].e) - 1); j--; p += t)
          memcpy(p, op, t);
      }
    }
  }
  *(--p) = '\n';
  p[1] = 0;
  fputs(buf, outfile);
}

/* Callback function called by filter_rels */

static void *
thread_callback (void * context_data, earlyparsed_relation_ptr rel)
{
  FILE *outfile = (FILE *) context_data;

  nrels_read++;
  int ret = process_one_relation (rel);
  int is_printable;

  if (ret == 0)
  {
    nrels_ok++;
    is_printable = 1;
  }
  else if ((ret & REL_FULLY_FIXED) && !(ret & FACTOR_ABOVE_LPB))
  {
    nrels_fullyfixed++;
    is_printable = 1;
  }
  else
  {
    nrels_err++;
    is_printable = 0;
  }

  if (ret & FACTOR_DO_NOT_DIVIDE)
    nrels_donotdivide++;
  if (ret & FACTOR_NOT_PRIME)
  {
    if (ret & REL_FULLY_FIXED && is_printable)
      nrels_fixednotprime++;
    else if (!(ret & REL_FULLY_FIXED))
      nrels_notprime++;
  }
  if (ret & FACTORIZATION_NOT_COMPLETE)
  {
    if (ret & REL_FULLY_FIXED && is_printable)
      nrels_fullycompleted++;
    else if (!(ret & REL_FULLY_FIXED))
      nrels_notcomplete++;
  }
  if (ret & FACTOR_ABOVE_LPB)
    nrels_abovelpb++;

  if (outfile && is_printable)
    print_relation (outfile, rel);

  return NULL;
}

static void
declare_usage (param_list pl)
{
  param_list_decl_usage(pl, "filelist", "file containing a list of input files");
  param_list_decl_usage(pl, "basepath", "path added to all file in filelist");
  param_list_decl_usage (pl, "poly", "polynomials file (mandatory)");
  param_list_decl_usage (pl, "abhexa", "(switch) read and write a and b as hexa "
                                        "(instead of decimal)");
  param_list_decl_usage (pl, "fixit",  "(switch) Try to fix wrong relations");
  param_list_decl_usage (pl, "check_primality", "(switch) check primality of "
                                                "primes (default, no checking)");
  param_list_decl_usage (pl, "out", "optional output file for correct or fixed "
                                    "relations");
  param_list_decl_usage (pl, "lpb0", "large prime bound on side 0");
  param_list_decl_usage (pl, "lpb1", "large prime bound on side 1");
  param_list_decl_usage (pl, "lpbs", "large primes bounds (comma-separated list) "
                                     "(for MNFS)");
  param_list_decl_usage (pl, "v", "(switch) more verbose output");
  param_list_decl_usage(pl, "force-posix-threads", "(switch)");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
}

static void
usage (param_list pl, char *argv0)
{
  param_list_print_usage (pl, argv0, stderr);
  exit(EXIT_FAILURE);
}

int
main (int argc, char * argv[])
{
    char * argv0 = argv[0];
    int lpb_arg[NB_POLYS_MAX] = { 0 };
    unsigned int lpb0_arg = 0;
    unsigned int lpb1_arg = 0;
    FILE *outfile = NULL;

    param_list pl;
    param_list_init(pl);
    declare_usage(pl);

    param_list_configure_switch(pl, "abhexa", &abhexa);
    param_list_configure_switch(pl, "fixit", &fix_it);
    param_list_configure_switch(pl, "v", &verbose);
    param_list_configure_switch(pl, "check_primality", &check_primality);
    param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    argv++, argc--;
    if (argc == 0)
      usage (pl, argv0);


    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
    }

    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    param_list_parse_uint (pl, "lpb0", &lpb0_arg);
    param_list_parse_uint (pl, "lpb1", &lpb1_arg);
    int narg = param_list_parse_int_list (pl, "lpbs", lpb_arg, NB_POLYS_MAX, ",");

    const char * polyfilename = param_list_lookup_string(pl, "poly");
    const char *outfilename = param_list_lookup_string(pl, "out");
    const char * filelist = param_list_lookup_string(pl, "filelist");
    const char * basepath = param_list_lookup_string(pl, "basepath");
    const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");
    set_antebuffer_path (argv0, path_antebuffer);

    if (param_list_warn_unused(pl))
      usage (pl, argv0);

    if (basepath && !filelist)
    {
      fprintf(stderr, "-basepath only valid with -filelist\n");
      usage (pl, argv0);
    }

    if ((filelist != NULL) + (argc != 0) != 1)
    {
      fprintf(stderr, "Provide either -filelist or freeform file names\n");
      usage (pl, argv0);
    }

    if (!polyfilename)
    {
      fprintf (stderr, "Error, missing -poly command line argument\n");
      usage (pl, argv0);
    }

    cado_poly_init (cpoly);
    if (!cado_poly_read (cpoly, polyfilename))
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }

    if (narg == 0) /* lpbs were given as -lpb0 and -lpb1 */
    {
      if (cpoly->nb_polys > 2) /* With more than 2 polys, must use -lpbs. */
      {
        fprintf (stderr, "Error, missing -lpbs command line argument\n");
        usage (pl, argv0);
      }
      if (lpb0_arg == 0 || lpb1_arg == 0)
      {
        fprintf (stderr, "Error, missing -lpb0 and/or -lpb1 command line "
                         "argument\n");
        usage (pl, argv0);
      }
      lpb[0] = 1UL << lpb0_arg;
      lpb[1] = 1UL << lpb1_arg;
    }
    else /* lpbs were given as -lpbs x,x,x,... */
    {
      if (narg != cpoly->nb_polys)
      {
        fprintf (stderr, "Error, the number of values given in -lpbs does not "
                         "correspond to the number of polynomials\n");
        usage (pl, argv0);
      }
      if (lpb0_arg != 0)
        fprintf (stderr, "Warning, the value given by -lpb0 will be ignored, "
                         "the one given by -lpbs will be used\n");
      if (lpb1_arg != 0)
        fprintf (stderr, "Warning, the value given by -lpb1 will be ignored, "
                         "the one given by -lpbs will be used\n");
      for (int i = 0; i < cpoly->nb_polys; i++)
      {
        if (lpb_arg[i] <= 0)
        {
          fprintf (stderr, "Error, -lpbs command line argument cannot contain "
                           "non-positive values\n");
          usage (pl, argv0);
        }
        lpb[i] = 1UL << lpb_arg[i];
      }
    }

    lpb_max = lpb[0];
    for (int i = 1; i < cpoly->nb_polys; i++)
      lpb_max = MAX(lpb_max, lpb[i]);

    char ** files = filelist ? filelist_from_file(basepath, filelist, 0) : argv;

    if (outfilename)
    {
      outfile = fopen_maybe_compressed (outfilename, "w");
      printf("# Correct relations will be written in %s\n", outfilename);
    }

    printf ("# Verbose output: %s\n", verbose ? "yes" : "no");
    printf ("# Correct wrong relations if possible: %s\n",
                                              fix_it ? "yes" : "no");
    printf ("# Check primality of ideal: %s\n", check_primality ? "yes" : "no");
    for (int i = 0; i < cpoly->nb_polys; i++)
    {
      printf ("# On side %d, ", i);
      if (lpb[i] > 0)
        printf ("will check that norms of ideals are below %lu\n", lpb[i]);
      else
        printf ("will not check the size of the norms of the ideals\n");
    }

    timingstats_dict_init(stats);
    filter_rels(files, (filter_rels_callback_t) &thread_callback, (void*)outfile,
                EARLYPARSE_NEED_PRIMES |
                (abhexa ? EARLYPARSE_NEED_AB_HEXA : EARLYPARSE_NEED_AB_DECIMAL),
                NULL, stats);

    printf("Number of read relations: %" PRIu64 "\n", nrels_read);
    printf("Number of correct relations: %" PRIu64 "\n", nrels_ok);
    if (fix_it)
    {
      printf("Number of fixed relations: %" PRIu64 "\n"
             "   among which %" PRIu64 " were fully completed\n"
             "           and %" PRIu64 " contained at least 1 non-prime factor\n",
             nrels_fullyfixed, nrels_fullycompleted, nrels_fixednotprime);
      printf("Number of wrong relations: %" PRIu64 "\n", nrels_err);
      printf("   among which %" PRIu64 " had a factor not dividing the norm\n"
             "           and %" PRIu64 " could not be fully completed\n"
             "           and %" PRIu64 " contained an ideal larger than a lpb\n",
             nrels_donotdivide, nrels_notcomplete, nrels_abovelpb);
    }
    else
    {
      printf("Number of wrong relations: %" PRIu64 "\n", nrels_err);
      printf("   among which %" PRIu64 " had a factor not dividing the norm\n"
             "           and %" PRIu64 " contained at least 1 non-prime factor\n"
             "           and %" PRIu64 " were not complete\n"
             "           and %" PRIu64 " contained an ideal larger than a lpb\n",
             nrels_donotdivide, nrels_notprime, nrels_notcomplete,
             nrels_abovelpb);
    }

    if (filelist)
      filelist_clear(files);

    if (outfilename)
      fclose_maybe_compressed (outfile, outfilename);

    param_list_clear(pl);
    cado_poly_clear (cpoly);

    timingstats_dict_add_mythread(stats, "main");
    timingstats_dict_disp(stats);
    timingstats_dict_clear(stats);

    return (nrels_err) ? EXIT_FAILURE : EXIT_SUCCESS ;
}
